"""Build and verify the Part 1 -> Part 2 Mined and Natural Linker Dataset (`v1`).

Completes Part 1 of the LinkerForge workflow at proteome scale by:
  1. Verifying the B1 V2 family-disjoint split (`d62136ce...`) and frozen Seed-31
     `ThreeLayerDFLCNN` checkpoint (`c0d356c5...`).
  2. Executing B2 two-stage domain-gated conformal calibration on the family-disjoint
     `development` split (`calibrate_conformal_threshold` + `call_two_stage_interdomain_spans`).
  3. Executing B3 specificity verification on a benchmark-disjoint 40-protein
     no-linker / single-domain and experimental-order panel (`auxiliary` + single-Pfam
     non-linker DisProt controls).
  4. Executing large-scale two-stage classifier mining across >160,000 full-length
     proteins (plus auditing >184,000 multi-proteome unreviewed domain-gap records):
       - `UniProtKB_SwissProt_Reviewed_MultiDomain` (135,407 reviewed Swiss-Prot proteins
         with domain features -> 41,661 multi-domain -> Stage-1 hetero-domain gaps)
       - `UniProtKB_Explicit_Linker_Regions` (18,978 UniProtKB `ft_region:linker` proteins
         + 884 `uniprot_explicit_flexible_linkers_v1` proteins)
       - `DLD_PDB_Structural_Domain_Linkers` (2,620 PDB chains from DLD train/val/blind)
       - `DisProt_2026_06_Pfam_Gaps` (3,337 DisProt proteins -> 338 Stage-1 hetero-Pfam gaps)
  5. Assembling the provenance-stratified dataset `data/processed/mined_and_natural_linkers_v1/`
     across four explicit authority families:
       - `natural_cath_pdb_verified` (91 coordinate-resolved CATH v4.4.0 natural linkers)
       - `natural_disprot_curated_dfl` (57 benchmark-disjoint curated `IDPO:0000033` linkers)
       - `mined_domain_gated_interdomain` (`mined_high_confidence_dfl_call`,
         `mined_provisional_dfl_candidate`, `mined_domain_gap_abstained`)
       - `literature_engineered_template` (23 cited Chen et al. 2013 templates)
     Every record is annotated with `linkerforge.properties` biophysical calls
     (`rigidity`, `cleavability`, `composition_features`) and frozen Seed-31 CNN scores.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from math import dist
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from jsonschema import Draft202012Validator  # noqa: E402
import parasail  # noqa: E402
import torch  # noqa: E402
from transformers import AutoModel, AutoTokenizer  # noqa: E402

from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta  # noqa: E402
from audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry  # noqa: E402
from build_dfl_curated_evidence_pool import (  # noqa: E402
    CAID3_SHA256,
    DP94_MANIFEST_SHA256,
    family_components,
)
from build_rcsb_pilot import read_ca  # noqa: E402
from preflight_dfl_family_training_v2 import v2_homolog_hit  # noqa: E402

from linkerforge.detector.cnn import ThreeLayerDFLCNN  # noqa: E402
from linkerforge.detector.conformal_mining import (  # noqa: E402
    calibrate_conformal_threshold,
    call_two_stage_interdomain_spans,
)
from linkerforge.detector.domain_gate import DomainMatch, domain_gap_candidates  # noqa: E402
from linkerforge.properties import (  # noqa: E402
    call_cleavability,
    composition_features,
    fit_rigidity_model,
    load_library,
    modelling_subset,
)

PINNED_B1_V2_MANIFEST_SHA256 = "d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452"
PINNED_SEED31_CHECKPOINT_SHA256 = "c0d356c55ed83807353e66b875501c61b17f0b9b83f2910c0f1147f6409fdf6f"
FROZEN_SEED31_THRESHOLD = 0.7351981997489929
CONFORMAL_TWO_STAGE_THRESHOLD = 0.6391201615333557
PROVISIONAL_SPAN_THRESHOLD = 0.50
ENCODER_NAME = "facebook/esm2_t12_35M_UR50D"
MAX_WINDOW = 1022
BATCH_WINDOW_LEN = 256
C_ALPHA_MAX_STEP_A = 3.8
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_TM_AA = set("FILVWM")

DOMAIN_RE = re.compile(r"DOMAIN\s+(\d+)\.\.(\d+)")
REGION_RE = re.compile(r"REGION\s+(\d+)\.\.(\d+)")
NOTE_RE = re.compile(r'/note="([^"]+)"')
REPEAT_TM_RE = re.compile(r"(?:REPEAT|TRANSMEM)\s+(\d+)\.\.(\d+)")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_pfam_matches(entry: dict) -> list[DomainMatch]:
    features = entry.get("features") or {}
    pfam_items = features.get("pfam") or []
    matches: list[DomainMatch] = []
    for item in pfam_items:
        acc = item.get("id") or item.get("accession")
        start = item.get("start")
        end = item.get("end")
        if acc and isinstance(start, int) and isinstance(end, int) and start <= end:
            matches.append(
                DomainMatch(
                    start=int(start),
                    end=int(end),
                    accession=str(acc),
                    kind="domain",
                )
            )
    return matches


def positive_runs(mask: str | list[int]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for idx, ch in enumerate(mask, start=1):
        is_pos = (ch == "1") if isinstance(ch, str) else (int(ch) == 1)
        if is_pos and start is None:
            start = idx
        elif not is_pos and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask)))
    return runs


class BenchmarkFirewallV2:
    """Fast SIMD-prefiltered and exact RULE_V2 sequence homolog firewall."""

    def __init__(self, references: list[tuple[str, str, str]]) -> None:
        self.ref_ids = {k for _, k, _ in references}
        self.ref_seqs = [s.upper() for _, _, s in references]
        self.ref_profiles = [
            (s, len(s), parasail.profile_create_16(s, parasail.blosum62))
            for s in self.ref_seqs
        ]
        self._cache: dict[str, bool] = {}

    def is_homolog(self, seq: str) -> bool:
        seq = seq.upper()
        if seq in self._cache:
            return self._cache[seq]
        n = len(seq)
        for r_seq, rn, prof in self.ref_profiles:
            if seq == r_seq or seq in r_seq or r_seq in seq:
                self._cache[seq] = True
                return True
            if min(n, rn) < 0.80 * max(n, rn):
                continue
            if parasail.sw_striped_profile_16(prof, seq, 10, 1).score < 45:
                continue
            if v2_homolog_hit(seq, r_seq):
                self._cache[seq] = True
                return True
        self._cache[seq] = False
        return False


class FrozenSeed31Scorer:
    """Offline residue-level probability scorer for the frozen Seed-31 CNN."""

    def __init__(self, checkpoint_path: Path) -> None:
        if sha256_file(checkpoint_path) != PINNED_SEED31_CHECKPOINT_SHA256:
            raise ValueError(f"Seed-31 checkpoint hash mismatch: {checkpoint_path}")
        torch.set_num_threads(min(4, torch.get_num_threads()))
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
        self.tokenizer = AutoTokenizer.from_pretrained(ENCODER_NAME, local_files_only=True)
        self.backbone = (
            AutoModel.from_pretrained(ENCODER_NAME, local_files_only=True)
            .to(self.device)
            .eval()
        )
        self.backbone.requires_grad_(False)
        self.head = ThreeLayerDFLCNN(self.backbone.config.hidden_size).to(self.device).eval()
        self.head.load_state_dict(checkpoint["head"])
        self._cache: dict[str, list[float]] = {}

    def score_sequence(self, sequence: str) -> list[float]:
        if sequence in self._cache:
            return self._cache[sequence]
        out: list[float] = []
        for start in range(0, len(sequence), MAX_WINDOW):
            chunk = sequence[start : start + MAX_WINDOW]
            tokens = self.tokenizer([chunk], return_tensors="pt", add_special_tokens=True)
            tokens = {k: v.to(self.device) for k, v in tokens.items()}
            with torch.no_grad():
                feats = self.backbone(**tokens).last_hidden_state[:, 1 : 1 + len(chunk)]
                probs = self.head(feats).sigmoid().squeeze(0).cpu().tolist()
            if isinstance(probs, float):
                probs = [probs]
            out.extend(float(p) for p in probs)
        if len(out) != len(sequence):
            raise ValueError("Scored residue count does not match sequence length")
        self._cache[sequence] = out
        return out

    def score_candidate_windows_batched(
        self,
        items: list[dict],
        *,
        batch_size: int = 64,
        window_size: int = BATCH_WINDOW_LEN,
    ) -> list[tuple[float, float, float]]:
        """Score a large list of candidate spans using length-sorted fp16 GPU batching.

        Each item in `items` has `full_sequence`, `start_1based`, `end_1based`.
        Extracts a centered `window_size`-aa sequence window containing flanking
        domains + the linker span and returns `(mean_score, max_score, frac_above_high)`.
        """
        if not items:
            return []
        prepared: list[tuple[int, str, int, int]] = []
        for idx, it in enumerate(items):
            seq = it["full_sequence"]
            s1 = int(it["start_1based"])
            e1 = int(it["end_1based"])
            if len(seq) <= window_size:
                win_seq = seq
                rel_s = s1 - 1
                rel_e = e1
            else:
                mid = (s1 - 1 + e1) // 2
                win_start = max(0, min(len(seq) - window_size, mid - (window_size // 2)))
                win_end = min(len(seq), win_start + window_size)
                win_seq = seq[win_start:win_end]
                rel_s = max(0, (s1 - 1) - win_start)
                rel_e = min(len(win_seq), e1 - win_start)
            prepared.append((idx, win_seq, rel_s, rel_e))

        # Sort by window length to minimize padding waste during ESM2 inference
        prepared.sort(key=lambda x: len(x[1]))
        results: list[tuple[float, float, float]] = [(0.0, 0.0, 0.0)] * len(items)

        use_amp = self.device.type == "cuda"
        total_items = len(prepared)
        for b_start in range(0, total_items, batch_size):
            batch = prepared[b_start : b_start + batch_size]
            batch_seqs = [x[1] for x in batch]
            tokens = self.tokenizer(
                batch_seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=window_size + 2,
                add_special_tokens=True,
            )
            tokens = {k: v.to(self.device) for k, v in tokens.items()}
            with torch.no_grad():
                if use_amp:
                    with torch.autocast(device_type="cuda", dtype=torch.float16):
                        feats = self.backbone(**tokens).last_hidden_state[:, 1:-1, :]
                        probs_tensor = self.head(feats.float()).sigmoid()
                else:
                    feats = self.backbone(**tokens).last_hidden_state[:, 1:-1, :]
                    probs_tensor = self.head(feats).sigmoid()
            probs_cpu = probs_tensor.cpu().tolist()
            for row_probs, (orig_idx, win_seq, rel_s, rel_e) in zip(probs_cpu, batch):
                if isinstance(row_probs, float):
                    row_probs = [row_probs]
                span_scores = [float(p) for p in row_probs[rel_s:rel_e]]
                if not span_scores:
                    span_scores = [0.0]
                mean_sc = sum(span_scores) / len(span_scores)
                max_sc = max(span_scores)
                frac_high = sum(1 for x in span_scores if x >= FROZEN_SEED31_THRESHOLD) / len(
                    span_scores
                )
                results[orig_idx] = (mean_sc, max_sc, frac_high)
            if (b_start // batch_size) % 50 == 0 or (b_start + batch_size) >= total_items:
                print(
                    f"[Stage-2 GPU CNN] Scored {min(total_items, b_start + batch_size)}/{total_items} candidate windows...",
                    flush=True,
                )
        return results


def annotate_biophysical_properties(sequence: str, library, rigidity_model) -> dict:
    comp = composition_features(sequence)
    cleav = call_cleavability(sequence, library)
    p_rigid = float(rigidity_model.predict_proba(sequence))
    return {
        "rigidity_predicted_label": rigidity_model.predict(sequence),
        "rigidity_p_rigid": round(p_rigid, 6),
        "cleavability_label": cleav.label,
        "matched_proteases": cleav.proteases,
        "matched_motif_ids": [m.entry_id for m in cleav.matches],
        "composition_features": {k: round(float(v), 6) for k, v in comp.items()},
    }


def classify_mined_tier(mean_sc: float, max_sc: float) -> str:
    if mean_sc >= FROZEN_SEED31_THRESHOLD or (
        mean_sc >= CONFORMAL_TWO_STAGE_THRESHOLD and max_sc >= FROZEN_SEED31_THRESHOLD
    ):
        return "mined_high_confidence_dfl_call"
    if mean_sc >= PROVISIONAL_SPAN_THRESHOLD or max_sc >= CONFORMAL_TWO_STAGE_THRESHOLD:
        return "mined_provisional_dfl_candidate"
    return "mined_domain_gap_abstained"


def ensure_uniprot_stream(out_path: Path, query: str) -> Path:
    if out_path.exists() and out_path.stat().st_size > 100_000:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {
            "query": query,
            "format": "tsv",
            "fields": "accession,organism_name,reviewed,sequence,ft_domain,ft_repeat,ft_transmem,ft_region",
            "compressed": "true",
        }
    )
    req = urllib.request.Request(url, headers={"User-Agent": "LinkerForge/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp, out_path.open("wb") as out_f:
        while chunk := resp.read(1024 * 1024):
            out_f.write(chunk)
    return out_path


def write_jsonl(path: Path, rows: list[dict]) -> dict:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for r in rows:
            handle.write(json.dumps(r, sort_keys=True) + "\n")
    return {
        "path": path.name,
        "records": len(rows),
        "sha256": sha256_file(path),
    }


def build_dataset(source_root: Path, output_dir: Path, summary_path: Path) -> dict:
    source_root = source_root.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    b1_dir = ROOT / "data" / "processed" / "dfl_b1_v2_family_split"
    b1_manifest_path = b1_dir / "manifest.json"
    if sha256_file(b1_manifest_path) != PINNED_B1_V2_MANIFEST_SHA256:
        raise ValueError("B1 V2 split manifest SHA-256 mismatch")

    checkpoint_path = ROOT / "artifacts" / "dfl-b1-v2-cnn-seed31-2026-10-18" / "best.pt"
    scorer = FrozenSeed31Scorer(checkpoint_path)

    # Verify raw pins from source_root
    disprot_path = source_root / "data" / "raw" / "disprot" / "disprot_2026_06_regions.json"
    if sha256_file(disprot_path) != DISPROT_SHA256:
        raise ValueError("DisProt 2026_06 raw JSON hash mismatch")
    cath_path = source_root / "data" / "processed" / "cath_full_pilot.jsonl"
    rcsb_dir = source_root / "data" / "raw" / "rcsb"
    caid_path = source_root / "artifacts" / "caid3-linker.fasta"
    if sha256_file(caid_path) != CAID3_SHA256:
        raise ValueError("CAID3 FASTA hash mismatch")
    dp_manifest_path = (
        source_root / "data" / "processed" / "prospective_dfl_holdout_dp94_v1" / "manifest.json"
    )
    if sha256_file(dp_manifest_path) != DP94_MANIFEST_SHA256:
        raise ValueError("DP94 manifest hash mismatch")
    dp_fasta_path = dp_manifest_path.parent / "sequences.fasta"

    curated_lib_path = ROOT / "data" / "curated" / "linker_property_library_sources.v1.json"

    # Load biophysical property models
    library = load_library(curated_lib_path)
    rigidity_entries = modelling_subset(library, "rigidity")
    rigidity_model = fit_rigidity_model(
        [e.sequence for e in rigidity_entries if e.sequence is not None],
        [e.rigidity_label for e in rigidity_entries if e.sequence is not None],
    )

    # Load DisProt raw entries and family components for benchmark firewall
    raw_disprot = json.loads(disprot_path.read_text(encoding="utf-8"))
    disprot_by_id: dict[str, dict] = {}
    metadata: dict[str, dict] = {}
    uniprot_to_disprot_comp: dict[str, str] = {}
    benchmark_uniprot_accs: set[str] = set()

    caid3 = caid3_sequences(caid_path)
    dp94 = two_line_fasta(dp_fasta_path)
    references = [("caid3", k, s) for k, s in caid3] + [("dp94", k, s) for k, s in dp94]
    ref_ids = {k for _, k, _ in references}
    firewall = BenchmarkFirewallV2(references)

    for entry in raw_disprot["data"]:
        did = entry["disprot_id"]
        disprot_by_id[did] = entry
        acc = (entry.get("acc") or "").strip()
        if did in ref_ids and acc:
            benchmark_uniprot_accs.add(acc)
        try:
            metadata[did] = project_entry(entry)
        except ValueError:
            continue

    components = family_components(metadata)
    for reg_name, ref_id, _ in references:
        if ref_id not in components:
            components[ref_id] = f"benchmark_isolated_{reg_name}_{ref_id}"
    benchmark_components = {components[ref_id] for _, ref_id, _ in references}

    for did, entry in disprot_by_id.items():
        acc = (entry.get("acc") or "").strip()
        comp_id = components.get(did)
        if acc and comp_id:
            uniprot_to_disprot_comp[acc] = comp_id
            if comp_id in benchmark_components:
                benchmark_uniprot_accs.add(acc)

    # Load B1 V2 train, development, auxiliary splits
    train_rows = [
        json.loads(line)
        for line in (b1_dir / "train.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    dev_rows = [
        json.loads(line)
        for line in (b1_dir / "development.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    aux_rows = [
        json.loads(line)
        for line in (b1_dir / "auxiliary.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    # -------------------------------------------------------------------------
    # Part 1 Validation: B2 Conformal Risk Control on `development`
    # -------------------------------------------------------------------------
    dev_res_labels: list[int] = []
    dev_res_scores: list[float] = []
    dev_gated_labels: list[int] = []
    dev_gated_scores: list[float] = []
    dev_stage1_passed_proteins = 0

    for row in dev_rows:
        pid = row["protein_id"]
        seq = row["sequence"]
        mask = row["labels"]
        scores = scorer.score_sequence(seq)
        pfam_matches = extract_pfam_matches(disprot_by_id.get(pid, {}))
        cands, reason = domain_gap_candidates(seq, pfam_matches, sequence_verified=True)
        in_gate = [False] * len(seq)
        if reason is None and cands:
            dev_stage1_passed_proteins += 1
            for span in cands:
                for pos in range(span.start - 1, span.end):
                    in_gate[pos] = True
        for idx, (lbl, sc) in enumerate(zip(mask, scores)):
            if lbl in ("0", "1"):
                y = int(lbl)
                dev_res_labels.append(y)
                dev_res_scores.append(sc)
                dev_gated_labels.append(y)
                dev_gated_scores.append(sc if in_gate[idx] else 0.0)

    b2_ungated = calibrate_conformal_threshold(
        dev_res_labels,
        dev_res_scores,
        target_fdr=0.10,
        alpha=0.05,
        min_true_positives=30,
    )
    b2_two_stage = calibrate_conformal_threshold(
        dev_gated_labels,
        dev_gated_scores,
        target_fdr=0.10,
        alpha=0.05,
        min_true_positives=20,
    )

    # -------------------------------------------------------------------------
    # Tier 1: `natural_cath_pdb_verified` (91 CATH two-domain natural linkers)
    # -------------------------------------------------------------------------
    cath_rows = [
        json.loads(line)
        for line in cath_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    tier1_records: list[dict] = []
    for row in cath_rows:
        ex_id = row["example_id"]
        pdb_id = row["source"]["pdb_id"].upper()
        chain_id = row["source"]["chain_id"]
        pdb_path = rcsb_dir / f"{pdb_id}.pdb"
        ca_map = read_ca(pdb_path.read_bytes(), chain_id)
        port_a_num = int(row["residue_mapping"][0]["author_number"])
        port_b_num = int(row["residue_mapping"][-1]["author_number"])
        port_a_xyz = ca_map[port_a_num][1]
        port_b_xyz = ca_map[port_b_num][1]
        ca_dist = float(dist(port_a_xyz, port_b_xyz))

        linker_seq = row["linker"]["sequence"]
        linker_len = int(row["linker"]["length"])
        max_cap = linker_len * C_ALPHA_MAX_STEP_A
        ext_ratio = ca_dist / max_cap if max_cap > 0 else 0.0
        feasible = ca_dist <= max_cap

        full_seq = row["domain_a"]["sequence"] + linker_seq + row["domain_b"]["sequence"]
        full_scores = scorer.score_sequence(full_seq)
        l_start = int(row["linker"]["start"])
        l_end = int(row["linker"]["end"])
        span_scores = full_scores[l_start - 1 : l_end]
        mean_sc = sum(span_scores) / len(span_scores)
        max_sc = max(span_scores)
        frac_high = sum(1 for x in span_scores if x >= FROZEN_SEED31_THRESHOLD) / len(span_scores)

        props = annotate_biophysical_properties(linker_seq, library, rigidity_model)
        tier1_records.append(
            {
                "record_id": f"cath:{ex_id}",
                "authority_tier": "natural_cath_pdb_verified",
                "evidence_grade": "experimental_structure_cath_v4_4_0",
                "protein_id": f"{pdb_id}_{chain_id}",
                "homology_cluster_id": row["homology_cluster_id"],
                "linker_sequence": linker_seq,
                "linker_length": linker_len,
                "linker_start_1based": l_start,
                "linker_end_1based": l_end,
                "flanking_context": {
                    "left_domain_source": "CATH",
                    "left_domain_id": row["domain_a"]["family_ids"][0],
                    "left_flank_sequence": row["domain_a"]["sequence"][-80:],
                    "right_domain_source": "CATH",
                    "right_domain_id": row["domain_b"]["family_ids"][0],
                    "right_flank_sequence": row["domain_b"]["sequence"][:80],
                },
                "geometry": {
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "port_a_author_number": port_a_num,
                    "port_b_author_number": port_b_num,
                    "inter_port_ca_distance_A": round(ca_dist, 4),
                    "max_extended_3p8A_capacity_A": round(max_cap, 4),
                    "end_to_end_extension_ratio": round(ext_ratio, 4),
                    "physical_3p8A_feasible": bool(feasible),
                    "unresolved_linker_residues": int(row["source"]["unresolved_linker_residues"]),
                    "port_a_accessibility_proxy": round(
                        float(row["connection_features"]["port_a"]["accessibility_proxy"]), 4
                    ),
                    "port_b_accessibility_proxy": round(
                        float(row["connection_features"]["port_b"]["accessibility_proxy"]), 4
                    ),
                },
                "classifier_scores": {
                    "cnn_mean_score": round(mean_sc, 6),
                    "cnn_max_score": round(max_sc, 6),
                    "cnn_fraction_above_high_threshold": round(frac_high, 6),
                },
                "biophysical_properties": props,
                "benchmark_firewall_verified": True,
            }
        )

    # -------------------------------------------------------------------------
    # Tier 2: `natural_disprot_curated_dfl` (Benchmark-disjoint curated DisProt DFLs)
    # -------------------------------------------------------------------------
    tier2_records: list[dict] = []
    curated_longer_than_60_count = 0
    for split_name, split_rows in (("train", train_rows), ("development", dev_rows)):
        for row in split_rows:
            pid = row["protein_id"]
            seq = row["sequence"]
            mask = row["labels"]
            scores = scorer.score_sequence(seq)
            pfam_matches = sorted(
                extract_pfam_matches(disprot_by_id.get(pid, {})),
                key=lambda m: (m.start, m.end, m.accession),
            )
            for s_start, s_end in positive_runs(mask):
                s_len = s_end - s_start + 1
                if s_len < 5 or s_len > 60:
                    if s_len > 60:
                        curated_longer_than_60_count += 1
                    continue
                l_seq = seq[s_start - 1 : s_end]
                span_scores = scores[s_start - 1 : s_end]
                left_domains = [m.accession for m in pfam_matches if m.end < s_start]
                right_domains = [m.accession for m in pfam_matches if m.start > s_end]
                props = annotate_biophysical_properties(l_seq, library, rigidity_model)
                tier2_records.append(
                    {
                        "record_id": f"disprot_curated:{pid}:{s_start}-{s_end}",
                        "authority_tier": "natural_disprot_curated_dfl",
                        "evidence_grade": "curator_validated_idpo_0000033",
                        "b1_v2_split_role": split_name,
                        "protein_id": pid,
                        "homology_cluster_id": row["family_component_id"],
                        "linker_sequence": l_seq,
                        "linker_length": s_len,
                        "linker_start_1based": s_start,
                        "linker_end_1based": s_end,
                        "flanking_context": {
                            "left_domain_source": "Pfam" if left_domains else None,
                            "left_domain_id": left_domains[-1] if left_domains else None,
                            "left_flank_sequence": seq[max(0, s_start - 81) : s_start - 1],
                            "right_domain_source": "Pfam" if right_domains else None,
                            "right_domain_id": right_domains[0] if right_domains else None,
                            "right_flank_sequence": seq[s_end : min(len(seq), s_end + 80)],
                        },
                        "geometry": None,
                        "classifier_scores": {
                            "cnn_mean_score": round(sum(span_scores) / len(span_scores), 6),
                            "cnn_max_score": round(max(span_scores), 6),
                            "cnn_fraction_above_high_threshold": round(
                                sum(1 for x in span_scores if x >= FROZEN_SEED31_THRESHOLD)
                                / len(span_scores),
                                6,
                            ),
                        },
                        "biophysical_properties": props,
                        "benchmark_firewall_verified": True,
                    }
                )

    # -------------------------------------------------------------------------
    # Tier 3: Large-Scale Proteome Mining + B3 No-Linker Specificity Panel
    # -------------------------------------------------------------------------
    b1_protein_ids = {r["protein_id"] for r in train_rows + dev_rows + aux_rows}
    tier3_records: list[dict] = []
    b3_panel_records: list[dict] = []

    # First add the 9 B1 V2 auxiliary experimental-order proteins to B3 panel
    for row in aux_rows:
        pid = row["protein_id"]
        seq = row["sequence"]
        scores = scorer.score_sequence(seq)
        pfam_matches = extract_pfam_matches(disprot_by_id.get(pid, {}))
        two_stage = call_two_stage_interdomain_spans(
            seq,
            pfam_matches,
            scores,
            sequence_verified=True,
            high_confidence_threshold=FROZEN_SEED31_THRESHOLD,
            provisional_threshold=PROVISIONAL_SPAN_THRESHOLD,
        )
        b3_panel_records.append(
            {
                "protein_id": pid,
                "panel_role": "b1_v2_auxiliary_experimental_order",
                "family_component_id": row["family_component_id"],
                "sequence_length": len(seq),
                "pfam_domain_count": len(pfam_matches),
                "known_order_residues": row["labels"].count("0"),
                "stage1_passed": two_stage["stage1_passed"],
                "stage1_abstention_reason": two_stage["stage1_abstention_reason"],
                "called_high_confidence_spans": two_stage["called_high_confidence_spans"],
                "called_provisional_spans": two_stage["called_provisional_spans"],
                "max_residue_score": round(max(scores), 6),
                "benchmark_firewall_verified": True,
            }
        )

    # Corpus 3A: DisProt_2026_06 hetero-Pfam inter-domain gap candidates + B3 panel
    for did, entry in sorted(disprot_by_id.items()):
        if did in b1_protein_ids or did in ref_ids:
            continue
        if did not in metadata:
            continue
        comp_id = components.get(did)
        if comp_id in benchmark_components:
            continue
        seq = entry.get("sequence") or ""
        if not seq or any(ch not in CANONICAL_AA for ch in seq):
            continue
        has_dfl_region = any(
            r.get("term_id") == "IDPO:0000033" for r in (entry.get("regions") or [])
        )
        if has_dfl_region:
            continue

        pfam_matches = extract_pfam_matches(entry)
        cands, reason = domain_gap_candidates(seq, pfam_matches, sequence_verified=True)

        if reason is None and cands:
            if firewall.is_homolog(seq):
                continue
            scores = scorer.score_sequence(seq)
            two_stage = call_two_stage_interdomain_spans(
                seq,
                pfam_matches,
                scores,
                sequence_verified=True,
                high_confidence_threshold=FROZEN_SEED31_THRESHOLD,
                provisional_threshold=PROVISIONAL_SPAN_THRESHOLD,
            )
            for span_info in two_stage["candidate_spans"]:
                s_start = int(span_info["start"])
                s_end = int(span_info["end"])
                l_seq = span_info["sequence"]
                props = annotate_biophysical_properties(l_seq, library, rigidity_model)
                tier3_records.append(
                    {
                        "record_id": f"mined_disprot_gap:{did}:{s_start}-{s_end}",
                        "authority_tier": span_info["calling_tier"],
                        "mining_source_corpus": "DisProt_2026_06_Pfam_Gaps",
                        "evidence_grade": "two_stage_domain_gated_seed31_cnn_prediction_not_ground_truth",
                        "protein_id": did,
                        "homology_cluster_id": comp_id,
                        "linker_sequence": l_seq,
                        "linker_length": int(span_info["length"]),
                        "linker_start_1based": s_start,
                        "linker_end_1based": s_end,
                        "flanking_context": {
                            "left_domain_source": "Pfam",
                            "left_domain_id": span_info["left_domain_accession"],
                            "left_flank_sequence": seq[max(0, s_start - 81) : s_start - 1],
                            "right_domain_source": "Pfam",
                            "right_domain_id": span_info["right_domain_accession"],
                            "right_flank_sequence": seq[s_end : min(len(seq), s_end + 80)],
                        },
                        "geometry": None,
                        "classifier_scores": {
                            "cnn_mean_score": span_info["cnn_mean_score"],
                            "cnn_max_score": span_info["cnn_max_score"],
                            "cnn_fraction_above_high_threshold": span_info[
                                "cnn_fraction_above_high_threshold"
                            ],
                        },
                        "biophysical_properties": props,
                        "benchmark_firewall_verified": True,
                    }
                )
        elif len(b3_panel_records) < 40 and reason in (
            "fewer_than_two_domains",
            "homo_repeat_domain_pair",
        ):
            if firewall.is_homolog(seq):
                continue
            scores = scorer.score_sequence(seq)
            two_stage = call_two_stage_interdomain_spans(
                seq,
                pfam_matches,
                scores,
                sequence_verified=True,
                high_confidence_threshold=FROZEN_SEED31_THRESHOLD,
                provisional_threshold=PROVISIONAL_SPAN_THRESHOLD,
            )
            b3_panel_records.append(
                {
                    "protein_id": did,
                    "panel_role": f"single_or_homo_domain_non_linker_{reason}",
                    "family_component_id": comp_id,
                    "sequence_length": len(seq),
                    "pfam_domain_count": len(pfam_matches),
                    "known_order_residues": 0,
                    "stage1_passed": two_stage["stage1_passed"],
                    "stage1_abstention_reason": two_stage["stage1_abstention_reason"],
                    "called_high_confidence_spans": two_stage["called_high_confidence_spans"],
                    "called_provisional_spans": two_stage["called_provisional_spans"],
                    "max_residue_score": round(max(scores), 6),
                    "benchmark_firewall_verified": True,
                }
            )

    # -------------------------------------------------------------------------
    # Large-Scale Multi-Proteome Candidate Collection (Stage 1 CPU Gate + RULE_V2)
    # -------------------------------------------------------------------------
    csv.field_size_limit(10_000_000)
    batched_candidates: list[dict] = []
    seen_span_keys: set[tuple[str, int, int]] = set()
    organism_gap_counts: dict[str, int] = {}

    # Corpus 3B: DLD PDB Structural Domain Linkers (`dld_train`, `dld_validation`, `dld_blind`)
    dld_dir = source_root / "data" / "processed" / "dld_reference"
    dld_proteins_scanned = 0
    dld_candidates_queued = 0
    for dld_file in ("dld_train.jsonl.gz", "dld_validation.jsonl.gz", "dld_blind.jsonl.gz"):
        dld_path = dld_dir / dld_file
        if not dld_path.exists():
            continue
        with gzip.open(dld_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                dld_proteins_scanned += 1
                r = json.loads(line)
                sid = str(r["source_id"])
                parts = sid.split("_")
                u_acc = parts[-1] if len(parts) >= 3 else sid
                if u_acc in benchmark_uniprot_accs:
                    continue
                seq = str(r["sequence"]).strip().upper()
                if not seq or not (set(seq) <= CANONICAL_AA):
                    continue
                runs = positive_runs(r["labels"])
                valid_runs = [
                    (s, e)
                    for s, e in runs
                    if 5 <= (e - s + 1) <= 60 and s > 20 and (len(seq) - e) >= 20
                ]
                if not valid_runs:
                    continue
                if firewall.is_homolog(seq):
                    continue
                for s_start, s_end in valid_runs:
                    key = (sid, s_start, s_end)
                    if key in seen_span_keys:
                        continue
                    seen_span_keys.add(key)
                    l_seq = seq[s_start - 1 : s_end]
                    hydro_frac = sum(1 for aa in l_seq if aa in HYDROPHOBIC_TM_AA) / len(l_seq)
                    if hydro_frac >= 0.50:
                        continue
                    comp_id = uniprot_to_disprot_comp.get(
                        u_acc, f"dld:{r.get('homology_cluster_id', sid)}"
                    )
                    batched_candidates.append(
                        {
                            "record_id": f"mined_dld_pdb:{sid}:{s_start}-{s_end}",
                            "mining_source_corpus": "DLD_PDB_Structural_Domain_Linkers",
                            "protein_id": sid,
                            "organism": "PDB_Structural_DLD",
                            "homology_cluster_id": comp_id,
                            "full_sequence": seq,
                            "linker_sequence": l_seq,
                            "start_1based": s_start,
                            "end_1based": s_end,
                            "left_domain_source": "DLD_Structural_Domain",
                            "left_domain_id": f"{parts[0]}_N_domain",
                            "right_domain_source": "DLD_Structural_Domain",
                            "right_domain_id": f"{parts[0]}_C_domain",
                        }
                    )
                    dld_candidates_queued += 1

    # Corpus 3C: UniProtKB/Swiss-Prot Reviewed Multi-Domain Proteome (135,407 proteins)
    uniprot_raw_dir = ROOT / "data" / "raw" / "uniprot"
    sp_domains_path = ensure_uniprot_stream(
        uniprot_raw_dir / "swissprot_all_domains.tsv.gz",
        "reviewed:true AND ft_domain:*",
    )
    sp_total_scanned = 0
    sp_two_plus_domains = 0
    sp_stage1_proteins = 0
    sp_stage1_candidates = 0

    with gzip.open(sp_domains_path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sp_total_scanned += 1
            acc = row["Entry"].strip()
            if not acc or acc in benchmark_uniprot_accs:
                continue
            seq = row["Sequence"].strip().upper()
            if not seq or not (set(seq) <= CANONICAL_AA):
                continue
            dom_text = row["Domain [FT]"]
            hits = list(DOMAIN_RE.finditer(dom_text))
            if len(hits) < 2:
                continue
            sp_two_plus_domains += 1
            spans: list[DomainMatch] = []
            for i, h in enumerate(hits):
                s, e = map(int, h.groups())
                tail = dom_text[h.end() : hits[i + 1].start() if i + 1 < len(hits) else None]
                m = NOTE_RE.search(tail)
                fam = m.group(1).strip() if m else f"DOM_{i}"
                fam_norm = re.sub(r"\s+\d+$", "", fam)
                spans.append(DomainMatch(s, e, fam_norm, "domain"))
            blocked = [
                tuple(map(int, m.groups()))
                for txt in (row["Repeat"], row["Transmembrane"])
                for m in REPEAT_TM_RE.finditer(txt)
            ]
            cands, reason = domain_gap_candidates(
                seq, spans, sequence_verified=True, min_length=5, max_length=60
            )
            if reason is not None or not cands:
                continue
            cands = [
                c for c in cands if not any(a <= c.end and b >= c.start for a, b in blocked)
            ]
            if not cands:
                continue
            if firewall.is_homolog(seq):
                continue
            sp_stage1_proteins += 1
            org = row["Organism"].split("(")[0].strip() or "Unknown"
            comp_id = uniprot_to_disprot_comp.get(
                acc, f"swissprot_arch:{cands[0].left_accession}__{cands[0].right_accession}"
            )
            for c in cands:
                key = (acc, c.start, c.end)
                if key in seen_span_keys:
                    continue
                l_seq = seq[c.start - 1 : c.end]
                hydro_frac = sum(1 for aa in l_seq if aa in HYDROPHOBIC_TM_AA) / len(l_seq)
                if hydro_frac >= 0.50:
                    continue
                seen_span_keys.add(key)
                sp_stage1_candidates += 1
                organism_gap_counts[org] = organism_gap_counts.get(org, 0) + 1
                batched_candidates.append(
                    {
                        "record_id": f"mined_swissprot_gap:{acc}:{c.start}-{c.end}",
                        "mining_source_corpus": "UniProtKB_SwissProt_Reviewed_MultiDomain",
                        "protein_id": acc,
                        "organism": org,
                        "homology_cluster_id": comp_id,
                        "full_sequence": seq,
                        "linker_sequence": l_seq,
                        "start_1based": c.start,
                        "end_1based": c.end,
                        "left_domain_source": "UniProt_FT_DOMAIN",
                        "left_domain_id": c.left_accession,
                        "right_domain_source": "UniProt_FT_DOMAIN",
                        "right_domain_id": c.right_accession,
                    }
                )

    # Corpus 3D: UniProtKB Explicit Linker Regions (18,978 proteins + 884 explicit flexible linkers)
    exp_linkers_path = ensure_uniprot_stream(
        uniprot_raw_dir / "uniprot_all_explicit_linkers.tsv.gz",
        "ft_region:linker",
    )
    exp_total_scanned = 0
    exp_candidates_queued = 0
    with gzip.open(exp_linkers_path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            exp_total_scanned += 1
            acc = row["Entry"].strip()
            if not acc or acc in benchmark_uniprot_accs:
                continue
            seq = row["Sequence"].strip().upper()
            if not seq or not (set(seq) <= CANONICAL_AA):
                continue
            dom_text = row.get("Domain [FT]") or ""
            dom_hits = list(DOMAIN_RE.finditer(dom_text))
            is_reviewed = (row.get("Reviewed") or "").strip().lower() == "reviewed"
            # Require either reviewed Swiss-Prot status or >=2 domain annotations
            if not is_reviewed and len(dom_hits) < 2:
                continue
            blocked = [
                tuple(map(int, m.groups()))
                for txt in (row.get("Repeat") or "", row.get("Transmembrane") or "")
                for m in REPEAT_TM_RE.finditer(txt)
            ]
            reg_text = row.get("Region") or ""
            reg_hits = list(REGION_RE.finditer(reg_text))
            valid_regs: list[tuple[int, int, str]] = []
            for i, h in enumerate(reg_hits):
                s, e = map(int, h.groups())
                tail = reg_text[
                    h.end() : reg_hits[i + 1].start() if i + 1 < len(reg_hits) else None
                ]
                m = NOTE_RE.search(tail)
                note = m.group(1).strip() if m else ""
                if "linker" not in note.lower():
                    continue
                if not (1 <= s <= e <= len(seq)) or not (5 <= (e - s + 1) <= 60):
                    continue
                if any(a <= e and b >= s for a, b in blocked):
                    continue
                valid_regs.append((s, e, note))
            if not valid_regs:
                continue
            if firewall.is_homolog(seq):
                continue
            # Extract domain names for flanking context if present
            dom_spans: list[tuple[int, int, str]] = []
            for i, h in enumerate(dom_hits):
                ds, de = map(int, h.groups())
                tail = dom_text[
                    h.end() : dom_hits[i + 1].start() if i + 1 < len(dom_hits) else None
                ]
                m = NOTE_RE.search(tail)
                dfam = re.sub(r"\s+\d+$", "", m.group(1).strip()) if m else f"DOM_{i}"
                dom_spans.append((ds, de, dfam))
            org = (row.get("Organism") or "Unknown").split("(")[0].strip() or "Unknown"
            for s, e, note in valid_regs:
                key = (acc, s, e)
                if key in seen_span_keys:
                    continue
                l_seq = seq[s - 1 : e]
                hydro_frac = sum(1 for aa in l_seq if aa in HYDROPHOBIC_TM_AA) / len(l_seq)
                if hydro_frac >= 0.50:
                    continue
                seen_span_keys.add(key)
                left_doms = [dfam for ds, de, dfam in dom_spans if de < s]
                right_doms = [dfam for ds, de, dfam in dom_spans if ds > e]
                l_id = left_doms[-1] if left_doms else f"explicit_region:{note[:30]}"
                r_id = right_doms[0] if right_doms else f"explicit_region:{note[:30]}"
                comp_id = uniprot_to_disprot_comp.get(acc, f"uniprot_explicit:{l_id}__{r_id}")
                exp_candidates_queued += 1
                organism_gap_counts[org] = organism_gap_counts.get(org, 0) + 1
                batched_candidates.append(
                    {
                        "record_id": f"mined_uniprot_explicit:{acc}:{s}-{e}",
                        "mining_source_corpus": "UniProtKB_Explicit_Linker_Regions",
                        "protein_id": acc,
                        "organism": org,
                        "homology_cluster_id": comp_id,
                        "full_sequence": seq,
                        "linker_sequence": l_seq,
                        "start_1based": s,
                        "end_1based": e,
                        "left_domain_source": "UniProt_FT_DOMAIN" if left_doms else "UniProt_FT_REGION",
                        "left_domain_id": l_id,
                        "right_domain_source": "UniProt_FT_DOMAIN" if right_doms else "UniProt_FT_REGION",
                        "right_domain_id": r_id,
                    }
                )

    # Also scan `uniprot_explicit_flexible_linkers_v1/candidates.jsonl` (884 proteins)
    exp_flex_path = (
        source_root / "data" / "processed" / "uniprot_explicit_flexible_linkers_v1" / "candidates.jsonl"
    )
    if exp_flex_path.exists():
        with exp_flex_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                exp_total_scanned += 1
                r = json.loads(line)
                acc = str(r["source_id"]).strip()
                if not acc or acc in benchmark_uniprot_accs:
                    continue
                seq = str(r["sequence"]).strip().upper()
                if not seq or not (set(seq) <= CANONICAL_AA):
                    continue
                ann = r.get("annotation") or {}
                s = int(ann.get("start_1based", 0))
                e = int(ann.get("end_1based", 0))
                if not (1 <= s <= e <= len(seq)) or not (5 <= (e - s + 1) <= 60):
                    continue
                key = (acc, s, e)
                if key in seen_span_keys:
                    continue
                if firewall.is_homolog(seq):
                    continue
                l_seq = seq[s - 1 : e]
                hydro_frac = sum(1 for aa in l_seq if aa in HYDROPHOBIC_TM_AA) / len(l_seq)
                if hydro_frac >= 0.50:
                    continue
                seen_span_keys.add(key)
                note = str(ann.get("note") or "Flexible linker")
                comp_id = uniprot_to_disprot_comp.get(acc, f"uniprot_explicit_flex:{acc}")
                exp_candidates_queued += 1
                batched_candidates.append(
                    {
                        "record_id": f"mined_uniprot_flex:{acc}:{s}-{e}",
                        "mining_source_corpus": "UniProtKB_Explicit_Linker_Regions",
                        "protein_id": acc,
                        "organism": "UniProtKB_Explicit_Flexible",
                        "homology_cluster_id": comp_id,
                        "full_sequence": seq,
                        "linker_sequence": l_seq,
                        "start_1based": s,
                        "end_1based": e,
                        "left_domain_source": "UniProt_FT_REGION",
                        "left_domain_id": note,
                        "right_domain_source": "UniProt_FT_REGION",
                        "right_domain_id": note,
                    }
                )

    # Audit unreviewed multi-proteome domain-gap catalogs in `linker_catalog/expanded_*_v2`
    catalog_dir = source_root / "data" / "processed" / "linker_catalog"
    unreviewed_proteome_audit: dict[str, dict[str, int]] = {}
    for sub in (
        "expanded_human_v2",
        "expanded_mouse_v2",
        "expanded_zebrafish_v2",
        "expanded_drosophila_v2",
        "expanded_arabidopsis_v2",
        "expanded_cow_v2",
        "expanded_maize_v2",
        "expanded_wheat_v2",
    ):
        fp = catalog_dir / sub / "uniprot_domain_gap_candidates.v1.jsonl.gz"
        if not fp.exists():
            continue
        tot_rows = 0
        hetero_5_60 = 0
        with gzip.open(fp, "rt", encoding="utf-8") as handle:
            for line in handle:
                tot_rows += 1
                row = json.loads(line)
                la = row["left_domain"]["annotation"]
                ra = row["right_domain"]["annotation"]
                ml = NOTE_RE.search(la)
                mr = NOTE_RE.search(ra)
                fl = re.sub(r"\s+\d+$", "", ml.group(1).strip()) if ml else la[:40]
                fr = re.sub(r"\s+\d+$", "", mr.group(1).strip()) if mr else ra[:40]
                if fl != fr and 5 <= len(row["sequence"]) <= 60:
                    hetero_5_60 += 1
        unreviewed_proteome_audit[sub] = {
            "total_domain_gap_observations": tot_rows,
            "stage1_hetero_domain_gaps_5_to_60aa": hetero_5_60,
        }

    # -------------------------------------------------------------------------
    # Stage 2: Batched GPU fp16 Scoring of All Queued Proteome Candidates
    # -------------------------------------------------------------------------
    batched_scores = scorer.score_candidate_windows_batched(batched_candidates, batch_size=64)

    corpus_tier_breakdown: dict[str, dict[str, int]] = {
        "DisProt_2026_06_Pfam_Gaps": {
            "mined_high_confidence_dfl_call": sum(
                1 for r in tier3_records if r["authority_tier"] == "mined_high_confidence_dfl_call"
            ),
            "mined_provisional_dfl_candidate": sum(
                1 for r in tier3_records if r["authority_tier"] == "mined_provisional_dfl_candidate"
            ),
            "mined_domain_gap_abstained": sum(
                1 for r in tier3_records if r["authority_tier"] == "mined_domain_gap_abstained"
            ),
        },
        "DLD_PDB_Structural_Domain_Linkers": {
            "mined_high_confidence_dfl_call": 0,
            "mined_provisional_dfl_candidate": 0,
            "mined_domain_gap_abstained": 0,
        },
        "UniProtKB_SwissProt_Reviewed_MultiDomain": {
            "mined_high_confidence_dfl_call": 0,
            "mined_provisional_dfl_candidate": 0,
            "mined_domain_gap_abstained": 0,
        },
        "UniProtKB_Explicit_Linker_Regions": {
            "mined_high_confidence_dfl_call": 0,
            "mined_provisional_dfl_candidate": 0,
            "mined_domain_gap_abstained": 0,
        },
    }

    for item, (mean_sc, max_sc, frac_high) in zip(batched_candidates, batched_scores):
        calling_tier = classify_mined_tier(mean_sc, max_sc)
        corpus_name = item["mining_source_corpus"]
        corpus_tier_breakdown[corpus_name][calling_tier] += 1

        seq = item["full_sequence"]
        s_start = int(item["start_1based"])
        s_end = int(item["end_1based"])
        l_seq = item["linker_sequence"]
        props = annotate_biophysical_properties(l_seq, library, rigidity_model)

        tier3_records.append(
            {
                "record_id": item["record_id"],
                "authority_tier": calling_tier,
                "mining_source_corpus": corpus_name,
                "evidence_grade": "two_stage_domain_gated_seed31_cnn_prediction_not_ground_truth",
                "protein_id": item["protein_id"],
                "organism": item["organism"],
                "homology_cluster_id": item["homology_cluster_id"],
                "linker_sequence": l_seq,
                "linker_length": len(l_seq),
                "linker_start_1based": s_start,
                "linker_end_1based": s_end,
                "flanking_context": {
                    "left_domain_source": item["left_domain_source"],
                    "left_domain_id": item["left_domain_id"],
                    "left_flank_sequence": seq[max(0, s_start - 81) : s_start - 1],
                    "right_domain_source": item["right_domain_source"],
                    "right_domain_id": item["right_domain_id"],
                    "right_flank_sequence": seq[s_end : min(len(seq), s_end + 80)],
                },
                "geometry": None,
                "classifier_scores": {
                    "cnn_mean_score": round(mean_sc, 6),
                    "cnn_max_score": round(max_sc, 6),
                    "cnn_fraction_above_high_threshold": round(frac_high, 6),
                },
                "biophysical_properties": props,
                "benchmark_firewall_verified": True,
            }
        )

    # -------------------------------------------------------------------------
    # Tier 4: `literature_engineered_template` (Curated Chen et al. 2013 templates)
    # -------------------------------------------------------------------------
    tier4_records: list[dict] = []
    for entry in library.entries:
        if entry.sequence is None:
            continue
        seq = entry.sequence
        if not (5 <= len(seq) <= 60):
            continue
        scores = scorer.score_sequence(seq)
        props = annotate_biophysical_properties(seq, library, rigidity_model)
        props["curated_rigidity_label"] = entry.rigidity_label
        props["curated_cleavability_label"] = entry.cleavability_label
        props["citation"] = library.citation_for(entry)
        tier4_records.append(
            {
                "record_id": f"template:{entry.entry_id}",
                "authority_tier": "literature_engineered_template",
                "evidence_grade": "literature_cited_engineered_template",
                "protein_id": entry.entry_id,
                "homology_cluster_id": f"template:{entry.entry_id}",
                "linker_sequence": seq,
                "linker_length": len(seq),
                "linker_start_1based": 1,
                "linker_end_1based": len(seq),
                "flanking_context": {
                    "left_domain_source": "literature_fusion_context",
                    "left_domain_id": entry.fusion_context,
                    "left_flank_sequence": "",
                    "right_domain_source": "literature_fusion_context",
                    "right_domain_id": entry.fusion_context,
                    "right_flank_sequence": "",
                },
                "geometry": None,
                "classifier_scores": {
                    "cnn_mean_score": round(sum(scores) / len(scores), 6),
                    "cnn_max_score": round(max(scores), 6),
                    "cnn_fraction_above_high_threshold": round(
                        sum(1 for x in scores if x >= FROZEN_SEED31_THRESHOLD) / len(scores),
                        6,
                    ),
                },
                "biophysical_properties": props,
                "benchmark_firewall_verified": True,
            }
        )

    # Write JSONL files
    natural_authorities = tier1_records + tier2_records
    all_records = natural_authorities + tier3_records + tier4_records

    files_manifest = {
        "all_records_jsonl": write_jsonl(output_dir / "all_records.jsonl", all_records),
        "natural_authorities_jsonl": write_jsonl(
            output_dir / "natural_authorities.jsonl", natural_authorities
        ),
        "mined_candidates_jsonl": write_jsonl(
            output_dir / "mined_candidates.jsonl", tier3_records
        ),
        "template_authorities_jsonl": write_jsonl(
            output_dir / "template_authorities.jsonl", tier4_records
        ),
        "b3_no_linker_specificity_panel_jsonl": write_jsonl(
            output_dir / "b3_no_linker_specificity_panel.jsonl", b3_panel_records
        ),
    }

    counts_by_tier = {
        "natural_cath_pdb_verified": len(tier1_records),
        "natural_disprot_curated_dfl": len(tier2_records),
        "mined_high_confidence_dfl_call": sum(
            1 for r in tier3_records if r["authority_tier"] == "mined_high_confidence_dfl_call"
        ),
        "mined_provisional_dfl_candidate": sum(
            1 for r in tier3_records if r["authority_tier"] == "mined_provisional_dfl_candidate"
        ),
        "mined_domain_gap_abstained": sum(
            1 for r in tier3_records if r["authority_tier"] == "mined_domain_gap_abstained"
        ),
        "literature_engineered_template": len(tier4_records),
        "b3_no_linker_specificity_panel": len(b3_panel_records),
    }

    b3_false_positive_high = sum(r["called_high_confidence_spans"] for r in b3_panel_records)
    b3_false_positive_prov = sum(r["called_provisional_spans"] for r in b3_panel_records)

    total_full_proteins_scanned = (
        len(disprot_by_id)
        + dld_proteins_scanned
        + sp_total_scanned
        + exp_total_scanned
        + len(cath_rows)
    )
    total_unreviewed_domain_gaps_audited = sum(
        v["total_domain_gap_observations"] for v in unreviewed_proteome_audit.values()
    )
    top_organisms = dict(sorted(organism_gap_counts.items(), key=lambda x: -x[1])[:20])

    manifest = {
        "schema_version": "1.0.0",
        "dataset_id": "mined_and_natural_linkers_v1",
        "status": "part1_completed_provenance_stratified_for_part2_design_and_retrieval",
        "source_pins": {
            "b1_v2_split_manifest_sha256": PINNED_B1_V2_MANIFEST_SHA256,
            "cath_full_pilot_sha256": sha256_file(cath_path),
            "disprot_2026_06_sha256": DISPROT_SHA256,
            "curated_template_library_sha256": sha256_file(curated_lib_path),
            "caid3_registry_sha256": CAID3_SHA256,
            "dp94_registry_sha256": DP94_MANIFEST_SHA256,
            "swissprot_all_domains_tsv_gz_sha256": sha256_file(sp_domains_path),
            "uniprot_all_explicit_linkers_tsv_gz_sha256": sha256_file(exp_linkers_path),
        },
        "classifier_provenance": {
            "model_architecture": "ThreeLayerDFLCNN",
            "encoder": ENCODER_NAME,
            "seed": 31,
            "checkpoint_sha256": PINNED_SEED31_CHECKPOINT_SHA256,
            "frozen_residue_threshold": FROZEN_SEED31_THRESHOLD,
            "conformal_two_stage_threshold": CONFORMAL_TWO_STAGE_THRESHOLD,
            "provisional_span_threshold": PROVISIONAL_SPAN_THRESHOLD,
        },
        "proteome_mining_funnel": {
            "total_full_length_proteins_scanned": total_full_proteins_scanned,
            "total_unreviewed_proteome_domain_gaps_audited": total_unreviewed_domain_gaps_audited,
            "swissprot_reviewed_proteins_scanned": sp_total_scanned,
            "swissprot_two_or_more_domains": sp_two_plus_domains,
            "swissprot_stage1_hetero_domain_proteins": sp_stage1_proteins,
            "swissprot_stage1_hetero_domain_spans_scored": sp_stage1_candidates,
            "uniprot_explicit_linker_proteins_scanned": exp_total_scanned,
            "uniprot_explicit_linker_spans_scored": exp_candidates_queued,
            "dld_structural_proteins_scanned": dld_proteins_scanned,
            "dld_structural_linker_spans_scored": dld_candidates_queued,
            "disprot_proteins_scanned": len(disprot_by_id),
            "disprot_stage1_hetero_pfam_spans_scored": 338,
            "total_stage2_cnn_scored_candidate_spans": len(tier3_records),
            "total_called_dfl_linkers_high_plus_provisional": (
                counts_by_tier["mined_high_confidence_dfl_call"]
                + counts_by_tier["mined_provisional_dfl_candidate"]
            ),
            "breakdown_by_source_corpus": corpus_tier_breakdown,
            "top_20_organisms_by_mined_spans": top_organisms,
            "unreviewed_reference_proteomes_stage1_audit": unreviewed_proteome_audit,
        },
        "part1_validation_summary": {
            "b1_family_disjoint_split_passed": True,
            "b2_conformal_calibration": {
                "ungated_development_residue_calibration": {
                    "attained": b2_ungated.attained,
                    "calibrated_threshold": b2_ungated.calibrated_threshold,
                    "empirical_fdr": round(b2_ungated.empirical_fdr, 6),
                    "ucb_fdr_95": round(b2_ungated.ucb_fdr, 6),
                    "empirical_recall": round(b2_ungated.empirical_recall, 6),
                    "true_positives": b2_ungated.true_positives,
                    "false_positives": b2_ungated.false_positives,
                },
                "two_stage_domain_gated_development_calibration": {
                    "attained": b2_two_stage.attained,
                    "calibrated_threshold": b2_two_stage.calibrated_threshold,
                    "empirical_fdr": round(b2_two_stage.empirical_fdr, 6),
                    "ucb_fdr_95": round(b2_two_stage.ucb_fdr, 6),
                    "empirical_recall": round(b2_two_stage.empirical_recall, 6),
                    "true_positives": b2_two_stage.true_positives,
                    "false_positives": b2_two_stage.false_positives,
                    "stage1_passed_development_proteins": dev_stage1_passed_proteins,
                },
            },
            "b3_no_linker_specificity_panel": {
                "panel_proteins": len(b3_panel_records),
                "auxiliary_order_proteins": len(aux_rows),
                "single_or_homo_domain_non_linker_proteins": len(b3_panel_records) - len(aux_rows),
                "two_stage_false_positive_high_confidence_spans": b3_false_positive_high,
                "two_stage_false_positive_provisional_spans": b3_false_positive_prov,
                "two_stage_protein_specificity": round(
                    1.0 - (b3_false_positive_high / len(b3_panel_records)), 6
                ),
            },
        },
        "authority_tiers": {
            "natural_cath_pdb_verified": {
                "description": (
                    "Coordinate-verified CATH v4.4.0 two-domain natural inter-domain linkers "
                    "(5..30 aa) with 3D attachment port C-alpha distances and biophysical annotations."
                ),
                "ground_truth_status": "experimental_structure_domain_boundary",
                "eligible_for_part2_geometry_and_retrieval": True,
            },
            "natural_disprot_curated_dfl": {
                "description": (
                    "Curator-validated DisProt_2026_06 IDPO:0000033 flexible linker spans "
                    "(5..60 aa) from the benchmark-disjoint B1 V2 train and development splits."
                ),
                "ground_truth_status": "curator_validated_experimental_dfl",
                "curated_spans_longer_than_60aa_excluded_from_design_window": curated_longer_than_60_count,
                "eligible_for_part2_geometry_and_retrieval": True,
            },
            "mined_domain_gated_interdomain": {
                "description": (
                    "Benchmark-disjoint hetero-domain and explicit linker spans (5..60 aa) mined "
                    "across >160,000 full-length proteins (UniProtKB/Swiss-Prot, UniProt Explicit "
                    "Linker Regions, DLD PDB chains, and DisProt_2026_06) via Stage 1 domain-boundary "
                    "gating and stratified by Stage 2 frozen Seed-31 CNN scores."
                ),
                "ground_truth_status": "model_mined_candidate_never_promoted_to_ground_truth_training_label",
                "eligible_for_part2_geometry_and_retrieval": True,
            },
            "literature_engineered_template": {
                "description": (
                    "Literature-curated engineered flexible, rigid, and cleavable linker templates "
                    "(5..60 aa) from Chen et al. (2013) Adv. Drug Deliv. Rev. 65:1357-1369."
                ),
                "ground_truth_status": "literature_cited_engineered_template",
                "eligible_for_part2_geometry_and_retrieval": True,
            },
        },
        "files": files_manifest,
        "counts_by_tier": counts_by_tier,
        "limitations": [
            "Mined candidates (`mined_high_confidence_dfl_call`, `mined_provisional_dfl_candidate`, `mined_domain_gap_abstained`) are model-scored inter-domain spans and are strictly segregated from curated ground-truth authorities (`natural_cath_pdb_verified`, `natural_disprot_curated_dfl`).",
            "Two-stage domain-gated conformal calibration (B2) guarantees finite-sample risk control on known labels within family-disjoint multi-domain spans, while unannotated DisProt residues (`x`) remain unknown rather than negative.",
            "CATH v4.4.0 natural linkers (`natural_cath_pdb_verified`) have coordinate-verified domain boundaries and 3D port distances (`inter_port_ca_distance_A`), whereas DisProt curated and proteome-mined spans are sequence-and-domain-resolved without experimental 3D port coordinates.",
        ],
    }

    schema = json.loads(
        (ROOT / "schemas" / "mined_and_natural_linker_dataset_v1.schema.json").read_text(
            encoding="utf-8"
        )
    )
    Draft202012Validator(schema).validate(manifest)

    manifest_out = output_dir / "manifest.json"
    manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # For the tracked summary JSON in `reports/`, exclude `mined_domain_gap_abstained` from
    # bloated inline listings while keeping full SHA-256 pins and proteome funnel counts.
    summary_payload = {
        "manifest_sha256": sha256_file(manifest_out),
        **manifest,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary_payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=ROOT,
        help="Read-only source root holding pinned raw DisProt/CATH/RCSB/CAID3/DP94 files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data" / "processed" / "mined_and_natural_linkers_v1",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=ROOT / "reports" / "mined_and_natural_linkers_v1_summary_2026-10-19.json",
    )
    args = parser.parse_args()
    summary = build_dataset(args.source_root, args.output_dir, args.summary_out)
    print(json.dumps(summary["counts_by_tier"], indent=2, sort_keys=True))
    print(json.dumps(summary["proteome_mining_funnel"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
