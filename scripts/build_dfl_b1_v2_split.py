"""Build and verify the preregistered V2 family-disjoint DFL split (`B1` resolution).

Implements `literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md` using the
pinned local `DisProt_2026_06` release (`1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140`),
`CAID3` reference FASTA (`4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`),
and `DP94` sequence-only registry (`fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`).
Benchmark labels are never read or parsed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from scripts.audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from scripts.audit_dfl_firewall_null import shuffled
    from scripts.build_dfl_curated_evidence_pool import CAID3_SHA256, DP94_MANIFEST_SHA256, family_components, mask_entry
    from scripts.preflight_dfl_family_training_v2 import RULE_V2, verify, v2_homolog_hit
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from audit_dfl_firewall_null import shuffled
    from build_dfl_curated_evidence_pool import CAID3_SHA256, DP94_MANIFEST_SHA256, family_components, mask_entry
    from preflight_dfl_family_training_v2 import RULE_V2, verify, v2_homolog_hit


SOURCE_POLICIES = {
    "DisProt_2026_06": {
        "task": "DFL_curated_positive_and_experimental_order",
        "positive_term": "IDPO:0000033",
        "negative_term": "IDPO:0000006",
        "zero_semantics": "explicit_negative_curated",
        "unknown_semantics": "all_unannotated_and_other_function_residues_are_x",
        "rights": "CC BY 4.0",
    },
    "DisProt_2026_06_order_panel": {
        "task": "DFL_experimental_order_negative_panel",
        "positive_term": None,
        "negative_term": "IDPO:0000006",
        "zero_semantics": "explicit_negative_curated",
        "unknown_semantics": "all_unannotated_and_other_function_residues_are_x",
        "rights": "CC BY 4.0",
    },
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def run_shuffle_null(candidates: list[tuple[str, str]], dp94: list[tuple[str, str]],
                     seeds: tuple[int, ...] = (173, 271, 811), sample_size: int = 40) -> dict:
    obs_hits = 0
    res_hits = 0
    blk_hits = 0
    total_pairs = 0
    for seed in seeds:
        rng = random.Random(seed)
        c_sub = rng.sample(candidates, sample_size)
        d_sub = rng.sample(dp94, sample_size)
        total_pairs += sample_size * sample_size
        for _, c_seq in c_sub:
            res_seq = shuffled(c_seq, rng, 1)
            blk_seq = shuffled(c_seq, rng, 12)
            for _, b_seq in d_sub:
                obs_hits += bool(v2_homolog_hit(c_seq, b_seq))
                res_hits += bool(v2_homolog_hit(res_seq, b_seq))
                blk_hits += bool(v2_homolog_hit(blk_seq, b_seq))
    return {
        "seeds": list(seeds),
        "sample_size_each_side": sample_size,
        "total_pairs": total_pairs,
        "observed_hits": obs_hits,
        "residue_shuffle_hits": res_hits,
        "block12_shuffle_hits": blk_hits,
        "residue_shuffle_rate": res_hits / total_pairs,
        "block12_shuffle_rate": blk_hits / total_pairs,
    }


def allocate_components(rows_by_comp: dict[str, list[dict]], seed: int,
                        fractions: tuple[float, ...], role_names: tuple[str, ...]) -> dict[str, list[dict]]:
    total = sum(len(rows) for rows in rows_by_comp.values())
    ordered_comps = sorted(
        rows_by_comp.keys(),
        key=lambda comp: (-len(rows_by_comp[comp]), digest_text(f"{seed}:{comp}")),
    )
    assigned: dict[str, list[dict]] = {name: [] for name in role_names}
    for comp in ordered_comps:
        target_role = max(
            role_names,
            key=lambda name: (
                fractions[role_names.index(name)] * total - len(assigned[name]),
                -role_names.index(name),
            ),
        )
        assigned[target_role].extend(rows_by_comp[comp])
    return assigned


def build(source_root: Path, output_dir: Path, *, seed: int = 17) -> dict:
    source_root, output_dir = source_root.resolve(), output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite V2 split directory: {output_dir}")

    raw_path = source_root / "data/raw/disprot/disprot_2026_06_regions.json"
    if sha(raw_path) != DISPROT_SHA256:
        raise ValueError("Pinned DisProt 2026_06 raw JSON hash mismatch")
    caid_path = source_root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != CAID3_SHA256:
        raise ValueError("Pinned CAID3 FASTA hash mismatch")
    dp_manifest_path = source_root / "data/processed/prospective_dfl_holdout_dp94_v1/manifest.json"
    if sha(dp_manifest_path) != DP94_MANIFEST_SHA256:
        raise ValueError("Pinned DP94 manifest hash mismatch")
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    dp_fasta_path = dp_manifest_path.parent / "sequences.fasta"
    dp_canonical = dp_fasta_path.read_text(encoding="utf-8").replace("\r\n", "\n")
    if digest_text(dp_canonical) != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("Pinned DP94 FASTA hash mismatch")

    caid3 = caid3_sequences(caid_path)
    dp94 = two_line_fasta(dp_fasta_path)
    references = [("caid3", k, s) for k, s in caid3] + [("dp94", k, s) for k, s in dp94]
    ref_ids = {k for _, k, _ in references}

    cand_path = source_root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl"
    cand_rows = [(r["protein_id"], r["sequence"]) for r in map(json.loads, cand_path.read_text("utf-8").splitlines())]
    null_verification = run_shuffle_null(cand_rows, dp94)

    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    metadata: dict[str, dict] = {}
    by_id: dict[str, dict] = {}
    for entry in raw["data"]:
        key = entry["disprot_id"]
        by_id[key] = entry
        try:
            metadata[key] = project_entry(entry)
        except ValueError:
            continue

    components = family_components(metadata)
    for reg_name, ref_id, ref_seq in references:
        if ref_id not in components:
            components[ref_id] = f"benchmark_isolated_{reg_name}_{ref_id}"
    benchmark_components = {components[ref_id] for _, ref_id, _ in references}

    pipeline_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, cwd=Path(__file__).resolve().parents[1]).strip()
    config_hash = digest_text(json.dumps({"rule_v2": RULE_V2, "policies": SOURCE_POLICIES, "seed": seed}, sort_keys=True))

    pos_by_comp: dict[str, list[dict]] = defaultdict(list)
    order_by_comp: dict[str, list[dict]] = defaultdict(list)
    quarantine_counts = Counter()

    for key, entry in sorted(by_id.items()):
        # Include curator-validated, non-obsolete, released, PMID/DOI-cited regions
        row = mask_entry(entry, include_unpublished_flag=True)
        if row is None:
            continue
        if key not in metadata or metadata[key]["sequence"] != row["sequence"]:
            quarantine_counts["unresolved_family_metadata"] += 1
            continue
        comp_id = components[key]
        if key in ref_ids or comp_id in benchmark_components:
            quarantine_counts["benchmark_id_or_family_component"] += 1
            continue
        if any(v2_homolog_hit(row["sequence"], ref_seq) for _, _, ref_seq in references):
            quarantine_counts["benchmark_v2_sequence_homolog"] += 1
            continue
        row["family_component_id"] = comp_id
        row["provenance"] = {
            "source_id": key,
            "input_sha256": DISPROT_SHA256,
            "pipeline_commit": pipeline_commit,
            "config_hash": config_hash,
            "python_version": platform.python_version(),
        }
        if "1" in row["labels"]:
            row["source"] = "DisProt_2026_06"
            row["task"] = "DFL_curated_positive_and_experimental_order"
            pos_by_comp[comp_id].append(row)
        else:
            row["source"] = "DisProt_2026_06_order_panel"
            row["task"] = "DFL_experimental_order_negative_panel"
            order_by_comp[comp_id].append(row)

    # Any order-only record in the same family component as a DFL-positive record
    # stays in the same component so whole components never cross splits.
    disjoint_order_by_comp: dict[str, list[dict]] = {}
    for comp_id, rows in order_by_comp.items():
        if comp_id in pos_by_comp:
            for row in rows:
                row["source"] = "DisProt_2026_06"
                row["task"] = "DFL_curated_positive_and_experimental_order"
            pos_by_comp[comp_id].extend(rows)
        else:
            disjoint_order_by_comp[comp_id] = rows

    # Also enforce V2 pairwise sequence disjunction across components
    all_kept_comps = list(pos_by_comp.keys()) + list(disjoint_order_by_comp.keys())
    comp_parent = {c: c for c in all_kept_comps}

    def find_c(c: str) -> str:
        while comp_parent[c] != c:
            comp_parent[c] = comp_parent[comp_parent[c]]
            c = comp_parent[c]
        return c

    all_rows_by_comp = {**pos_by_comp, **disjoint_order_by_comp}
    for i, c1 in enumerate(all_kept_comps):
        for c2 in all_kept_comps[i + 1:]:
            if find_c(c1) == find_c(c2):
                continue
            if any(v2_homolog_hit(r1["sequence"], r2["sequence"])
                   for r1 in all_rows_by_comp[c1] for r2 in all_rows_by_comp[c2]):
                comp_parent[find_c(c1)] = find_c(c2)

    merged_pos: dict[str, list[dict]] = defaultdict(list)
    merged_order: dict[str, list[dict]] = defaultdict(list)
    grouped_rows: dict[str, list[dict]] = defaultdict(list)
    for comp_id, rows in all_rows_by_comp.items():
        root_comp = find_c(comp_id)
        grouped_rows[root_comp].extend(rows)

    for root_comp, rows in grouped_rows.items():
        for row in rows:
            row["family_component_id"] = root_comp
            components[row["protein_id"]] = root_comp
        if any("1" in r["labels"] for r in rows):
            for r in rows:
                r["source"] = "DisProt_2026_06"
                r["task"] = "DFL_curated_positive_and_experimental_order"
            merged_pos[root_comp] = sorted(rows, key=lambda r: r["protein_id"])
        else:
            merged_order[root_comp] = sorted(rows, key=lambda r: r["protein_id"])

    pos_alloc = allocate_components(merged_pos, seed=seed, fractions=(0.70, 0.30), role_names=("train", "development"))
    order_alloc = allocate_components(
        merged_order,
        seed=seed + 100,
        fractions=(0.50, 0.25, 0.25),
        role_names=("train", "development", "auxiliary"),
    )

    for role in ("train", "development"):
        for r in order_alloc[role]:
            r["source"] = "DisProt_2026_06"
            r["task"] = "DFL_curated_positive_and_experimental_order"

    split_rows = {
        "train": sorted(pos_alloc["train"] + order_alloc["train"], key=lambda r: r["protein_id"]),
        "development": sorted(pos_alloc["development"] + order_alloc["development"], key=lambda r: r["protein_id"]),
        "auxiliary": sorted(order_alloc["auxiliary"], key=lambda r: r["protein_id"]),
    }

    output_dir.mkdir(parents=True)
    reg_dir = output_dir / "registries"
    reg_dir.mkdir()

    def write_registry(name: str, status: str, seal: str, release: str, records: list[tuple[str, str]]) -> dict:
        fasta_rel = f"{name}.fasta"
        fasta_text = "".join(f">{k}\n{s}\n" for k, s in records)
        (reg_dir / fasta_rel).write_text(fasta_text, encoding="utf-8", newline="\n")
        reg_doc = {
            "schema_version": "1.0.0",
            "name": name,
            "status": status,
            "label_seal_state": seal,
            "source_release": release,
            "fasta_path": fasta_rel,
            "fasta_sha256_canonical_lf": digest_text(fasta_text),
            "ids": [k for k, _ in records],
        }
        reg_path = reg_dir / f"{name}.registry.json"
        reg_path.write_text(json.dumps(reg_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        return {"name": name, "path": f"registries/{reg_path.name}", "sha256": sha(reg_path)}

    reg_entries = [
        write_registry("caid3", "sequence_only_protected", "retrospective_spent_exclusion_only", "CAID3_linker_reference", caid3),
        write_registry("dp94", "sequence_only_labels_sealed", "unopened", "DisoFLAG_DP94_sequence_only", dp94),
    ]

    membership = {}
    for _, ref_id, _ in references:
        membership[ref_id] = components[ref_id]
    for role_rows in split_rows.values():
        for r in role_rows:
            membership[r["protein_id"]] = r["family_component_id"]

    family_doc = {
        "schema_version": "2.0.0",
        "status": "family_clustered_v2",
        "method": "transitive_Pfam_Gene3D_UniRef50_plus_RULE_V2_reciprocal_80pct_coverage",
        "version": "DisProt_2026_06+parasail_sw_trace_blosum62",
        "source_sha256": DISPROT_SHA256,
        "membership": dict(sorted(membership.items())),
    }
    family_path = output_dir / "family_manifest.json"
    family_path.write_text(json.dumps(family_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

    files_info = {}
    counts = {}
    mask_counts = {}
    for role, rows in split_rows.items():
        rel_name = f"{role}.jsonl"
        file_path = output_dir / rel_name
        file_path.write_text(
            "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows),
            encoding="utf-8",
            newline="\n",
        )
        files_info[role] = {"path": rel_name, "records": len(rows), "sha256": sha(file_path)}
        counts[role] = len(rows)
        mask_counts[role] = {label: sum(r["labels"].count(label) for r in rows) for label in "10x"}

    manifest = {
        "schema_version": "2.0.0",
        "status": "preregistered_v2_family_disjoint_split",
        "preregistration_doc": "literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md",
        "similarity_rule": RULE_V2,
        "shuffle_null_verification": null_verification,
        "source_policies": SOURCE_POLICIES,
        "benchmark_registries": reg_entries,
        "family_manifest": {"path": family_path.name, "sha256": sha(family_path)},
        "files": files_info,
        "counts": counts,
        "mask_counts": mask_counts,
        "quarantine_counts": dict(quarantine_counts),
        "limitations": [
            "Curator-validated released DisProt_2026_06 IDPO:0000033 (1) and IDPO:0000006 (0) regions only; all other residues remain x.",
            "Includes publicly released (2022_03-2026_06), PMID/DOI-cited, curator-validated regions that retained DisProt's internal unpublished=true flag.",
            "Auxiliary split is a family-disjoint IDPO:0000006 experimental-order negative panel; single-domain PDB crystallographic panel remains recommended for full proteome mining calibration.",
            "No benchmark labels were opened or parsed.",
        ],
    }
    alias_doc = {
        "schema_version": "2.0.0",
        "purpose": "Auxiliary experimental-order negative panel alias for bounded DFL launcher preflight",
        "files": {
            "pretrain": {
                "path": files_info["auxiliary"]["path"],
                "records": files_info["auxiliary"]["records"],
                "sha256": files_info["auxiliary"]["sha256"],
            }
        },
    }
    alias_path = output_dir / "auxiliary_alias_manifest.json"
    alias_path.write_text(json.dumps(alias_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    verification = verify(manifest_path)
    return {
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha(manifest_path),
        "auxiliary_alias_manifest_sha256": sha(alias_path),
        "verification": verification,
        "quarantine_counts": dict(quarantine_counts),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    print(json.dumps(build(args.source_root, args.output_dir, seed=args.seed), indent=2, sort_keys=True))
