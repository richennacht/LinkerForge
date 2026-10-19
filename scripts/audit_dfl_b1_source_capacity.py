"""CPU-only B1 source census; benchmark labels are never read or scored.

This does not build a split. It estimates an upper bound on the available
DFL-positive proteins after the declared sequence firewall. Family filtering
can only lower that upper bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from collections import Counter
from pathlib import Path

if __package__:
    from scripts.preflight_dfl_family_training import RULE, local_hit
else:
    from preflight_dfl_family_training import RULE, local_hit


POOL_FILES = ("published_dfl_train", "disprot_positive_only")
IDPFUNNET_TR552_URL = ("https://raw.githubusercontent.com/IDRIDP/IDPFunNet/"
                      "a9830001da05b6c74c914c0283b0d4b4d3c69df7/dataset/TR552.fasta")
IDPFUNNET_TR552_SHA256 = "89e374e30924a277fb4bdcf1b7a82776017e4ab70397757b185fa9d0b16954b2"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def two_line_fasta(path: Path) -> list[tuple[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or len(lines) % 2:
        raise ValueError(f"Malformed two-line sequence registry: {path}")
    records = []
    for at in range(0, len(lines), 2):
        if not lines[at].startswith(">") or not lines[at + 1]:
            raise ValueError(f"Malformed sequence registry: {path}")
        records.append((lines[at][1:].split()[0], lines[at + 1].upper()))
    return records


def caid3_sequences(path: Path) -> list[tuple[str, str]]:
    """Read header/sequence pairs; consume but never retain the third line."""
    records = []
    with path.open(encoding="utf-8") as stream:
        while header := stream.readline():
            sequence = stream.readline().strip()
            unused_third_line = stream.readline()
            if not header.startswith(">") or not sequence or not unused_third_line:
                raise ValueError("Malformed official CAID3 three-line FASTA")
            records.append((header[1:].split()[0], sequence.upper()))
    if len(records) != 31:
        raise ValueError("Unexpected CAID3 reference count")
    return records


def audit(root: Path, *, include_idpfunnet: bool = False,
          weak_candidates: Path | None = None) -> dict:
    root = root.resolve()
    pool = root / "data/processed/linker_pred2_style_pool_v1"
    manifest_path = pool / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "prepared_not_linker_pred2_exact_split_not_training_cleared":
        raise ValueError("Unexpected DFL pool status")
    extra_path = root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl"
    extra_manifest_path = root / "data/processed/dfl_positive_candidates_2026_10_12.manifest.json"
    extra_manifest = json.loads(extra_manifest_path.read_text(encoding="utf-8"))
    if extra_manifest.get("status") != "positive_unknown_candidates_not_trainable":
        raise ValueError("Unexpected DisProt candidate status")
    dp_root = root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest.get("status") != "candidate_not_certified_independent":
        raise ValueError("DP94 candidate registry changed")
    dp_fasta = dp_root / "sequences.fasta"
    canonical = dp_fasta.read_text(encoding="utf-8").replace("\r\n", "\n")
    if hashlib.sha256(canonical.encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 sequence registry hash changed")
    dp94 = two_line_fasta(dp_fasta)
    if len(dp94) != dp_manifest["holdout_proteins"] or [key for key, _ in dp94] != dp_manifest["holdout_ids"]:
        raise ValueError("DP94 roster changed")
    caid_path = root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != manifest["source_sha256"]["caid3"]:
        raise ValueError("Official CAID3 source hash changed")
    caid3 = caid3_sequences(caid_path)
    references = (("CAID3", caid3), ("DP94", dp94))
    output = {"status": "sequence_screen_only_not_B1", "rule": RULE,
              "source_manifest_sha256": sha(manifest_path),
              "benchmark_sequences": {"CAID3": {"count": len(caid3), "sha256": sha(caid_path)},
                                      "DP94": {"count": len(dp94), "manifest_sha256": sha(dp_manifest_path),
                                               "fasta_sha256_canonical_lf": dp_manifest["sequence_fasta_sha256"]}},
              "sources": {}, "limitations": [
                  "Counts are an upper bound before family/domain clustering and negative-label adjudication.",
                  "DP94 is a candidate, not a certified DFL benchmark.",
                  "Only benchmark IDs and sequences were used; no benchmark labels or scores were inspected."]}
    all_survivors = []
    all_positive = []
    sources = [(name, pool / f"{name}.jsonl", manifest["files"][name]["sha256"],
                manifest["files"][name]["records"]) for name in POOL_FILES]
    sources.append(("disprot_candidate_catalog", extra_path, extra_manifest["output_sha256"],
                    extra_manifest["records"]))
    output["candidate_catalog_manifest_sha256"] = sha(extra_manifest_path)
    for name, path, expected_sha, expected_records in sources:
        if sha(path) != expected_sha:
            raise ValueError(f"Source hash changed: {name}")
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
        if len(rows) != expected_records:
            raise ValueError(f"Source row count changed: {name}")
        counts = Counter()
        survivors = []
        for row in rows:
            key, sequence, labels = row["protein_id"], row["sequence"].upper(), row["labels"]
            if not sequence or len(sequence) != len(labels) or set(labels) - set("01x"):
                raise ValueError(f"Invalid DFL mask: {name}/{key}")
            if name == "disprot_candidate_catalog":
                if row.get("record_type") != "disprot_dfl_positive_candidate" or "0" in labels:
                    raise ValueError(f"Invalid positive-only candidate: {name}/{key}")
            elif row["task"] not in {"DFL", "DFL_positive_unlabeled"}:
                raise ValueError(f"Non-DFL source substituted: {name}/{key}")
            if row.get("task") == "DFL_positive_unlabeled" and "0" in labels:
                raise ValueError(f"Positive-only unknown relabelled zero: {name}/{key}")
            counts["records"] += 1
            if "1" not in labels:
                counts["no_observed_DFL_positive"] += 1
                continue
            counts["observed_DFL_positive"] += 1
            all_positive.append((key, sequence))
            matched = []
            for registry_name, registry in references:
                if any(key == ref_key or local_hit(sequence, ref_seq) for ref_key, ref_seq in registry):
                    matched.append(registry_name)
            for registry_name in matched:
                counts[f"flagged_{registry_name}"] += 1
            if not matched:
                survivors.append((key, sequence))
        counts["positive_not_flagged_by_either_registry"] = len(survivors)
        counts["unique_surviving_sequences"] = len({sequence for _, sequence in survivors})
        output["sources"][name] = {"sha256": sha(path), **counts}
        all_survivors.extend(survivors)
    if include_idpfunnet:
        request = urllib.request.Request(IDPFUNNET_TR552_URL, headers={"User-Agent": "LinkerForge-B1-source-audit"})
        with urllib.request.urlopen(request, timeout=25) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != IDPFUNNET_TR552_SHA256:
            raise ValueError("Pinned IDPFunNet TR552 hash changed")
        lines = raw.decode("utf-8-sig").splitlines()
        if len(lines) != 552 * 9:
            raise ValueError("Unexpected IDPFunNet TR552 record layout")
        counts = Counter()
        survivors = []
        for at in range(0, len(lines), 9):
            record = lines[at:at + 9]
            key, sequence, labels = record[0][1:], record[1].upper(), record[8]
            if (not record[0].startswith(">") or not sequence or len(sequence) != len(labels)
                    or set(labels) - set("01xX-?")):
                raise ValueError(f"Malformed IDPFunNet DFL training record: {at // 9}")
            counts["records"] += 1
            if "1" not in labels:
                counts["no_observed_DFL_positive"] += 1
                continue
            counts["observed_DFL_positive"] += 1
            all_positive.append((key, sequence))
            matched = []
            for registry_name, registry in references:
                if any(key == ref_key or local_hit(sequence, ref_seq) for ref_key, ref_seq in registry):
                    matched.append(registry_name)
            for registry_name in matched:
                counts[f"flagged_{registry_name}"] += 1
            if not matched:
                survivors.append((key, sequence))
        counts["positive_not_flagged_by_either_registry"] = len(survivors)
        counts["unique_surviving_sequences"] = len({sequence for _, sequence in survivors})
        output["sources"]["idpfunnet_TR552_training_only"] = {
            "sha256": IDPFUNNET_TR552_SHA256, "source_url": IDPFUNNET_TR552_URL, **counts,
            "rights_note": "No redistribution licence verified; author training labels are not certified negatives."}
        all_survivors.extend(survivors)
    output["combined_positive_not_flagged"] = len(all_survivors)
    output["combined_unique_surviving_sequences"] = len({sequence for _, sequence in all_survivors})
    # A direct non-hit can still join a benchmark through a third sequence.
    # This is a conservative local-similarity graph, not a domain-family proof.
    nodes = {}
    ids = []
    protected = []
    for _, registry in references:
        for key, sequence in registry:
            nodes.setdefault(sequence, len(nodes))
            ids.append((key, sequence))
            protected.append(sequence)
    for key, sequence in all_positive:
        nodes.setdefault(sequence, len(nodes))
        ids.append((key, sequence))
    if len(nodes) > 1000:
        raise ValueError("Local graph cap exceeded; use a versioned external clustering tool")
    sequences = list(nodes)
    parent = list(range(len(sequences)))

    def find(at):
        while parent[at] != at:
            parent[at] = parent[parent[at]]
            at = parent[at]
        return at

    def union(left, right):
        parent[find(left)] = find(right)

    by_id = {}
    for key, sequence in ids:
        if key in by_id:
            union(nodes[sequence], by_id[key])
        else:
            by_id[key] = nodes[sequence]
    edges = 0
    for left in range(len(sequences)):
        for right in range(left + 1, len(sequences)):
            if local_hit(sequences[left], sequences[right]):
                union(left, right)
                edges += 1
    protected_components = {find(nodes[sequence]) for sequence in protected}
    eligible = {sequence for _, sequence in all_positive
                if find(nodes[sequence]) not in protected_components}
    component_sizes = Counter(find(nodes[sequence]) for sequence in eligible)
    output["declared_local_graph"] = {
        "unique_sequences_including_benchmarks": len(nodes),
        "local_similarity_edges": edges,
        "positive_sequences_outside_benchmark_components": len(eligible),
        "unprotected_positive_component_sizes": sorted(component_sizes.values(), reverse=True),
        "limitation": "Transitive local-alignment components are not validated biological/domain families."}
    if weak_candidates is not None:
        weak_candidates = weak_candidates.resolve()
        weak_manifest_path = weak_candidates / "manifest.json"
        weak_manifest = json.loads(weak_manifest_path.read_text(encoding="utf-8"))
        if weak_manifest.get("status") != "explicit_flexible_linker_candidates_not_DFL_gold":
            raise ValueError("Weak UniProt source status changed")
        weak_path = weak_candidates / "candidates.jsonl"
        if sha(weak_path) != weak_manifest["output_sha256"][weak_path.name]:
            raise ValueError("Weak UniProt source hash changed")
        weak_rows = [json.loads(line) for line in weak_path.read_text(encoding="utf-8").splitlines()]
        if len(weak_rows) != weak_manifest["validated"]:
            raise ValueError("Weak UniProt source count changed")
        weak_counts = Counter()
        unflagged = set()
        for row in weak_rows:
            key, sequence, labels = row["source_id"], row["sequence"], row["labels"]
            if (row.get("record_type") != "uniprot_explicit_flexible_linker_candidate"
                    or row.get("training_eligibility") != "candidate_not_DFL_gold_or_negative_adjudicated"
                    or len(sequence) != len(labels) or set(labels) - set("1x") or "1" not in labels):
                raise ValueError(f"Weak annotation mask/task invalid: {key}")
            tier = row["annotation"]["evidence_tier"]
            weak_counts[f"{tier}_total"] += 1
            matched = [name for name, registry in references
                       if any(key == ref_key or local_hit(sequence, ref_seq)
                              for ref_key, ref_seq in registry)]
            for name in matched:
                weak_counts[f"{tier}_flagged_{name}"] += 1
            if not matched:
                weak_counts[f"{tier}_direct_unflagged"] += 1
                unflagged.add(sequence)
        output["weak_uniprot_candidates"] = {
            "manifest_sha256": sha(weak_manifest_path), "records_sha256": sha(weak_path),
            "counts": dict(weak_counts), "unique_direct_unflagged_sequences": len(unflagged),
            "limitation": "Explicit flexible-linker notes are not certified DisProt DFL labels or adjudicated negatives."}
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--include-idpfunnet", action="store_true",
                        help="Fetch only pinned author TR552 training file; never any test labels")
    parser.add_argument("--weak-candidates", type=Path,
                        help="Optional versioned full-protein UniProt flexible-linker candidate directory")
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root, include_idpfunnet=args.include_idpfunnet,
                           weak_candidates=args.weak_candidates), indent=2, sort_keys=True))
