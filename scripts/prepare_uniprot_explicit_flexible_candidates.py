"""Reconstruct full proteins for explicitly annotated flexible-linker regions.

This creates an evidence-tiered candidate source, NOT DFL gold labels. It never
queries benchmark labels. Current UniProt sequences must exactly contain the
previously pinned region peptide at its stated coordinates, else quarantine.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import urllib.parse
import urllib.request

import jsonschema


SOURCE_URL = "https://rest.uniprot.org/uniprotkb/stream"
CATALOG_REL = Path("data/processed/linker_catalog/annotated_only_v1/annotated_linker_sequences.v1.jsonl.gz")
CATALOG_MANIFEST_REL = CATALOG_REL.with_name("manifest.v1.json")
SCHEMA_REL = Path("schemas/uniprot_explicit_flexible_linker_candidate.v1.schema.json")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def catalog_candidates(root: Path) -> tuple[list[dict], str]:
    source = root / CATALOG_REL
    manifest = json.loads((root / CATALOG_MANIFEST_REL).read_text(encoding="utf-8"))
    raw = source.read_bytes()
    if sha(raw) != manifest["output_sha256"][source.name]:
        raise ValueError("Pinned annotated-linker catalog hash changed")
    selected = []
    for line in gzip.decompress(raw).decode("utf-8").splitlines():
        row = json.loads(line)
        annotation = row.get("qualifying_sources", {}).get("uniprot_explicit_linker_region")
        if not annotation:
            continue
        note = annotation["note"].lower()
        if "flexible linker" not in note and "disordered linker" not in note:
            continue
        selected.append({"accession": annotation["accession"], "peptide": row["sequence"].upper(),
                         "start": annotation["start_1based"], "end": annotation["end_1based"],
                         "note": annotation["note"], "evidence": annotation["evidence"]})
    if len(selected) != 884 or len({row["accession"] for row in selected}) != 884:
        raise ValueError("Explicit flexible/disordered annotation roster changed")
    return selected, sha(raw)


def fetch_sequences(accessions: list[str]) -> tuple[dict[str, str], list[dict]]:
    results = {}
    responses = []
    for at in range(0, len(accessions), 40):
        batch = accessions[at:at + 40]
        query = " OR ".join(f"accession:{accession}" for accession in batch)
        url = SOURCE_URL + "?" + urllib.parse.urlencode({
            "format": "tsv", "fields": "accession,sequence", "query": f"({query})"})
        request = urllib.request.Request(url, headers={"User-Agent": "LinkerForge-evidence-tiered-source-audit"})
        with urllib.request.urlopen(request, timeout=40) as response:
            raw = response.read()
        parsed = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"))
        if not parsed or set(parsed[0]) != {"Entry", "Sequence"}:
            raise ValueError("Unexpected UniProt TSV schema")
        for row in parsed:
            key = row["Entry"]
            if key not in batch or key in results or not row["Sequence"]:
                raise ValueError(f"Unexpected or duplicate UniProt accession: {key}")
            results[key] = row["Sequence"].upper()
        responses.append({"batch_first": batch[0], "requested": len(batch), "returned": len(parsed),
                          "response_sha256": sha(raw)})
    return results, responses


def make_record(candidate: dict, sequence: str, provenance: dict) -> dict:
    start, end = candidate["start"], candidate["end"]
    if start < 1 or end < start or end > len(sequence):
        raise ValueError("Region coordinates outside current full sequence")
    if sequence[start - 1:end] != candidate["peptide"]:
        raise ValueError("Pinned linker peptide does not match current full sequence")
    if set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
        raise ValueError("Unexpected residue in full sequence")
    tier = ("pubmed_supported_region_not_DFL_adjudicated"
            if "ECO:0000269" in candidate["evidence"] and "PubMed:" in candidate["evidence"]
            else "explicit_region_weak")
    return {"schema_version": "1.0.0", "record_type": "uniprot_explicit_flexible_linker_candidate",
            "source_id": candidate["accession"], "sequence": sequence,
            "labels": "x" * (start - 1) + "1" * (end - start + 1) + "x" * (len(sequence) - end),
            "annotation": {"start_1based": start, "end_1based": end, "note": candidate["note"],
                           "evidence": candidate["evidence"], "evidence_tier": tier},
            "provenance": {**provenance, "uniprot_sequence_sha256": sha(sequence.encode())},
            "training_eligibility": "candidate_not_DFL_gold_or_negative_adjudicated"}


def build(root: Path, output: Path) -> dict:
    root, output = root.resolve(), output.resolve()
    code_root = Path(__file__).resolve().parents[1]
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    candidates, catalog_hash = catalog_candidates(root)
    sequences, responses = fetch_sequences([row["accession"] for row in candidates])
    schema_path = code_root / SCHEMA_REL
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema)
    timestamp = datetime.now(timezone.utc).isoformat()
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=code_root, text=True).strip()
    base_provenance = {"catalog_sha256": catalog_hash, "source_url": SOURCE_URL,
                       "retrieved_at_utc": timestamp, "pipeline_commit": commit,
                       "script_sha256": sha(Path(__file__).read_bytes()),
                       "schema_sha256": sha(schema_path.read_bytes()),
                       "software": {"python": sys.version.split()[0], "jsonschema": jsonschema.__version__}}
    valid, quarantine = [], []
    for candidate in candidates:
        accession = candidate["accession"]
        sequence = sequences.get(accession)
        if not sequence:
            quarantine.append({"accession": accession, "reason": "current_sequence_missing"})
            continue
        try:
            record = make_record(candidate, sequence, base_provenance)
            validator.validate(record)
        except (ValueError, jsonschema.ValidationError) as error:
            quarantine.append({"accession": accession, "reason": str(error)})
            continue
        valid.append(record)
    output.mkdir(parents=True)
    records_path = output / "candidates.jsonl"
    quarantine_path = output / "quarantine.jsonl"
    records_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in valid), encoding="utf-8")
    quarantine_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in quarantine), encoding="utf-8")
    manifest = {"schema_version": "1.0.0", "status": "explicit_flexible_linker_candidates_not_DFL_gold",
                "created_utc": timestamp, "candidate_catalog_sha256": catalog_hash,
                "source_manifest_sha256": sha((root / CATALOG_MANIFEST_REL).read_bytes()),
                "script_sha256": base_provenance["script_sha256"],
                "schema_sha256": base_provenance["schema_sha256"], "pipeline_commit": commit,
                "retrieval": responses, "candidate_roster": len(candidates), "validated": len(valid),
                "quarantined": len(quarantine),
                "pubmed_tier": sum(row["annotation"]["evidence_tier"].startswith("pubmed") for row in valid),
                "output_sha256": {records_path.name: sha(records_path.read_bytes()),
                                  quarantine_path.name: sha(quarantine_path.read_bytes())},
                "limitations": ["UniProt flexible-linker annotation does not certify the DisProt DFL target.",
                                "All non-annotated residues are x unknown, not negative.",
                                "No benchmark sequence/family screen or licence review has been applied.",
                                "Do not train or redistribute this candidate source without source-level review."]}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.source_root, args.output), indent=2, sort_keys=True))
