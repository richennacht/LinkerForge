"""Development-only audit of a conservative InterPro/Pfam domain gate.

Live UniProt and InterPro responses are checked against the frozen B1 V2
development sequences. This script does not train, open benchmark labels, or
infer DFL labels from domain gaps. Output is an isolated diagnostic snapshot.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from linkerforge.detector.domain_gate import DomainMatch, domain_gap_candidates

DEV = ROOT / "data/processed/dfl_b1_v2_family_split/development.jsonl"
API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{}/?page_size=200"
FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str) -> tuple[int, bytes]:
    request = Request(url, headers={"Accept": "application/json, text/plain", "User-Agent": "LinkerForge-domain-gate-audit/1.0"})
    with urlopen(request, timeout=20) as response:
        return response.status, response.read()


def parse_fasta(raw: bytes) -> str:
    lines = raw.decode("utf-8").splitlines()
    if not lines or not lines[0].startswith(">") or len(lines) < 2:
        raise ValueError("Invalid FASTA response")
    return "".join(lines[1:]).upper()


def parse_pfam(raw: bytes, accession: str, length: int) -> list[DomainMatch]:
    document = json.loads(raw)
    if document.get("next") is not None:
        raise ValueError("Paginated Pfam response requires review")
    matches = []
    for entry in document["results"]:
        metadata = entry["metadata"]
        if metadata["source_database"].lower() != "pfam":
            raise ValueError("Unexpected non-Pfam entry")
        for protein in entry["proteins"]:
            if protein["accession"].upper() != accession.upper() or protein["protein_length"] != length:
                raise ValueError("InterPro protein identity/length mismatch")
            for location in protein["entry_protein_locations"]:
                fragments = location["fragments"]
                if len(fragments) != 1 or fragments[0].get("dc-status") != "CONTINUOUS":
                    raise ValueError("Discontinuous or multi-fragment domain annotation")
                part = fragments[0]
                matches.append(DomainMatch(int(part["start"]), int(part["end"]),
                                           metadata["accession"], metadata["type"]))
    return matches


def audit(rows: list[dict]) -> dict:
    counts = Counter()
    details = []
    for row in rows:
        accession, sequence, labels = row["uniprot_accession"], row["sequence"], row["labels"]
        if len(sequence) != len(labels):
            raise ValueError("Frozen development sequence/mask length mismatch")
        item = {"protein_id": row["protein_id"], "uniprot_accession": accession}
        counts["proteins"] += 1
        counts["known_positive_residues"] += labels.count("1")
        counts["known_negative_residues"] += labels.count("0")
        counts["unknown_residues"] += labels.count("x")
        try:
            fasta_status, fasta_raw = fetch(FASTA.format(accession))
            if fasta_status != 200:
                raise ValueError("Unexpected UniProt FASTA status")
            item["uniprot_fasta_sha256"] = sha(fasta_raw)
            exact = parse_fasta(fasta_raw) == sequence
            if not exact:
                item["abstention"] = "uniprot_sequence_mismatch"
                counts[item["abstention"]] += 1
                details.append(item)
                continue
            pfam_status, pfam_raw = fetch(API.format(accession))
            item["interpro_http_status"] = pfam_status
            item["interpro_json_sha256"] = sha(pfam_raw)
            if pfam_status == 204 and not pfam_raw:
                matches = []  # InterPro's documented no-content/no-match response.
            elif pfam_status == 200:
                matches = parse_pfam(pfam_raw, accession, len(sequence))
            else:
                raise ValueError("Unexpected InterPro HTTP status")
            candidates, reason = domain_gap_candidates(sequence, matches, sequence_verified=True)
            if reason:
                item["abstention"] = reason
                counts[reason] += 1
            else:
                counts["proteins_with_candidate"] += 1
            item["domain_matches"] = len(matches)
            item["candidates"] = [{"start": span.start, "end": span.end,
                                   "left_accession": span.left_accession,
                                   "right_accession": span.right_accession} for span in candidates]
            for span in candidates:
                counts["candidate_spans"] += 1
                segment = labels[span.start - 1:span.end]
                counts["candidate_known_positive_residues"] += segment.count("1")
                counts["candidate_known_negative_residues"] += segment.count("0")
                counts["candidate_unknown_residues"] += segment.count("x")
        except (HTTPError, URLError, TimeoutError, ValueError, KeyError, json.JSONDecodeError) as error:
            item["abstention"] = "source_fetch_or_schema_failure"
            item["error_type"] = type(error).__name__
            counts["source_fetch_or_schema_failure"] += 1
        details.append(item)
    return {
        "schema_version": "1.0.0",
        "claim_status": "development_gate_coverage_only_not_precision_validation",
        "development_sha256": sha(DEV.read_bytes()),
        "counts": dict(counts),
        "details": details,
        "limitations": ["No CNN score or held-out precision computed", "Unknown residues are not negatives",
                        "Missing, ambiguous, or mismatched annotations abstain", "Live database responses are hashed but may change"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rows = [json.loads(line) for line in DEV.read_text(encoding="utf-8").splitlines()]
    result = audit(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result["counts"], indent=2, sort_keys=True))
