"""Reproducibly screen 100 Europe PMC abstracts for DFL dataset leads.

This is discovery, never training-data ingestion. The output is a ranked,
machine-assisted literature triage table; datasets require manual source,
licence, label-semantics and benchmark-lineage review before admission.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ENDPOINT = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
QUERIES = [
    "TITLE_ABS:(disordered AND flexible AND linker)",
    "TITLE:(APOD)",
    "TITLE:(TransDFL)",
    "TITLE:(DFLpred)",
    "TITLE:(LINKER-Pred)",
    "TITLE_ABS:(protein AND linker AND prediction)",
    "TITLE_ABS:(protein AND linker AND domain)",
    "TITLE_ABS:(interdomain AND linker AND protein)",
    "TITLE_ABS:(intrinsically AND disordered AND linker)",
    "TITLE_ABS:(protein AND linker AND sequence)",
    "TITLE_ABS:(intrinsically AND disordered AND function AND dataset)",
    "TITLE_ABS:(protein AND linker AND annotation)",
    "TITLE_ABS:(protein AND disorder AND benchmark)",
]
MAX_PER_QUERY = 80
TARGET = 100


def clean(value: str | None) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]*>", " ", value or ""))).strip()


def key_of(item: dict) -> str:
    title = re.sub(r"\W+", "", clean(item.get("title")).lower())
    if title:
        return "title:" + title
    if item.get("doi"):
        return "doi:" + item["doi"].lower().strip()
    if item.get("pmid"):
        return "pmid:" + item["pmid"]
    return "source:" + str(item.get("source")) + ":" + str(item.get("id"))


def classify(item: dict) -> tuple[int, str, str]:
    title = clean(item.get("title")).lower()
    abstract = clean(item.get("abstractText")).lower()
    body = title + " " + abstract
    score = 0
    for term, weight in (("disordered flexible linker", 40), ("dflpred", 35),
                         ("transdfl", 35), ("linker-pred", 30),
                         ("disoflag", 25), ("interdomain linker", 20),
                         ("domain linker", 15), ("protein linker", 15),
                         ("flexible linker", 12), ("disprot", 12),
                         ("dataset", 7), ("benchmark", 7),
                         ("intrinsic disorder", 4), ("sequence", 2)):
        if term in title:
            score += weight * 2
        elif term in body:
            score += weight
    if "linker" in body and ("protein" in body or "domain" in body):
        score += 15
    if "linker" in body and re.search(r"\bapod\b", title):
        score += 30
    if "dataset" in body or "data set" in body:
        score += 5
    exclusion = ("neuroanatomical sites of apod", "linker histone", "cross-linker",
                 "linkerability of protein ligands", "linker-containing signature peptide",
                 "linker domain in a dimeric compound")
    if any(phrase in title for phrase in exclusion):
        lane = "other_or_false_positive"
    elif "disordered flexible linker" in body or "dflpred" in body or "transdfl" in body:
        lane = "DFL_direct_candidate"
    elif "linker" in body and ("protein" in body or "domain" in body):
        lane = "protein_linker_related_not_verified_DFL"
    elif "disorder" in body or "disprot" in body:
        lane = "disorder_auxiliary_only"
    else:
        lane = "other_or_false_positive"
    reason = "abstract_mention_only_manual_full_text_data_and_rights_check_required"
    return score, lane, reason


def fetch(query: str) -> tuple[dict, str]:
    url = ENDPOINT + "?" + urllib.parse.urlencode(
        {"query": query, "pageSize": MAX_PER_QUERY, "format": "json", "resultType": "core"})
    req = urllib.request.Request(url, headers={"User-Agent": "LinkerForge-literature-audit/1.0"})
    with urllib.request.urlopen(req, timeout=45) as response:
        raw = response.read()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def survey(output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite survey: {output_dir}")
    candidates: dict[str, dict] = {}
    query_audit = []
    for query in QUERIES:
        result, payload_sha = fetch(query)
        rows = result.get("resultList", {}).get("result", [])
        query_audit.append({"query": query, "reported_hits": int(result["hitCount"]),
                            "retrieved": len(rows), "response_sha256": payload_sha})
        for item in rows:
            key = key_of(item)
            if key in candidates:
                candidates[key]["queries"].append(query)
            else:
                candidates[key] = {"source": item.get("source"), "id": item.get("id"),
                                   "pmid": item.get("pmid"), "pmcid": item.get("pmcid"),
                                   "doi": item.get("doi"), "year": item.get("pubYear"),
                                   "title": clean(item.get("title")),
                                   "abstract": clean(item.get("abstractText")),
                                   "queries": [query]}
    ranked = []
    for key, item in candidates.items():
        score, lane, reason = classify(item)
        ranked.append({"key": key, "score": score, "lane": lane,
                       "reason": reason, **item})
    ranked.sort(key=lambda row: (-row["score"], -len(row["queries"]), row["key"]))
    if len(ranked) < TARGET:
        raise ValueError(f"Only {len(ranked)} distinct papers; no fabricated 100-paper survey")
    related = [row for row in ranked if row["lane"] != "other_or_false_positive"]
    if len(related) < TARGET:
        raise ValueError(f"Only {len(related)} potentially relevant abstracts; no false positives added to fill 100")
    chosen = related[:TARGET]
    output_dir.mkdir(parents=True)
    jsonl = output_dir / "screened_100.jsonl"
    jsonl.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in chosen),
                     encoding="utf-8", newline="\n")
    tsv = output_dir / "screened_100.tsv"
    with tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
            "rank", "score", "lane", "year", "title", "doi", "pmid", "pmcid", "key", "reason"])
        writer.writeheader()
        for rank, row in enumerate(chosen, start=1):
            writer.writerow({col: rank if col == "rank" else row.get(col, "")
                             for col in writer.fieldnames})
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"schema_version": "1.0.0", "status": "abstract_triage_not_dataset_admission",
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "endpoint": ENDPOINT, "query_audit": query_audit,
                "distinct_candidates": len(ranked), "screened_count": len(chosen),
                "class_counts": {lane: sum(row["lane"] == lane for row in chosen)
                                 for lane in sorted({row["lane"] for row in chosen})},
                "files": {name: sha(output_dir / name) for name in ("screened_100.jsonl", "screened_100.tsv")},
                "limitations": ["Title/abstract ranking is machine-assisted, not full-text review of all 100.",
                                "No paper datasets or annotations are admitted by this survey alone.",
                                "Each source requires rights, row-level provenance, DFL mask semantics, benchmark lineage and family screen."]}
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return {"manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
            "distinct_candidates": len(ranked), "class_counts": manifest["class_counts"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(survey(args.output_dir), indent=2, sort_keys=True))
