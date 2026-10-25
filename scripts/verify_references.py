"""Verify every bibliography entry against Crossref / Europe PMC and regenerate the .bib and .bbl.

Nothing in `paper/references.bib` is typed by hand: author lists, volumes, pages and years come from
the registry record for the DOI (or PMCID).  An entry is only marked `verified` when the registry
title matches the title we expected to cite (difflib ratio >= 0.80).  Entries that do not resolve
or do not match are written to the report and are NOT emitted to the .bib.

Usage:  python scripts/verify_references.py [--tex paper/linkerforge_arxiv_2026.tex]
"""

from __future__ import annotations

import argparse
import difflib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UA = "LinkerForge-reference-check/1.0 (https://github.com/richennacht/LinkerForge)"

# key -> (identifier kind, identifier, expected title)
REFERENCES: dict[str, tuple[str, str, str]] = {
    "crasto2000linker": ("doi", "10.1093/protein/13.5.309", "LINKER: a program to generate linker sequences for fusion proteins"),
    "george2002analysis": ("doi", "10.1093/protein/15.11.871", "An analysis of protein domain linkers: their classification and role in protein folding"),
    "bae2005hmm": ("doi", "10.1093/bioinformatics/bti363", "Prediction of protein interdomain linker regions by a hidden Markov model"),
    "miyazaki2006neural": ("doi", "10.1186/1471-2105-7-323", "Identification of putative domain linkers by a neural network"),
    "ebina2011drop": ("doi", "10.1093/bioinformatics/btq700", "DROP: an SVM domain linker predictor trained with optimal features selected by random forest"),
    "shatnawi2014rf": ("doi", "10.1186/1471-2105-15-S16-S8", "Protein inter-domain linker prediction using Random Forest and amino acid physiochemical properties"),
    "meng2016dflpred": ("doi", "10.1093/bioinformatics/btw280", "DFLpred: High-throughput prediction of disordered flexible linker regions in protein sequences"),
    "peng2020apod": ("doi", "10.1093/bioinformatics/btaa808", "APOD: accurate sequence-based predictor of disordered flexible linkers"),
    "wu2023transdfl": ("doi", "10.1016/j.gpb.2022.10.004", "TransDFL: Identification of Disordered Flexible Linkers in Proteins by Transfer Learning"),
    "idpelm2024": ("doi", "10.1021/acs.jcim.3c01202", "Accurate and Fast Prediction of Intrinsically Disordered Protein by Multiple Protein Language Models and Ensemble Learning"),
    "dld2026": ("doi", "10.1016/j.jmb.2026.169994", "Expanding the Landscape of Disordered Flexible Linkers: A Structural and Computational Framework for DLD Dataset Assembly"),
    "linkerpred2025": ("doi", "10.1101/2025.10.29.685455", "LINKER-Pred: A Deep Learning Method and Web-Server for the Prediction of Disordered Flexible Linkers in Proteins"),
    "liu2015synlinker": ("doi", "10.1093/bioinformatics/btv447", "SynLinker: an integrated system for designing linkers and synthetic fusion proteins"),
    "chen2013fusion": ("doi", "10.1016/j.addr.2012.09.039", "Fusion protein linkers: property, design and functionality"),
    "arai2001helical": ("doi", "10.1093/protein/14.8.529", "Design of the linkers which effectively separate domains of a bifunctional fusion protein"),
    "vanrosmalen2017gs": ("doi", "10.1021/acs.biochem.7b00902", "Tuning the flexibility of glycine-serine linkers to allow rational design of multidomain proteins"),
    "hubbard1994proteolysis": ("doi", "10.1002/pro.5560030505", "Modeling studies of the change in conformation required for cleavage of limited proteolytic sites"),
    "pauling1951sheet": ("doi", "10.1073/pnas.37.5.251", "The pleated sheet, a new layer configuration of polypeptide chains"),
    "pauling1951helix": ("doi", "10.1073/pnas.37.4.205", "The structure of proteins: two hydrogen-bonded helical configurations of the polypeptide chain"),
    "kratky1949wlc": ("doi", "10.1002/recl.19490681203", "Röntgenuntersuchung gelöster Fadenmoleküle"),
    "shrake1973sasa": ("doi", "10.1016/0022-2836(73)90011-9", "Environment and exposure to solvent of protein atoms. Lysozyme and insulin"),
    "lee1971sasa": ("doi", "10.1016/0022-2836(71)90324-X", "The interpretation of protein structures: Estimation of static accessibility"),
    "bondi1964vdw": ("doi", "10.1021/j100785a001", "van der Waals Volumes and Radii"),
    "tien2013maxasa": ("doi", "10.1371/journal.pone.0080635", "Maximum Allowed Solvent Accessibilites of Residues in Proteins"),
    "rost1994rsa": ("doi", "10.1002/prot.340200303", "Conservation and prediction of solvent accessibility in protein families"),
    "chen2021cirpred": ("doi", "10.1186/s12859-021-04403-1", "CirPred, the first structure modeling and linker design system for circularly permuted proteins"),
    "berrondo2018loodo": ("doi", "10.1002/prot.25445", "A computational method for the design of nested proteins by loop-directed domain insertion"),
    "prodomino2025": ("doi", "10.1038/s41592-025-02741-z", "Rational engineering of allosteric protein switches by in silico prediction of domain insertion sites"),
    "splitseekpro2026": ("doi", "10.1038/s41467-026-74059-z", "SplitSeek-Pro: accurate prediction of splittable sites on protein structures"),
    "multimodal2026": ("doi", "10.1016/j.synbio.2026.02.003", "Functional customization of peptide linkers in fusion proteins through multimodal deep learning approach"),
    "mcbride2026pose": ("doi", "10.1002/prot.70062", "Predicting Pose Distribution of Protein Domains Connected by Flexible Linkers Is an Unsolved Problem"),
    "insertionfamilies2026": ("doi", "10.1002/pro.70586", "Family-level specialization in protein domain insertion architectures"),
    "clopper1934": ("doi", "10.1093/biomet/26.4.404", "The use of confidence or fiducial limits illustrated in the case of the binomial"),
    "ltt2021": ("doi", "10.1214/24-AOAS1998", "Learn then Test: Calibrating Predictive Algorithms to Achieve Risk Control"),
    "bekker2020pu": ("doi", "10.1007/s10994-020-05877-5", "Learning from positive and unlabeled data: a survey"),
    "esm2_2023": ("doi", "10.1126/science.ade2574", "Evolutionary-scale prediction of atomic-level protein structure with a language model"),
    "pedelacq2006sfgfp": ("doi", "10.1038/nbt1172", "Engineering and characterization of a superfolder green fluorescent protein"),
    "shu2006mcherry": ("doi", "10.1021/bi060773l", "Novel chromophores and buried charges control color in mFruits"),
    "sillitoe2021cath": ("doi", "10.1093/nar/gkaa1079", "CATH: increased structural coverage of functional space"),
    "uniprot2023": ("doi", "10.1093/nar/gkac1052", "UniProt: the Universal Protein Knowledgebase in 2023"),
    "quaglia2022disprot": ("doi", "10.1093/nar/gkab1082", "DisProt in 2022: improved quality and accessibility of protein intrinsic disorder annotation"),
    "paysan2023interpro": ("doi", "10.1093/nar/gkac993", "InterPro in 2022"),
    "berman2000pdb": ("doi", "10.1093/nar/28.1.235", "The Protein Data Bank"),
    "caid2linker2024": ("pmcid", "PMC10968229", "Assessment of Disordered Linker Predictions in the CAID2 Experiment"),
    "caid3_2025": ("pmcid", "PMC12750029", "Critical Assessment of Protein Intrinsic Disorder Round 3 - Predicting Disorder in the era of Protein Language Models"),
}


def _get_json(url: str, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            time.sleep(1.5 * (attempt + 1))
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    return None


def _norm(title: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", re.sub(r"<[^>]+>", "", title).lower()).strip()


def fetch_crossref(doi: str) -> dict | None:
    data = _get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi))
    if not data:
        return None
    m = data["message"]
    issued = (m.get("issued") or {}).get("date-parts") or [[None]]
    return {
        "title": (m.get("title") or [""])[0],
        "authors": [(a.get("family", ""), a.get("given", "")) for a in m.get("author", []) if a.get("family")],
        "container": (m.get("container-title") or [""])[0],
        "volume": m.get("volume"),
        "issue": m.get("issue"),
        "pages": m.get("page") or m.get("article-number"),
        "year": issued[0][0],
        "doi": doi,
        "type": m.get("type"),
        "publisher": m.get("publisher"),
    }


def fetch_europepmc(pmcid: str) -> dict | None:
    q = urllib.parse.quote(f"PMCID:{pmcid}")
    data = _get_json(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&resultType=core&format=json"
    )
    results = ((data or {}).get("resultList") or {}).get("result") or []
    if not results:
        return None
    r = results[0]
    authors = []
    for a in (r.get("authorList") or {}).get("author", []):
        if a.get("lastName"):
            authors.append((a["lastName"], a.get("firstName", "")))
    j = (r.get("journalInfo") or {}).get("journal") or {}
    return {
        "title": (r.get("title") or "").rstrip("."),
        "authors": authors,
        "container": j.get("title") or "",
        "volume": (r.get("journalInfo") or {}).get("volume"),
        "issue": (r.get("journalInfo") or {}).get("issue"),
        "pages": r.get("pageInfo"),
        "year": int(r["pubYear"]) if r.get("pubYear") else None,
        "doi": r.get("doi"),
        "type": "article",
        "publisher": None,
        "pmcid": pmcid,
    }


def _tex(s: str) -> str:
    s = re.sub(r"<[^>]+>", "", s or "")
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_")):
        s = s.replace(a, b)
    return s


def to_bibtex(key: str, rec: dict) -> str:
    authors = " and ".join(f"{fam}, {giv}".strip(", ") for fam, giv in rec["authors"]) or "Anonymous"
    kind = "article" if rec.get("container") else "misc"
    fields = [("author", _tex(authors)), ("title", "{" + _tex(rec["title"]) + "}")]
    if rec.get("container"):
        fields.append(("journal", _tex(rec["container"])))
    for fld in ("volume", "issue", "pages"):
        if rec.get(fld):
            fields.append(("number" if fld == "issue" else fld, _tex(str(rec[fld]).replace("-", "--") if fld == "pages" else str(rec[fld]))))
    if rec.get("year"):
        fields.append(("year", str(rec["year"])))
    if rec.get("doi"):
        fields.append(("doi", rec["doi"]))
    if rec.get("pmcid"):
        fields.append(("note", f"PMCID: {rec['pmcid']}"))
    if kind == "misc" and rec.get("publisher"):
        fields.append(("howpublished", _tex(rec["publisher"] + " preprint")))
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields)
    return f"@{kind}{{{key},\n{body}\n}}\n"


def bbl_item(key: str, rec: dict) -> str:
    def initials(given: str) -> str:
        return " ".join(p[0] + "." for p in re.split(r"[ \-]+", given) if p)

    names = [f"{initials(g)} {_tex(f)}".strip() for f, g in rec["authors"]]
    if len(names) > 6:
        names = names[:6] + ["et~al."]
    auth = ", ".join(names[:-1]) + (", and " if len(names) > 2 else " and ") + names[-1] if len(names) > 1 else (names[0] if names else "")
    parts = [f"\\bibitem{{{key}}}", f"{auth}.", f"\\newblock {_tex(rec['title'])}."]
    venue = _tex(rec.get("container") or "preprint")
    vol = f"{rec['volume']}" if rec.get("volume") else ""
    if rec.get("issue"):
        vol += f"({rec['issue']})"
    pages = f":{_tex(str(rec['pages']).replace('-', '--'))}" if rec.get("pages") else ""
    parts.append(f"\\newblock {{\\em {venue}}}, {vol}{pages}, {rec.get('year', '')}.")
    if rec.get("doi"):
        parts.append(f"\\newblock \\texttt{{doi:{_tex(rec['doi'])}}}.")
    return "\n".join(parts) + "\n"


def cited_keys_in_order(tex_path: Path) -> list[str]:
    text = tex_path.read_text(encoding="utf-8")
    keys: list[str] = []
    for grp in re.findall(r"\\cite[a-z]*\{([^}]*)\}", text):
        for k in grp.split(","):
            k = k.strip()
            if k and k not in keys:
                keys.append(k)
    return keys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", type=Path, default=ROOT / "paper" / "linkerforge_arxiv_2026.tex")
    ap.add_argument("--bib", type=Path, default=ROOT / "paper" / "references.bib")
    ap.add_argument("--bbl", type=Path, default=ROOT / "paper" / "linkerforge_arxiv_2026.bbl")
    ap.add_argument("--report", type=Path, default=ROOT / "reports" / "reference_verification_2026-10-19.json")
    args = ap.parse_args()

    sys.stdout.reconfigure(encoding="utf-8")
    results: dict[str, dict] = {}
    for key, (kind, ident, expected) in REFERENCES.items():
        rec = fetch_crossref(ident) if kind == "doi" else fetch_europepmc(ident)
        time.sleep(0.2)
        if rec is None:
            results[key] = {"status": "not_found_in_registry", "identifier": ident, "expected_title": expected}
            continue
        ratio = difflib.SequenceMatcher(None, _norm(rec["title"]), _norm(expected)).ratio()
        contains = _norm(expected) in _norm(rec["title"]) or _norm(rec["title"]) in _norm(expected)
        ok = ratio >= 0.80 or contains
        if not ok and kind == "doi":
            # Some publisher records carry a placeholder title; retry via Europe PMC by DOI.
            alt = _get_json(
                "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
                + urllib.parse.quote(f"DOI:{ident}")
                + "&resultType=core&format=json"
            )
            hits = ((alt or {}).get("resultList") or {}).get("result") or []
            if hits:
                alt_title = (hits[0].get("title") or "").rstrip(".")
                alt_ratio = difflib.SequenceMatcher(None, _norm(alt_title), _norm(expected)).ratio()
                if alt_ratio >= 0.80:
                    rec = dict(rec, title=alt_title)
                    ratio, ok = alt_ratio, True
        results[key] = {
            "status": "verified" if ok else "title_mismatch",
            "identifier": ident,
            "title_similarity": round(ratio, 3),
            "expected_title": expected,
            "registry": rec,
        }
        print(f"{key:26s} {results[key]['status']:22s} sim={ratio:.2f}  {rec['title'][:70]}")

    verified = {k: v for k, v in results.items() if v["status"] == "verified"}
    args.bib.write_text(
        "% Generated by scripts/verify_references.py from Crossref / Europe PMC records.\n"
        "% Do not edit by hand.\n\n"
        + "\n".join(to_bibtex(k, v["registry"]) for k, v in verified.items()),
        encoding="utf-8",
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(results, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")

    bad = [k for k, v in results.items() if v["status"] != "verified"]
    print(f"\nverified={len(verified)} unverified={len(bad)}: {bad}")

    if args.tex.exists():
        order = cited_keys_in_order(args.tex)
        missing = [k for k in order if k not in verified]
        if missing:
            print("CITED BUT NOT VERIFIED:", missing)
        items = [bbl_item(k, verified[k]["registry"]) for k in order if k in verified]
        args.bbl.write_text(
            "\\begin{thebibliography}{%d}\n\n%s\n\\end{thebibliography}\n" % (len(items), "\n".join(items)),
            encoding="utf-8",
        )
        print(f"wrote {args.bbl.name} with {len(items)} entries in citation order")
        return 1 if missing else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())

