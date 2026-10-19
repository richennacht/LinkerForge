# External P-Linker archive intake audit (2026-10-19)

Status: **quarantined as candidate intervals, not admitted as DFL training positives or a validation set**. This audit does not assert that every candidate is biologically non-linker. The source ZIP was read without extraction or executing bundled code. Reproduce with:

```powershell
& '.venv-gpu\Scripts\python.exe' scripts\audit_external_linker_zip.py --archive 'C:\Users\Perthro\Downloads\Linker.zip' --split-manifest data\processed\dfl_b1_v2_family_split\manifest.json
```

Archive SHA-256: `572e2b672b5b7bda39ac7b0e8466a76de8ffdbfa8d893e1c2b3742f4cda6632f` (648,952,221 bytes; 9,166 ZIP members). Many members are bundled `.venv` and `.git` contents; neither is required or safe to execute for intake. JSONL member hashes: gold `afebab5b13b9a2fda043f12cf72f5c4e03f0902da44cdf6f61bc77b332c1be3f`, silver `dcdf62f3febac83c92648545193325f8a9496d0c4ebc6445d7cdb99e958133c7`, review `b04739a118134a51f94e081568ff2e8686564ee228094ac69fd6eada1f837e24`, pre-sequence candidates `f98b03d551fb9cb261e9e5f3164751caa9563417ccb443a2595fd1407e5e1979`.

## What the records actually establish

| File | Rows | Supported interpretation |
| --- | ---: | --- |
| `dataset_v1.jsonl` (called gold) | 10 | PDB-linked engineered-fusion *candidates* requiring independent paper, entity, chain and junction review. |
| `dataset_v1_needs_review.jsonl` | 12 | Explicit review queue; not training labels. |
| `dataset_v2.jsonl` (called silver) | 7,892 | UniProt sequences with a 3–40-residue gap between the only two qualifying Pfam **matches**. These are geometrically extracted intervals, not annotated disordered flexible linkers. |
| `mined_linkers.jsonl` | 7,899 | Pre-sequence-fetch candidates, not an additional label tier. |

The archive's data card says 11 gold and 7,903 total. The actual gold and silver files contain 10 + 7,892 = **7,902** rows. The silver generator (`pipeline/io/mine_linkers.py`) filters two Pfam matches of length ≥30, protein length 100–1,000, and gap 3–40. It has no disorder, flexibility, biological-junction, or explicit-negative verification. The sequence-fetcher labels the gap `silver`; that is a heuristic provenance tag, not a DFL observation.

The silver intervals pass basic sequence-slice and coordinate arithmetic checks, but the distribution is highly concentrated: 4,211 `PF00001/PF00001`, 2,333 `PF00004/PF00004`, 666 `PF00002/PF00002`, 502 `PF00005/PF00005`, 178 `PF00003/PF00003`, and two singletons—just **seven** ordered Pfam pairs. `PF00001` is the seven-transmembrane rhodopsin receptor family ([NCBI CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=pfam00001)); 4,212 rows involve it. A gap between two matches to the same family is not necessarily a linker between two domains. 1,724 candidate gap sequences are ≥50% hydrophobic by the audit's simple `AILMFWVY` count, and 80 contain a noncanonical amino-acid character. There are 7,786 unique full sequences among 7,892 rows (106 duplicate-sequence rows). These observations warrant biological review, not wholesale rejection or promotion.

The purported gold file has eight unique full sequences and six unique linker motifs across ten rows. Two motifs are two residues long, outside both the archive's 3–40 silver rule and LinkerForge's 5–30 target. Relative to the file's own domain coordinates, 7/10 domain-A ends are not immediately before the extracted linker, 10/10 domain-B starts are not immediately after it, and 5/10 domain-B ends exceed the supplied full sequence. The bundled builder hardcodes PDB IDs and coordinates but then uses the first `full_seq.find(linker_seq)` to locate a motif; this is not a chain/junction alignment. Moreover, its source string cites “Chen, Zaro & Shen 2013 — Table 1” for every row, whereas [that paper's Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC3726540/) reports properties of natural linkers, not these PDB-specific junction annotations. Its citation therefore does not independently verify these ten rows.

## Leakage and use decision

Against the B1 V2 manifest (`d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`), the 7,892 silver rows have zero exact accession and zero exact full-sequence matches. **This is not a family/homology or local-sequence firewall**, and no prospective benchmark audit was performed. None of these rows is currently eligible as a *verified DFL-positive training record* or independent precision validation record. Do not convert unannotated residues to negatives, treat unknowns as 0, use gap geometry as a DFL label, merge these rows into the existing B1 split, or cite the nominal 7,903 as a supervised DFL sample size.

Preserve the original ZIP as a source artifact. If useful, a separate versioned *candidate mining* table may hold accession, sequence and Pfam-gap coordinates with `annotation_status=unknown`, source-member hash and method version; it must stay outside positive/negative training masks and benchmark metrics. To promote any interval, independently verify the original protein sequence and annotation, disorder/flexibility and actual inter-domain junction, negative/unknown semantics, rights, duplicate/family clusters, and exclusion against all held-out targets. For engineered fusions, separately verify each cited primary experiment, PDB polymer entity/chain and residue-number mapping; those labels remain distinct from natural DFL.

No training or precision claim was made from this archive. In particular, these records cannot validate a 99% precision mining threshold because they do not supply adjudicated DFL negatives or a label-sealed test.
