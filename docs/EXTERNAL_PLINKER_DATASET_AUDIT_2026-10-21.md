# External P-Linker dataset audit

Date: 21 October 2026
Input: `Linker.zip`
SHA-256: `572e2b672b5b7da39ac7b0e8466a76de8ffdbfa8d893e1c2b3742f4cda6632f`

## Decision

Do **not** merge this archive into LinkerForge's training or curation data. It
is not the missing fusion-linker compatibility dataset and does not unblock the
compatibility ranker. Keep it quarantined as a possible source for reproducing
and auditing a domain-gap mining method after ownership, licensing and source
provenance are resolved.

Nothing in the archive was executed. The audit reads selected data files
directly from the ZIP and treats archived documentation and code as untrusted
source material.

## What the archive contains

- 9,166 ZIP entries and 981,130,558 uncompressed bytes.
- 8,986 entries, including directories, belong to an embedded Python virtual
  environment.
- 91 entries, including directories, belong to an embedded Git repository.
- 7,892 rule-mined `silver` domain-gap records.
- 10 `high` literature/PDB records, despite the data card claiming 11 gold
  records and 7,903 total records.
- 7,899 raw mined records, seven more than the sequence-enriched silver file.
- A 672,806,868-byte InterPro parquet shard and 24 PDB files.
- No licence file.

The embedded Git remote names `ramlovespaneer/p-linker`, at commit
`83b1a9a5b38b5e16a73ff397f6d6783191d38880`, but the repository was not
accessible through GitHub during this audit. The archive therefore is not a
substitute for a durable, licensed upstream source.

## Silver-record findings

| Check | Result |
| --- | ---: |
| Records / unique accessions | 7,892 / 7,892 |
| Unique linker sequences | 6,702 |
| Linker-length range | 3–40 aa |
| Within LinkerForge's frozen 5–30 aa range | 5,049 |
| Outside that range | 2,843 |
| Same-Pfam-family domain pairs | 7,890 |
| Cross-family pairs | 2 |
| Seven-transmembrane-family rows | 5,057 |
| Full sequences containing ambiguous residues | 241 |
| Linkers containing ambiguous residues | 80 |
| Coordinate/slice consistency errors | 0 |
| Record-level release/retrieval/licence/source metadata | 0 |

Only seven domain-family pair combinations occur. The five dominant Pfam
accessions are PF00001, PF00002, PF00003, PF00004 and PF00005. InterPro identifies
the first three as seven-transmembrane receptor families, PF00004 as an AAA
ATPase domain and PF00005 as an ABC-transporter domain. This narrow, primarily
same-family repeat collection is not representative of arbitrary soluble
two-domain fusion design and conflicts with LinkerForge's frozen soluble-domain
scope.

The silver coordinates are internally consistent, so the file can be useful
for reproducing its mining query. However, an annotation gap between two Pfam
matches is not itself a validated flexible linker and is not evidence that the
sequence is compatible with a different engineered domain pair.

## Literature/PDB-record findings

The ten records labelled `high` are not ready for import. All ten disagree with
their own declared domain boundaries under the natural invariant
`domain A end + 1 = linker start` and `linker end + 1 = domain B start`.
The archived builder locates the linker by calling `full_sequence.find()` on a
short motif. This chooses the first occurrence, which demonstrably mislocates
short or repeated sequences such as `GT`, `GGGGG` and `GGSGG`.

These records may be useful only as a manually reconstructed public comparison
panel after checking each construct against its primary paper, PDB entity,
chain, residue numbering and deposited construct sequence.

## Comparison with the existing 27,766-record corpus

The 7,892 silver records have zero exact protein-identifier overlap and zero
unique linker-sequence overlap with the current composite under the stored
identifiers. Novelty alone is not a quality signal. Adding them now would
increase quantity while weakening scope consistency and provenance.

## Permitted future use

1. Recover an accessible upstream source and explicit redistribution licence.
2. Pin the exact InterPro/Pfam and UniProt releases and retrieval dates.
3. Re-run the mining query from clean dependencies; do not reuse the archived
   `.venv` or embedded `.git` directory.
4. Remove transmembrane/out-of-scope proteins, ambiguous sequences and linkers
   outside 5–30 aa.
5. Reconcile overlapping domain annotations and cluster by protein/domain
   family before sampling.
6. Send only a stratified subset to independent linker-identity review.
7. Keep any accepted records as natural/weak observations. Do not label them as
   engineered-fusion success or compatibility outcomes.

Even after those steps, a separate engineered-fusion outcome dataset or a
carefully defined within-context contrastive ranking dataset remains necessary
for scientific compatibility-ranker training.

Authoritative family metadata: [PF00001](https://www.ebi.ac.uk/interpro/entry/pfam/PF00001/),
[PF00002](https://www.ebi.ac.uk/interpro/entry/pfam/PF00002/),
[PF00003](https://www.ebi.ac.uk/interpro/entry/pfam/PF00003/),
[PF00004](https://www.ebi.ac.uk/interpro/entry/pfam/PF00004/), and
[PF00005](https://www.ebi.ac.uk/interpro/entry/pfam/PF00005/).
