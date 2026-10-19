# B1 classifier split feasibility audit — 18 October 2026

Status: **B1 remains open.** This is a CPU-only, sequence-only source-capacity
audit, not a split, training run, benchmark score or approval to mine proteins.
The executable is `scripts/audit_dfl_b1_source_capacity.py`; synthetic tests
are `tests/test_audit_dfl_b1_source_capacity.py`. Reproduce from a checkout
with the locally pinned ignored source files:

Executable SHA-256: `55d59dc61fb08ccffa4e4edb733440ee2d252fe141cf35927810de4afdad044d`.

```powershell
.venv-gpu\Scripts\python.exe scripts\audit_dfl_b1_source_capacity.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --include-idpfunnet
.venv-gpu\Scripts\python.exe -m unittest tests.test_audit_dfl_b1_source_capacity tests.test_preflight_dfl_family_training -q
```

The optional weak-source reconstruction is reproducible from the pinned
annotated-subsequence catalog (input data and generated full-protein records
remain gitignored):

```powershell
.venv-gpu\Scripts\python.exe scripts\prepare_uniprot_explicit_flexible_candidates.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --output C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\data\processed\uniprot_explicit_flexible_linkers_v1
.venv-gpu\Scripts\python.exe scripts\audit_dfl_b1_source_capacity.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion --weak-candidates C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\data\processed\uniprot_explicit_flexible_linkers_v1
```

The preparer refuses to overwrite an existing output directory. Its source
code SHA-256 is `da171da3c3717a1cff57ecb1bb073b0d4b7613a3d44f44ed040153332de8dff2`;
schema SHA-256 is `3acaf7778c61e8e125d68f91037e832760b67fd42910aaf87d211476bde305e8`.

The script verifies source SHA-256s and counts before analysis. CAID3's official
source SHA-256 is `4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`
(31 sequences; the third annotation line is consumed without parsing). The
DP94 sequence-only manifest SHA-256 is
`fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`
(98 sequences). No CAID3/DP94 labels, predictions, thresholds or scores were
used. The screen is BLOSUM62 Smith–Waterman, gap 10/1, at least 50 aligned
residues, 25% identity and 50% shorter-sequence coverage, plus exact ID,
sequence and containment. The official CAID3/DP94 benchmark annotation tracks
were not opened. Candidate *source* annotations may refer to proteins later
removed by the sequence firewall; they are not used for model selection or
admitted to a training split.

| DFL-labelled training candidate | Source SHA-256 | Observed-positive proteins | Positive proteins with no direct hit to either registry |
| --- | --- | ---: | ---: |
| Published DFL train | `a1a5f3b3b30548ad65f6f378535c7c732272b8c9b60048b38f6500fca9bb1aed` | 74 of 158 | 9 |
| DisProt positive-only pool | `0265b062ca9ff88cde2576789b613860355f33842fa9714788ba1105a19df42a` | 159 of 159 | 23 |
| DisProt positive-candidate catalog | `2bc6726545f47256bcec70376bd51d689db0b5c80a64f45a221b9d88bb0c4c8d` | 205 of 205 | 27 |
| IDPFunNet author TR552 training file (pinned revision) | `89e374e30924a277fb4bdcf1b7a82776017e4ab70397757b185fa9d0b16954b2` | 124 of 552 | 17 |

Those rows overlap. Together they yield **40 unique DFL-positive full-protein
sequences** with no *direct* hit to CAID3 or DP94. A conservative transitive
graph across 404 unique source-plus-benchmark sequences contains 1,805 local
edges and leaves only **5** positive sequences outside benchmark-touching
components. The 5 are in five singleton components. This graph is a sensitivity
analysis, **not** a validated biological or Pfam-domain family assignment:
the permissive local rule can connect unrelated low-complexity proteins.
Conversely, the 40 direct non-hits are not certified family-independent. A
qualified B1 split cannot simply be declared from either number.

## Other routes checked

1. The 24,763-record annotated-linker catalog has **zero** full-protein DFL
   segmentation-training-ready records according to its local readiness
   manifest (SHA-256 `d383255e2d84214317ee1ba2507df6e9cdaee56a817ac6b65007bd30d3e77903`).
   It contains labelled linker subsequences, including generic/structural
   classes, not automatically `1/0/x` masks for the DFL task. No inferred
   domain-gap candidates were promoted to ground truth.
2. Within that catalog, 884 distinct UniProt accessions have a region note
   containing `Flexible linker` or `Disordered linker`; just 10 are marked
   with a direct ECO:0000269/PubMed evidence tier. We implemented
   `scripts/prepare_uniprot_explicit_flexible_candidates.py` with formal
   `schemas/uniprot_explicit_flexible_linker_candidate.v1.schema.json`, fetched
   current full sequences from UniProt in 23 hash-recorded batches, and
   verified **all 884** annotated peptides at their exact 1-based coordinates
   (zero missing/mismatched records). The local ignored candidate manifest is
   SHA-256 `6b0c4f704ac9502d2bf94bf95ae0f64d67483d7d04e8421183414adbed6de614`;
   output JSONL SHA-256 `83aaa2e25eabebd6dd5b29f132c307b44d10b06aad7721462ca36f7f8c7256a3`.
   Every other residue remains `x` unknown. Direct CAID3/DP94 screening leaves
   **11** unique full-protein sequences, only **1** from the PubMed tier.
   These annotations are not independently DFL-adjudicated, and most are
   computational/homology-derived. They cannot be silently merged into the
   gold DFL train/development set; rights and family review remain open.
3. DLD's thousands of structural IDL/DDL labels remain a separate task and
   were not converted to DFL positives or negatives. Published `0=non-DFL`
   tracks also lack evidence that every zero is an adjudicated negative; the
   audit uses their observed positives only.

## Decision and next gate

Do not mark B1 complete or launch a GPU job on the old split. The current
`check_detector_autoresearch.py` still returns `data_integrity_audit_only` for
a DP94 sequence collision. Completing B1 requires a source-qualified pool of
full-protein DFL masks, an independently justified family/domain clustering
method and manifest, a new versioned train/development/auxiliary split, and
both required preflights passing on that split. The most promising next
investigation is source-level reconstruction of the explicit UniProt
flexible-linker records, retaining evidence tier and verifying full sequences,
spans, rights and family disjointness. Those weak annotations must not be
silently treated as experimentally adjudicated DFL ground truth.

DP94 is still a candidate, not a certified untouched DFL test. CAID3 remains
a spent retrospective. No benchmark labels may enter this decision.
