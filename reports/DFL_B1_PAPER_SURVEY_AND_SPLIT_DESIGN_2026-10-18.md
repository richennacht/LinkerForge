# B1 literature survey and split-design checkpoint — 18 October 2026

Status: **candidate split design, not a family-disjoint split or training
authorization**. The current binding preflights still reject the old split.
No CAID3/DP94 labels, predictor scores or per-protein errors were used here.

## What the primary literature says

| Source | Direct finding relevant to this decision | Consequence for LinkerForge |
| --- | --- | --- |
| [DFLpred (Meng & Kurgan, 2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4908364/) | DFLs can be intra- as well as inter-domain; the paper built 204 functionally annotated proteins, kept unannotated residues out of assessment, and allocated BLASTClust clusters rather than individual sequences across folds. | The target is curated *disordered flexible linker*, not every domain gap. Keep `x` masked and split whole families. Do not require a two-domain protein for every positive. |
| [APOD (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7773485/) | The authors assigned 194 BLASTClust clusters of 248 proteins as units to TR166 and TE82. | Component-level allocation has direct precedent for DFL prediction, but the paper's 25%-identity BLASTClust procedure is not evidence that our current BLOSUM62 Smith–Waterman 25%/50%-shorter rule measures the same relationship. |
| [MMseqs2 user guide](https://mmseqs.com/latest/userguide.pdf) | For full-length proteins, its bidirectional coverage mode counts aligned overlap relative to the longer sequence; the guide recommends high coverage to preserve multidomain architecture. | Use a reciprocal full-length screen as one *candidate* homology criterion, while separately testing shared domains. The chosen version, mask policy and thresholds must be frozen from sequence-only controls. |
| [DataSAIL (2025)](https://doi.org/10.1038/s41467-025-58606-8) | Biomolecular ML splits should reflect the similarity structure relevant to intended inference; metadata-only or naive random partitions can leak. | Sequence clusters are necessary but not sufficient. Audit domain families, accession/isoform lineage, near duplicates and the intended unseen-protein use case. |
| [DisProt curation report (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767923/) and [curation practices (2024)](https://pubmed.ncbi.nlm.nih.gov/38507044/) | DisProt functional regions are literature-curated with explicit experimental evidence and provenance; coverage is incomplete. | Treat curated DFL regions as positives; do not turn unannotated residues into negatives. Maintain reference/evidence IDs and create a separately adjudicated no-linker set before making a precision claim. |
| [CAID3 methods (2025)](https://www.research.unipd.it/retrieve/260dc782-65d7-4f89-809d-e77d87ed983a/2025_caid3.pdf) | The Linker reference derives from DisProt flexible-linker annotations; the test is a temporal DisProt delta. | CAID3 is task-relevant but has already been viewed here, so it remains retrospective and cannot select another architecture or threshold. |

This survey does **not** endorse an exact numerical threshold for our split.
Protein-language-model pretraining exposure and published-competitor training
lineage are additional limits on any eventual independent comparison; see
[Hermann et al. (2024)](https://proceedings.mlr.press/v261/hermann24a.html).

## Executable sequence-only sensitivity check

`scripts/audit_dfl_full_length_screen.py` verifies the source and registry
hashes, reads only 279 unique observed-DFL-positive candidate sequences from
three pinned local sources, and compares two rules against CAID3 (31) and DP94
(98) *sequences*. The proposed diagnostic rule uses a gapped BLOSUM62 local
alignment but requires at least 30% identity, 50 aligned columns and at least
80% of **each** full protein's residues to align. Unlike the legacy
alignment-length/shorter-length ratio, each coverage fraction is bounded by
1. It does not yet include domain-family evidence and is **not** a preflight
replacement.

| Screen | Unique positive sequences flagged against either registry | Unflagged |
| --- | ---: | ---: |
| Binding local rule (25% identity; alignment length / shorter length ≥0.5) | 243 / 279 | 36 / 279 |
| Proposed full-length sensitivity rule (30% identity; both residue coverages ≥0.8) | 37 / 279 | 242 / 279 |

All 37 proposed-rule flags also occur under the binding rule. The 206 changed
decisions are **not** automatically independent families: a shared Pfam domain
may be real homology even when two complete multidomain proteins fail the 80%
screen. Conversely, the prior composition-preserving shuffle null audit found
many legacy-rule hits on rearranged sequences; see
`DFL_B1_FIREWALL_REPAIR_PROPOSAL_2026-10-18.md`. These two controls establish
that family adjudication has high leverage. They do not authorize dropping
the 206 exclusions.

### Independent MMseqs2 sensitivity run

The official Windows MMseqs2 release 18-8cc5c (binary commit
`8cc5ce367b5638c4306c2d7cfc652dd099a4643f`; downloaded ZIP SHA-256
`9752087ca4e25e93b830775832a75044dc94ae348d66d4681f9ff0355eecfcbf`)
was run on the same **408** exported sequences: 279 unique DFL-positive
candidates plus 31 CAID3 and 98 DP94 sequence-only records. The export
manifest SHA-256 is
`a543631aebf7684afdea7e2196aa40f3bff0c08394607a4c5ed8e0ce1bc12075`;
FASTA SHA-256 is
`2c50b594341c92a048989ec272db852bc92eb250b1e9440329d0e3889dcb8ed9`.
The official `easy-cluster` workflow failed on this Windows host because it
attempted a `/bin/sh` script; the lower-level `createdb`, `prefilter`, `align`,
`clust --cluster-mode 1`, and `createtsv` modules ran successfully. The
pre-filter used default compositional-bias correction and `--mask 1`; align
used `--min-seq-id 0.3 -c 0.8 --cov-mode 0` and the default E-value cutoff
0.001. Clusters were connected components. The resulting TSV SHA-256 is
`374c3eb97059122a7de5ccd93b7aa56a65f5a13dad661bcea77682b63a1d11e6`.

The hash-checked `scripts/summarize_dfl_mmseqs_audit.py` found **333
components**, 110 benchmark-touching and 223 candidate-only; **37 candidate
sequences** occur in benchmark-touching components and **242** do not. This
agrees in count with the separate full-length parasail diagnostic. It is
still not a domain-family certificate: at high reciprocal coverage, two
multidomain proteins can share a genuine domain without clustering together.
The old protected-data preflight and its similarity threshold were not changed.

The [DisoFLAG authors' dataset documentation](https://bliulab.net/DisoFLAG/download)
explicitly defines line 9 as DFL 1/0 annotations for its DP93/DP94/CAID2
files, and the [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10762911/)
reports DP94 linker evaluation. Therefore the earlier shorthand “DP94 is
not a DFL test” is incorrect. DP94 *is DFL-labelled*, but its per-residue
negative/unknown semantics, suitability as a common-target comparison,
published-model exposure and independence from our sources still need audit.
Its labels remained unopened in this work, and it remains protected under the
binding contract; we did not retire it to enlarge training capacity.

### Domain-family and transitive-component sensitivity

The local DisProt 2026_06 release also includes Pfam, Gene3D and UniRef50
metadata. `scripts/audit_dfl_disprot_family_metadata.py` projects only entry
ID, accession, exact sequence and those family identifiers from the pinned
source JSON (SHA-256
`1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140`).
It does **not** read region or consensus labels. Nevertheless the raw source
file contains annotations, so only its family projection—not that source
file—may enter any future label-sealed evaluation process. The audit requires
an exact ID-and-sequence match and verifies the MMseqs input/cluster hashes.

Of 408 input records, 397 have exact-source metadata; 6 candidate and 5
benchmark records are unresolved. Among candidates, 97 directly share a Pfam
ID with a benchmark sequence, 37 share a Gene3D ID and 16 share a UniRef50
ID. Taking the union with full-length MMseqs clusters directly flags 111
candidates; **163 mapped candidates** are not directly flagged. Crucially,
whole connected components over *all* MMseqs/Pfam/Gene3D/UniRef50 links are
more conservative: **134 candidate sequences in 112 candidate-only
components** remain outside every benchmark-touching component. The largest
remaining component contains 5 candidate sequences. The 6 unresolved
candidate sequences are ineligible for automatic split assignment. This
demonstrates why a full-length screen alone would overstate available
family-disjoint training data.

These 134 are a *candidate pool*, not a certified split: shared domain IDs
can overconnect unrelated proteins and absent metadata can miss remote
homology. The legacy local-alignment firewall still rejects the current
split. Label-mask and source-rights checks, a preregistered rule review,
train/development component allocation, and both binding preflights remain
before B1 can close. No training or benchmark scoring followed this audit.

### Binding-rule ceiling and mask semantics

After the transitive family screen, applying the *unchanged binding direct
local-alignment exclusion* to **both** protected registries leaves just **24
unique observed-positive sequences in 24 components**. The source
appearances overlap: 19 catalog, 17 DisProt-positive-only, and 5 published
DFL-train records. A hash-pinned, sequence-hash-only ignored capacity manifest
has SHA-256
`bca74cae3932d8b086d5d292d18adb41f0c9ce39e9744cc4b05d63fea9912368`.
`scripts/audit_dfl_binding_survivor_masks.py` verifies its source IDs and
sequences against all three pinned input files. Among the 24 are 1,106
union-positive residues, **200 published zero residues whose negative meaning
is not independently verified**, and no cross-source 1/0 conflict. All other
DisProt residues remain `x`. This is a small positive/unknown pilot ceiling,
not a credible precision-calibration set or an adjudicated no-linker set.

For diagnosis *only*, protecting CAID3 alone would leave 108 mapped
sequences under the same family-plus-binding screens; protecting DP94 alone
would leave 34. Those counts show DP94 drives much of the current capacity
loss. Because DP94 is documented as DFL-labelled and remains nominated, the
108-sequence scenario is **not** a permissible replacement split. Removing a
benchmark after seeing these counts would compromise the claimed comparison.

There was also a registry-schema issue: the new split preflight previously
required `unopened` for *every* registry even though CAID3 had already been
viewed retrospectively. This audit repaired that metadata gate: CAID3 must
be `retrospective_spent_exclusion_only`, while DP94 remains `unopened` in the
local project. Both sequence firewalls remain mandatory. A false unopened
CAID3 registry now fails the preflight test.

To reproduce the family summary after exporting the inputs and running the
MMseqs2 modules above (paths are relative to the audit worktree):

```powershell
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' scripts/audit_dfl_disprot_family_metadata.py --source-root 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion' --input-dir artifacts/b1-mmseqs-input-2026-10-18 --cluster-tsv artifacts/b1-mmseqs-modules-2026-10-18/clusters.tsv --cluster-tsv-sha256 374c3eb97059122a7de5ccd93b7aa56a65f5a13dad661bcea77682b63a1d11e6
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' scripts/audit_dfl_binding_survivor_masks.py --source-root 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion' --capacity-manifest artifacts/b1-mmseqs-input-2026-10-18/binding-family-survivors-v2.json --capacity-manifest-sha256 bca74cae3932d8b086d5d292d18adb41f0c9ce39e9744cc4b05d63fea9912368
```

The generated inputs and Windows MMseqs executable live in ignored artifact
directories, not the Git repository; the checked-in scripts and hashes specify
how to rebuild and verify them on a machine with the pinned raw inputs.

## B1 decision from current assets

**No-go for a high-precision mining classifier on the current data.** The
binding two-registry firewall plus exploratory family graph yields an upper
bound of 24 positive proteins, 24 components, and zero independently
adjudicated no-linker proteins. Even if a tiny train/development partition
were mechanically possible, it would not estimate proteome-scale precision.
This 24-count is for the three pinned DFL candidate sources in this audit,
not a census of every future source. The separately audited IDPFunNet TR552
file adds at most 17 direct-registry non-hits *before* family review, with
redistribution rights and negative semantics still unverified; the 884
explicit UniProt flexible-linker candidates add only 11 direct non-hits,
mostly from weaker evidence tiers. See
`reports/DFL_B1_FEASIBILITY_AUDIT_2026-10-18.md`. Neither alternative
supplies an adjudicated no-linker panel.
For intuition, with zero false positives on 299 independent negative units,
the one-sided 95% binomial upper bound on false-positive probability is still
`1 - 0.05^(1/299) = 0.00997` (roughly 1%); correlations or prevalence shift
would worsen the interpretation. Our current audited set has **no such
negative units**, so there is no empirical bound at all.

The scientifically defensible way to make B1 actionable is to retain DP94's
sequence protection, obtain a materially larger independently sourced DFL
positive pool, and commission a separately adjudicated, family-disjoint
no-linker panel with explicit unknown regions. Every new source must be
screened against both protected registries and assigned by component before
training. If the team instead wants to use DP94 proteins for development,
that is a *new, less independent evaluation contract* and cannot be presented
as a blind comparison. Likewise, the reciprocal full-length rule cannot
replace the binding local rule merely because it retains more records.

The B1 infrastructure repair is narrower: the new split preflight now
requires CAID3 to be declared `retrospective_spent_exclusion_only` and DP94
`unopened` in their respective sequence-only registries. This prevents a
false claim that CAID3's labels were never viewed; it does not relax any
sequence or family exclusion and does not make the 24-protein pool train-ready.

Reproduction commands after running
`scripts/prepare_dfl_mmseqs_audit_input.py` with the pinned ignored sources:

```powershell
mmseqs createdb sequences.fasta input
mmseqs prefilter input input prefilter --mask 1 --threads 2
mmseqs align input input prefilter align --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads 2
mmseqs clust input align cluster --cluster-mode 1 --threads 2
mmseqs createtsv input input cluster clusters.tsv
```

Pinned source SHA-256s: published DFL train
`a1a5f3b3b30548ad65f6f378535c7c732272b8c9b60048b38f6500fca9bb1aed`,
DisProt positive-only
`0265b062ca9ff88cde2576789b613860355f33842fa9714788ba1105a19df42a`,
and DisProt candidate catalog
`2bc6726545f47256bcec70376bd51d689db0b5c80a64f45a221b9d88bb0c4c8d`.
CAID3 reference SHA-256
`4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`;
DP94 sequence-only manifest SHA-256
`fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`.
The run used the local parasail 1.3.4 implementation. Three synthetic unit
tests verify exact/containment hard hits and that a short shared segment does
not alone constitute full-length homology.

```powershell
C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe scripts\audit_dfl_full_length_screen.py --source-root C:\Users\Perthro\Documents\ChatGPT\Folddiffusion
C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe -m unittest tests.test_audit_dfl_full_length_screen -v
```

## Split protocol to complete B1

1. Freeze the source releases and both benchmark **sequence-only** registries.
   Keep exact accession, sequence and containment collisions excluded.
2. Obtain a reproducible family/domain annotation layer: pinned MMseqs2
   full-length clusters with low-complexity handling, plus Pfam/InterPro domain
   IDs or profiles. Verify sequence-to-accession matches and source rights.
   Same-domain or ambiguous remote-homology pairs must be reviewed or excluded,
   not waved through because their full proteins differ.
3. Use the shuffle-null and known homolog controls to preregister the
   sequence-only edge rule. Save all old-versus-new edge changes and review
   decisions. A contract change requires explicit scientific review; until
   then the old rule remains binding.
4. Build connected components over source and benchmark sequences with exact,
   validated full-length and domain-family edges. Remove benchmark-touching
   components. Allocate *whole remaining components* to train/development
   (and auxiliary only for its actual source task). Keep source-specific
   `1/0/x` masks; no inferred domain gap or unknown-as-negative relabeling.
5. Create a versioned split manifest, repeat both B1 preflights with pinned
   hashes, and wire the same checks into the launcher before GPU allocation.
   If insufficient positive clusters or no defensible negative semantics
   remain, B1 stays blocked and more curated annotations are required.

An independently adjudicated no-linker set is also required before a
high-precision mining claim. The prospective DP94 DFL track still needs
fitness, coverage and competitor-lineage validation. No model has been
trained or evaluated from this report.
