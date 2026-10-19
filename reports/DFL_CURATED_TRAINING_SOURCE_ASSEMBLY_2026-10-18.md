# DFL curated training-source assembly — 18 October 2026

Decision: **a versioned, provenance-bearing candidate pool is assembled, but no defensible B1 training split exists yet.** No GPU training, held-out-label reading, or benchmark scoring occurred. The candidate records are not a licence to train. The dataset is local and gitignored; the builder, schema, tests, source hashes and counts are tracked here.

## What can be borrowed from the papers

| Published source | Relevant dataset/preparation | What we use; what we do not infer |
| --- | --- | --- |
| [DFLpred](https://pmc.ncbi.nlm.nih.gov/articles/PMC4908364/) | DisProt 6.0.2 functional annotations; 82 DFL and 122 fully annotated non-DFL proteins, cluster-aware partitioning. | Use its explicit DFL-versus-other/unknown distinction as a preparation precedent. Its non-DFL convention is not independent proof that every residue lacking a DFL annotation is negative. |
| [APOD](https://pmc.ncbi.nlm.nih.gov/articles/PMC7773485/) | DisProt 8.0, 248 proteins, including 117 DFL and 131 other functionally annotated proteins; TR166/TE82 grouped by BLASTClust. | Use family-unit splitting and exclusion of unannotated residues. The published data overlap later DisProt/benchmarks and cannot be copied into an ostensibly independent split. |
| [TransDFL](https://pmc.ncbi.nlm.nih.gov/articles/PMC10626177/) | General-disorder pretraining on thousands of IDR examples, followed by DFL-specific fine-tuning on TR166; TE82/TE64 tests. | Keep disorder auxiliary data separate. The thousands are not thousands of verified DFL masks. |
| [LINKER-Pred2 release](https://github.com/deemeng/punch_linker2) and [DLD download](https://pcrgwd.ucd.ie/linker/) | More than 2,000 mixed DLD/DisProt linker examples are described; DLD provides structural IDL/DDL annotation. | Structural DLD records remain an auxiliary task, not DFL positive or negative ground truth. |
| [DisProt 2026_06](https://disprot.org/) | Current local JSON includes curated region IDs, IDPO terms, reference IDs, validation, release and source accession. DisProt states CC BY 4.0. | Reconstruct explicit full-protein `1/0/x` evidence masks with attribution, preserving the source flags and retaining unknown residues. |

The previously assembled APOD/DFLpred-style local pool (`data/processed/linker_pred2_style_pool_v1`) has 158 TR166-derived records, of which 74 carry observed DFL positives; see [B1 feasibility audit](DFL_B1_FEASIBILITY_AUDIT_2026-10-18.md). It is not a ready-made independent train set. The 24,763 annotated-linker catalog is likewise not full-protein DFL segmentation ground truth. The separate UniProt explicit-region candidate pool has 884 matched full sequences but only 11 direct benchmark non-hits, one in its direct PubMed-evidence tier. Neither was silently promoted here.

## Assembled DisProt evidence pool

Run `scripts/build_dfl_curated_evidence_pool.py` with the pinned DisProt 2026_06 raw JSON, CAID3 **sequences** and DP94 sequence-only registry. The official CAID3 annotation line is discarded without parsing. The source JSON SHA-256 is `1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140`; CAID3 FASTA is `4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`; DP94 manifest is `fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`. The locally generated `data/processed/dfl_curated_evidence_pool_v6/manifest.json` has SHA-256 `9cdfd62a63d99cf20a434e9b46fafe9582cff5a20131343f149057aa75d3019c`.

The builder accepts only validated, cited, released region records. `IDPO:0000033` flexible-linker spans are `1`; `IDPO:0000006` experimental-order spans are `0`; all unannotated and other-functional residues are `x`; conflicting 1/0 positions revert to `x`. An explicit ordered *region* is negative evidence at those residues, not certification of a whole no-linker protein. Coordinates, region/reference/ECO IDs, flags, source ID, input/config/schema hashes and producing commit are retained per record. `schemas/dfl_curated_evidence_candidate.v1.schema.json` is enforced during generation. Output paths are immutable: rerunning refuses to overwrite the directory.

| Tier | Before sequence/family screening | Surviving CAID3+DP94 screen | Interpretation |
| --- | ---: | ---: | --- |
| Strict: validated, cited, released, `unpublished=false` | 57 proteins; 1,718 positive and 290 order-negative residues | **2 proteins; 36 positive, 0 order-negative, 1,065 unknown residues** | Source-qualified *candidates*, not a B1 train/development split or specificity set. |
| Separate flag-review tier: otherwise qualifying regions with `unpublished=true` | 186 changed-mask proteins; 6,654 positive and 1,639 order-negative residues | 15 proteins; 443 positive, 83 order-negative, 6,874 unknown residues | Not trainable until DisProt flag semantics and source evidence are reviewed; may overlap strict proteins. |

The strict two survivors are DP00034 (22 positive residues) and DP02785 (14), each with zero explicit ordered residues. The broad local screen quarantined 50 strict candidates by ID/direct sequence relation and a further 5 by transitive shared Pfam/Gene3D/UniRef50 component. It is deliberately conservative and can overconnect, but it cannot be weakened merely to increase sample size. The family graph is limited to local DisProt metadata and is **not** a certified remote-homology screen.

`unpublished=true` is not interpreted as proof that the cited paper itself is unpublished: public DisProt pages can display a region carrying a PubMed reference and release date. Its precise field meaning has not been verified, so the review tier remains separate rather than being silently discarded or admitted. The counts above are **not additive**: a review row can be a more extensively masked/annotated version of a strict row.

Reproduce with the pinned ignored source files:

```powershell
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' -m unittest tests.test_build_dfl_curated_evidence_pool -v
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' scripts/build_dfl_curated_evidence_pool.py --source-root 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion' --output-dir 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\data\processed\dfl_curated_evidence_pool_v7'
```

The `v6` build and four synthetic tests passed on 18 October. A rebuilt directory name or repository commit will change per-record provenance and output hashes; verify its manifest rather than expecting byte identity across commits.

## B1 gate and high-leverage next action

The present protected benchmark contract excludes CAID3 and DP94 families and cannot be altered on the basis of available sample counts. With two strict positives and no known negatives, it is impossible to form a useful family-disjoint train/development split or calibrate high-precision mining. **Do not launch a three-layer CNN on this pool or claim B1 complete.**

To make a defensible dataset rather than a larger relabelled one, the next source effort must (1) resolve the `unpublished` flag with DisProt documentation or curator confirmation; (2) seek *new* independently curated DFL positives outside both protected benchmark families, recording reference and sequence lineage; (3) create a separately adjudicated family-disjoint no-linker/ordered-region panel with explicit unknown masks; and (4) freeze a validated family-domain graph and whole-component train/development assignment, then pass both preflights. General-disorder IDRs and DLD structural linkers can be auxiliary inputs only. If a prospective benchmark is later changed, that must be a documented new evaluation contract, not a covert relaxation of this firewall.
