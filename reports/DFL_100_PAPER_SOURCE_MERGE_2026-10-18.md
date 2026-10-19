# Hundred-paper DFL source screen and conservative merge — 18 October 2026

Status: **100 title/abstract records screened; four pinned local DFL evidence streams deduplicated into a 24-protein positive/unknown candidate pool; B1 is still blocked.** This is not a 100-dataset download, full-text review of every paper, an approved training split, or a precision estimate. No benchmark labels or scores were opened, and no GPU job was launched.

## Search and triage

`scripts/survey_dfl_literature_sources.py` queried the [Europe PMC article API](https://www.ebi.ac.uk/europepmc/webservices/rest/search) with twelve declared title/abstract queries on 18 October 2026. It retrieved metadata/abstracts for 549 distinct title-normalized candidates, ranked them with disclosed task terms, excluded clear non-protein or cross-linker false positives, and retained 100 for the [tracked screening table](DFL_100_PAPER_ABSTRACT_SCREEN_2026-10-18.tsv). The table includes DOI/PMID/PMCID where available, title, year, tentative lane and an explicit manual-verification warning. The complete abstract JSONL and per-query response hashes remain in the local ignored `artifacts/dfl-literature-screen-2026-10-18-v3/`; its manifest SHA-256 is `92eec5275744068bff11a075d598063202d2891bfc4ef2054ccc43814c2f0ae1` and table SHA-256 is `d74f8411605a50b075de6b48f8878578541e5dd5e55f22f8eac4022ca266dda2`.

Triage counts are **5 directly DFL-related papers**, 59 protein-linker-related papers without verified DFL row labels, and 36 general-disorder/auxiliary papers. Those classifications are *abstract-level leads*, not an assertion that all 100 contain a downloadable DFL dataset. Title normalization collapses a paper/preprint pair with the same title; a deliberate false-positive filter removes unrelated apoD neurology, linker histones and chemical cross-linkers. A manual spot check of the top and bottom titles confirmed that many remaining entries are single-protein biology or disorder-method papers, useful for candidate curation or background, **not** training data. The search is a reproducible prioritized sample, not a claim that these are the globally most cited or exhaustive 100 papers.

## Dataset lineage decisions

| Direct source or useful lead | Published evidence | Disposition |
| --- | --- | --- |
| [DFLpred](https://pmc.ncbi.nlm.nih.gov/articles/PMC4908364/) and [APOD](https://pmc.ncbi.nlm.nih.gov/articles/PMC7773485/) | DisProt-derived DFL regions; DFLpred describes 82 DFL plus 122 fully annotated other proteins; APOD describes 117 DFL plus 131 other-function proteins and TR166/TE82 clustering. Both exclude unannotated residues in evaluation. | Only locally pinned training-side source rows with observed `1` are candidate positives. Historical paper/test partitions are not pulled across the protected benchmark firewall; their `0` labels remain unverified for our precision objective. |
| [TransDFL](https://pmc.ncbi.nlm.nih.gov/articles/PMC10626177/) | General-disorder pretraining on 2,645 IDR train and 1,077 validation proteins, followed by DFL fine-tuning on TR166; TE82/TE64 are tests. | The thousands are **not** DFL labels. Disorder data stay in a separately labelled auxiliary tier and are not merged into the DFL masks. |
| [LINKER-Pred2 documentation](https://pcrgwd.ucd.ie/linker_pred/help/) and [DLD portal](https://pcrgwd.ucd.ie/linker/) | CNNs trained on DisProt plus a structural DLD collection; DLD IDL/DDL classes derive from structural domain annotation. | No DLD structural `1` is promoted to disordered-flexible-linker gold. Preserve as a separate structural task until exact common-target semantics and rights are established. |
| [DisoFLAG author download](https://bliulab.net/DisoFLAG/download) / [author repository](https://github.com/YihePang/DisoFLAG) | Multi-function labels include a DFL track; author TR552 data have been separately pinned and capacity-audited. | Not merged here: redistribution licence, source zero semantics and benchmark lineage need verification. DP94 remains a protected, label-sealed candidate and is never mined as training data. |
| [DisProt 2026_06](https://disprot.org/) | Curated IDPO flexible-linker and experimental-order spans with references and accession lineage; CC BY 4.0 requires attribution. | Strict validated/cited/released, nonflagged source gives two benchmark-screened proteins. Both already occur among the 24 merged sequences; source observations are retained. Flagged annotations stay in a separate review tier. |

Other papers in the 100-row screen may supply individual biological examples, but a paper title or abstract mentioning a flexible linker is **not** a residue-labelled dataset. No inferred domain gap, AlphaFold low-confidence segment, generic disorder label, or engineered GS linker was turned into a DFL positive. No third-party dataset was rehosted in Git.

## Exact-sequence merge and rejection rules

`scripts/merge_dfl_published_positive_candidates.py` takes four hash-pinned local streams: published TR166-style train (`a1a5f3b3b30548ad65f6f378535c7c732272b8c9b60048b38f6500fca9bb1aed`), DisProt positive-only (`0265b062ca9ff88cde2576789b613860355f33842fa9714788ba1105a19df42a`), DisProt candidate catalog (`2bc6726545f47256bcec70376bd51d689db0b5c80a64f45a221b9d88bb0c4c8d`), and the separately built strict DisProt region tier (`07e054bbf8d3be9d17438cd843d0610bda84cdd7baa0d1e5b6083e4ae2b98845`). It requires the pre-existing **24-sequence** binding-rule-plus-family capacity allowlist (SHA-256 `bca74cae3932d8b086d5d292d18adb41f0c9ce39e9744cc4b05d63fea9912368`) and independently rechecks those records against CAID3 and DP94 **IDs/sequences**. It refuses changed hashes, source-task substitution, malformed masks, missing allowlist members and output overwrites. The earlier report's published-train hash contained a one-character transcription error; the hash above matches the file and its source manifest, and the two earlier reports were corrected.

The local immutable candidate output is `data/processed/dfl_published_positive_merged_v3/manifest.json`, SHA-256 `9e0def03d4d3d95f9cb4ddea6841b9e8973070868e5e052cc594f278069c1299`. Every JSONL row validates against `schemas/dfl_merged_positive_candidate.v1.schema.json` and retains source IDs, input hashes, family component, producing commit, config/schema hashes and discarded-zero counts.

| Merge result | Count |
| --- | ---: |
| Pinned source observations that survived the prior screens | 43: 5 published train, 17 DisProt positive-only, 19 DisProt catalog, 2 strict-curated |
| Exact-sequence duplicates removed | 19 |
| Unique positive/unknown full-protein candidates | **24** |
| Union of observed DFL-positive residues | **1,123** |
| Unknown residues | **12,721** |
| Source `0` observations discarded as unverified negatives | **200** |
| Cross-source `1`/`0` disagreements masked to `x` | **0** |
| Independently adjudicated no-linker proteins | **0** |

The source `0` observations are counted, never accepted as whole-protein no-linker ground truth. Positives are unioned only when the exact full sequence is identical; any future source `1`/`0` disagreement is masked to `x`. Apparent duplicates with altered sequence are not forcibly merged by accession. Structural DLD and disorder-only data are not in this output. These are conservative **candidate** records, not a frozen family-disjoint train/development allocation or precision-calibration panel.

## Reproduction and remaining gate

The original raw datasets and benchmark sequence registries are ignored locally, so a clean clone needs licensed/downloaded source inputs with the pinned hashes. From the audit worktree, pointing to the source checkout:

```powershell
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' -m unittest tests.test_merge_dfl_published_positive_candidates tests.test_build_dfl_curated_evidence_pool tests.test_preflight_dfl_family_training -q
& 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\.venv-gpu\Scripts\python.exe' scripts/merge_dfl_published_positive_candidates.py --source-root 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion' --capacity-manifest 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion-pr-merge-review\artifacts\b1-mmseqs-input-2026-10-18\binding-family-survivors-v2.json' --output-dir 'C:\Users\Perthro\Documents\ChatGPT\Folddiffusion\data\processed\dfl_published_positive_merged_rebuild'
```

The focused 21-test suite passed. The evidence at hand **does not close B1**: 24 positive families and zero independently adjudicated no-linker proteins cannot calibrate proteome-scale high precision, and the exploratory family graph is not a certified remote-homology split. The highest-leverage next data acquisition is curator-verified positive/negative DFL evidence **outside both protected registries' families**, with an explicit no-linker panel and source rights. A new versioned split and both preflights must pass before any new CNN training. Neither CAID3's spent labels nor DP94's sealed labels may be used to enlarge the train set.
