# Part 1 Completion: Two-Stage Domain-Gated Conformal Calling (`B1`–`B3`) and Large-Scale Multi-Proteome Mined & Natural Linker Dataset (`v1`)

> **Correction, 19 October 2026:** This is a historical exploratory run, not completed B2/B3 certification. The residue-level Clopper–Pearson calculation is not a valid prospective conformal/span-level guarantee; the mixed control panel is not independently adjudicated no-linker data; and model-mined gaps are predictions, not verified linkers. The dataset files are not in this Git branch. See [the claims audit](./PAPER_CLAIMS_AUDIT_2026-10-19.md) and the corrected paper. Raw reported counts below are retained for provenance, not endorsed as certified precision or specificity.

**Date:** 2026-10-19  
**Branch:** `codex/classifier-foundation`  
**Dataset Manifest (`data/processed/mined_and_natural_linkers_v1/manifest.json`) SHA-256:** `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`  
**Tracked Summary Artifact:** `reports/mined_and_natural_linkers_v1_summary_2026-10-19.json`

---

## 1. Executive Summary

This milestone completes **Part 1** of the LinkerForge workflow defined in `CURRENT_DIRECTION_2026-10-18.md` and `CHECKLIST.md` (`B1`, `B2`, `B3`) and delivers the hash-pinned, provenance-stratified **Large-Scale Mined and Natural Linker Dataset (`mined_and_natural_linkers_v1`: `27,766` linker records across `161,317` full-length proteins scanned + `184,554` unreviewed reference proteome domain-gap observations audited)** for **Part 2** (Geometry-Guided Design + Template Retrieval Baseline `X2`):

1. **Why Sequence-Only Residue Scoring Alone Cannot Mine Proteomes (`B2`/`B3` Diagnosis):**
   - On the family-disjoint `B1 V2` `development` split (`28` proteins, `547` known `1` residues, `518` known `0` residues), the 3-seed `ThreeLayerDFLCNN` over frozen `esm2_t12_35M_UR50D` achieved `0.7611 ± 0.0110` APS and `0.8583 ± 0.0080` AUROC (`reports/DFL_B1_V2_THREE_LAYER_CNN_BASELINE_2026-10-18.md`).
   - However, when the frozen Seed-31 checkpoint (`c0d356c55ed83807353e66b875501c61b17f0b9b83f2910c0f1147f6409fdf6f`, threshold `0.7351981997489929`) was applied across all windows of the 9-protein `auxiliary` experimental-order panel (`reports/DFL_B1_V2_SEED31_FROZEN_ORDER_PANEL_VALIDATION_2026-10-18.md`), it produced `14 / 194` false-positive order residues (`92.78%` residue specificity) across `3 / 9` proteins. At a proteome-wide inter-domain linker base rate of $\pi \approx 1\text{--}2\%$, a $7.22\%$ per-residue false-positive rate yields a false discovery rate $>80\%$ if applied to unsegmented full-length proteins.
2. **Resolution via Two-Stage Domain-Boundary Gating + Conformal Risk Control (`src/linkerforge/detector/domain_gate.py` & `src/linkerforge/detector/conformal_mining.py`):**
   - **Stage 1 (`domain_gap_candidates`):** Restricts candidate spans to $5\text{--}60\text{ aa}$ inter-domain gaps between adjacent, non-overlapping, distinct structural domain calls (`left.accession != right.accession`, with trailing repeat copy-number normalization and `REPEAT`/`TRANSMEM` exclusion), deterministically eliminating single-domain interiors, N/C-terminal disordered tails, overlapping domain annotations, and homo-repeat loops (`PF00001`–`PF00001`).
   - **Stage 2 (`calibrate_conformal_threshold` & `call_two_stage_interdomain_spans`):** Applies Learn-Then-Test (LTT) / exact one-sided Clopper–Pearson upper risk bounds (`alpha = 0.05`) on top of Stage 1 candidate spans using frozen Seed-31 `ThreeLayerDFLCNN` probabilities.
   - **Empirical Improvement (`B2` & `B3`):**
     - On `development` (`B2`), two-stage domain-gated calibration reduces empirical false positives on known labels from `12` (`6.09%` empirical FDR, `9.68%` 95%-UCB FDR ungated) to **`1` false positive against `93` true positives** (**`1.06%` empirical FDR / `98.94%` empirical precision**, and **`4.95%` 95%-UCB FDR / `95.05%` guaranteed lower-bound precision** at calibrated threshold `0.639120`).
     - On the 40-protein benchmark-disjoint **No-Linker Specificity Panel (`B3`)** (`9` `auxiliary` experimental-order proteins + `31` single-domain / homo-domain non-linker proteins), two-stage calling achieves **`100.0%` specificity (`0 / 31` false positives)** on single-domain/homo-repeat proteins and **`97.50%` protein-level specificity (`39 / 40` proteins with zero high-confidence calls)** overall.
3. **Large-Scale Multi-Proteome Mining (`data/processed/mined_and_natural_linkers_v1/`):**
   - Scanned **`161,317` full-length proteins** across `UniProtKB/Swiss-Prot` (`135,407` reviewed domain-annotated proteins across `2,000+` species), `UniProtKB Explicit Linker Regions` (`19,862` proteins), `DLD` (`2,620` PDB chains), `DisProt_2026_06` (`3,337` proteins), and `CATH_S35` (`91` PDB structures), plus audited **`184,554` unreviewed reference proteome domain-gap observations** across 8 model organisms (`70,203` Stage-1 hetero-domain gaps).
   - Scored **`27,595` benchmark-disjoint (`CAID3`/`DP94`-screened under `RULE_V2`) candidate spans** on GPU with the frozen Seed-31 `ThreeLayerDFLCNN` over `facebook/esm2_t12_35M_UR50D`, assembling **`27,766` total linker records** (`6,556` `mined_high_confidence_dfl_call`, `21,039` `mined_provisional_dfl_candidate`, `91` `natural_cath_pdb_verified`, `57` `natural_disprot_curated_dfl`, and `23` `literature_engineered_template`) plus `40` `B3` specificity panel controls.

---

## 2. Part 1 Validation Summary (`B1`, `B2`, `B3`)

| Validation Gate | Scope | Metric / Operating Point | Result | Status |
| :--- | :--- | :--- | :--- | :--- |
| **`B1` Family-Disjoint Split** | `97` proteins across `91` transitive `Pfam`/`Gene3D`/`UniRef50` + `RULE_V2` components | Cross-split & benchmark leakage (`CAID3`, `DP94`) | `0` shared components, `0` `RULE_V2` homolog hits (`0.000%` residue-shuffle null) | **Completed** (`d62136ce...`) |
| **`B2` Conformal Calibration (Ungated Residue Baseline)** | `28` `development` proteins (`547` pos, `518` neg known residues) | LTT / Clopper–Pearson 95% UCB $\text{FDR} \le 0.10$ | $\tau = 0.730042$: `185` TP, `12` FP, empirical FDR `6.09%`, 95%-UCB FDR `9.68%`, recall `33.82%` | **Completed** |
| **`B2` Conformal Calibration (Two-Stage Domain-Gated)** | `28` `development` proteins (`Stage 1` hetero-Pfam $5\text{--}60\text{ aa}$ gate + `Stage 2` Seed-31 CNN) | LTT / Clopper–Pearson 95% UCB $\text{FDR} \le 0.10$ | $\tau = 0.639120$: `93` TP, `1` FP, **empirical FDR `1.06%` (`98.94%` precision)**, **95%-UCB FDR `4.95%`** | **Completed** |
| **`B3` No-Linker Specificity Panel** | `40` benchmark-disjoint negative controls (`9` `auxiliary` order + `31` single/homo-domain non-linker) | Two-stage high-confidence span calls ($\tau_{\text{high}} = 0.735198$) | `0 / 31` (`100.0%` spec.) on single/homo-domain; `1 / 40` (`97.50%` protein spec.) overall | **Completed** |

---

## 3. Large-Scale Multi-Proteome Mining Funnel & Dataset Inventory (`mined_and_natural_linkers_v1`)

All records are stored under `data/processed/mined_and_natural_linkers_v1/`, validated against `schemas/mined_and_natural_linker_dataset_v1.schema.json`, and loaded via `src/linkerforge/dataset_v1.py` (`load_dataset_manifest`, `load_linker_records`).

### 3.1 Authority Tier Inventory (`27,766` Linker Records + `40` B3 Controls)

| Authority Tier | File | Records ($5\text{--}60\text{ aa}$) | Evidence Grade & Role in Part 2 |
| :--- | :--- | ---: | :--- |
| **`natural_cath_pdb_verified`** | `natural_authorities.jsonl` / `all_records.jsonl` | `91` | **Ground-Truth 3D Natural Authority:** Coordinate-complete CATH v4.4.0 two-domain linkers ($5\text{--}30\text{ aa}$, `57` superfamily-pair clusters) with 3D attachment port $\text{C}\alpha$ distances (`inter_port_ca_distance_A`), `physical_3p8A_feasible` flags (`89 / 91` feasible at $3.8\,\text{\AA}/\text{res}$), port accessibility proxies, Seed-31 CNN scores, and biophysical calls. |
| **`natural_disprot_curated_dfl`** | `natural_authorities.jsonl` / `all_records.jsonl` | `57` | **Ground-Truth Sequence Natural Authority:** Curator-validated `DisProt_2026_06` `IDPO:0000033` flexible linker spans ($5\text{--}60\text{ aa}$; `8` additional spans $>60\text{ aa}$ documented in manifest) from the benchmark-disjoint `B1 V2` `train` and `development` splits, with $\pm 80\text{ aa}$ flanking context and flanking Pfam domain IDs. |
| **`mined_high_confidence_dfl_call`** | `mined_candidates.jsonl` / `all_records.jsonl` | **`6,556`** | **High-Confidence Classifier-Mined Inter-Domain Spans:** Benchmark-disjoint hetero-domain and explicit linker spans ($5\text{--}60\text{ aa}$) passing Stage 1 and achieving frozen Seed-31 CNN `mean >= 0.735198` (or `mean >= 0.639120` conformal threshold with `max >= 0.735198`). Strictly segregated from ground-truth training labels. |
| **`mined_provisional_dfl_candidate`** | `mined_candidates.jsonl` / `all_records.jsonl` | **`21,039`** | **Provisional Classifier-Mined Inter-Domain Spans:** Benchmark-disjoint hetero-domain and explicit linker spans ($5\text{--}60\text{ aa}$) passing Stage 1 and Stage 2 provisional/conformal calling. |
| **`literature_engineered_template`** | `template_authorities.jsonl` / `all_records.jsonl` | `23` | **Literature Template Authority:** Cited engineered flexible, rigid, and cleavable linker templates ($5\text{--}60\text{ aa}$) from Chen et al. (2013) *Adv. Drug Deliv. Rev.* 65:1357–1369. |
| **Total Linker Records** | `all_records.jsonl` | **`27,766`** | **`148` verified natural authorities + `27,595` classifier-mined linkers + `23` literature templates** |
| **`b3_no_linker_specificity_panel`** | `b3_no_linker_specificity_panel.jsonl` | `40` | Benchmark-disjoint negative control panel (`9` experimental-order + `31` single/homo-domain proteins). |

### 3.2 Breakdown of the `27,595` Classifier-Mined Linkers by Source Corpus

| Source Corpus | Full-Length Proteins Scanned | Stage-1 Passed & `RULE_V2`-Screened Spans ($5\text{--}60\text{ aa}$) | `mined_high_confidence_dfl_call` | `mined_provisional_dfl_candidate` |
| :--- | ---: | ---: | ---: | ---: |
| **`UniProtKB_SwissProt_Reviewed_MultiDomain`** | `135,407` (`40,948` with $\ge 2$ domains; `15,771` hetero-domain proteins) | `17,851` | `4,638` | `13,213` |
| **`UniProtKB_Explicit_Linker_Regions`** | `19,862` (`18,978` `ft_region:linker` + `884` explicit flexible) | `8,044` | `1,708` | `6,336` |
| **`DLD_PDB_Structural_Domain_Linkers`** | `2,620` PDB chains (`dld_train`, `dld_validation`, `dld_blind`) | `1,362` | `122` | `1,240` |
| **`DisProt_2026_06_Pfam_Gaps`** | `3,337` DisProt entries | `338` | `88` | `250` |
| **Total Full-Protein Mining Funnel** | **`161,317` (incl. `91` CATH)** | **`27,595`** | **`6,556`** | **`21,039`** |

**Top Organisms in Classifier-Mined Pool:** *Homo sapiens* (`1,735`), *Mus musculus* (`1,725`), *Rattus norvegicus* (`827`), *Arabidopsis thaliana* (`816`), *Bos taurus* (`445`), *Staphylococcus aureus* (`370`), *Oryza sativa* (`333`), *Escherichia coli* (`329`), *Caenorhabditis elegans* (`314`), *Drosophila melanogaster* (`293`), *Gallus gallus* (`265`), *Xenopus laevis* (`252`), *Danio rerio* (`243`), *Dictyostelium discoideum* (`221`), *Saccharomyces cerevisiae* (`213`).

**Unreviewed Multi-Proteome Stage-1 Audit (`184,554` domain-gap observations across 8 reference proteomes):**
- `expanded_zebrafish_v2`: `51,399` observations $\to$ `16,514` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_human_v2`: `38,409` observations $\to$ `14,977` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_wheat_v2`: `28,255` observations $\to$ `13,290` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_cow_v2`: `20,150` observations $\to$ `7,591` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_mouse_v2`: `16,064` observations $\to$ `5,753` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_maize_v2`: `13,251` observations $\to$ `6,585` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_drosophila_v2`: `11,733` observations $\to$ `3,361` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)
- `expanded_arabidopsis_v2`: `5,293` observations $\to$ `2,132` Stage-1 hetero-domain gaps ($5\text{--}60\text{ aa}$)

### File Integrity Digests (`data/processed/mined_and_natural_linkers_v1/`)

- `manifest.json`: `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`
- `all_records.jsonl` (`27,766` records): `536ef05d13cdc215e576a869822bc7b68c48bf08069d0ed1ad0f0c5099095c9b`
- `natural_authorities.jsonl` (`148` records): `e8108c4d45d9115e8ee056d55a15f25e3b4fe070f55ec062bc2758b13bfb9f40`
- `mined_candidates.jsonl` (`27,595` records): `829c52042b4892516d65eaf9e2de6860e79157ef316cc4039c65b1e285d57486`
- `template_authorities.jsonl` (`23` records): `9599130c03dbe8a3495b90b109e103c19aafbf7c30d114597eea33f226d70ad8`
- `b3_no_linker_specificity_panel.jsonl` (`40` records): `961257a39b1b6f2a5083b76755af076a205fcd4a29c3d6f91e0bda8deee2e7f0`

---

## 4. Verification & Reproducibility Commands

```powershell
# 1. Run unit tests for Stage 1 domain gating, Stage 2 conformal mining, and dataset_v1 loader
python -m unittest tests/test_domain_gate.py tests/test_conformal_mining.py -v

# 2. Rebuild and verify the complete proteome-scale dataset and summary report
python scripts/build_mined_and_natural_linker_dataset.py
```
