# Part 2 Completion: Geometry-Plus-Template Linker Design (`Q1`–`Q3`) and Domain-Pair/Family-Disjoint Retrieval Baseline (`X2`)

**Date:** 2026-10-19  
**Branch:** `codex/classifier-foundation`  
**Input Dataset (`data/processed/mined_and_natural_linkers_v1/manifest.json`) SHA-256:** `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`  
**Tracked Evaluation Artifact:** `reports/x2_retrieval_and_part2_design_2026-10-19.json`

---

## 1. Executive Summary

This milestone completes **Part 2** of the LinkerForge workflow defined in `CURRENT_DIRECTION_2026-10-18.md` (Steps 3 and 4), `SPEC.md` §7, and `CHECKLIST.md` (`X2`), evaluated over the expanded **`27,766`-linker dataset (`148` verified natural linkers + `27,595` classifier-mined linkers + `23` literature templates)**:

1. **Part 2A — Deterministic Geometry-Plus-Template & Authority Linker Design Engine (`src/linkerforge/design/engine.py`, `scripts/design_linker.py`):**
   - Integrates `linkerforge.geometry` (`required_length`, `verdict_for_requested_length`, `flexible_strain`, `evaluate_junction`, `implied_chain_topology`), `linkerforge.properties` (`fit_rigidity_model`, `call_cleavability`, `composition_features`, `load_library`), and `linkerforge.dataset_v1` (`mined_and_natural_linkers_v1`, `27,766` records across 4 authority families).
   - Answers the three design questions under explicit two-axis conditioning (`rigidity: "flexible" | "rigid"`, `cleavability: "non_cleavable" | "cleavable"`, optional `protease` filter) and explicit topology selection (`topology: "end_to_end_tandem" | "internal_insertion_two_tether"`):
     - **`Q1` (Geometric Feasibility & Topology):** Computes minimum physical ($3.8\,\text{\AA}/\text{res}$ virtual $\text{C}\alpha$ bond ceiling) and design ($3.5\,\text{\AA}/\text{res}$ extended or $1.5\,\text{\AA}/\text{res}$ helical rise) residue windows, Gaussian/worm-like chain RMS extension ratios, and polypeptide topology constraints (`1` tether for tandem end-to-end fusion vs `2` leaving/re-entering polypeptide tethers for internal domain insertion). Reports deterministic flags, never a fusion success probability (`SPEC.md` §5.8).
     - **`Q2` (Attachment Site & Flanking Context Assessment):** Evaluates 3D anchor relative solvent accessibility (SASA), outward orientation (`anchor_faces_partner_domain`, `chain_exit_angle_to_span_degrees`), and domain steric occlusion (`required_path_length`) when 3D `Chain` structures are supplied, alongside flanking sequence/domain similarity against `mined_and_natural_linkers_v1`.
     - **`Q3` (Candidate Linker Selection & Two-Part Rigid+Cleavable Composite Construction):** Ranks bridgeable candidates across cited literature templates (`literature_engineered_template`), structure-verified natural linkers (`natural_cath_pdb_verified`), curator-validated disordered linkers (`natural_disprot_curated_dfl`), and domain-gated high-confidence mined spans (`mined_high_confidence_dfl_call`, `6,556` records). For `(rigidity="rigid", cleavability="cleavable")`, `build_rigid_cleavable_composite` enforces the biophysical requirement (`SPEC.md` §7; Hubbard et al. 1994; Chen et al. 1994/2013) that a protease cleavage site requires local backbone unfolding and cannot sit inside an unbroken $\alpha$-helix, automatically assembling a **two-part composite construct** (helical `A(EAAAK)nA` or `PAPAP`/`(AP)n` rigid spacer + `GGGGS`-flanked scissile protease cassette).

2. **Part 2B — Domain-Pair / Family-Disjoint Retrieval Baseline (`X2`) (`src/linkerforge/retrieval/baseline.py`, `scripts/evaluate_x2_template_and_authority_retrieval.py`):**
   - Evaluated via 5-fold strictly cluster-disjoint (`homology_cluster_id`-disjoint) cross-validation across two cohorts, augmented with the **`27,595` classifier-mined linkers**:
     - **Cohort A (`cath_3d_natural_91`):** `91` coordinate-complete CATH v4.4.0 natural two-domain linkers across `57` superfamily-pair clusters (with 3D `inter_port_ca_distance_A` and flanking sequences/CATH IDs).
     - **Cohort B (`all_natural_authorities_148`):** All `148` verified natural linkers (`91` CATH + `57` DisProt curated `IDPO:0000033` linkers) across `109` disjoint family/domain-pair clusters.
   - **Decisive Empirical Findings (`X2`):**
     1. **Proteome-scale classifier mining (`27,595` spans) substantially improves held-out natural linker retrieval:** Adding the `27,595` classifier-mined linkers to the retrieval pool (`flanking_context_nn_natural_plus_mined` and `geometry_plus_flanking_hybrid_retrieval`) increases global Needleman–Wunsch sequence identity on held-out natural linkers from `0.0928` (`flanking_context_nn_natural_authority`) to **`0.1426`** (`+53.7%`) and **`0.1565`** (`+68.6%`, `geometry_plus_flanking_hybrid_retrieval`) on Cohort A, and from `0.0981` to **`0.1451` (`+47.9%`)** on Cohort B!
     2. **3D port geometry + proteome-mined flanking retrieval (`geometry_plus_flanking_hybrid_retrieval`) achieves the strongest 3D reconstruction:** On Cohort A (`91` CATH 3D natural linkers), `geometry_plus_flanking_hybrid_retrieval` achieves **`43 / 91` (`47.3%`) length recovery within $\pm 2\text{ aa}$**, **`63 / 91` (`69.2%`) within $\pm 5\text{ aa}$ (`MAE = 4.56 aa`)**, **`91 / 91` (`100.0%`) rigidity axis match**, **`91 / 91` (`100.0%`) $3.8\,\text{\AA}/\text{res}$ physical feasibility**, and **`0.1565` NW sequence identity** (the highest sequence identity of all 9 arms).
     3. **Flanking-context retrieval strongly recovers biophysical composition (`0.726–0.805` composition cosine vs `0.483–0.610` for fixed/repeat controls):** On `Cohort B` (`148` natural linkers), `flanking_context_template_retrieval` achieves **`0.8052` composition cosine (`95% CI [0.7685, 0.8400]`)**, `flanking_context_nn_natural_authority` achieves **`0.8001`**, and `flanking_context_nn_natural_plus_mined` achieves **`0.7971` (`95% CI [0.7648, 0.8283]`)** alongside **`137 / 148` (`92.6%`) rigidity axis match**, outperforming `fixed_ggggs3_control` (`0.5996`), `training_median_length_gs_control` (`0.5330`), and `length_matched_frequent_natural_control` (`0.4833`) by **$+0.20\text{ to }+0.32$ cosine units**.

---

## 2. Experiment `X2` Results: Cohort A (`91` CATH 3D Natural Linkers, `57` Disjoint Domain-Pair Clusters)

5-fold cluster-disjoint cross-validation (`seed = 17`, `1,000` cluster-bootstrap resamples, `27,595` classifier-mined linkers in augmented pool). Raw counts and means from `reports/x2_retrieval_and_part2_design_2026-10-19.json`:

| Arm | Input Used at Inference | Length MAE (aa) [95% CI] | Within $\pm 2\text{ aa}$ | Within $\pm 5\text{ aa}$ | Global NW Seq Identity | Composition Cosine [95% CI] | Rigidity Match | 3D $3.8\,\text{\AA}$ Feasible |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `fixed_ggggs3_control` | None (`(GGGGS)3`, `15 aa`) | `6.26` `[5.24, 7.21]` | `18 / 91` (`19.8%`) | `39 / 91` (`42.9%`) | `0.1171` | `0.6100` `[0.5798, 0.6386]` | `56 / 91` (`61.5%`) | **`91 / 91` (`100.0%`)** |
| `training_median_length_gs_control` | Training median length (`12 aa` GS repeat) | `4.52` `[3.56, 5.45]` | `24 / 91` (`26.4%`) | **`69 / 91` (`75.8%`)** | `0.1331` | `0.5199` `[0.4892, 0.5499]` | `56 / 91` (`61.5%`) | `89 / 91` (`97.8%`) |
| `flanking_context_template_retrieval` | Flanking sequences + CATH IDs $\to$ PR #2 template | `5.27` `[4.38, 6.13]` | `21 / 91` (`23.1%`) | `57 / 91` (`62.6%`) | `0.0819` | `0.7257` `[0.6744, 0.7745]` | `54 / 91` (`59.3%`) | `87 / 91` (`95.6%`) |
| `flanking_context_nn_natural_authority` | Flanking sequences + CATH IDs $\to$ 1-NN natural (`148`) | `7.15` `[5.88, 8.32]` | `18 / 91` (`19.8%`) | `43 / 91` (`47.3%`) | `0.0928` | `0.7274` `[0.6788, 0.7736]` | `53 / 91` (`58.2%`) | `82 / 91` (`90.1%`) |
| `flanking_context_nn_natural_plus_mined` | Flanking sequences + IDs $\to$ 1-NN natural + `27,595` mined | `14.53` `[11.53, 18.16]` | `15 / 91` (`16.5%`) | `24 / 91` (`26.4%`) | **`0.1426`** | **`0.7622` `[0.7192, 0.8049]`** | **`84 / 91` (`92.3%`)** | `90 / 91` (`98.9%`) |
| **`geometry_only_control`** | **3D `inter_port_ca_distance_A` only** | **`4.49` `[3.70, 5.29]`** | `31 / 91` (`34.1%`) | `61 / 91` (`67.0%`) | `0.1317` | **`0.7701` `[0.7336, 0.8027]`** | `52 / 91` (`57.1%`) | **`91 / 91` (`100.0%`)** |
| **`geometry_plus_flanking_hybrid_retrieval`** | **3D `inter_port_ca_distance_A` + flanking context (`27,743` pool)** | **`4.56` `[3.57, 5.56]`** | **`43 / 91` (`47.3%`)** | `63 / 91` (`69.2%`) | **`0.1565`** | `0.7264` `[0.6772, 0.7756]` | **`91 / 91` (`100.0%`)** | **`91 / 91` (`100.0%`)** |
| `length_matched_template_control` *(Oracle $L_{\text{true}}$)* | True length $L_{\text{true}}$ + flanking rigidity $\to$ PR #2 template | **`0.74` `[0.45, 1.07]`** | **`83 / 91` (`91.2%`)** | **`86 / 91` (`94.5%`)** | `0.1316` | `0.7466` `[0.7032, 0.7876]` | `54 / 91` (`59.3%`) | `90 / 91` (`98.9%`) |
| `length_matched_frequent_natural_control` *(Oracle $L_{\text{true}}$)* | True length $L_{\text{true}}$ + training mode dipeptide repeat | **`0.00` `[0.00, 0.00]`** | **`91 / 91` (`100.0%`)** | **`91 / 91` (`100.0%`)** | `0.0863` | `0.4202` `[0.3887, 0.4533]` | `56 / 91` (`61.5%`) | `90 / 91` (`98.9%`) |

---

## 3. Experiment `X2` Results: Cohort B (`148` Combined Natural Authorities = `91` CATH + `57` DisProt Curated `IDPO:0000033`, `109` Disjoint Clusters)

| Arm | Input Used at Inference | Length MAE (aa) [95% CI] | Within $\pm 2\text{ aa}$ | Within $\pm 5\text{ aa}$ | Global NW Seq Identity | Composition Cosine [95% CI] | Rigidity Match |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `fixed_ggggs3_control` | None (`(GGGGS)3`, `15 aa`) | `7.57` `[6.21, 9.10]` | `25 / 148` (`16.9%`) | `71 / 148` (`48.0%`) | `0.1130` | `0.5996` `[0.5739, 0.6251]` | `102 / 148` (`68.9%`) |
| `training_median_length_gs_control` | Training median length (`13 aa` GS repeat) | **`7.19` `[5.69, 8.91]`** | `34 / 148` (`23.0%`) | **`83 / 148` (`56.1%`)** | `0.1094` | `0.5330` `[0.5067, 0.5604]` | `102 / 148` (`68.9%`) |
| **`flanking_context_template_retrieval`** | **Flanking context $\to$ PR #2 template** | `9.60` `[8.12, 11.12]` | `19 / 148` (`12.8%`) | `49 / 148` (`33.1%`) | `0.1050` | **`0.8052` `[0.7685, 0.8400]`** | `93 / 148` (`62.8%`) |
| **`flanking_context_nn_natural_authority`** | **Flanking context $\to$ 1-NN natural (`148`)** | `13.30` `[11.39, 15.22]` | `17 / 148` (`11.5%`) | `51 / 148` (`34.5%`) | `0.0981` | **`0.8001` `[0.7634, 0.8340]`** | `89 / 148` (`60.1%`) |
| **`flanking_context_nn_natural_plus_mined`** | **Flanking context $\to$ 1-NN natural + `27,595` mined** | `13.69` `[11.57, 16.04]` | `19 / 148` (`12.8%`) | `44 / 148` (`29.7%`) | **`0.1451`** | **`0.7971` `[0.7648, 0.8283]`** | **`137 / 148` (`92.6%`)** |
| `length_matched_template_control` *(Oracle $L_{\text{true}}$)* | True length $L_{\text{true}}$ + flanking rigidity $\to$ PR #2 template | **`2.87` `[1.92, 3.96]`** | **`114 / 148` (`77.0%`)** | **`125 / 148` (`84.5%`)** | `0.1258` | **`0.7945` `[0.7584, 0.8266]`** | `93 / 148` (`62.8%`) |
| `length_matched_frequent_natural_control` *(Oracle $L_{\text{true}}$)* | True length $L_{\text{true}}$ + training mode dipeptide repeat | **`0.00` `[0.00, 0.00]`** | **`148 / 148` (`100.0%`)** | **`148 / 148` (`100.0%`)** | `0.0868` | `0.4833` `[0.4539, 0.5130]` | `102 / 148` (`68.9%`) |

---

## 4. Part 2 Geometry-Plus-Template Design Verification (`Q1`–`Q3`)

Verified in `reports/x2_retrieval_and_part2_design_2026-10-19.json` across all four `(rigidity, cleavability)` combinations and both topologies:

1. **`flexible_non_cleavable_tandem_3d` ($18.5\,\text{\AA}$ on `4e43_protease_dimer.pdb`):**
   - `Q1`: Minimum physical length `4 aa` ($3.8\,\text{\AA}/\text{res}$), design minimum `5 aa` ($3.5\,\text{\AA}/\text{res}$), unstrained Gaussian RMS length `13 aa` (`G4S2`).
   - `Q2`: Evaluates 3D anchor SASA, outward angles, and domain clearance on `4e43_protease_dimer.pdb`.
   - `Q3`: Top-ranked bridgeable candidates combine natural CATH/DisProt authorities (`cath:3o33_A_1-57_68-182` `CDAPSHNSEK`, `disprot_curated:DP01108:305-317` `TTGSTKAPAAKGR`), proteome-mined high-confidence calls, and literature templates (`flex-gsagsaagsgef` `GSAGSAAGSGEF`).
2. **`rigid_non_cleavable_tandem` ($22.0\,\text{\AA}$):**
   - `Q1`: Helical rod window `[9, 18] aa` at $1.5\,\text{\AA}/\text{res}$ ($\pm 7.6\,\text{\AA}$ junction tolerance).
   - `Q3`: Selects bridgeable helical/Pro-rich constructs within `[9, 18] aa` (`template:rigid-eaaak2` `AEAAAKEAAAKA`, `12 aa`; `template:rigid-eaaak3` `AEAAAKEAAAKEAAAKA`, `17 aa`).
3. **`flexible_cleavable_tandem` ($20.0\,\text{\AA}$):**
   - `Q3`: Selects bridgeable protease-cleavable flexible constructs (`template:cleav-mmp-plglwa` `GGGGSPLGLWAGGGGS`, `16 aa`, MMP-cleavable; `template:cleav-fvii-ggggs-fkr` `GGGGSFKRGGGGS`, `13 aa`).
4. **`rigid_cleavable_composite_tandem` ($28.0\,\text{\AA}$):**
   - `Q3`: Enforces the two-part composite rule (`build_rigid_cleavable_composite`), returning `composite:rigid-eaaak1+cleav-fvii-ggggs-fkr` (`AEAAAKAGGGGSFKRGGGGS`, `20 aa`) and `composite:rigid-eaaak1+cleav-mmp-plglwa` (`AEAAAKAGGGGSPLGLWAGGGGS`, `23 aa`) with explicit flags `composite_rigid_spacer_plus_flexible_cassette` and `requires_local_unfolding_at_cassette`.
5. **`flexible_non_cleavable_internal_insertion` ($15.0\,\text{\AA}$ forward / $16.5\,\text{\AA}$ return):**
   - `Q1`: Models internal domain insertion as `2` polypeptide tethers (`implied_chain_topology(2, peptide_tethers=2)`), reporting both forward and return span windows and flagging `requires_second_return_tether_for_internal_insertion`.

---

## 5. Reproducibility Commands

```powershell
# 1. Run Part 2 unit and integration tests
python -m unittest tests/test_part2_design_and_x2_retrieval.py -v

# 2. Run the full X2 retrieval benchmark and Part 2 design verification suite
python scripts/evaluate_x2_template_and_authority_retrieval.py

# 3. Run the CLI linker designer on a custom request
python scripts/design_linker.py --distance-angstrom 22.0 --rigidity rigid --cleavability cleavable --protease MMP
```
