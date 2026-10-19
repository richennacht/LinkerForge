# Part 4 Completion Report: Pairwise Compatibility Ranker, Hard Structural Rejection Scorecard, 6-Baseline Benchmark, and `sfGFP (2B3P) — linker — mCherry (2H5Q)` Fusion Demonstration

> **Correction, 19 October 2026:** The ranking results below compare against generated synthetic corruptions and internal heuristics, not published predictors or experimentally measured fusion outcomes. Exported PDBs are geometric illustrations, not validated complete-fusion structures. See [the claims audit](./PAPER_CLAIMS_AUDIT_2026-10-19.md).

**Date:** 2026-10-19  
**Dataset Authority:** `mined_and_natural_linkers_v1` (`27,766` total linker records; `6,727` verified natural + literature + high-confidence proteome-mined authorities indexed; `manifest.json` SHA-256 `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`)  
**Module:** `src/linkerforge/pipeline/ranker.py`  
**Evaluation Script:** `scripts/evaluate_part4_end_to_end_pipeline.py`  
**JSON Artifact:** `reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json`  
**3D Structure & PyMOL Artifacts:** `reports/structures/sfgfp_mcherry_flexible.{pdb,pml}`, `reports/structures/sfgfp_mcherry_rigid.{pdb,pml}`, `reports/structures/sfgfp_mcherry_tev_cleavable.{pdb,pml}`

---

## 1. Architectural Compliance (`ARCHITECTURE_HANDOFF.md` §7 & §8)

1. **Deterministic Hard Structural Rejection Scorecard (`evaluate_hard_rejection_scorecard`):**
   - Evaluates 6 explicit physical & sequence rules before ranking:
     1. `physical_span_closure_pass`: Checks $L \times 3.8\text{ Å} \ge d_{\text{path}}$ (physical virtual $C_\alpha$ ceiling).
     2. `domain_clash_clearance_pass`: Checks steric detour ratio $\le 2.2$, forbids $0\text{-aa}$ direct fusion when $d_{\text{path}} > 3.8\text{ Å}$, and flags $>2.2\times$ entropic loop excess.
     3. `anchor_accessibility_pass`: Requires both flanking domain attachment ports to have relative solvent accessibility $\text{rSASA} \ge 0.15$.
     4. `hydrophobic_stretch_pass`: Requires hydrophobic fraction (`FILVWM`) $< 0.50$ and maximum contiguous hydrophobic run $< 5$ residues.
     5. `prohibited_cleavage_pass`: Enforces absence of unintended protease sites when `cleavable=False` and presence of requested protease recognition site when `cleavable=True`.
     6. `composite_rigid_cleavable_pass`: Enforces that protease sites are flanked by flexible `GS` spacers rather than buried inside a rigid `EAAAK` $\alpha$-helix (`Hubbard et al. 1994`; `Chen et al. 2013`).
   - **Strict Scope Compliance (`SPEC.md` §5.8 / `ARCHITECTURE_HANDOFF.md` §8):** Never reports a calibrated "fusion success probability".

2. **Six Synthetic Negative Corruptions (`generate_synthetic_negatives`):**
   - `too_short_by_50pct`: Truncates linker length by $50\%$ so physical span closure fails.
   - `too_long_by_200pct`: Triples linker length ($+200\%$) inducing excessive entropic loop penalty.
   - `rigid_helical_where_flexible_needed`: Substitutes rigid `(EAAAK)_n` helical rod when flexible coil is required (or vice versa).
   - `all_hydrophobic_core`: Substitutes `FILVWM`-rich hydrophobic core (`LVLFLILVWMF...`) inducing collapse/aggregation.
   - `charge_clash_with_junction_anchors`: Constructs an uncompensated homopolymeric charge block (`KKKK...` or `EEEE...`) clashing with flanking port electrostatics.
   - `protease_site_buried_in_rigid_helix`: Embeds `ENLYFQS` directly inside an unbroken `EAAAK` $\alpha$-helix (`EAAAKEAAAKENLYFQSEAAAKEAAAK`) without flexible `GS` spacers.

3. **Learned Pairwise Compatibility Ranker (`PairwiseCompatibilityRanker`):**
   - Trained with margin-regularized pairwise logistic ranking loss over 12 interpretable geometric, biophysical, and junction-compatibility features on cluster-disjoint (`homology_cluster_id`) splits from `mined_and_natural_linkers_v1`.

---

## 2. Cluster-Disjoint 5-Fold Cross-Validation Results

### 2.1 Natural / Mined Linkers vs. 6 Synthetic Negative Corruptions

| Benchmark Cohort | Evaluated Records | Top-1 Natural Rank Accuracy | Mean Reciprocal Rank (MRR) | vs. `too_short_by_50pct` | vs. `too_long_by_200pct` | vs. `rigid_helical_wrong` | vs. `all_hydrophobic_core` | vs. `charge_clash` | vs. `buried_protease` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Verified Natural Authorities (`148`)** | `148` | **`99.32%`** | **`0.9966`** | `100.00%` | `100.00%` | `99.32%` | `100.00%` | `100.00%` | `100.00%` |
| **Proteome-Mined High-Confidence (`600`)** | `600` | **`98.67%`** | **`0.9931`** | `100.00%` | `100.00%` | `98.67%` | `100.00%` | `99.83%` | `100.00%` |

### 2.2 Six Required Baselines Comparison (`148` Verified Natural Authorities)

| Baseline / Method | Hard-Rejection Pass Rate | Mean Compatibility Score | Length Within $\pm 2\text{ aa}$ | Rigidity Axis Alignment | Win Rate vs. 6 Synthetic Negatives |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`linkerforge_ranker` (Full Pipeline)** | **`100.00%`** | **`10.4340` (Rank #1)** | **`99.32%`** | **`0.9803`** | **`100.00%`** |
| `random_length_matched_gs` | `100.00%` | `10.1562` | `100.00%` | `0.9865` | `99.77%` |
| `nearest_neighbour_retrieval` (`X2` over `6,727` authorities) | `100.00%` | `9.7191` | `95.95%` | `0.8185` | `99.66%` |
| `geometry_only` | `100.00%` | `9.3986` | `100.00%` | `0.5092` | `100.00%` |
| `fixed_ggggs3` (`15 aa`) | `85.14%` | `8.8044` | `17.57%` | `0.9865` | `94.03%` |
| `direct_0aa_fusion` (`0 aa`) | `0.00%` | `-2.1164` | `0.00%` | `0.0000` | `0.00%` |

### 2.3 Six Required Baselines Comparison (`600` High-Confidence Proteome-Mined Linkers)

| Baseline / Method | Hard-Rejection Pass Rate | Mean Compatibility Score | Length Within $\pm 2\text{ aa}$ | Rigidity Axis Alignment | Win Rate vs. 6 Synthetic Negatives |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`linkerforge_ranker` (Full Pipeline)** | **`100.00%`** | **`10.2848` (Rank #1)** | **`100.00%`** | **`0.9966`** | **`100.00%`** |
| `random_length_matched_gs` | `100.00%` | `10.2215` | `100.00%` | `0.9850` | `99.81%` |
| `nearest_neighbour_retrieval` (`X2` over `6,727` authorities) | `100.00%` | `9.6496` | `98.00%` | `0.8027` | `99.75%` |
| `geometry_only` | `100.00%` | `9.3471` | `100.00%` | `0.5299` | `100.00%` |
| `fixed_ggggs3` (`15 aa`) | `46.00%` | `6.7658` | `10.83%` | `0.9850` | `79.08%` |
| `direct_0aa_fusion` (`0 aa`) | `0.00%` | `-2.1647` | `0.00%` | `0.0000` | `0.00%` |

---

## 3. End-to-End `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)` Fusion Demonstration

- **Domain A (`sfGFP`, PDB `2B3P`, `238 aa`):** Selected C-terminal port `Lys238` (`...HGMDELYK`, $\text{rSASA} = 0.62$).
- **Domain B (`mCherry`, PDB `2H5Q`, `236 aa`):** Selected N-terminal port `Met1` (`MVSKGEED...`, $\text{rSASA} = 0.58$).
- **3D Junction Geometry:** Straight-line $C_\alpha$ distance $d_{CA} = 24.8\text{ Å}$, collision-free surface detour path $d_{\text{path}} = 28.4\text{ Å}$ (detour ratio $1.145$).
- **Top Ranked Candidates by User Intent Mode:**
  1. **`flexible_15aa` Mode:**
     - **Rank #1 (`linkerforge_ranker`):** `GGGGSGGGGSGGGGS` (`15 aa`, Score `10.3379`, `HardPass=True`, contour `57.0 Å` $\ge 28.4\text{ Å}$).
     - **Top Proteome-Mined Authority (`X2`):** `ARGGQYSNQGVVAST` (`15 aa`, `mined_swissprot_gap:Q9XTD6:272-286`, Score `10.1476`, `HardPass=True`).
     - **Waldo GFP-Compatible Template:** `GSAGSAAGSGEFGSA` (`15 aa`, Score `10.1155`, `HardPass=True`).
     - **All 6 Synthetic Negatives & `direct_0aa_fusion`:** Deterministically rejected (`HardPass=False`) with explicit rejection reasons (`physical_span_closure_failed`, `domain_clash_or_entropic_excess`, `hydrophobic_collapse_risk`, `unintended_protease_cleavage_motif_present`, `protease_site_buried_inside_rigid_alpha_helix_without_gs_spacers`).
  2. **`rigid_18aa` Mode:**
     - **Rank #1 (`linkerforge_ranker`):** `AEAAAKEAAAKEAAAKA` (`17 aa`, Chen et al. 2013 helical rod, `HardPass=True`).
  3. **`flexible_tev_cleavable_19aa` Mode:**
     - **Rank #1 (`linkerforge_ranker`):** `GSAGSAENLYFQSGSAGSA` (`19 aa`, TEV recognition site flanked by flexible `GSAGSA` spacers, `HardPass=True`).
- **3D Coordinate & PyMOL Artifacts:** Complete `489–493` residue (`1,956–1,972` backbone atoms `N, CA, C, O`) 3D fusion `.pdb` files and `.pml` PyMOL visualization scripts exported to `reports/structures/`.
