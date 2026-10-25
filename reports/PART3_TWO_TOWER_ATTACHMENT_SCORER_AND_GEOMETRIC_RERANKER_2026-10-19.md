# Part 3 Diagnostic Report: Two-Tower Attachment-Site Scorer, Position-Matched Evaluation & 3D Geometric Site-Pair Re-Ranker

> **Audit & Remediation Note (2026-10-19):** The earlier Part 3 evaluation defined positives as the 3 residues adjacent to the linker at one end of the flanking window and negatives as interior residues of the same window while keeping window-coordinate features (`term_prox`, `rel_interior`) active and passing label-conditioned SASA/secondary-structure proxies in the structure arm. In this revision (`schema_version: 2.0.0`), positional weights are fixed to zero (`POSITIONAL_FEATURE_INDICES = (9, 11)`), label-derived structure proxies are removed, an explicit leakage diagnostic quantifies the original positional shortcut, and the primary benchmark is replaced with a **position-matched** 5-fold cluster-disjoint evaluation (`matched_offset_positions`).

**Date:** 2026-10-19  
**Dataset Authority:** `mined_and_natural_linkers_v1` (`manifest.json` SHA-256 `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`)  
**Module:** `src/linkerforge/attachment/scorer.py`  
**Evaluation Script:** `scripts/evaluate_part3_attachment_scorer.py`  
**JSON Artifact:** `reports/part3_attachment_scorer_2026-10-19.json`

---

## 1. Architectural Compliance (`ARCHITECTURE_HANDOFF.md` §6 / Component 4)

1. **Two-Tower Domain Attachment Encoder (`TwoTowerAttachmentScorer`):**
   - Scores residue positions $i \in \{1, \dots, N\}$ of a flanking protein domain using local biophysical and structural features.
   - **Strict Tick-Box Separation (`ARCHITECTURE_HANDOFF.md` §6):** User intent tick boxes (`flexible`, `rigid`, `cleavable`, `non_cleavable`) **never enter the two-tower encoder** (`tick_boxes_excluded_from_encoder: true`).
   - **Positional Shortcut Elimination:** Window-coordinate features (`term_prox`, `rel_interior`, indices `9` and `11`) have weights fixed to `0.0` during fitting (`positional_feature_weights_fixed_to_zero: true`) so the model cannot separate candidates by their index within a truncated flanking window.

2. **Closed-Form 3D Geometric Site-Pair Re-Ranker (`rank_attachment_site_pairs`):**
   - Given candidate attachment sites on Domain A and Domain B with 3D coordinates, calls `evaluate_junction` (`linkerforge.geometry.junction`) and `implied_chain_topology` (`linkerforge.geometry.tether`) to compute straight-line $C_\alpha$ distance ($d_{CA}$), obstacle-aware detour path length ($d_{\text{path}}$), outward exit-vector orientation, and required linker length windows.

---

## 2. Leakage Diagnostic on Original Unconstrained Task

In the original task design (`last/first 3 residues vs interior residues of the same window` across `109,226` candidate residues, `12,852` positives), single window-coordinate features alone separate the labels almost perfectly without reading amino-acid identity or 3D structure:
- **`term_prox` alone:** $\text{AUROC} = 1.0000$
- **`rel_interior` alone:** $\text{AUROC} = 0.9998$

This confirms that the earlier `1.0000` AUROC reflected window-edge position rather than sequence or structural preference.

---

## 3. Position-Matched Cluster-Disjoint 5-Fold Cross-Validation Results

To remove window-coordinate confounding, `matched_offset_positions` samples 3 positive residues at offsets $\{0, 1, 2\}$ from the linker-adjacent end of each flanking window ($\text{length} \ge 12\text{ aa}$) and 3 negative residues at the exact mirror offsets $\{0, 1, 2\}$ from the opposite end of the same window (`chance_auroc = 0.5000`), evaluated across 5 `homology_cluster_id`-disjoint folds with 1,000-replicate cluster bootstrap 95% confidence intervals:

### 3.1 Verified Natural Linker Authorities (`148` records; `293` windows; `109` homology clusters; `1,758` candidate residues)

| Scoring Arm | Pooled AUROC | 95% Cluster Bootstrap CI | Mean Within-Window AUROC |
| :--- | :---: | :---: | :---: |
| `single_residue_flexibility_propensity` | `0.5180` | `[0.4874, 0.5476]` | `0.5188` |
| `random_control` | `0.4953` | `[0.4723, 0.5187]` | `0.5009` |
| `sequence_context_logistic` (`TwoTowerAttachmentScorer`) | `0.4756` | `[0.4162, 0.5300]` | `0.5059` |

### 3.2 High-Confidence Proteome-Mined Linkers Subsample (`2,000` sampled of `6,556` records; `3,991` windows; `461` homology clusters; `23,946` candidate residues)

| Scoring Arm | Pooled AUROC | 95% Cluster Bootstrap CI | Mean Within-Window AUROC |
| :--- | :---: | :---: | :---: |
| `random_control` | `0.5059` | `[0.4998, 0.5115]` | `0.5076` |
| `single_residue_flexibility_propensity` | `0.4917` | `[0.4743, 0.5070]` | `0.4924` |
| `sequence_context_logistic` (`TwoTowerAttachmentScorer`) | `0.4589` | `[0.4211, 0.4958]` | `0.4646` |

**Scientific Conclusion:** Once distance from the window edge is controlled, local sequence-only amino-acid propensities within a 30-residue domain flank do not distinguish the linker-facing boundary from the opposite end of the window ($\text{AUROC} \approx 0.46\text{--}0.52$, indistinguishable from `0.50` chance). Meaningful attachment-site ranking therefore depends on **3D structural geometry, solvent accessibility (`rSASA`), secondary-structure core avoidance, and steric path clearance**, not 1D sequence propensities alone.

---

## 4. 3D Structure Site-Pair Re-Ranking (`tests/fixtures/geometry/4e43_protease_dimer.pdb`)

- **`end_to_end_tandem` Top-1 Site Pair:** Domain A surface loop `Ile50` ($\text{rSASA} = 0.9484$, curvature $74.02^\circ$) $\rightarrow$ Domain B surface loop `Pro81` ($\text{rSASA} = 0.8365$, curvature $90.77^\circ$), $d_{CA} = 6.713\text{ Å}$, $d_{\text{path}} = 6.713\text{ Å}$, `closed_form_geometric_score = 0.9038`, `combined_site_pair_score = 0.8465`, flags `["both_anchors_surface_accessible", "span_window_bridgeable"]`.
- **`internal_insertion_two_tether` Top-1 Site Pair:** Host Domain A surface loop `Ile50` $\leftrightarrow$ Inserted Domain B surface port `Pro81`, validated with `implied_chain_topology(2, peptide_tethers=2)` (`tether_count = 2`).
