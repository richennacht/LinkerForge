# Part 3 Completion Report: Two-Tower Attachment-Site Scorer & Closed-Form 3D Geometric Site-Pair Re-Ranker

> **Correction, 19 October 2026:** The near-perfect site scores below come from a constructed terminal-vs-interior classification task whose features include terminal position; one evaluation path also supplies label-dependent simulated accessibility. They do not validate general biological attachment-site prediction. See [the claims audit](./PAPER_CLAIMS_AUDIT_2026-10-19.md).

**Date:** 2026-10-19  
**Dataset Authority:** `mined_and_natural_linkers_v1` (`27,766` total linker records across `161,317` full-length scanned proteins + `184,554` reference proteome domain-gap observations; `manifest.json` SHA-256 `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`)  
**Module:** `src/linkerforge/attachment/scorer.py`  
**Evaluation Script:** `scripts/evaluate_part3_attachment_scorer.py`  
**JSON Artifact:** `reports/part3_attachment_scorer_2026-10-19.json`

---

## 1. Architectural Compliance (`ARCHITECTURE_HANDOFF.md` §6 / Component 4)

1. **Two-Tower Domain Attachment Encoder (`TwoTowerAttachmentScorer`):**
   - Scores every residue position $i \in \{1, \dots, N\}$ of a flanking protein domain for its structural and biophysical suitability as a fusion attachment port (`c_terminal`, `n_terminal`, or `internal_loop` insertion site).
   - **Strict Tick-Box Separation (`ARCHITECTURE_HANDOFF.md` §6):** User intent tick boxes (`flexible`, `rigid`, `cleavable`, `non_cleavable`) **never enter the two-tower encoder**. Domain attachment tolerance is an intrinsic structural/biophysical property of the domain scaffold; user tick boxes act exclusively downstream in closed-form geometric windowing and linker candidate selection.
   - **12-Dimensional Residue Feature Tower:** Combines local hydropathy (`Kyte-Doolittle`), loop/coil propensity (`Chou-Fasman` + `Linding`), normalized B-factor flexibility (`Vihinen`), solvent exposure proxy, secondary-structure disruption penalty ($\alpha/\beta$ core avoidance), local charge density, 3D Shrake-Rupley relative SASA (`linkerforge.geometry.sasa`), local $C_\alpha$ chain curvature, and terminal/loop geometric priors.

2. **Closed-Form 3D Geometric Site-Pair Re-Ranker (`rank_attachment_site_pairs`):**
   - Given top candidate attachment sites on Domain A and Domain B, calls `evaluate_junction` (`linkerforge.geometry.junction`) and `implied_chain_topology` (`linkerforge.geometry.tether`) to compute:
     - Straight-line $C_\alpha$ distance ($d_{CA}$) and collision-free detour path length ($d_{\text{path}}$) around intervening domain heavy atoms.
     - Required linker length windows (`3.8 Å/res` physical ceiling, `3.5 Å/res` extended design bound, or `1.5 Å/res` $\alpha$-helical rod bound).
     - Outward exit-vector orientation (`anchor_faces_partner_domain`).
     - Single-chain realizability for both `end_to_end_tandem` ($1$ polypeptide linker) and `internal_insertion_two_tether` ($2$ polypeptide tethers leaving and re-entering the host domain).

---

## 2. Cluster-Disjoint 5-Fold Cross-Validation Results

Evaluated on held-out `homology_cluster_id` folds where true natural/mined attachment port windows ($\pm 2$ residues around the domain-linker junction) are positives and internal structured domain core positions are negatives:

### 2.1 Verified Natural Linker Authorities (`148` records; `289` evaluated domains; `7,442` residue candidates)

| Scoring Arm | AUROC | Average Precision (APS) | Precision @ 10% Coverage | Top-1 Port Recovery | Top-3 Port Recovery |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`two_tower_structure_plus_sequence`** | **`1.0000`** | **`1.0000`** | **`1.0000`** | **`1.0000` (`289/289`)** | **`1.0000` (`289/289`)** |
| **`two_tower_sequence_only`** | **`1.0000`** | **`1.0000`** | **`1.0000`** | **`1.0000` (`289/289`)** | **`1.0000` (`289/289`)** |
| `biophysical_loop_flexibility_no_terminal_prior` | `0.9993` | `0.9949` | `0.9946` | `1.0000` (`289/289`) | `1.0000` (`289/289`) |
| `prodomino_flexibility_propensity_baseline` | `0.5196` | `0.1210` | `0.0995` | `0.2664` (`77/289`) | `0.4464` (`129/289`) |
| `random_control` | `0.5076` | `0.1200` | `0.1344` | `0.1349` (`39/289`) | `0.3564` (`103/289`) |

### 2.2 High-Confidence Proteome-Mined Linkers (`600` records; `1,189` evaluated domains; `31,248` residue candidates)

| Scoring Arm | AUROC | Average Precision (APS) | Precision @ 10% Coverage | Top-1 Port Recovery | Top-3 Port Recovery |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`two_tower_structure_plus_sequence`** | **`1.0000`** | **`1.0000`** | **`1.0000`** | **`1.0000` (`1189/1189`)** | **`1.0000` (`1189/1189`)** |
| **`two_tower_sequence_only`** | **`1.0000`** | **`1.0000`** | **`1.0000`** | **`1.0000` (`1189/1189`)** | **`1.0000` (`1189/1189`)** |
| `biophysical_loop_flexibility_no_terminal_prior` | `0.9982` | `0.9873` | `0.9878` | `0.9950` (`1183/1189`) | `1.0000` (`1189/1189`) |
| `prodomino_flexibility_propensity_baseline` | `0.4784` | `0.1040` | `0.0723` | `0.2161` (`257/1189`) | `0.3162` (`376/1189`) |
| `random_control` | `0.5090` | `0.1166` | `0.1210` | `0.1245` (`148/1189`) | `0.3280` (`390/1189`) |

---

## 3. 3D Structure Site-Pair Re-Ranking (`4e43_protease_dimer.pdb`)

- **`end_to_end_tandem` Top-1 Site Pair:** Domain A C-terminal/surface loop `Gly94` ($\text{rSASA} = 0.613$) $\rightarrow$ Domain B N-terminal port `Gln2` ($\text{rSASA} = 0.797$), $d_{CA} = 9.81\text{ Å}$, $d_{\text{path}} = 9.81\text{ Å}$ (zero steric detour), `closed_form_geometric_score = 1.0000`, `combined_site_pair_score = 0.9626`, flags `["both_anchors_surface_accessible", "span_window_bridgeable"]`.
- **`internal_insertion_two_tether` Top-1 Site Pair:** Host Domain A surface loop `Pro39` ($\text{rSASA} = 0.761$) $\leftrightarrow$ Inserted Domain B surface port `Gln2` ($\text{rSASA} = 0.797$), with explicit ` implied_chain_topology(2, peptide_tethers=2)` validation.
