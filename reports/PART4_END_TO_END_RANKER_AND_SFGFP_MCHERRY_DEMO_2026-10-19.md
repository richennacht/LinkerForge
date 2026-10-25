# Part 4 Diagnostic Report: Hard Rejection Scorecard, Cross-Cluster Real-Decoy Ranking, Non-Oracle Generator Comparison, and Real-Coordinate `sfGFP (2B3P) — linker — mCherry (2H5Q)` Schematic Assembly

> **Audit & Remediation Note (2026-10-19):** This revision (`schema_version: 2.0.0`) replaces circular baseline comparisons and procedural PDB coordinates with four non-circular evaluations (`src/linkerforge/pipeline/evaluation.py`) and real experimental crystal coordinates downloaded from RCSB (`data/raw/pdb/2B3P.pdb` and `data/raw/pdb/2H5Q.pdb` in `src/linkerforge/pipeline/fusion_demo.py`). Missing geometry or SASA fields are reported in `not_evaluated` rather than filled with invented defaults, and `PairwiseCompatibilityRanker` is trained from zero initialization (`init="zero"`) without hand-seeded weights.

**Date:** 2026-10-19  
**Dataset Authority:** `mined_and_natural_linkers_v1` (`manifest.json` SHA-256 `b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`)  
**Modules:** `src/linkerforge/pipeline/ranker.py`, `src/linkerforge/pipeline/evaluation.py`, `src/linkerforge/pipeline/fusion_demo.py`  
**Evaluation Script:** `scripts/evaluate_part4_end_to_end_pipeline.py`  
**JSON Artifact:** `reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json`  
**Schematic Structure Artifacts:** `reports/structures/sfgfp_mcherry_flexible_non_cleavable.{pdb,pml}`, `reports/structures/sfgfp_mcherry_rigid_non_cleavable.{pdb,pml}`, `reports/structures/sfgfp_mcherry_flexible_cleavable_factor_xa.{pdb,pml}`

---

## 1. Deterministic Hard Structural Rejection Scorecard on Real Linkers (`scorecard_on_real_linkers`)

Evaluates the 6 deterministic necessary-condition rules (`physical_span_closure_pass` using $(L+1)\times 3.8\text{ Å} \ge d_{\text{path}}$, `domain_clash_clearance_pass`, `anchor_accessibility_pass`, `hydrophobic_stretch_pass`, `prohibited_cleavage_pass`, `composite_rigid_cleavable_pass`) across all `27,766` records in `mined_and_natural_linkers_v1`, evaluating geometry rules only when real 3D coordinates exist (`n_with_real_geometry = 91` on `natural_cath_pdb_verified`):

| Dataset Tier | Total Records | Real 3D Geometry | All Evaluated Rules Pass | Rule Failures |
| :--- | :---: | :---: | :---: | :--- |
| `natural_cath_pdb_verified` | `91` | `91` | **`90 / 91` (`98.90%`)** | `physical_span_closure_pass: 1` |
| `natural_disprot_curated_dfl` | `57` | `0` | **`57 / 57` (`100.00%`)** | None |
| `literature_engineered_template` | `23` | `0` | **`21 / 23` (`91.30%`)** | `hydrophobic_stretch_pass: 1`, `prohibited_cleavage_pass: 1` |
| `mined_high_confidence_dfl_call` | `6,556` | `0` | **`6,546 / 6,556` (`99.85%`)** | `hydrophobic_stretch_pass: 10` |
| `mined_provisional_dfl_candidate` | `21,039` | `0` | **`20,981 / 21,039` (`99.72%`)** | `hydrophobic_stretch_pass: 54`, `prohibited_cleavage_pass: 4` |

---

## 2. Non-Circular Cross-Cluster Real-Decoy Ranking Benchmark (`decoy_ranking`)

To test whether `PairwiseCompatibilityRanker` (fit from `init="zero"` on 5 `homology_cluster_id`-disjoint folds) captures genuine sequence–flank compatibility beyond length and rigidity, each held-out true linker is ranked against **19 real linkers from other homology clusters** matched on exact length ($\pm 1\text{ aa}$ fallback) and rigidity class (`mean_candidates_per_query = 20.0`, chance `Hit@1 = 0.0500`, `Hit@5 = 0.2500`, `MRR = 0.1799`):

### 2.1 Verified Natural Linker Queries vs. Cross-Cluster Real Decoys (`148` queries)

| Ranking Method | Hit@1 (95% Cluster CI) | Hit@5 (95% Cluster CI) | MRR (95% Cluster CI) |
| :--- | :---: | :---: | :---: |
| **`linkerforge_ranker` (`init="zero"`, 5-fold OOF)** | **`0.2635` (`[0.1875, 0.3355]`)** | **`0.5811` (`[0.4902, 0.6624]`)** | **`0.4260` (`[0.3515, 0.4933]`)** |
| `flank_conditioned_nn` (flank Kidera/charge NN) | `0.0878` (`[0.0376, 0.1384]`) | `0.2568` (`[0.1793, 0.3241]`) | `0.2039` (`[0.1553, 0.2483]`) |
| `random` | `0.0608` (`[0.0268, 0.1027]`) | `0.3108` (`[0.2414, 0.3741]`) | `0.1994` (`[0.1628, 0.2398]`) |
| `training_composition_prior` | `0.0541` (`[0.0178, 0.0942]`) | `0.2568` (`[0.1688, 0.3423]`) | `0.1851` (`[0.1390, 0.2276]`) |
| Exact Uniform Chance (1 in 20) | `0.0500` | `0.2500` | `0.1799` |

### 2.2 High-Confidence Mined Queries vs. Unsampled Cross-Cluster Real Decoys (`598` evaluated queries)

| Ranking Method | Hit@1 (95% Cluster CI) | Hit@5 (95% Cluster CI) | MRR (95% Cluster CI) |
| :--- | :---: | :---: | :---: |
| **`linkerforge_ranker` (`init="zero"`, 5-fold OOF)** | **`0.1405` (`[0.1041, 0.1844]`)** | **`0.4398` (`[0.3724, 0.5092]`)** | **`0.2996` (`[0.2597, 0.3435]`)** |
| `flank_conditioned_nn` | `0.0719` (`[0.0446, 0.0946]`) | `0.3110` (`[0.2573, 0.3582]`) | `0.2141` (`[0.1836, 0.2413]`) |
| `random` | `0.0569` (`[0.0391, 0.0738]`) | `0.2876` (`[0.2523, 0.3199]`) | `0.1990` (`[0.1795, 0.2174]`) |
| `training_composition_prior` | `0.0552` (`[0.0331, 0.0803]`) | `0.2542` (`[0.2083, 0.3008]`) | `0.1845` (`[0.1595, 0.2107]`) |
| Exact Uniform Chance (`19.96` candidates/query) | `0.0501` | `0.2506` | `0.1802` |

---

## 3. Non-Oracle Generator Comparison on Coordinate-Verified CATH Linkers (`91` queries)

All 6 generators receive **only** the inter-port 3D distance $d_{CA}$ and flanking sequences (never the held-out natural linker's true residue length), under two pre-declared geometry length policies (`unstrained_reach` and `design_minimum_plus_2`):

| Generator Arm (`design_minimum_plus_2` policy) | Span Closure Pass | All Scorecard Rules Pass | Length Within $\pm 2\text{ aa}$ of Natural | Mean Abs. Length Error (aa) | Composition Cosine to Natural |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `linkerforge_ranker` | `1.0000` | `1.0000` | `0.4066` | `7.09` | `0.4092` |
| `nearest_neighbour_retrieval` | `1.0000` | `1.0000` | `0.3956` | `7.04` | `0.4095` |
| `geometry_only` | `1.0000` | `1.0000` | `0.3956` | `7.60` | `0.3712` |
| `random_length_matched_gs` | `1.0000` | `1.0000` | `0.3956` | `7.60` | `0.3811` |
| `fixed_ggggs3` (`15 aa`) | `1.0000` | `1.0000` | `0.1978` | `6.14` | `0.3712` |
| `direct_0aa_fusion` (`0 aa`) | `0.1099` | `0.1099` | `0.0000` | `11.15` | `0.0000` |

*(Note: The synthetic corruption sanity check `corruption_check` is retained separately in the JSON artifact as a by-construction training-objective verification, achieving `95.27%` top-1 / `0.9764` MRR on `natural_148` and `97.17%` top-1 / `0.9854` MRR on `mined_high_confidence_sample_600`.)*

---

## 4. Real-Coordinate `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)` Schematic Assembly (`fusion_demo`)

- **Domain A (`sfGFP`, experimental PDB `2B3P` chain A):** `228` resolved residues, last resolved C-terminal residue `Gly232` ($\text{rSASA} = 1.0461$), `14` unresolved C-terminal tail residues in SEQRES.
- **Domain B (`mCherry`, experimental PDB `2H5Q` chain A):** `217` resolved residues, first resolved N-terminal residue `Asn4` ($\text{rSASA} = 0.3476$), `8` unresolved N-terminal tail residues in SEQRES.
- **Seeded Rigid-Body Placement (`seed=7`, `400` rotations):** Minimum inter-domain heavy-atom distance `4.005 Å` ($\ge 4.0\text{ Å}$ steric clearance threshold), straight-line port $C_\alpha$ distance $d_{CA} = 5.062\text{ Å}$, detour path length $d_{\text{path}} = 5.062\text{ Å}$.
- **Evaluated Design Scenarios (`status: schematic_assembly_on_real_domain_coordinates_not_a_structure_prediction`):**
  1. **`flexible_non_cleavable` (`target_length = 5 aa`):** Top candidate `IKDSN` (`x2:mined_uniprot_explicit:Q9ZDE6:141-145`, score `6.2607`, `hard_rejection_pass = True`), followed by `SKGQA` (`6.1827`) and CATH authority `GGADN` (`x2:cath:5xix_A_1-102_108-437`, `6.1376`).
  2. **`rigid_non_cleavable` (`target_length = 5 aa`):** Top candidate `APEKK` (`x2:mined_swissprot_gap:Q4VC17:1177-1181`, score `5.3959`, `hard_rejection_pass = True`), followed by literature template `PAPAP` (`library:rigid-papap`, Chen et al. 2013, score `5.3270`).
  3. **`flexible_cleavable_factor_xa` (`target_length = 5 aa`, Factor Xa motif `IEGR`):** Top candidate `GGIEGRGS` (`8 aa`, score `5.4406`, `hard_rejection_pass = True`), while uncleavable baselines (`GGSSG`, `GGGGSGGGGSGGGGS`) are deterministically rejected with `missing_requested_protease_motif`.
