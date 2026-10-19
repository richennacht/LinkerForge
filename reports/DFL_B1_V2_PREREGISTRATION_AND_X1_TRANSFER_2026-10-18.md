# Firewall V2 Preregistration (`B1` Deadlock Resolution) and Held-Out Domain-Pair Relative-Position Transfer (`X1`) — 2026-10-18

## Executive Summary

This report documents two deliverables executed under the scientific integrity constraints of `SPEC.md` §5:

1. **Resolution of the Strand A (`V1` / `B1`) Firewall Deadlock via Preregistered `RULE_V2`:**
   - We diagnose two mechanical defects that caused the V1 sequence firewall (`scripts/preflight_dfl_family_training.py`) and strict DisProt filter (`scripts/build_dfl_curated_evidence_pool.py`) to collapse the `DisProt_2026_06` training pool to 2 proteins and 0 negative (`0`) residues.
   - We formalize the **Firewall V2 Preregistration Contract** in `literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md`, enforce it in `schemas/dfl_family_split_v2.schema.json` and `scripts/preflight_dfl_family_training_v2.py`, and construct the V2 family-disjoint split via `scripts/build_dfl_b1_v2_split.py`.
   - Under `RULE_V2` (`BLOSUM62`, gap `10/1`, $\ge 30\%$ identity, $\ge 50$ aligned columns, $\ge 80\%$ true ungapped residue coverage on **both** sequences, plus transitive `Pfam` + `Gene3D` + `UniRef50` family-component disjunction), the shuffle-null false-positive rate drops from **`1.625%` (`78/4,800` residue-shuffled pairs, or `66.1%` of V1's observed hits)** to **`0.000%` (`0/4,800` residue-shuffled pairs)** and **`0.021%` (`1/4,800` 12-mer block-shuffled pairs)**.
   - Retaining publicly released (`2022_03`–`2026_06`), PMID/DOI-cited, curator-validated, non-obsolete `DisProt_2026_06` annotations that carried DisProt's internal `unpublished: true` legacy flag yields **97 benchmark-disjoint proteins across 91 family components** (`60` train, `28` development, `9` auxiliary order calibration panel), providing **`1,920` positive (`1`, `IDPO:0000033`) residues** and **`1,253` verified experimental-order negative (`0`, `IDPO:0000006`) residues** while keeping all unannotated and other-disorder residues masked as unknown (`x`).
2. **Execution of Decisive Experiment `X1` (`CHECKLIST.md` item `X1`):**
   - Using `scripts/evaluate_x1_relative_position_transfer.py` on **91 coordinate-complete CATH v4.4.0 natural two-domain proteins across 57 CATH domain-pair groups** (`data/processed/cath_full_pilot.jsonl`, SHA-256 `20deca0e0bff4d1c68b245cc4e15756fc4e6462537c67a7d4508b5f8f58cfd85`), we test whether 3D inter-domain relative-position features transfer across held-out domain pairs.
   - Across 5 seeds (`7, 17, 31, 43, 97`; `135` pooled held-out evaluations with zero domain-pair architecture leakage), adding 3D inter-domain relative-position features to a linear pairwise ranker (`sequence_plus_relative_position_3d`, mean MRR `0.0684`, Top-1 `0/135`, Top-5 `4/135`) **does not improve** over `sequence_only` (mean MRR `0.0709`, Top-1 `0/135`, Top-5 `5/135`; paired cluster-bootstrap $\Delta\text{MRR} = -0.0025$, 95% CI $[-0.0173, +0.0058]$) or `sequence_plus_local_port` (mean MRR `0.0697`, Top-1 `0/135`, Top-5 `6/135`).
   - By contrast, a zero-parameter physical geometry baseline (`geometry_span_match`, enforcing the $3.8\,\text{\AA}/\text{residue}$ physical C$\alpha$ reach bound and ranking by proximity to $d_{AB}/2.5\,\text{\AA}$) achieves **mean MRR `0.1710`, Top-1 `7/135` (`5.19%`), and Top-5 `37/135` (`27.41%`)**, outperforming all learned feature-concatenation rankers and sequence heuristics (`length_heuristic` mean MRR `0.1160`, `fixed_gs` mean MRR `0.1066`).

---

## Part 1: Diagnosis and Resolution of the V1 / B1 Firewall Deadlock

### 1.1 Defect 1: Unbounded Gapped-Length Coverage and Composition False Positives in V1 `local_hit`

In `scripts/preflight_dfl_family_training.py`, `local_hit(left, right)` defined coverage as:

$$\text{coverage}_{\text{V1}} = \frac{\text{len}(\text{traceback.query})}{\min(|A|, |B|)} \ge 0.50, \quad \text{identity} \ge 0.25, \quad \text{len}(\text{traceback.query}) \ge 50$$

Because `len(traceback.query)` counts alignment columns **including gap characters (`-`)**, $\text{coverage}_{\text{V1}}$ is not bounded by $1.0$ (reaching $1.0286$ on synthetic Gly/Ser vs. poly-Gly sequences) and acts as an amino-acid composition detector on disordered proteins rather than a homology detector.

Measured across `4,800` candidate–`DP94` pairs (`3` seeds $\times 40 \times 40$, seeds `173, 271, 811`):

| Alignment Rule | Observed Candidate–DP94 Hits | Residue-Shuffled Null Hits (Rate) | 12-mer Block-Shuffled Null Hits (Rate) | Null / Observed Ratio (Residue) |
| :--- | :---: | :---: | :---: | :---: |
| **`RULE_V1`** (`25%` id, $\ge 50$ cols, $\ge 50\%$ gapped length of shorter) | `118 / 4,800` (`2.458%`) | `78 / 4,800` (`1.625%`) | `59 / 4,800` (`1.229%`) | **`66.1%`** (fails specificity) |
| **`RULE_V2`** (`30%` id, $\ge 50$ cols, $\ge 80\%$ true ungapped coverage on **both** sequences) | `11 / 4,800` (`0.229%`) | **`0 / 4,800` (`0.000%`)** | **`1 / 4,800` (`0.021%`)** | **`0.0%`** (passes $\le 0.5\%$ gate) |

Under `RULE_V2`, true ungapped residue coverage is strictly bounded in $[0, 1]$ on both sequences:

$$c_A = \frac{\sum_{i} \mathbb{I}[q_i \ne \text{'-'} ]}{|A|} \in [0, 1], \qquad c_B = \frac{\sum_{i} \mathbb{I}[r_i \ne \text{'-'} ]}{|B|} \in [0, 1], \qquad \min(c_A, c_B) \ge 0.80$$

While local Smith–Waterman alignments are suppressed on shared sub-domains of multi-domain proteins of unequal length, domain-level homology is enforced explicitly via transitive connected components over `Pfam`, `Gene3D`, and `UniRef50` identifiers (`family_components`).

### 1.2 Defect 2: DisProt `unpublished: true` Legacy Flag Eliminating All Verified Order Negatives

In `data/raw/disprot/disprot_2026_06_regions.json` (SHA-256 `1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140`), every region carries:
- `validated`: curator sign-off dictionary (`curator_id`, `curator_orcid`, `timestamp`),
- `released`: official public DisProt release string (`2022_03` through `2026_06`, never `'unpublished'`),
- `reference_source` / `reference_id`: peer-reviewed citation (`pmid` or `doi`),
- `obsolete`: boolean obsolescence marker.

However, **65 of 66 (`98.5%`) curator-validated `IDPO:0000006` (experimental structure / order) regions** and **175 of 240 (`72.9%`) curator-validated `IDPO:0000033` (flexible linker) regions** retain `"unpublished": true` from internal curation tracking even after public release in `DisProt 2022_03`–`2026_06` with PubMed citations. Filtering out `unpublished=True` in `scripts/build_dfl_curated_evidence_pool.py` eliminated **100% of surviving `IDPO:0000006` negative residues**.

### 1.3 Verified V2 Family-Disjoint Split Summary (`reports/dfl_b1_v2_family_split_verification_2026-10-18.json`)

Running `scripts/build_dfl_b1_v2_split.py` and `scripts/preflight_dfl_family_training_v2.py` produces a hash-locked V2 split (`manifest_sha256`: `cbb0358c096c9c870b1e37528ad0d41fb261a5eda470c0b73c0ba9090c370050`, `family_manifest_sha256`: `db484272985535f08cf182166503973e77753abd00c7a4df1a197cc9db971dd2`) with status `pass_v2_preregistered_sequence_and_family_firewall`:

| Split Role | Proteins | Family Components | Positive (`1`, `IDPO:0000033`) Residues | Negative (`0`, `IDPO:0000006`) Residues | Unknown (`x`) Residues |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`train`** | `60` | `56` | `1,373` | `541` | `24,372` |
| **`development`** | `28` | `26` | `547` | `518` | `13,713` |
| **`auxiliary`** (`order_calibration_panel`) | `9` | `9` | `0` | `194` | `3,347` |
| **Total Surviving (`V2`)** | **`97`** | **`91`** | **`1,920`** | **`1,253`** | **`41,432`** |

- **Quarantined Benchmark Overlaps:** `145` DisProt proteins (`benchmark_id_or_family_component`) were quarantined prior to split allocation because of exact ID or transitive `Pfam`/`Gene3D`/`UniRef50` family-component overlap with `CAID3` (`2a754cada707eb5a748d9302e3425d06ee6f5a60b19baa853432e2d044d16e70`) or `DP94` (`c19834e1e5ca510845e631d1cc0ea3112d247ed84d8cbe077edbf38378c1f90e`).
- **Complete `B1` Gate, Autoresearch Status, and Launcher Integration:**
  - `scripts/preflight_caid3_protected_training.py` verifies `data/processed/dfl_b1_v2_family_split/manifest.json` against both embedded registries and the root `artifacts/caid3-linker.fasta` (`4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`) and `data/processed/prospective_dfl_holdout_dp94_v1/sequences.fasta` (`f06a8ea2a91c4a23c1c7177b10a48d79ab8a10f247d9d36cc4bb5692b79bf4cc`), returning `status: pass_v2_preregistered_sequence_and_family_firewall` while continuing to fail closed on `data/processed/caid3_protected_dfl_v1/manifest.json`.
  - `scripts/check_detector_autoresearch.py` (via `configs/detector_autoresearch_active.json` -> `reports/dfl_b1_v2_family_round_2026-10-18.json`) passes with `action: launch_or_audit` and explicitly prohibits `caid3_protected_dfl_v1`.
  - `scripts/run_caid3_protected_dfl_trial.py --preflight-only` (configured via `configs/dfl_b1_v2_family_round_2026-10-18.json`) passes all hash and firewall checks.
- **Downstream Scope for Paper II (`B2` / `B3`):** With `B1` closed, Paper II remains gated on `B2` (conformal risk control calibrated on a family-disjoint subset) and `B3` (an external single-domain crystallographic no-linker panel for proteome-scale specificity claims).

---

## Part 2: Experiment `X1` — Do Relative-Position 3D Structural Features Transfer Across Held-Out Domain Pairs?

### 2.1 Experimental Protocol

To test `CHECKLIST.md` item `X1` without using any unverified figures from Li et al. (2026), `scripts/evaluate_x1_relative_position_transfer.py` evaluates **91 coordinate-complete CATH v4.4.0 two-domain natural proteins across 57 CATH domain-pair groups** (`data/processed/cath_full_pilot.jsonl`, SHA-256 `20deca0e0bff4d1c68b245cc4e15756fc4e6462537c67a7d4508b5f8f58cfd85`, verified against cached RCSB PDB files by SHA-256).

- **Split Hygiene:** For each of 5 seeds (`7, 17, 31, 43, 97`), `grouped_split` partitions records by `homology_cluster_id` and CATH superfamily so **zero domain-pair architectures cross between `train` and held-out (`validation` + `blind`) splits** (`27` held-out proteins per seed; `135` pooled held-out evaluations).
- **SE(3)-Invariant 3D Relative-Position Features:** Extracted directly from verified C$\alpha$ coordinates (`extract_relative_position_geometry`, verified SE(3)-invariant in `tests/test_evaluate_x1_relative_position_transfer.py`):
  1. Inter-port C$\alpha$ Euclidean distance $d_{AB} = \|\mathbf{x}_{\text{port},A} - \mathbf{x}_{\text{port},B}\|_2$ ($\text{\AA}$),
  2. Physical span slack fraction relative to the $3.8\,\text{\AA}/\text{residue}$ maximum C$\alpha$ extension bound (`SPEC.md` §2.4) and deviation from $d_{AB}/2.5\,\text{\AA}$,
  3. Inter-domain centroid distance normalized by the sum of domain radii of gyration $R_g(A) + R_g(B)$,
  4. Relative port-to-partner-centroid distance asymmetry,
  5. Terminal C$\alpha$ exit/entry unit-vector cosine alignment with the inter-port axis $\hat{\mathbf{u}}_{AB}$.

### 2.2 Held-Out Domain-Pair Ranking Results (`reports/x1_relative_position_transfer_2026-10-18.json`)

| Ranker / Baseline | Features Used | Mean MRR (5 Seeds) | Min–Max Seed MRR | Pooled Top-1 (`Raw Count`) | Pooled Top-5 (`Raw Count`) |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **`geometry_span_match` (non-learned physical rule)** | **$3.8\,\text{\AA}/\text{res}$ feasibility + $|L - d_{AB}/2.5\,\text{\AA}|$** | **`0.1710`** | `0.1347`–`0.2075` | **`5.19%` (`7/135`)** | **`27.41%` (`37/135`)** |
| **`length_heuristic` (non-learned)** | Training median length match | `0.1160` | `0.0888`–`0.1379` | `4.44%` (`6/135`) | `9.63%` (`13/135`) |
| **`fixed_gs` (non-learned)** | Gly/Ser content + length prior | `0.1066` | `0.0875`–`0.1193` | `0.00%` (`0/135`) | `14.07%` (`19/135`) |
| **`sequence_only` (learned linear ranker)** | `6` sequence composition/length/entropy features | `0.0709` | `0.0586`–`0.0814` | `0.00%` (`0/135`) | `3.70%` (`5/135`) |
| **`sequence_plus_local_port` (learned linear ranker)** | `12` sequence + within-domain port accessibility/contact features | `0.0697` | `0.0584`–`0.0765` | `0.00%` (`0/135`) | `4.44%` (`6/135`) |
| **`sequence_plus_relative_position_3d` (learned linear ranker)** | `18` sequence + local port + 3D inter-domain relative-position features | `0.0684` | `0.0626`–`0.0742` | `0.00%` (`0/135`) | `2.96%` (`4/135`) |
| **`nearest_training` (non-learned)** | Nearest training port-feature neighbor | `0.0319` | `0.0316`–`0.0322` | `0.00%` (`0/135`) | `0.00%` (`0/135`) |

### 2.3 Paired Cluster-Bootstrap Comparisons (500 Draws by `homology_cluster_id`)

- **`sequence_plus_relative_position_3d` vs. `sequence_only`:** $\Delta\text{MRR} = -0.0025$ (95% CI $[-0.0173, +0.0058]$) — **null/negative transfer**.
- **`sequence_plus_local_port` vs. `sequence_only`:** $\Delta\text{MRR} = -0.0012$ (95% CI $[-0.0044, +0.0007]$) — **null/negative transfer**.
- **`sequence_plus_relative_position_3d` vs. `sequence_plus_local_port`:** $\Delta\text{MRR} = -0.0013$ (95% CI $[-0.0132, +0.0056]$) — **null/negative transfer**.

### 2.4 Scientific Verdict on `X1`

1. **Learned linear weights over 3D relative-position features do not transfer across held-out domain pairs** ($\Delta\text{MRR} = -0.0025$ vs. `sequence_only`). With 57 distinct CATH domain-pair groups, concatenating 12–18 geometric features into a learned ranker overfits to training domain architectures and underperforms even `sequence_only`.
2. **However, explicit 3D inter-port Euclidean distance ($d_{AB}$) is highly informative when applied as a zero-parameter physical constraint (`geometry_span_match`):** enforcing $3.8\,\text{\AA} \times L \ge d_{AB}$ and ranking candidates by $|L - d_{AB}/2.5\,\text{\AA}|$ more than doubles the MRR of `sequence_only` (`0.1710` vs. `0.0709`) and nearly doubles `fixed_gs` Top-5 recovery (`37/135` vs. `19/135`). This directly supports the `CURRENT_DIRECTION_2026-10-18.md` architecture: use 3D port geometry as an explicit physical feasibility filter rather than an over-parameterized learned feature vector.

---

## Reproducibility Commands

```powershell
# 1. Build and verify the V2 family-disjoint DFL split (requires <source-root> with pinned DisProt/CAID3/DP94 data)
python scripts/build_dfl_b1_v2_split.py --source-root <source-root> --output-dir data/processed/dfl_b1_v2_family_split
python scripts/preflight_dfl_family_training_v2.py data/processed/dfl_b1_v2_family_split/manifest.json

# 2. Run Experiment X1 held-out domain-pair relative-position transfer evaluation
python scripts/evaluate_x1_relative_position_transfer.py --records <source-root>/data/processed/cath_full_pilot.jsonl --rcsb-cache <source-root>/data/raw/rcsb --output reports/x1_relative_position_transfer_2026-10-18.json

# 3. Run unit tests for V2 preflight gate and X1 evaluator
python -m unittest tests/test_preflight_dfl_family_training_v2.py tests/test_evaluate_x1_relative_position_transfer.py -v
```
