# Preregistration: DFL Sequence & Domain-Family Firewall V2 (`B1` Deadlock Resolution)

- **Date:** 2026-10-18
- **Status:** Preregistered V2 benchmark-isolation and family-split contract (`schema_version: "2.0.0"`).
- **Predecessors:** `literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md`, `reports/DFL_B1_FIREWALL_REPAIR_PROPOSAL_2026-10-18.md`, `reports/DFL_B1_PAPER_SURVEY_AND_SPLIT_DESIGN_2026-10-18.md`, `reports/DFL_B1_TRAINING_START_GATE_2026-10-18.md`.
- **Integrity boundary:** This preregistration uses **benchmark IDs and sequences only**. No CAID3 or DP94 residue labels, prediction scores, or per-protein errors were opened, parsed, or inspected in designing this rule.

---

## 1. Diagnosis of the V1 Deadlock

Under the V1 preflight (`scripts/preflight_dfl_family_training.py`) and curated pool builder (`scripts/build_dfl_curated_evidence_pool.py`), `B1` reached a methodological deadlock (`24` surviving positive sequences and `0` verified negative residues in `reports/DFL_B1_TRAINING_START_GATE_2026-10-18.md`) due to two distinct software defects rather than genuine biological overlap:

### Defect 1: Unbounded Gapped-Alignment Ratio & Unmasked Composition Matching in V1 `local_hit`
The V1 pairwise screen (`RULE`) rejected any candidate pair satisfying:
```python
hit = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
hit.length >= 50 and hit.matches / hit.length >= 0.25 and hit.length / min(len(left), len(right)) >= 0.5
```
1. **Unbounded denominator bug:** `hit.length` is the Smith–Waterman alignment length *including gap columns*, not the number of aligned residues in the shorter sequence. When gaps extend a local alignment across low-complexity regions, `hit.length / min(len(left), len(right))` exceeds `1.0` and does not measure sequence coverage.
2. **False-positive saturation on composition-shuffled controls:** Because disordered proteins share strong compositional bias (Gly/Ser/Pro/Glu/Lys enrichment), `scripts/audit_dfl_firewall_null.py` demonstrated that across 4,800 sampled candidate–DP94 pairs (seeds `173`, `271`, `811`), V1 flagged **118 observed pairs**, **78 whole-residue-shuffled pairs ($66.1\%$ of observed)**, and **59 12-residue-block-shuffled pairs ($50.0\%$ of observed)**.
3. **Destruction of family-clear training capacity:** In `scripts/audit_dfl_disprot_family_metadata.py`, 134 candidate sequences in 112 connected components were completely disjoint from both `CAID3` and `DP94` across `MMseqs2` (`--mask 1 --min-seq-id 0.3 -c 0.8 --cov-mode 0`), `Pfam`, `Gene3D`, and `UniRef50`, yet V1's unmasked local alignment rejected 110 of those 134 family-clear sequences.

### Defect 2: Misinterpretation of DisProt's Internal `unpublished` Boolean on Released, PMID-Cited Regions
In `data/raw/disprot/disprot_2026_06_regions.json`, every region annotation carries:
- `validated`: a curator sign-off object (`{"curator_id": ..., "timestamp": ...}`),
- `released`: the official public DisProt release identifier (`"2022_03"` through `"2026_06"`),
- `reference_source` and `reference_id`: peer-reviewed literature citation (`"pmid"` or `"doi"`),
- `unpublished`: a legacy internal curation-batch flag that remained `True` on **65 of 66** curator-validated `IDPO:0000006` (experimental structure / order) regions and **175 of 240** curator-validated `IDPO:0000033` (flexible linker) regions even after their public release in DisProt `2022_03`–`2026_06` with PubMed IDs.

Filtering out `region.get("unpublished") == True` discarded 98.5% of curator-validated experimental-order (`IDPO:0000006`) negative regions, leaving `0` surviving negative residues and triggering the fail-closed stop in `scripts/check_dfl_b1_training_readiness.py`.

---

## 2. Preregistered V2 Sequence & Family Firewall (`RULE_V2`)

To replace V1 without relaxing benchmark protection, a split manifest under `schema_version: "2.0.0"` must pass **all four** layers below in `scripts/preflight_dfl_family_training_v2.py`:

### Layer 1: Sealed Benchmark Registries
Both benchmark sequence registries remain mandatory and hash-pinned:
- `caid3`: `status = "sequence_only_protected"`, `label_seal_state = "retrospective_spent_exclusion_only"` (31 sequences, SHA-256 `4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8`).
- `dp94`: `status = "sequence_only_labels_sealed"`, `label_seal_state = "unopened"` (98 sequences, SHA-256 `fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051`).

### Layer 2: Exact, Containment, and True Ungapped Reciprocal-Coverage Homology (`RULE_V2`)
Any training, development, or auxiliary/calibration sequence $q$ is rejected if, against any benchmark sequence $r \in \text{CAID3} \cup \text{DP94}$:
1. **Identifier match:** $\text{id}(q) = \text{id}(r)$, or any source observation of $q$ carries $\text{id}(r)$.
2. **Exact or substring containment:** $q = r$, $q \subset r$, or $r \subset q$.
3. **True ungapped reciprocal residue coverage (`RULE_V2`):** Under Smith–Waterman traceback (`BLOSUM62`, `gap_open = 10`, `gap_extend = 1`):
   - Alignment columns $L_{\text{align}} \ge 50$,
   - Exact amino-acid identity across alignment columns $\ge 0.30$,
   - **True ungapped query residue coverage** $c_q = \frac{|\{i : q_{\text{align}}[i] \neq \text{'-'}\}|}{|q|} \ge 0.80 \in [0, 1]$, AND
   - **True ungapped reference residue coverage** $c_r = \frac{|\{i : r_{\text{align}}[i] \neq \text{'-'}\}|}{|r|} \ge 0.80 \in [0, 1]$.

```json
{
  "rule_version": "2.0.0",
  "matrix": "BLOSUM62",
  "gap_open": 10,
  "gap_extend": 1,
  "minimum_identity": 0.30,
  "minimum_aligned_residues": 50,
  "minimum_each_sequence_residue_coverage": 0.80,
  "coverage_definition": "ungapped_aligned_residues_divided_by_sequence_length",
  "maximum_shuffle_null_rate": 0.005
}
```

### Layer 3: Mandatory Shuffle-Null Specificity Gate
Before any split is accepted, `RULE_V2` is evaluated on the 3-seed composition-preserving null control (`seeds = [173, 271, 811]`, $40 \times 40 = 1{,}600$ pairs per seed, $4{,}800$ total pairs per control):

| Metric | V1 (`local_hit`) | V2 (`RULE_V2`) | Preregistered V2 Gate |
| --- | ---: | ---: | --- |
| Observed candidate–DP94 hits / 4,800 | 118 (`2.458%`) | 11 (`0.229%`) | Diagnostic record |
| Whole-residue shuffle hits / 4,800 | 78 (`1.625%`) | **0 (`0.000%`)** | $\le 0.005$ (`0.50%`) |
| 12-residue block-shuffle hits / 4,800 | 59 (`1.229%`) | **1 (`0.021%`)** | $\le 0.005$ (`0.50%`) |

### Layer 4: Transitive Domain-Family & Composition-Masked MMseqs2 Component Disjunction ($\mathcal{G}_{\text{V2}}$)
Because full-length reciprocal coverage (`0.80`) alone does not detect single-domain sharing between multidomain proteins of different total lengths, `preflight_dfl_family_training_v2.py` enforces **strict component disjunction** over the union graph $\mathcal{G}_{\text{V2}}$:
- **Nodes:** All candidate proteins and all sealed benchmark proteins with verified exact-sequence metadata in `DisProt_2026_06` (plus optional MMseqs2 `--mask 1 --min-seq-id 0.3 -c 0.8 --cov-mode 0` cluster representatives when provided).
- **Edges:** Two proteins are connected if they share **any**:
  1. `Pfam` domain identifier (`features.pfam`),
  2. `Gene3D` superfamily identifier (`features.gene3D`),
  3. `UniRef50` cluster identifier (`uniref50`),
  4. `MMseqs2` composition-masked cluster representative (`cluster_tsv`, when present), or
  5. Direct `RULE_V2` sequence match.
- **Fail-closed component rules:**
  1. Any protein lacking an exact sequence match in the family metadata source is **quarantined** (`unresolved_family_metadata`) and cannot enter `train`, `development`, or `auxiliary`.
  2. Any protein whose connected component $C \in \pi_0(\mathcal{G}_{\text{V2}})$ contains **any** `CAID3` or `DP94` member is **quarantined** (`benchmark_connected_family_component`).
  3. Every remaining connected component $C$ is assigned **atomically as a whole unit** to exactly one split role (`train`, `development`, or `auxiliary`). Cross-split component overlap or cross-split `RULE_V2` sequence overlap immediately aborts preflight.

---

## 3. Preregistered Label & Negative-Supervision Semantics (`1/0/x`)

1. **Curator-validated released evidence criterion:** A region in `DisProt_2026_06` (`disprot_2026_06_regions.json`, SHA-256 `1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140`) is eligible iff:
   - `bool(region.get("validated")) is True` (signed off by a named DisProt curator),
   - `not region.get("obsolete")`,
   - `bool(region.get("released")) is True` (publicly released in a numbered DisProt release), and
   - `region.get("reference_source") in {"pmid", "doi"}` with non-empty `reference_id`.
   Each record retains `unpublished_flag_count` in its provenance metadata so strict (`unpublished=False`) vs. legacy-flagged (`unpublished=True`, `released` in `2022_03`–`2026_06`) regions remain auditable.
2. **Residue mask encoding (`1/0/x`):**
   - `1` (**DFL positive**): residue covered by at least one qualifying `IDPO:0000033` (flexible linker) region and zero qualifying `IDPO:0000006` regions.
   - `0` (**Verified experimental-order negative**): residue covered by at least one qualifying `IDPO:0000006` (experimental structure / order) region and zero qualifying `IDPO:0000033` regions.
   - `x` (**Unknown / excluded from loss**): all unannotated residues, residues annotated only with other disorder/binding functions (`IDPO:0000002`, etc.), and any residue where `IDPO:0000033` and `IDPO:0000006` overlap.
   - **Hard prohibition:** `x` is never converted to `0` during training, threshold selection, or evaluation.
3. **Role requirements across splits:**
   - `train`: allocated from whole family components; must contain both observed `1` (`IDPO:0000033`) and verified `0` (`IDPO:0000006`) residues.
   - `development`: allocated from whole family components disjoint from `train`; must contain both observed `1` and verified `0` residues.
   - `auxiliary` (`order_calibration_panel`): allocated from whole family components of `IDPO:0000006`-only proteins disjoint from `train`, `development`, `CAID3`, and `DP94`, providing an untouched family-disjoint negative/order specificity panel.
