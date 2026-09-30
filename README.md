# LinkerForge

LinkerForge is a **computational research prototype** for recommending a peptide linker between two independently folding soluble protein domains. It learns compatibility patterns from naturally occurring multidomain proteins, ranks constrained linker candidates, and evaluates complete fusion-protein predictions through a reproducible structural-validation pipeline.

This README is the canonical, self-contained project specification. A new team member, coding agent, cloud notebook, or clean laptop should be able to recover the project scope and continue without the original planning conversation.

> **Current state:** the React demonstration interface, CI/CD, and an embedded Mol* reference-structure viewer are implemented. The viewer displays experimental PDB entries and is not a prediction service. The scientific Python pipeline, dataset, trained model, and cloud jobs below are the MVP specification and must not be represented as already implemented.

## One-minute description

LinkerForge learns from natural proteins of the form `Domain A — natural linker — Domain B`. It represents the domains, attachment-point geometry, and candidate linker sequence, then learns to rank the observed natural linker above corrupted or transplanted alternatives. A frozen model is tested by hiding linkers from completely unseen protein families. Finally, it recommends a flexible linker for an `sfGFP — linker — mCherry` demonstration, predicts complete fused structures with AlphaFold/ColabFold, visualises them in PyMOL, and rejects candidates with closure failures, clashes, domain deformation, protected-site obstruction, or unacceptable uncertainty.

The output is a **highest-ranked computational candidate**, not proof that a protein will express or function in a laboratory.

---

## 1. Frozen scientific scope

### Initial target

```text
sfGFP — designed peptide linker — mCherry
```

- The two protein domains remain fixed; LinkerForge designs only the linker.
- Initial orientation: C-terminus of sfGFP to N-terminus of mCherry.
- Linker class: flexible, non-cleavable, primarily glycine/serine-rich.
- Length: 5–30 amino acids.
- Initially permitted auxiliary residues: `A`, `T`, `N`, `Q`, `E`, and `K`.
- Reject configured cleavage motifs and long hydrophobic stretches.
- Standard control: `(GGGGS)3 = GGGGSGGGGSGGGGS`.
- The demonstration is computational; fluorescence testing is future work.

### Required baselines

Every final experiment compares LinkerForge with:

1. Direct fusion with no linker.
2. Fixed `(GGGGS)3`.
3. Random length-matched GS-rich linkers.
4. Linker length selected only from terminal geometry.
5. Nearest-neighbour retrieval from training data.
6. A feature-based XGBoost or equivalent non-neural ranker.

### MVP non-goals

- No wet-lab work or experimental-success claim.
- No membrane-protein or protease-cleavable-linker design.
- No unrestricted whole-protein generation.
- No custom molecular-dynamics engine.
- No large diffusion model trained from scratch.
- No assumption that low linker pLDDT alone means failure.

---

## 2. Terminology

- **Natural multidomain protein:** one continuous chain containing two or more domains.
- **Interdomain linker:** residues between two same-chain domain assignments.
- **Engineered fusion:** a deliberately constructed sequence joining originally separate domains.
- **Candidate:** one complete `Domain A + linker + Domain B` sequence.
- **Natural compatibility example:** an observed natural combination; not automatically an optimal engineered linker.
- **Synthetic negative:** a deliberately corrupted or transplanted candidate; not an experimentally demonstrated failure.
- **Blind masked-linker reconstruction:** evaluation in which a held-out natural linker is hidden until predictions are frozen.
- **Hard rejection:** a critical failure that cannot be compensated by a high average score.

Do not call two domains in the same chain a protein complex. Reserve *complex* for associated chains unless a source uses another definition.

---

## 3. Research gap and evidence boundary

The defensible gap is a synthesis, not a claim that one paper requested LinkerForge verbatim:

> Linker selection remains context-dependent and often empirical; engineered fusion structures and interdomain orientations remain difficult to predict; and existing tools do not provide one explainable workflow that recommends a linker from two-domain structural context and evaluates the complete fusion under a frozen, reproducible contract.

Core sources:

- Chen, Zaro, and Shen, *Fusion Protein Linkers: Property, Design and Functionality*. [Open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC3726540/) · [DOI](https://doi.org/10.1016/j.addr.2012.09.039)
- Sapsford et al., *Biomolecular Engineering for Nanobio/Bionanotechnology*. [Article](https://link.springer.com/article/10.1186/s40580-017-0103-4)
- Kumar and Kim, *Artificial Intelligence in Fusion Protein Three-Dimensional Structure Prediction*. [Open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294035/) · [DOI](https://doi.org/10.1002/ctm2.1789)
- *Recent Progress of Protein Tertiary Structure Prediction*. [Open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC10893003/)
- *Fusion of Two Unrelated Protein Domains in a Chimera Protein and Its 3D Prediction*. [Open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9796088/) · [DOI](https://doi.org/10.1002/prot.26398)

Literature supports the problem. The hypothesis that the proposed representation and ranker improve linker selection must be tested by this project.

---

## 4. MVP architecture

```text
Natural structures and annotations
        ↓
Domain/linker extraction and evidence grading
        ↓
Feature and connection-graph construction
        ↓
Frozen train/validation/blind splits
        ↓
Natural positives + labelled synthetic negatives
        ↓
Structure-conditioned compatibility ranker
        ↓
Constrained candidate generation and shortlist
        ↓
Rosetta closure/conformation sampling
        ↓
AlphaFold/ColabFold complete-fusion prediction
        ↓
Hard structural rejection and scorecard
        ↓
OpenMM minimisation; optional MD for 1–3 finalists
        ↓
React results and provenance interface
```

### Model decision

The MVP is a multimodal compatibility ranker, not a diffusion model:

1. **Domain sequence encoder:** frozen ESM-2 embeddings for Domains A and B.
2. **Connection encoder:** 3D graph neural network over residues surrounding both attachment points.
3. **Linker encoder:** small Transformer, BiLSTM, or 1D CNN over the candidate linker.
4. **Fusion network:** combines those representations.
5. **Prediction heads:** compatibility rank, class/length suitability, clash risk, domain-deformation risk, protected-site obstruction risk, and confidence.

Training objective:

```text
score(observed compatible linker) > score(corrupted or transplanted linker)
```

Use pairwise or contrastive ranking. Do not interpret output as an experimentally calibrated probability of success.

### Candidate-generation decision

The MVP uses constrained, auditable generation:

- enumerate standard `(GGGGS)n` variants;
- mutate GS-rich sequences within the permitted alphabet;
- retrieve structurally compatible natural linkers;
- vary length within 5–30 residues;
- remove candidates violating deterministic rules.

A learned autoregressive generator is optional after the ranker works. Diffusion is future scope.

---

## 5. End-to-end workflow

### Phase A — Prove structural modelling

1. Select one high-confidence natural two-domain development protein.
2. Verify domain boundaries and linker residues from source annotations.
3. Inspect it in PyMOL; colour both domains, linker, and protected residues separately.
4. Remove linker coordinates while preserving its known sequence.
5. Reconstruct 100–500 conformations using Rosetta loop modelling/KIC.
6. Verify closure, geometry, clash removal, and recovery of plausible conformations.

**Gate A:** do not scale collection or claim generative ability until known-linker reconstruction works reproducibly.

### Phase B — Build the natural-linker dataset

1. Collect experimentally resolved same-chain two-domain PDB structures.
2. Obtain UniProt sequences and InterPro/Pfam/CATH domain assignments.
3. Retain soluble, sufficiently complete, unambiguous examples.
4. Extract an intervening segment only after reconciling domain annotations.
5. Calculate secondary structure, solvent exposure, contacts, flexibility indicators, composition, and attachment geometry.
6. Grade evidence as `high`, `medium`, `uncertain`, or `rejected`.
7. Keep predicted structures as a separately labelled supplementary tier.
8. Record source identifiers, retrieval date, transformations, and software versions.

**Gate B:** every training example must be traceable and schema-valid.

### Phase C — Freeze leakage-resistant splits

1. Deduplicate PDB entries and alternative structures of the same protein.
2. Cluster by sequence similarity, domain family, domain-pair architecture, and linker similarity.
3. Assign whole clusters to training, validation, or blind test sets.
4. Publish the split manifest before final training.
5. Never tune features, thresholds, or model choice on the blind set.

Start near 70% training, 15% validation, and 15% blind testing, but cluster integrity is more important than exact percentages.

### Phase D — Create positives, negatives, and baselines

Keep each observed natural combination as a natural compatibility example. Generate labelled synthetic negatives by shortening, lengthening, shuffling, reversing, transplanting, hydrophobising, or otherwise corrupting its linker. Where feasible, include candidates that cause closure failures, clashes, or protected-site obstruction. Implement all required baselines before the final neural model.

### Phase E — Train and freeze the ranker

1. Precompute frozen ESM-2 embeddings.
2. Train the XGBoost baseline.
3. Train the connection GNN, linker encoder, and prediction heads.
4. Use grouped batches so comparisons occur within a domain pair.
5. Tune only on validation data.
6. Save configuration, checkpoint, seeds, dependency lock, and data hashes.
7. Freeze the model and hard-rejection thresholds before blind evaluation.

### Phase F — Blind masked-linker reconstruction

For every blind example:

1. Hide the natural sequence and, in the hardest track, its length.
2. Supply only the domains, allowed attachment points, and permitted context.
3. Save top-1, top-5, and top-10 predictions before revealing the reference.
4. Compare class, length, identity, edit distance, physicochemical properties, and natural-linker rank.
5. Model complete structures for shortlisted predictions.
6. Compare preservation, orientation, closure, clashes, and protected-site clearance.

Exact sequence recovery is secondary because several sequences may satisfy the same constraints.

### Phase G — Generate the sfGFP–mCherry demonstration

1. Fix the sfGFP C-terminus and mCherry N-terminus as initial ports.
2. Generate constrained 5–30-residue candidates.
3. Score candidates with the frozen ranker.
4. Retain the top 10–20.
5. Build complete FASTA sequences.
6. Model linker conformations with Rosetta where appropriate.
7. Predict complete fusion structures with AlphaFold/ColabFold.
8. Compare LinkerForge with every required baseline.
9. Energy-minimise finalists with OpenMM.
10. Run replicated OpenMM or GROMACS simulations only for the best 1–3 candidates if time and compute permit.

### Phase H — Report without overclaiming

Return the best passing candidate plus alternatives, scores, structures, rejection reasons, configuration, model version, and provenance.

Permitted conclusion:

> This is the highest-ranked computational linker under the predefined criteria.

Prohibited conclusion:

> This linker is experimentally proven to work.

---

## 6. Validation contract

### Primary endpoint

Percentage of unseen domain pairs for which at least one top-five LinkerForge candidate passes every frozen hard structural requirement.

### Secondary metrics

- Natural-linker top-1/top-5 rank and mean reciprocal rank.
- Linker-class accuracy.
- Length mean/median absolute error and percentage within ±2/±5 residues.
- Exact match, identity, and edit distance.
- Physicochemical-property similarity.
- Domain A/B structural preservation.
- Severe-clash count and chain-closure success.
- Protected-site clearance.
- Domain-specific pLDDT and interdomain PAE.
- Runtime and failure rate.

### Hard rejection rules

A candidate fails if any frozen rule is violated:

- chain cannot be closed;
- severe atomic clashes remain;
- either domain is substantially deformed;
- a protected functional region is obstructed;
- linker length or deterministic sequence constraints are violated;
- prediction or scoring output is invalid;
- interdomain uncertainty exceeds the declared threshold for the claimed analysis.

Choose thresholds on development/validation data, store them in configuration, and never change them silently after seeing blind results.

### AlphaFold interpretation

AlphaFold/ColabFold is a structure predictor and visualisation input, not experimental validation. Low linker pLDDT may reflect flexibility. Evaluate domain confidence, interdomain PAE, clashes, preservation, and geometry separately. If a blind reference predates or may appear in a predictor's training corpus, do not call that predictor an independent blind validator.

---

## 7. Data contracts

All components communicate through versioned JSON. These are minimum fields; extensions must be backward compatible or increment the schema.

### Natural example

```json
{
  "schema_version": "0.1.0",
  "example_id": "pdbid_chain_domainA_domainB",
  "source": {
    "pdb_id": "XXXX",
    "chain_id": "A",
    "uniprot_id": "P00000",
    "retrieved_at": "YYYY-MM-DD",
    "structure_tier": "experimental"
  },
  "domain_a": {"start": 1, "end": 120, "sequence": "...", "family_ids": []},
  "linker": {
    "start": 121,
    "end": 132,
    "sequence": "GSGTGGSSGSGS",
    "length": 12,
    "class": "flexible",
    "evidence_quality": "high"
  },
  "domain_b": {"start": 133, "end": 260, "sequence": "...", "family_ids": []},
  "connection_features": {},
  "provenance": {
    "pipeline_commit": "git-sha",
    "config_hash": "sha256",
    "software": {}
  }
}
```

### Candidate recommendation

```json
{
  "schema_version": "0.1.0",
  "run_id": "LF-YYYYMMDD-0001",
  "input_example_id": "...",
  "candidate_id": "cand-0001",
  "linker_sequence": "GGGGSGGGGSGGGGS",
  "linker_length": 15,
  "generation_method": "fixed_baseline",
  "scores": {
    "compatibility": null,
    "clash_risk": null,
    "domain_deformation_risk": null,
    "protected_site_obstruction_risk": null,
    "confidence": null
  },
  "rank": null,
  "model": {
    "name": "linkerforge-ranker",
    "version": "0.1.0",
    "checkpoint_hash": "sha256"
  },
  "provenance": {}
}
```

### Final result

```json
{
  "schema_version": "0.1.0",
  "run_id": "LF-YYYYMMDD-0001",
  "candidate_id": "cand-0001",
  "status": "pass",
  "hard_rejections": [],
  "structure_prediction": {
    "tool": "colabfold",
    "tool_version": "record-exact-version",
    "domain_a_plddt": null,
    "linker_plddt": null,
    "domain_b_plddt": null,
    "interdomain_pae": null,
    "artifact_paths": []
  },
  "structural_metrics": {
    "domain_a_rmsd": null,
    "domain_b_rmsd": null,
    "severe_clash_count": null,
    "chain_closed": null,
    "protected_sites_clear": null
  },
  "limitations": [],
  "provenance": {}
}
```

Never use fabricated values in scientific output. Frontend demonstration values must remain visibly labelled as examples.

---

## 8. Software stack

| Purpose | Preferred tool |
| --- | --- |
| Inspection, measurements, and figures | PyMOL |
| PDB and sequence processing | Biopython, PDBFixer |
| Domain annotation | InterPro, Pfam, CATH |
| Secondary structure | DSSP |
| Solvent accessibility | FreeSASA |
| Linker reconstruction/sampling | Rosetta/PyRosetta |
| Protein embeddings | ESM-2 |
| Machine learning | PyTorch |
| Geometric graph model | PyTorch Geometric |
| Classical baseline | scikit-learn/XGBoost |
| Complete-fusion prediction | AlphaFold/ColabFold |
| Minimisation/lightweight MD | OpenMM |
| Optional high-performance MD | GROMACS |
| Trajectory analysis | MDAnalysis |
| Interface | React, TypeScript, Vite |
| Long-running execution | College cloud/HPC |

Build custom orchestration, data, scoring, and provenance around validated tools. Do not build a custom force field or MD engine.

---

## 9. Proposed repository layout

Only `frontend/` exists today. Add the remaining directories as work begins.

```text
LinkerForge/
├── README.md
├── AGENTS.md
├── CONTRIBUTING.md
├── frontend/                 # React demonstration and experimental-reference Mol* viewer
├── schemas/                  # versioned JSON Schema files
├── configs/                  # immutable experiment configurations
├── src/linkerforge/
│   ├── data/                 # ingestion, boundaries, evidence grading
│   ├── features/             # sequence and structural features
│   ├── negatives/            # synthetic-negative generation
│   ├── models/               # baselines and neural ranker
│   ├── generation/           # constrained candidates
│   ├── structure/            # Rosetta/ColabFold/OpenMM adapters
│   ├── scoring/              # metrics and hard rejection
│   └── provenance/           # manifests, hashes, versions
├── tests/
├── scripts/                  # thin reproducible entry points
├── data/
│   ├── raw/                  # ignored; never hand-edited
│   ├── interim/              # ignored
│   ├── processed/            # ignored or external store
│   └── manifests/            # small tracked manifests
├── artifacts/                # ignored model/structure outputs
└── reports/                  # tracked summaries, not huge binaries
```

Large datasets, checkpoints, AlphaFold databases, trajectories, and structures belong in the college cloud or an artefact store. Track hashes and manifests in Git.

---

## 10. Plug-and-play setup

### Clone

```bash
git clone https://github.com/richennacht/LinkerForge.git
cd LinkerForge
```

The repository name and project branding are both `LinkerForge`.

### Run the implemented frontend

Requirements: Node.js 22 and npm.

```bash
cd frontend
npm ci
npm run lint
npm run test
npm run build
npm run dev
```

The current UI displays labelled demonstration data and is not connected to a trained scientific service.

### Scientific environment to create

The backend environment is not committed yet. The first backend change must add a reproducible `pyproject.toml` plus lock file, or an equivalent Conda definition, and test it on the college-cloud image.

Expected Python packages include:

```text
numpy pandas scipy biopython torch torch-geometric
scikit-learn xgboost mdanalysis openmm
FreeSASA bindings or CLI integration
```

Rosetta/PyRosetta, DSSP, PyMOL, ColabFold/AlphaFold, and GROMACS may require separate licences, containers, databases, or HPC modules. Document the institutional installation instead of claiming a universal one-command setup.

### Cloud execution contract

Long jobs run on the college cloud, not GitHub Pages or Vercel. Every job receives an immutable config, input-manifest hash, Git SHA, seeds, environment identifier, resource request, output directory, and completion marker. Every job writes logs, machine-readable status, software versions, checksums, and a failure reason if incomplete. Never commit credentials.

---

## 11. Twelve-week MVP plan

| Week | Work | Exit evidence |
| --- | --- | --- |
| 1 | Pin environment; choose development proteins; define schemas | Reproducible setup and verified examples |
| 2 | Build extraction and PyMOL inspection flow | Validated natural-example records |
| 3 | Reconstruct known linkers with Rosetta | Closure/reconstruction report |
| 4 | Collect, grade, and deduplicate examples | Source manifest and dataset report |
| 5 | Calculate sequence, geometry, and structural features | Schema-valid feature records |
| 6 | Freeze splits; create negatives; implement baselines | Leakage audit and baseline report |
| 7 | Build frozen-ESM, GNN, and linker encoders | Training smoke test and unit tests |
| 8 | Train, tune, calibrate rules, and freeze model | Checkpoint, config, seeds, validation report |
| 9 | Run blind reconstruction once | Immutable blind predictions and metrics |
| 10 | Run Rosetta/ColabFold scorecards | Structures and rejection report |
| 11 | Run sfGFP–mCherry comparison and finalist minimisation/MD | Final comparative demonstration |
| 12 | Reproduce cleanly; finish UI and paper artefacts | Reproduction log and evidence package |

---

## 12. Atomic workstreams and branches

Branch names describe work, not people. Any member may contribute, but each output has one owner at a time.

| Branch | Responsibility | Fixed output |
| --- | --- | --- |
| `work/build-port-records` | ingestion, boundaries, grading, deduplication, splits, provenance | schema-valid natural-example and split manifests |
| `work/train-linker-ranker` | baselines, embeddings, models, training, evaluation | recommendation records, checkpoint, evaluation report |
| `work/predict-and-score-fusions` | Rosetta, ColabFold, PyMOL metrics, OpenMM/GROMACS | final-result records and structural artefacts |
| `work/build-results-interface` | inputs, run status, comparisons, uncertainty, provenance | tested UI consuming shared schemas |

All work enters `integration` through pull requests. Reviewed releases move from `integration` to `main`.

---

## 13. CI/CD and deployment

- GitHub Actions runs frontend lint, tests, and build on `main`, `integration`, and `work/**`.
- GitHub Pages deploys the frontend from `main`.
- Vercel may provide previews after connection; it is not the scientific compute backend.
- Add backend unit, schema, fixture, and deterministic smoke tests as Python is introduced.
- Do not run production AlphaFold, Rosetta sampling, or MD in ordinary pull-request CI.

---

## 14. MVP definition of done

The MVP is complete only when:

- a versioned natural-linker dataset regenerates from a source manifest;
- every example has evidence quality and provenance;
- homology-aware splits are frozen and audited;
- all required baselines use the same evaluation contract;
- the ranker trains from a committed config and reproducible seeds;
- blind predictions are saved before references are revealed;
- top-five blind recovery is reported with raw counts;
- sfGFP–mCherry compares every baseline;
- complete fusion structures and confidence outputs are retained;
- hard rejections trace to named rules and raw measurements;
- the frontend reads real schema-valid results or clearly labelled examples;
- a clean environment reproduces the small end-to-end demonstration;
- limitations distinguish plausibility from experimental proof.

---

## 15. Immediate next actions

1. Add formal JSON Schema files for Section 7.
2. Pin the Python environment on the college-cloud image.
3. Select and document 10 development PDB chains with unambiguous two-domain annotations.
4. Implement one extraction record and inspect it in PyMOL.
5. Remove and reconstruct one known linker with Rosetta.
6. Add a tiny fixture dataset and CI test before scaling collection.

Do not begin full model training until these actions pass review.

## 16. Scientific and engineering rules

- Preserve raw data; transformations create new versioned outputs.
- Record seeds, configs, versions, hashes, and raw counts.
- Separate measured, predicted, bounded, and inferred claims.
- Never relabel synthetic negatives as experimental failures.
- Never tune on the blind set.
- Never report frontend mock data as a scientific result.
- Keep changes reviewable and one output contract per workstream.
- Do not change shared schemas without team review.
- Do not commit credentials, restricted data, or large cloud artefacts.

See [CONTRIBUTING.md](./CONTRIBUTING.md) for the contribution workflow.
