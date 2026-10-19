# LinkerForge agent handoff

Read `README.md` completely before changing this repository. It is the canonical scientific and MVP specification.

For the **active DFL detector** research, then read `RESEARCH_HANDOFF.md` and
`literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md`. The fusion-ranker
decisions below are historical/future scope, not the active detector training
objective. The current protected DFL split fails the prospective benchmark
firewall; never bypass preflight or use public CAID3 scores to select a model.

## Current state

- The React/Vite frontend exists and uses labelled demonstration data. Its
  two-sequence template workflow exports a complete fusion FASTA and accepts a
  sequence-matched predicted PDB; it does not produce a trained recommendation.
- An optional Tauri desktop shell can launch a separately installed local
  ColabFold command and display its full-fusion PDB. ColabFold/weights are not
  bundled, and a successful fold on the 6 GB laptop has not been verified.
- CI and GitHub Pages deployment exist.
- Some detector and design-pipeline code exists; the qualified dataset, trained
  fusion ranker, completed structural funnel, and cloud jobs do not.
- Do not imply that proposed architecture or example scores are completed results.

Start with `START_HERE.md` to choose the active detector, the untrained design
pipeline, the frontend, or the desktop shell. Do not treat code-only readiness
or a synthetic smoke test as biological validation.

## Frozen decisions

- Demonstration: `sfGFP — linker — mCherry`.
- Linker: flexible, non-cleavable, GS-rich, 5–30 residues.
- Fixed control: `(GGGGS)3 = GGGGSGGGGSGGGGS`.
- Data: natural two-domain proteins for pretraining; engineered fusions for stronger evidence where available.
- Learning task: pairwise/contrastive compatibility ranking, not experimental-success classification.
- Model: frozen ESM-2 embeddings + attachment-point geometric GNN + small linker encoder + multitask heads.
- Generation: constrained enumeration/mutation/retrieval for MVP; learned generation and diffusion are future work.
- Blind test: homology-aware held-out proteins with masked natural linkers; save predictions before revealing references.
- Structural funnel: Rosetta → AlphaFold/ColabFold → hard rejection → OpenMM; optional GROMACS for 1–3 finalists.
- PyMOL is for inspection, measurements, and figures, not the simulation engine.
- Result claim: highest-ranked computational candidate, never experimental proof.

## Development rules

- Preserve user work and unrelated files.
- Use the versioned JSON contracts in README Section 7; add formal JSON Schema before producing scientific records.
- Every derived record needs source IDs, commit SHA, config hash, versions, and evidence quality.
- Label natural observations, synthetic negatives, predictions, and experimental outcomes distinctly.
- Split by homology/domain-family clusters before model development; never tune on blind examples.
- Baselines and LinkerForge must use the same scoring and rejection rules.
- Keep long scientific jobs on the college cloud, not web hosting or ordinary CI.
- Frontend example values must remain visibly labelled until replaced by provenance-bearing outputs.

## Branch contract

- `main`: reviewed release.
- `integration`: combined testing.
- `work/build-port-records`: data and connection records.
- `work/train-linker-ranker`: baselines and ML.
- `work/predict-and-score-fusions`: structural modelling and scoring.
- `work/build-results-interface`: UI and provenance presentation.

Run frontend verification before handing off UI changes:

```bash
cd frontend
npm ci
npm run lint
npm run test
npm run build
```

If no backend environment is committed, do not invent a working installation claim. Pin and test it on the target cloud image first.
