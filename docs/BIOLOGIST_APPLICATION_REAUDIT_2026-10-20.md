# LinkerForge application re-audit: biologist workflow

> **Historical workflow snapshot.** This document records the earlier
> biologist-facing gap analysis and is retained for traceability; its status
> table is not the current completion ledger. For the implemented
> PAE/reference/export/history/protected-site/OpenMM/matched-panel state and the
> publication decision, use
> `PROFESSIONAL_PUBLICATION_READINESS_AUDIT_2026-10-20.md`.

Audit date: 20 October 2026 (Asia/Calcutta)
Scope: local React/Tauri workflow, deterministic linker templates, ColabFold handoff,
3Dmol.js/PyMOL visualization, result interpretation, and research traceability.

## Decision

The application is now a clearer **sequence-template and local structure-prediction
prototype**, but it is not yet a validated linker-recommendation product. A biologist
can document a construct, provide two canonical domain sequences, compare auditable
linker templates and the two immediately available controls, select the exact complete
fusion, run or import a sequence-matched structure prediction, inspect it, and see an
honest result record. The interface does not claim that a template is biologically
successful.

Professional recommendation remains blocked by the absent reviewed fusion dataset,
validated linker ranker, complete matched structural baseline suite, hard structural
rejection pipeline, cloud handoff, MD workflow, and wet-lab evidence.

## Why the Candidates page exists before Structure

ColabFold does not predict an abstract linker. It predicts one exact complete sequence:

`Domain A + selected linker + Domain B`

The Candidates page is therefore a pre-compute decision gate. It must let a researcher:

1. verify the biological context and domain order;
2. distinguish a design template from a control;
3. inspect the exact linker and complete fusion sequence;
4. understand why the candidate exists and what evidence supports it;
5. assess local-compute risk before spending GPU time; and
6. record which sequence advanced to prediction.

The page must not present deterministic templates as trained recommendations. A
professional comparison eventually needs every shortlisted design and baseline folded
and scored under identical settings.

## Biologist-facing requirements

| Need | Current application state | Audit result |
| --- | --- | --- |
| Construct/project identity | Project name and domain labels are captured and exportable as project.v1 JSON | Implemented; reference coordinates are deliberately separate |
| Biological objective | Free-text objective is captured and shown with results | Implemented in-session |
| Expression context | Intended expression system is captured | Implemented in-session |
| Functional/protected sites | Domain-relative residues are validated and mapped onto fusion numbering | Partial; no accessibility/obstruction scoring |
| Exact sequences and order | FASTA accepted; canonical-residue and length checks applied | Implemented |
| Domain-boundary evidence | Researcher is warned to document/verify it | Missing formal source fields |
| Tags/signal peptides/PTMs/cofactors | Input preflight explains required decisions and limitations | Guidance only |
| Linker class and length | Flexible/rigid, cleavable/non-cleavable and 5–30 aa templates | Implemented as deterministic templates |
| Candidate provenance | Template/control role, rationale and evidence quality shown | Implemented |
| Required controls | Direct fusion and fixed `(GGGGS)3` included | Partial baseline coverage |
| Learned recommendation | No approved fusion ranker exists | Correctly unavailable |
| Sequence privacy | Public ColabFold/MMseqs2 disclosure gates local prediction | Implemented per session |
| Compute suitability | Complete length and local 4 GB GPU risk are shown | Coarse preflight only |
| Exact prediction input | Complete fusion sequence can be reviewed/exported | Implemented |
| Prediction progress | Start/running/complete/failed state, elapsed time and run path | Implemented |
| Structural confidence | Mean pLDDT for whole fusion and each sequence region | Implemented from PDB B-factors |
| Interdomain uncertainty | Matching ColabFold score JSON is ingested; pTM and mean A↔B PAE are reported | Implemented; acceptance threshold not frozen |
| Closure/clashes | No Rosetta/scoring adapter is connected | Missing |
| Domain preservation | Each input domain accepts its own sequence-matched reference PDB; C-alpha RMSD is reported | Implemented; needs an additional metric/frozen thresholds |
| Protected-site clearance | Sites map to full-fusion numbering | Partial; no solvent/interface clearance |
| Visualization | Browser 3Dmol.js and sequence-checked PyMOL handoff | Implemented; PyMOL licence unresolved |
| Simulation | No OpenMM/GROMACS workflow exists | Missing; PyMOL is not simulation |
| Experimental evidence | No expression/activity/stability measurements | Not performed |
| Research export | FASTA, PDB, PyMOL script and versioned result JSON export exist | Partial; no one-click checksummed ZIP bundle |
| Reopen/reproduce | Project JSON save/open and timestamp/hash run directories exist | Partial; no signed run index/history UI |

## Changes made in this re-audit

- Reframed Candidates as a decision gate before GPU prediction.
- Added project name, domain names, biological objective, expression system, and
  protected-site documentation.
- Added a biologist input-preflight checklist.
- Added explicit public-MSA disclosure and consent before local ColabFold can start.
- Added direct-fusion and fixed `(GGGGS)3` controls to custom candidate sets.
- Added candidate role, rationale, evidence classification, exact full sequence,
  domain/linker lengths, and local compute-risk messaging.
- Replaced the misleading provenance-only Results page with the selected construct,
  prediction state, structure source, per-region pLDDT, run identifiers, claim boundary,
  and a required-evidence checklist.
- Added PAE/pTM ingestion, per-input sequence-matched domain references, C-alpha
  preservation RMSD, protected-site numbering, coordinate pre-screening, versioned
  project/result exports and formal schemas. Rosetta, matched structural controls,
  minimization/MD and experiments remain explicitly unavailable.
- Added direct-fusion support to PyMOL domain colouring.

## Remaining priority work

### P0: required before calling this a recommendation system

1. Build and review a leakage-resistant fusion-linker dataset and validated ranker.
2. Add a checksummed immutable run index/history UI, complete bundle archive and CI
   schema validation.
3. Add cancellation, stale-job recovery and complete tool/model/MSA provenance.
4. Calculate frozen closure, clash, domain-preservation and
   interdomain-uncertainty rejection rules.
5. Run all required baselines under the same prediction/scoring contract.

### P1: required for routine biological interpretation

1. Add UniProt/PDB/domain-boundary source identifiers and residue-number mapping.
2. Add domain reference structures and per-domain RMSD/TM-score comparison.
3. Add construct annotations for tags, localization signals, cleavage motifs, cofactors,
   disulfides, oligomeric state and known interfaces.
4. Add an exportable research bundle containing JSON, FASTA, PDB/mmCIF, PAE, plots,
   logs, manifests, citations, hashes and limitations.
5. Add an institutional/cloud execution handoff for sequences unsuitable for the local
   4 GB GPU.

### P2: later structural funnel

Add Rosetta linker sampling, OpenMM minimization, and—only for reviewed finalists—a
validated replicated MD protocol. Keep PyMOL limited to inspection and figures.

## Permitted interpretation

The current result may say:

> This is a sequence-matched ColabFold structure prediction for the selected deterministic
> template or control, with the displayed local-confidence measurements.

It must not say:

> This linker is predicted to work biologically or is experimentally validated.
