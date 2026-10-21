# Deferred structural-validation roadmap

Status: **future scope / execution on hold as of 21 October 2026**.

This roadmap applies to the fusion-design application. It does not alter the
active DFL detector objective, authorize compatibility-ranker training, relax
the protected-benchmark firewall, or replace independent curation. The five
analyses and production-length OpenMM simulation below are deliberately not
reported as completed evidence.

## Why this work is deferred

The current application can ingest sequence-matched AlphaFold2/ColabFold
coordinates and confidence records, compare reference-domain C-alpha RMSD,
screen continuity, inspect protected-site distance, report fixed-coordinate
surface burial, and run a pinned Reduce/Probe overlap audit. These are useful
measurements, but they do not yet establish reproducibility across prediction
models, reference-quality fold preservation, functional-site preservation,
secondary-structure preservation, or production-timescale molecular-dynamics
behaviour.

The following order keeps dependencies explicit:

> ensemble agreement -> lDDT/TM-score -> automatic functional-site mapping ->
> DSSP -> post-OpenMM change analysis

Each stage must consume immutable, sequence-checked inputs and write a
versioned result with hashes, tool versions, configuration, timestamps and a
claim boundary. Thresholds must be frozen on development data before any
pass/fail or comparative-performance claim.

## 1. AlphaFold2 ensemble agreement

### Why it is needed

A single rank-1 structure can hide disagreement among the five AlphaFold2
parameter models. Agreement does not prove correctness, but disagreement is an
important warning that an apparent linker pose or domain orientation is not
robust even within the AlphaFold2 model family.

### How it will be implemented

1. Require all expected ranked PDB and matching score JSON files for the exact
   full-fusion sequence; reject mixed sequences, partial panels and duplicate
   model identifiers.
2. Split every model into Domain A, linker and Domain B using the immutable
   submitted boundaries.
3. Superpose Domain A and Domain B separately, then calculate pairwise domain
   C-alpha RMSD, linker/junction displacement, interdomain centre/orientation
   variation, pLDDT, pTM and cross-domain PAE summaries.
4. Report the full pairwise matrix, median, range and outliers rather than
   selecting only the most visually appealing model.
5. Label the result `within_alphafold2_family_agreement`; it is not an
   independent-model comparison. ESMFold, Boltz/Chai or another validated model
   family would be needed for architecture-level sensitivity.

Expected output: a versioned ensemble-agreement JSON record plus CSV/figure
artifacts linked to every input hash.

## 2. Reference-domain lDDT and TM-score

### Why they are needed

The existing C-alpha RMSD is sensitive to outliers and protein length. lDDT
measures preservation of local interatomic distances without requiring a
global superposition, while TM-score provides a length-normalized global fold
similarity. Together with RMSD, they give a more interpretable view of whether
each input domain remains folded similarly after fusion.

### How they will be implemented

1. Require a sequence-matched, provenance-bearing reference structure for each
   domain; never guess a PDB from a label.
2. Create an explicit residue mapping and record unresolved residues,
   insertions, alternate locations and comparison coverage.
3. Calculate lDDT with a pinned OpenStructure implementation and TM-score with
   a pinned US-align binary for Domain A and Domain B separately.
4. Run identical calculations for all matched linker controls and all retained
   ensemble members.
5. Store raw tool output, parsed metrics, coverage and software hashes. Any
   later acceptance threshold must be selected on development examples and
   frozen before protected evaluation.

Expected output: per-domain lDDT/TM-score records and matched-control tables.
These scores measure structural similarity, not expression or activity.

Primary implementation references: [OpenStructure lDDT](https://openstructure.org/docs/2.10/mol/alg/lddt/)
and [US-align](https://zhanggroup.org/US-align/).

## 3. Automatic functional-site mapping

### Why it is needed

A fusion can retain the overall domain fold while burying an active site,
binding interface, ligand-contact residue or experimentally important
modification site. The current manually entered residue list and distance
screen cannot discover or validate those annotations automatically.

### How it will be implemented

1. Resolve user-supplied UniProt/PDB identifiers against pinned, dated RCSB PDB
   and UniProt records. Cache the exact response and its checksum; do not infer
   an identifier solely from sequence similarity.
2. Import residue-level annotations with source and evidence type, including
   catalytic/binding residues, ligand contacts, interfaces, engineered
   mutations and relevant sequence features. Preserve unknown and unmapped
   states instead of treating them as clear.
3. Map source numbering to the submitted domain sequence and then to the fusion
   coordinates with an auditable alignment. Ambiguous or low-coverage mappings
   fail closed.
4. For every mapped site, calculate linker/domain/ligand distances, contact
   changes and differential solvent-accessible surface area. If a biological
   assembly or ligand is required, identify that coordinate context explicitly.
5. Present annotations and geometric measurements separately. An annotation
   plus predicted proximity is a risk flag, not proof of functional loss.

Expected output: a provenance-bearing site-mapping JSON record, mapping table,
unmapped-residue report and spatial-risk visualization. RCSB's
[Data API](https://data.rcsb.org/) is the planned structure metadata source;
release-specific UniProt records must also be retained.

## 4. DSSP secondary-structure analysis

### Why it is needed

Linker insertion can extend or disrupt terminal helices and strands near a
junction. pLDDT and whole-domain RMSD do not directly describe those local
secondary-structure changes.

### How it will be implemented

1. Run a pinned `mkdssp` executable on each reference, predicted ensemble
   member and matched control.
2. Map DSSP assignments back to the immutable fusion boundaries and retain
   missing/unassigned residues explicitly.
3. Compare the last residues of Domain A, the linker and the first residues of
   Domain B using per-residue states, helix/strand persistence and junction
   transition counts.
4. Report ensemble consensus and disagreement rather than treating one model's
   assignment as observed biology.
5. Link the DSSP record to the same sequence, coordinate and configuration
   hashes used by the other structural analyses.

Expected output: per-residue DSSP JSON/CSV records and a junction-focused
summary. See the [DSSP documentation](https://github.com/PDB-REDO/dssp/blob/trunk/doc/mkdssp.html).

## 5. Post-OpenMM change and convergence analysis

### Why it is needed

A successful minimization or a visually stable short trajectory is not a
stability result. Replicate-aware analysis is needed to show whether the
predicted fold, domain arrangement, linker contacts and protected sites remain
consistent under the declared simulation protocol.

### How it will be implemented

1. Verify that all replicas used the same immutable starting structure, force
   field, solvent model, temperature/pressure protocol and production length,
   while retaining independent random seeds.
2. Compare the starting prediction, minimized structure and trajectory blocks
   for domain-aligned RMSD, RMSF, radius of gyration, linker end-to-end distance,
   junction geometry, interdomain orientation, contacts and protected-site
   accessibility.
3. Calculate within-replica block summaries and between-replica variation.
   Preserve time series and do not reduce the analysis to a best frame.
4. Re-run selected static checks—lDDT/TM-score, DSSP, SASA and overlap/contact
   audits—on representative, deterministically selected frames.
5. Report lack of convergence, replica disagreement and failures plainly. Even
   a well-behaved 3 x 50 ns panel is computational plausibility evidence, not
   proof of expression, activity or long-timescale stability.

Expected output: versioned convergence/change JSON, CSV time series, summary
figures and a complete manifest. The [OpenMM user guide](https://docs.openmm.org/development/userguide/)
is the implementation reference.

## 6. Production OpenMM execution — on hold

The storage-aware protocol is configured for three independent 50 ns replicas,
10 ps metric/state/checkpoint intervals and 50 ps full-system DCD intervals.
Its last successful preflight projected 1.232 GiB of trajectory data and 31.76
GPU hours, with a 6.849 GiB free-space safety requirement. That preflight is
operational readiness only.

Production execution is now **future scope and must not start automatically**.
Before a future launch, the operator must:

1. select and freeze the finalist using criteria that do not inspect future MD
   outcomes;
2. confirm ColabFold and all other GPU jobs are stopped and the shared lock is
   free;
3. repeat the disk, dependency, CUDA, input-hash and configuration preflights;
4. use a new unique durable run directory—never overwrite the smoke run or a
   previous production attempt;
5. run all replicas, convergence analysis, evidence linking and schema
   validation in one chat-independent pipeline; and
6. report the production result as descriptive computational evidence only.

The completed 3 x 50 ps run remains an execution smoke test and must not be
promoted to stability evidence.

## Operational hold and ColabFold continuation

The public multi-pair ColabFold panel was intentionally stopped after four
complete constructs and 20 of 70 planned rankable models. The completed
`.done.txt` queries, PDB files, score JSON, MSA files, logs, manifests and input
FASTA remain in
`runs/folds/public-multipair-af2ptm-20261021/`. The interrupted fifth query is
not counted as complete and may be recomputed.

The resume runner checks the original FASTA hash, output path, construct count,
model count and prior manifest before reusing the directory. ColabFold's
default completed-query policy skips the four retained `.done.txt` queries.
Resume on or after 22 October with a new durable job-control directory:

```powershell
& 'C:\Users\HARDIK\.codex\scripts\Start-CodexRunPipeline.ps1' `
  -PipelineFile 'M:\ChatGPT\LinkerForge 2\configs\detached-public-multipair-af2ptm-resume-20261022.json' `
  -RunDirectory 'M:\ChatGPT\LinkerForge 2\runs\job-control\public-multipair-af2ptm-resume-20261022'
```

This resumes only the ColabFold panel and its summary/presentation stages.
Production OpenMM remains future scope and must not be started merely because
the ColabFold panel later completes.

## Relationship to the unresolved blockers

These future analyses strengthen structural interpretation only. They do not:

- replace two independent curators or adjudicate the 971-record packet;
- create an engineered-fusion outcome dataset;
- authorize compatibility-ranker training;
- provide missing protected detector cross-check files or relax the firewall;
- reproduce Rosetta KIC/loop sampling; or
- establish a family-disjoint comparison against independent structure-model
  families unless those additional predictors and targets are separately run.
