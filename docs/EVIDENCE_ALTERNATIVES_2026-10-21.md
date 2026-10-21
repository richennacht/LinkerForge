# Evidence alternatives and future integrations

This addendum covers the fusion application and its external-method audit.
It changes neither the active DFL training contract nor the protected-data
firewall. An available alternative answers a specific question; it does not
turn an unmeasured outcome into a positive result. Wet-lab execution remains
outside this project's requested scope.

| Evidence question | Practical method | Implementation and remaining scope |
| --- | --- | --- |
| Does a linker bury a domain residue in the predicted pose? | Differential solvent-accessible surface area (SASA) with Biopython Shrake–Rupley | Implemented CPU analysis, sequence and coordinate hash checks, per-residue outputs and evidence-record linking. Requires researcher-supplied functional annotations to interpret protected sites. |
| Is an annotated residue near the linker? | Existing heavy-atom distance/SASA report plus the app's Cα screen | Implemented. Ligands, binding partners and conformational ensembles remain future context. |
| Does the structure contain unrealistic atomic contacts? | Pinned Reduce 4.15 / Probe 2.18 hydrogen-aware severe-overlap screen | Implemented and completed on seven public predictions. It reports transparent unique-pair counts and never calls them the canonical MolProbity clashscore. No LinkerForge acceptance threshold is frozen. |
| Can strained coordinates relax? | Existing OpenMM minimization | Available and executed as part of the prior smoke protocol. Energy decrease does not establish biological stability or replace Rosetta loop reconstruction. |
| Can a missing linker conformation be reconstructed? | Rosetta KIC or another separately validated loop-sampling method | Future integration. Surface area, clash checks and minimization do not replace this test. |
| Are results reproducible across public constructs? | Matched ColabFold controls across three public pairs | Existing single-pair panel complete; additional two-pair durable prediction job is separate and ongoing. This panel is descriptive, not a family-disjoint prospective evaluation. |
| Does an external detector recognize inserted linkers? | Published DFLpred implementation on the existing 28-construct public panel | Executed comparator. It assesses flexible-linker propensity. LINKER-Pred, TransDFL and DisoFLAG remain future pinned comparisons. |
| How did historical LinkerForge detector heads compare on identical public targets? | Frozen-score retrospective on the 27-protein CAID3 common set | Completed negative result: best historical LinkerForge APS/Fmax was 0.089/0.165 versus 0.482/0.491 for IPA-AF2-Linker, 0.429/0.473 for LINKER-Pred-Lite and 0.409/0.464 for DisoFLAG-DFL. This spent public result cannot guide tuning or selection, contains only linker-positive proteins, and does not assess fusion compatibility. Fresh pinned inference remains future work. |
| Are structural conclusions robust across architectures? | ESMFold; another validated independent model family; OpenFold for AlphaFold-family reproducibility | Future integration with the same sequences and analysis rules. Five AF2 parameter models are not five independent model families. |
| Is the compatibility ranker scientifically calibrated? | Independently evidenced engineered-fusion outcomes and an untouched grouped evaluation | No computational substitution establishes this. Future scientific development. The mined corpus remains useful for discovery/retrieval and reviewed weak supervision. |
| Are MD trajectories long enough to support stability claims? | Replicated production simulations with convergence checks | Storage-aware v2 protocol and launch preflight are implemented and pass: 10 ps metrics/checkpoints are retained while full-system DCD frames are stored every 50 ps. Execution remains future until the active ColabFold GPU panel finishes; the existing 3 × 50 ps run verifies software execution only. |
| Can the app run without institutional cloud access? | Existing local Tauri/WSL execution and durable run manifests | Implemented for bounded public inputs. Approved institutional/cloud execution remains future deployment. |
| Can research results be preserved without a signing service? | Existing checksummed research ZIP and SHA-256 history | Implemented integrity alternative. Publisher signing and institutional signatures remain future release/deployment work. |

## Surface-burial calculation

`scripts/assess_fusion_surface_burial.py` reads a completed structure-comparison
record and the public sequence manifest. Each input must match both its recorded
PDB hash, complete sequence, and domain/linker boundaries. It rejects ambiguous chains, models,
alternate locations, missing backbone atoms and unsupported residues/elements.
It records source identifiers, configuration, executable, arguments, code hash,
commit and software versions in `linkerforge.fusion-surface-burial.v1` JSON.

For every domain residue it calculates:

1. SASA in the full fusion.
2. SASA after removing linker atoms while retaining both domains in place.
3. SASA in its isolated domain at those same coordinates.

`SASA_without_linker − SASA_full` estimates surface additionally buried by the
linker in that pose. `SASA_isolated_domain − SASA_full` measures burial by the
linker and the other domain together. Both calculations use the same 1.4 Å
probe, heavy-atom radii, and 960 surface points per atom. They do not simulate
the physically separated domain, restore a natural terminus, or relax it.

All domain residues are retained. No residue is called functional or protected
without an annotation. No acceptance cutoff or winner is inferred from the
measurements. High PAE and low local confidence must remain visible when
interpreting burial. Missing sidechains, excluded ligands/partners and alternate
poses are limitations. Domain totals sum lost per-residue surface; they are not
the half-sum convention sometimes used for interface area.

The research command and detached dependency chain are available in
`configs/detached-surface-evidence-20261021.json`. They use the seven completed
1UBQ–2GB1 structures and leave the running GPU prediction job untouched.
The analysis is currently a research CLI/report workflow; automatic per-design
UI execution/import is future interface integration.

Output directory: `runs/structure-audits/surface-burial-20261021/`.
The durable job's `status.json` is authoritative for completion. Existing run
directories and evidence records are never overwritten.

## Reduce/Probe severe-overlap screen

The locked CPU environment in `requirements/molprobity-tools/` contains Reduce
4.15 and Probe 2.18. `scripts/assess_reduce_probe_overlaps.py` first verifies
the complete sequence, PDB hash and domain/linker boundaries, uses Reduce to
place and optimize hydrogens, then runs Probe with main-chain contacts and
standard-bond topology exclusions. Reciprocal atom contacts are deduplicated.
Raw contacts, hydrogenated structures, exact arguments, versions and hashes are
retained under `runs/structure-audits/reduce-probe-20261021/`.

The durable three-stage run completed successfully. Counts below describe the
fixed unrelaxed predicted coordinates; there is no preregistered pass/fail
threshold and lower values do not independently establish a better linker.

| Construct | Unique severe pairs | Pairs / 1000 explicit atoms | Worst overlap (Å) |
| --- | ---: | ---: | ---: |
| bird18 | 35 | 14.89 | 1.285 |
| direct0 | 19 | 9.13 | 1.117 |
| EAAAK15 | 26 | 11.39 | 1.234 |
| GGGGS10 acceptance | 48 | 22.22 | 1.473 |
| GGGGS15 | 29 | 13.19 | 1.259 |
| GGGGS5 | 29 | 13.67 | 1.179 |
| WALDO12 | 60 | 27.17 | 1.424 |

This rate is deliberately named `unique severe-overlap pairs per 1000 explicit
atoms`; it is **not** the canonical MolProbity clashscore. The screen excludes
ligands, partners and waters and cannot replace Rosetta linker reconstruction,
conformational sampling, energetic ranking or experiment.

## Storage-aware production MD alternative

The original full-system trajectory policy stored coordinates every 10 ps and
failed its conservative storage gate. The new, separately versioned
`configs/openmm-publication-protocol-v2-storage-aware.json` keeps structural
metrics, state data and checkpoints at 10 ps but stores full-system DCD frames
every 50 ps. It does not reduce the explicit-solvent atom count, simulation
length, three independent replicas or force-field calculations.

The same-construct measured preflight now projects 3,000 DCD frames, 1.232 GiB
of trajectory data and 31.76 GPU hours. With the unchanged 1.5× allowance plus
5 GiB recovery reserve, 6.849 GiB is required and the current 10.545 GiB free
space passes. Artifact:
`reports/openmm-production-preflight-storage-aware-2026-10-21.json`.

This only removes the storage blocker. A durable six-stage production pipeline
is now running in a wait state at
`runs/job-control/openmm-production-v2-after-colabfold-20261021/`. It waits for
the ColabFold comparison artifact, repeats the disk preflight, acquires the
shared GPU lock, runs all three replicas, produces descriptive five-block and
between-replicate diagnostics, links the evidence, and validates schemas. A
completed run still would not prove biological stability.

## Sources checked

- [Biopython SASA API](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html):
  Shrake–Rupley implementation, probe radius and sphere-resolution controls.
- [MolProbity CLI documentation](https://github.com/rlabduke/MolProbity/blob/master/doc/command_line_tools.txt):
  hydrogen addition/Probe clash analysis and geometry tools.
- [MolProbity methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1933162/):
  why explicit hydrogens are important for all-atom contact validation.
- [OpenMM user guide](https://docs.openmm.org/latest/userguide/): minimization
  and simulation execution; interpret results according to the executed protocol.

The external-method rows in the main audit retain their measured status.
Unavailable methods belong to future scope. Invalid sequence matches, malformed
outputs and measured geometry failures remain actual failed checks.
