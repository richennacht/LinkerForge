# Verified public acceptance run

## Scope

The completed public panel contains three protein-domain pairs and six linker-containing constructs per pair. ColabFold generated five AlphaFold2-ptm models for each construct. The presentation reports the rank-1 model for all 18 linker-containing constructs, representing 90 rankable models. Direct fusions without a linker are excluded from the presentation comparison.

The run used public sequences, MMseqs2 UniRef/environmental MSAs, no templates, and seed 0. It did not read protected benchmark labels and did not train or select a LinkerForge compatibility model.

## Highest observed linker-containing result

The highest whole-fusion mean pLDDT among the 18 linker-containing constructs was the Crambin–BPTI fusion using `EAAAKEAAAKEAAAK`:

- Whole-fusion mean pLDDT: 93.81
- Linker mean pLDDT: 84.85
- pTM: 0.69
- Mean bidirectional cross-domain PAE: 9.75 Å
- Rank-1 model: `alphafold2_ptm_model_4_seed_000`
- Public references: PDB 1CRN chain A and PDB 4PTI chain A

This is a post-prediction confidence comparison. It supports inspection of the predicted structure, but it does not establish expression, activity, stability, compatibility, or biological superiority.

## Linker-only result matrix

| Protein pair | Linker | Sequence | Length | Fusion pLDDT | Linker pLDDT | pTM | Cross-domain PAE (Å) |
|---|---|---|---:|---:|---:|---:|---:|
| Crambin–BPTI | EAAAK15 | `EAAAKEAAAKEAAAK` | 15 | 93.81 | 84.85 | 0.69 | 9.75 |
| Protein G–repressor | GGGGS5 | `GGGGS` | 5 | 92.20 | 48.93 | 0.59 | 21.52 |
| Protein G–repressor | GGGGS10 | `GGGGSGGGGS` | 10 | 90.22 | 40.54 | 0.57 | 23.75 |
| Protein G–repressor | EAAAK15 | `EAAAKEAAAKEAAAK` | 15 | 89.73 | 52.77 | 0.56 | 24.28 |
| Protein G–repressor | BIRD18 | `KESGSVSSEQLAQFRSLD` | 18 | 89.72 | 64.91 | 0.54 | 26.02 |
| Protein G–repressor | GGGGS15 | `GGGGSGGGGSGGGGS` | 15 | 88.60 | 36.92 | 0.54 | 25.98 |
| Protein G–repressor | WALDO12 | `GSAGSAAGSGEF` | 12 | 87.94 | 41.61 | 0.55 | 25.44 |
| Ubiquitin–GB1 | EAAAK15 | `EAAAKEAAAKEAAAK` | 15 | 87.88 | 50.88 | 0.51 | 24.90 |
| Crambin–BPTI | GGGGS5 | `GGGGS` | 5 | 87.62 | 51.19 | 0.60 | 16.22 |
| Ubiquitin–GB1 | GGGGS5 | `GGGGS` | 5 | 87.03 | 48.57 | 0.56 | 21.14 |
| Crambin–BPTI | GGGGS10 | `GGGGSGGGGS` | 10 | 86.42 | 47.97 | 0.55 | 20.07 |
| Crambin–BPTI | BIRD18 | `KESGSVSSEQLAQFRSLD` | 18 | 85.00 | 60.91 | 0.57 | 18.80 |
| Crambin–BPTI | GGGGS15 | `GGGGSGGGGSGGGGS` | 15 | 84.20 | 52.32 | 0.50 | 24.08 |
| Ubiquitin–GB1 | WALDO12 | `GSAGSAAGSGEF` | 12 | 84.09 | 44.42 | 0.54 | 21.77 |
| Crambin–BPTI | WALDO12 | `GSAGSAAGSGEF` | 12 | 83.42 | 46.44 | 0.55 | 20.37 |
| Ubiquitin–GB1 | GGGGS15 | `GGGGSGGGGSGGGGS` | 15 | 82.75 | 38.99 | 0.53 | 23.43 |
| Ubiquitin–GB1 | GGGGS10 | `GGGGSGGGGS` | 10 | 81.74 | 37.73 | 0.53 | 22.77 |
| Ubiquitin–GB1 | BIRD18 | `KESGSVSSEQLAQFRSLD` | 18 | 80.21 | 44.87 | 0.52 | 24.34 |

## Metric interpretation

- Fusion pLDDT is AlphaFold local-coordinate confidence averaged across the full fusion.
- Linker pLDDT is the same confidence measure restricted to the inserted residues. A low value can be consistent with flexibility but does not prove it.
- pTM is AlphaFold's predicted global topology confidence and is not a linker-success probability.
- Cross-domain PAE is the mean bidirectional uncertainty between Domain A and Domain B. Lower values indicate more confident relative placement, not validated biological activity.

## Reproducibility record

The versioned summary schema is `linkerforge.public-multipair-structure-comparison.v1`. The local scientific record retains the immutable FASTA and panel manifest, exact PDB and score-JSON paths, per-file SHA-256 hashes, prediction settings, and all 105 rankable models including the three direct-fusion controls. The presentation deliberately uses the linker-only subset.

The source run is stored under `runs/folds/public-multipair-af2ptm-20261021`. Generated run artifacts remain outside version control because `/runs/` is excluded by repository policy. The committed presentation builder and JSON schemas reproduce the report from a provenance-bearing summary.

## Claim boundary and deferred work

This panel demonstrates repeatable structure prediction and confidence extraction across three short public constructs. It does not validate the proposed compatibility ranker or generalization to arbitrary proteins. Rosetta, production OpenMM, ensemble agreement, reference-domain lDDT/TM-score, automatic functional-site mapping, DSSP, and post-OpenMM change analysis remain future work. OpenMM was not started.
