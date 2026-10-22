# LinkerForge presentation test: highest-confidence linker-containing acceptance construct

## What this example is

This is the linker-containing construct with the highest whole-fusion mean pLDDT in the completed matched public panel of 21 constructs. It is a reproducible software and structure-prediction demonstration, not evidence that the linker is biologically superior.

The direct-fusion Protein G–repressor control had the highest whole-fusion mean pLDDT overall (94.9), but it contains no linker. For a linker demonstration, use the Crambin–BPTI construct below.

## Design input

| Application field | Value |
| --- | --- |
| Design name | `Crambin–BPTI · highest-confidence acceptance example` |
| Domain A name | `Crambin` |
| Domain A PDB accession | `1CRN` |
| Domain A chain | `A` |
| Domain A sequence | `TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN` |
| Domain A length | 46 aa |
| Domain B name | `BPTI` |
| Domain B PDB accession | `4PTI` |
| Domain B chain | `A` |
| Domain B sequence | `RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA` |
| Domain B length | 58 aa |
| Intended linker length | 15 aa |
| Linker behaviour | Rigid / non-cleavable control |
| Linker candidate to select | `EAAAKEAAAKEAAAK` (`eaaak_15`) |
| Full fusion length | 119 aa |
| Complete fusion sequence | `TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYANEAAAKEAAAKEAAAKRPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA` |

## Fast presentation path

1. Create a design with the two domain sequences above, target length 15, rigid behaviour, and non-cleavable linker constraints.
2. Generate candidates and select `EAAAKEAAAKEAAAK`.
3. On Prediction, select the Crambin individual-domain tab, enter `1CRN`, and fetch it from RCSB.
4. Select the BPTI individual-domain tab, enter `4PTI`, and fetch it from RCSB.
5. Return to Full fusion model and load the completed rank-1 PDB:
   `runs/folds/public-multipair-af2ptm-20261021/prediction/LinkerForge_1CRN_4PTI_eaaak_15_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb`
6. Load the matching confidence JSON:
   `runs/folds/public-multipair-af2ptm-20261021/prediction/LinkerForge_1CRN_4PTI_eaaak_15_scores_rank_001_alphafold2_ptm_model_4_seed_000.json`
7. In the embedded viewer, show Cartoon + Domains and linker first. Then switch to B-factor / pLDDT colouring, molecular surface, spin, reset view, and Save PNG.

Using the completed files is recommended for a timed presentation. Starting a new ColabFold job proves live execution but its duration and network-dependent MSA step are not deterministic.

## Verified values to explain

| Metric | Observed value | Correct interpretation |
| --- | ---: | --- |
| Whole-fusion mean pLDDT | 93.8 | High local-coordinate confidence averaged over the complete predicted fusion |
| Domain A mean pLDDT | 94.3 | High local-coordinate confidence for Crambin in this model |
| Linker mean pLDDT | 84.9 | Model confidence for the inserted residues; not a compatibility probability |
| Domain B mean pLDDT | 95.8 | High local-coordinate confidence for BPTI in this model |
| pTM | 0.69 | AlphaFold2-ptm global topology confidence, not biological success |
| Mean cross-domain PAE | 9.7 Å | Moderate model uncertainty in the relative placement of the two domains; lower is more confident |
| Descriptive panel rank | 2 of 21 | Post-hoc whole-fusion mean-pLDDT order; rank 1 is a no-linker direct-fusion control |

## One-sentence presentation claim

“Among the linker-containing constructs in our matched 21-construct public acceptance panel, this Crambin–EAAAKEAAAKEAAAK–BPTI prediction had the highest observed whole-fusion mean pLDDT; it demonstrates execution, confidence extraction, traceability, and interactive inspection, but not experimental linker superiority.”

## Provenance

- Summary: `runs/folds/public-multipair-af2ptm-20261021/three_pair_comparison.json`
- PDB SHA-256: `23fce4ce8acba965324ba608cab8b619662b2827dfdcd2cfc9192d30764562c4`
- Score JSON SHA-256: `10131927de17eca75bb9cb2ba233f0ca4408c54d76d29bc3c691cf883d9a1a52`
- Predictor/model: ColabFold, AlphaFold2-ptm model 4, seed 0, rank 1
