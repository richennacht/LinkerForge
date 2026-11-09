# Full-fusion local workflow

The design studio accepts two canonical amino-acid sequences (FASTA text is accepted), constructs explicit `Domain A + linker + Domain B` candidates, and exports the **entire fusion** as FASTA. Current custom candidates are deterministic GGGGS/EAAAK/TEV templates, not a learned compatibility ranking. A linker by itself must not be submitted to a structure predictor as a substitute for the fusion.

For a selected candidate:

1. Export its full-fusion FASTA from the viewer. Keep the domain order and selected linker unchanged.
2. On a machine with a separately installed, validated ColabFold `colabfold_batch` executable and its model/database resources, run:

   ```powershell
   python scripts/fold_fusion_colabfold.py --fasta linkerforge-fusion.fasta --output artifacts/fusion-fold-run-001 --timeout-minutes 120
   ```

   This script runs one bounded full-chain prediction. It records input hash, command, timestamps, completion/failure and selected PDB hash. It refuses to overwrite a nonempty output directory and requires a single-chain PDB whose C-alpha sequence matches the submitted FASTA. ColabFold is **not** installed or bundled by this repository, and this command has not been validated on the target cloud image.
3. Upload the selected PDB into LinkerForge. The browser checks the whole residue sequence against the selected candidate and colours Domain A green, the actual predicted linker residues magenta, and Domain B red. The view can expand to the viewport. Export the PyMOL colouring script if desktop PyMOL is desired; the in-browser renderer is 3Dmol.js, not PyMOL.

The uploaded file is user-supplied. Sequence matching cannot authenticate its provenance or establish biological function. The interface never fabricates PDB coordinates. A scientifically complete result still requires predictor confidence, interdomain PAE, closure/clash/domain-preservation checks, baselines, and experimental validation as appropriate.

## Individual domain references by PDB accession

The Prediction view accepts a four-character RCSB PDB accession, for example `1UBQ`, for each individual domain. LinkerForge downloads the public PDB coordinate file and searches every protein chain for the exact amino-acid sequence entered for that domain. A unique chain segment is reduced to a single-chain PDB before it is shown in 3Dmol.js or sent to the desktop PyMOL command. This makes the displayed reference depend on both the accession and the entered sequence rather than on a hard-coded example.

If no exact chain or chain segment matches, the app rejects the structure and reports the mismatch. It does not silently align gaps, substitute a predicted structure, or display a different chain. The local-file fallback remains available for private, unpublished, or predicted coordinates and applies the same exact sequence check. RCSB retrieval provides a public coordinate reference; it does not establish that the structure is the correct biological state, assembly, construct, or experimental comparator for a proposed fusion.
