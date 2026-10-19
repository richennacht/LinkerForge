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
