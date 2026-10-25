# Historical classifier and mined-linker data releases

These are the **actual byte-for-byte local outputs** used/reported in the October
2026 classifier work, published for inspection and reproduction. They are two
different datasets. Neither is a certified independent test or experimental
proof that a mined span is a functional fusion linker.

## 1. Original B1 V2 classifier split

[`dfl_b1_v2_family_split/`](./dfl_b1_v2_family_split/) contains the original
DisProt 2026_06-derived training split for the three-layer DFL CNN: 60 train,
28 development, and 9 auxiliary proteins, plus its family manifest and
**sequence-only** CAID3/DP94 exclusion registries. The registries do not
contain benchmark labels. The manifest SHA-256 is
`d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`.

Labels are residue-level: `1` is a curated DisProt flexible-linker region
(`IDPO:0000033`), `0` is explicitly curated experimental order
(`IDPO:0000006`), and `x` is unknown. An `x` must never be converted to a
negative. The auxiliary set is an experimental-order **proxy**, not a
fully adjudicated no-linker population. The data include some publicly
released, curator-validated regions carrying DisProt's internal
`unpublished=true` flag; the source policy is recorded in the manifest.

The original training command expects this folder at
`data/processed/dfl_b1_v2_family_split/`; copy the folder there only to
reproduce the historical run. The current prospective benchmark firewall
still governs **new** training, and publication of these files does not make
them eligible for a new blind comparison. See
[`CLASSIFIER_DATA_GATE.md`](../CLASSIFIER_DATA_GATE.md) and
[`RESEARCH_HANDOFF.md`](../RESEARCH_HANDOFF.md).

## 2. Mined-and-natural composite

[`mined_and_natural_linkers_v1.zip`](./mined_and_natural_linkers_v1.zip)
contains the original six-file output, including `all_records.jsonl` and
the source-tier partitions. The independently browsable
[`manifest.json`](./mined_and_natural_linkers_v1.manifest.json) has SHA-256
`b1b8c5b34eb9a98eb627005a810f367e4364be9a75447af9740b3bff6d3e3232`.
Uncompressed member checksums and row counts are pinned inside it.
The ZIP itself has SHA-256
`402071ae2e07e9060a36a0c6fa3a960edabbec545cbfa9e20cac6f9ea6e57da0`.

| Tier | Records | Meaning |
| --- | ---: | --- |
| High-score model-mined gaps | 6,556 | Predicted candidates, **not** verified DFL or training truth |
| Provisional model-mined gaps | 21,039 | Lower-confidence candidates, **not** verified DFL or training truth |
| CATH/PDB-derived natural segments | 91 | Structure/domain-boundary evidence; not functional fusion validation |
| DisProt-curated DFL segments | 57 | Curated DFL annotations, distinct from mined calls |
| Literature engineered templates | 23 | Cited sequences, distinct from natural observations |

There are 27,766 rows in `all_records.jsonl`; the three partition files
sum to exactly that count. `b3_no_linker_specificity_panel.jsonl` has 40
additional diagnostic controls and is **not** part of the 27,766. The archive
is compressed to keep Git manageable; `all_records.jsonl` duplicates the
partition content by design. To use the loader, extract these six files into
`data/processed/mined_and_natural_linkers_v1/` without changing bytes.

The archive retains its historical manifest unmodified for reproducibility.
Its statement that the B2 procedure *guarantees* finite-sample risk control
is **not supported**: the threshold search was adaptive on correlated
residues, and the production span caller uses a different decision rule.
Likewise, the 40-protein proxy does not establish proteome-wide specificity.
The corrected interpretation is in
[`reports/PAPER_CLAIMS_AUDIT_2026-10-19.md`](../reports/PAPER_CLAIMS_AUDIT_2026-10-19.md).
Do not turn the mined candidate tiers into positive labels for classifier
training or claim 99% precision from this release.

## Verify and attribute

Run `python scripts/validate_public_dataset_release.py` from the repository
root. It checks the archived member list, SHA-256 digests, row counts, B1
masks, tier counts and unique composite record IDs. For the stronger original
sequence/family preflight, run
`python scripts/preflight_dfl_family_training_v2.py datasets/dfl_b1_v2_family_split/manifest.json`
in the project detector environment.

Derived data cite [DisProt 2026_06](https://www.disprot.org/) ([CC BY 4.0](https://disprot.org/about)),
[UniProtKB](https://www.uniprot.org/api-documentation/support-data) (CC BY 4.0),
[CATH v4.4](https://www.cathdb.info/) ([CC BY 4.0](https://www.cathdb.info/wiki/doku/?id=start)),
[wwPDB](https://www.wwpdb.org/about/usage-policies) (CC0), and the
[DLD public download](https://pcrgwd.ucd.ie/linker/download/). The DLD
portal does not state an explicit redistribution licence; its contribution
here is short, source-identified **derived candidate spans**, not a mirror of
the DLD archive. Check the original source terms for downstream redistribution.
The 23 templates cite Chen, Zaro and Shen (2013) in their individual records.
The model checkpoint and raw third-party archives are not included.
