# PyMOL and Windows distribution policy

Date: 21 October 2026 (Asia/Calcutta)

## Release decision

LinkerForge installers do **not** bundle Incentive PyMOL, Open-Source PyMOL,
ColabFold, Rosetta, OpenMM, model parameters or sequence databases.

PyMOL remains an optional external visualization program. LinkerForge detects a
user- or institution-installed executable through documented paths or explicit
environment variables, verifies the selected PDB sequence, and then hands the
structure to that separate program. The in-application viewer remains 3Dmol.js.

This avoids representing PyMOL as the simulation engine and avoids silently
redistributing either a commercial Incentive build or a conda environment with
its own dependency licence set.

## PyMOL licence finding

The upstream Open-Source PyMOL repository publishes a BSD-like licence that
permits use, copying, modification and distribution provided the copyright and
permission notices accompany copies. Its trademark notice requires derivative
software containing PyMOL source to be clearly distinguished from Schrödinger's
products and restricts other endorsements using the PyMOL name.

Sources:

- <https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE>
- <https://github.com/schrodinger/pymol-open-source/blob/master/README.md>
- <https://www.pymol.org/>

Because LinkerForge does not contain or redistribute PyMOL source or binaries,
the Windows installer is not blocked by PyMOL redistribution. If a future
release bundles Open-Source PyMOL, it must separately inventory the exact
package and every bundled dependency, reproduce required notices, satisfy the
trademark terms, and undergo a new release review. The unlicensed Incentive
executable must never be bundled.

## Windows installer state

`src-tauri/tauri.bundle-preview.conf.json` enables an unsigned, current-user
NSIS installer for local packaging verification:

```powershell
cd frontend
npm run desktop:bundle:preview
```

This is a developer preview, not an official public release. The default
`tauri.conf.json` keeps bundling disabled so ordinary local builds cannot be
mistaken for signed distribution artifacts.

According to Tauri's Windows distribution documentation, NSIS produces a
Windows setup executable. Tauri's signing documentation explains that code
signing is needed for Microsoft Store listing and to avoid an untrusted
SmartScreen warning on browser-downloaded applications.

Sources:

- <https://v2.tauri.app/distribute/windows-installer/>
- <https://tauri.app/distribute/sign/windows/>

## Requirements before official publication

The release owner must supply information that cannot be inferred safely:

1. legal publisher name and distribution owner;
2. application identifier ownership and support URL;
3. Windows code-signing certificate or approved Azure Artifact Signing setup;
4. privacy, update and vulnerability-response policies;
5. complete third-party notices from the final dependency and installer bill of materials;
6. supported Windows versions and tested GPU/driver matrix;
7. signed installer hash, signature verification output and clean-machine installation test.

No secret, certificate or publisher identity belongs in repository configuration.
