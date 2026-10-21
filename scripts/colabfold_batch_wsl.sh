#!/usr/bin/env bash
set -euo pipefail

# LocalColabFold is installed in the dedicated M:-backed WSL distribution.
# Keep this wrapper on WSL's default PATH so the Tauri backend can invoke it
# with `wsl.exe --exec colabfold_batch ...` without activating Pixi first.
readonly LOCALCOLABFOLD_ROOT="${LINKERFORGE_LOCALCOLABFOLD_ROOT:-/home/linkerforge/localcolabfold}"
readonly LOCALCOLABFOLD_BIN="${LOCALCOLABFOLD_ROOT}/.pixi/envs/default/bin"

if [[ ! -x "${LOCALCOLABFOLD_BIN}/colabfold_batch" ]]; then
  printf 'LocalColabFold executable not found at %s\n' \
    "${LOCALCOLABFOLD_BIN}/colabfold_batch" >&2
  exit 127
fi

export PATH="${LOCALCOLABFOLD_BIN}:${PATH}"
export TF_FORCE_UNIFIED_MEMORY="${TF_FORCE_UNIFIED_MEMORY:-1}"
export TF_FORCE_GPU_ALLOW_GROWTH="${TF_FORCE_GPU_ALLOW_GROWTH:-true}"
export XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-4.0}"
export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"
export XLA_PYTHON_CLIENT_ALLOCATOR="${XLA_PYTHON_CLIENT_ALLOCATOR:-platform}"

exec "${LOCALCOLABFOLD_BIN}/colabfold_batch" "$@"
