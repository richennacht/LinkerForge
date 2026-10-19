"""Part 2 Geometry-Plus-Template and Authority Linker Design Package.

Integrates:
  - `linkerforge.geometry` (PR #3 deterministic 3D span, SASA, occlusion, and
    multi-tether satisfiability)
  - `linkerforge.properties` (PR #2 curated Chen et al. 2013 template library,
    cited protease motif lookup, and calibrated compositional rigidity caller)
  - `linkerforge.dataset_v1` (Part 1 provenance-stratified natural and domain-gated
    mined linker authorities)

Answers the three LinkerForge design questions (Q1 feasibility, Q2 attachment/flanking
context, Q3 candidate linker sequences) using deterministic flags rather than
uncalibrated fusion success probabilities (`SPEC.md` section 5.8).
"""

from .engine import (
    DesignCandidate,
    DesignReport,
    DesignRequest,
    build_rigid_cleavable_composite,
    design_linker_candidates,
)

__all__ = [
    "DesignCandidate",
    "DesignReport",
    "DesignRequest",
    "build_rigid_cleavable_composite",
    "design_linker_candidates",
]
