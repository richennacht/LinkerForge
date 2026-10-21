#!/usr/bin/env python3
"""Run a provenance-bearing replicated OpenMM minimization/MD protocol.

This is an execution workflow, not evidence that a short smoke trajectory has
converged. Publication interpretation requires the preregistered production
configuration plus convergence analysis across all completed replicates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform as py_platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import openmm
import numpy as np
from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, XmlSerializer, unit
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def energy_kj(state) -> float:
    return round(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 6)


def ca_indices(topology, protein_residue_count: int) -> list[int]:
    indices: list[int] = []
    for residue in list(topology.residues())[:protein_residue_count]:
        atoms = [atom for atom in residue.atoms() if atom.name == "CA"]
        if len(atoms) != 1:
            raise ValueError(f"Expected one CA atom in protein residue {residue.index + 1}; found {len(atoms)}")
        indices.append(atoms[0].index)
    return indices


def coordinates_angstrom(state, indices: list[int]) -> np.ndarray:
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    return np.asarray(positions[indices], dtype=float)


def aligned_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    mobile_centered = mobile - mobile.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    covariance = mobile_centered.T @ reference_centered
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    delta = mobile_centered @ rotation - reference_centered
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def frame_metrics(
    coordinates: np.ndarray,
    reference: np.ndarray,
    domain_a_length: int,
    linker_length: int,
    linker_contact_cutoff_angstrom: float,
    linker_contact_exclusion_residues: int,
) -> dict[str, float | int]:
    domain_b_start = domain_a_length + linker_length
    domain_a = coordinates[:domain_a_length]
    linker = coordinates[domain_a_length:domain_b_start]
    domain_b = coordinates[domain_b_start:]
    ref_domain_a = reference[:domain_a_length]
    ref_linker = reference[domain_a_length:domain_b_start]
    ref_domain_b = reference[domain_b_start:]
    centered = coordinates - coordinates.mean(axis=0)
    eligible_domain_a = domain_a[: max(0, len(domain_a) - linker_contact_exclusion_residues)]
    eligible_domain_b = domain_b[min(linker_contact_exclusion_residues, len(domain_b)) :]
    eligible_domains = np.concatenate([eligible_domain_a, eligible_domain_b], axis=0)
    distances = np.linalg.norm(linker[:, None, :] - eligible_domains[None, :, :], axis=2)
    return {
        "fusion_ca_rmsd_angstrom": aligned_rmsd(coordinates, reference),
        "domain_a_ca_rmsd_angstrom": aligned_rmsd(domain_a, ref_domain_a),
        "linker_ca_rmsd_angstrom": aligned_rmsd(linker, ref_linker),
        "domain_b_ca_rmsd_angstrom": aligned_rmsd(domain_b, ref_domain_b),
        "ca_radius_of_gyration_angstrom": float(np.sqrt(np.mean(np.sum(centered * centered, axis=1)))),
        "domain_centroid_distance_angstrom": float(np.linalg.norm(domain_a.mean(axis=0) - domain_b.mean(axis=0))),
        "minimum_linker_domain_ca_distance_angstrom": float(distances.min()),
        "linker_domain_ca_contacts": int(np.count_nonzero(distances <= linker_contact_cutoff_angstrom)),
    }


def summarize_metric_rows(rows: list[dict[str, float | int]]) -> dict[str, dict[str, float]]:
    keys = [key for key in rows[0] if key != "step"]
    return {
        key: {
            "mean": round(float(np.mean([float(row[key]) for row in rows])), 6),
            "minimum": round(float(np.min([float(row[key]) for row in rows])), 6),
            "maximum": round(float(np.max([float(row[key]) for row in rows])), 6),
            "final": round(float(rows[-1][key]), 6),
        }
        for key in keys
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit("Output directory must be absent or empty; refusing to overwrite a simulation record")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    started = datetime.now(timezone.utc).isoformat()

    pdb = PDBFile(str(args.pdb))
    protein_residue_count = int(config["domain_a_length"]) + int(config["linker_length"]) + int(config["domain_b_length"])
    if pdb.topology.getNumResidues() != protein_residue_count:
        raise ValueError(
            f"Configured fusion length is {protein_residue_count}, but input PDB has "
            f"{pdb.topology.getNumResidues()} residues"
        )
    original_atom_count = pdb.topology.getNumAtoms()
    fixer = PDBFixer(filename=str(args.pdb))
    fixer.findMissingResidues()
    if fixer.missingResidues:
        raise ValueError(
            "Input prediction has missing residues; refusing to model unresolved sequence during the MD smoke test"
        )
    fixer.findMissingAtoms()
    if fixer.missingAtoms:
        details = {
            f"{residue.chain.id}:{residue.index + 1}:{residue.name}": [atom.name for atom in atoms]
            for residue, atoms in fixer.missingAtoms.items()
        }
        raise ValueError(f"Input prediction has missing non-terminal heavy atoms: {details}")
    terminal_atoms = {
        f"{residue.chain.id}:{residue.index + 1}:{residue.name}": list(atom_names)
        for residue, atom_names in fixer.missingTerminals.items()
    }
    unexpected_terminal_atoms = {
        residue_id: atom_names
        for residue_id, atom_names in terminal_atoms.items()
        if atom_names != ["OXT"] or not residue_id.split(":")[1] == str(protein_residue_count)
    }
    if unexpected_terminal_atoms:
        raise ValueError(f"Unexpected terminal heavy-atom repair requested: {unexpected_terminal_atoms}")
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=float(config["ph"]))

    forcefield = ForceField(config["protein_forcefield"], config["water_forcefield"])
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        forcefield,
        model=config["water_model"],
        padding=float(config["solvent_padding_nm"]) * unit.nanometer,
        ionicStrength=float(config["ionic_strength_molar"]) * unit.molar,
        neutralize=True,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=float(config["nonbonded_cutoff_nm"]) * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=float(config["ewald_error_tolerance"]),
    )
    system.addForce(MonteCarloBarostat(
        float(config["pressure_bar"]) * unit.bar,
        float(config["temperature_kelvin"]) * unit.kelvin,
        int(config["barostat_interval_steps"]),
    ))
    platform_name = config["platform"]
    selected_platform = Platform.getPlatformByName(platform_name)
    properties = config.get("platform_properties", {})

    preparation_integrator = LangevinMiddleIntegrator(
        float(config["temperature_kelvin"]) * unit.kelvin,
        float(config["friction_per_ps"]) / unit.picosecond,
        float(config["timestep_fs"]) * unit.femtoseconds,
    )
    preparation = Simulation(modeller.topology, system, preparation_integrator, selected_platform, properties)
    preparation.context.setPositions(modeller.positions)
    initial = preparation.context.getState(getEnergy=True)
    preparation.minimizeEnergy(
        tolerance=float(config["minimization_tolerance_kj_mol_nm"]) * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=int(config["minimization_max_iterations"]),
    )
    minimized = preparation.context.getState(getPositions=True, getEnergy=True)
    minimized_positions = minimized.getPositions()
    protein_ca_indices = ca_indices(modeller.topology, protein_residue_count)
    minimized_reference = coordinates_angstrom(minimized, protein_ca_indices)
    prepared_pdb = args.output_dir / "prepared-minimized.pdb"
    with prepared_pdb.open("w", encoding="utf-8") as stream:
        PDBFile.writeFile(modeller.topology, minimized_positions, stream, keepIds=True)
    system_xml = args.output_dir / "system.xml"
    system_xml.write_text(XmlSerializer.serialize(system), encoding="utf-8")

    replicates = []
    for replicate in range(1, int(config["replicates"]) + 1):
        replicate_dir = args.output_dir / f"replicate-{replicate:03d}"
        replicate_dir.mkdir()
        replicate_system = XmlSerializer.deserialize(system_xml.read_text(encoding="utf-8"))
        seed = int(config["base_seed"]) + replicate - 1
        for force in replicate_system.getForces():
            if isinstance(force, MonteCarloBarostat):
                force.setRandomNumberSeed(seed)
        integrator = LangevinMiddleIntegrator(
            float(config["temperature_kelvin"]) * unit.kelvin,
            float(config["friction_per_ps"]) / unit.picosecond,
            float(config["timestep_fs"]) * unit.femtoseconds,
        )
        integrator.setRandomNumberSeed(seed)
        simulation = Simulation(modeller.topology, replicate_system, integrator, selected_platform, properties)
        simulation.context.setPositions(minimized_positions)
        simulation.context.setVelocitiesToTemperature(float(config["temperature_kelvin"]) * unit.kelvin, seed)
        simulation.step(int(config["equilibration_steps"]))
        report_interval = int(config["report_interval_steps"])
        trajectory_report_interval = int(config.get("trajectory_report_interval_steps", report_interval))
        production_steps = int(config["production_steps"])
        if production_steps % report_interval or production_steps % trajectory_report_interval:
            raise ValueError("Production steps must be divisible by metric and trajectory report intervals")
        simulation.reporters.append(DCDReporter(str(replicate_dir / "trajectory.dcd"), trajectory_report_interval))
        simulation.reporters.append(StateDataReporter(
            str(replicate_dir / "state.csv"), report_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            temperature=True, volume=True, density=True, speed=True,
            separator=",",
        ))
        simulation.reporters.append(CheckpointReporter(str(replicate_dir / "checkpoint.chk"), report_interval))
        production_started = datetime.now(timezone.utc).isoformat()
        metric_rows: list[dict[str, float | int]] = []
        completed_steps = 0
        while completed_steps < production_steps:
            step_count = min(report_interval, production_steps - completed_steps)
            simulation.step(step_count)
            completed_steps += step_count
            metric_state = simulation.context.getState(getPositions=True)
            metric_rows.append({
                "step": completed_steps,
                **frame_metrics(
                    coordinates_angstrom(metric_state, protein_ca_indices),
                    minimized_reference,
                    int(config["domain_a_length"]),
                    int(config["linker_length"]),
                    float(config["linker_contact_cutoff_angstrom"]),
                    int(config["linker_contact_exclusion_residues"]),
                ),
            })
        metrics_csv = replicate_dir / "structural-metrics.csv"
        with metrics_csv.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(metric_rows[0]))
            writer.writeheader()
            writer.writerows(metric_rows)
        final = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        final_pdb = replicate_dir / "final.pdb"
        with final_pdb.open("w", encoding="utf-8") as stream:
            PDBFile.writeFile(modeller.topology, final.getPositions(), stream, keepIds=True)
        simulation.saveState(str(replicate_dir / "final-state.xml"))
        replicates.append({
            "replicate": replicate,
            "seed": seed,
            "production_started_at": production_started,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "final_potential_energy_kj_mol": energy_kj(final),
            "trajectory_sha256": sha256(replicate_dir / "trajectory.dcd"),
            "trajectory_frames": production_steps // trajectory_report_interval,
            "trajectory_report_interval_steps": trajectory_report_interval,
            "state_csv_sha256": sha256(replicate_dir / "state.csv"),
            "structural_metrics_sha256": sha256(metrics_csv),
            "structural_metric_frames": len(metric_rows),
            "structural_metric_summary": summarize_metric_rows(metric_rows),
            "final_pdb_sha256": sha256(final_pdb),
        })
    manifest = {
        "schema_version": "linkerforge.openmm-run.v1",
        "status": "complete",
        "started_at": started,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_pdb": str(args.pdb),
        "input_pdb_sha256": sha256(args.pdb),
        "source_ids": [f"pdb-sha256:{sha256(args.pdb)}"],
        "config": config,
        "config_sha256": sha256(args.config),
        "prepared_pdb_sha256": sha256(prepared_pdb),
        "system_xml_sha256": sha256(system_xml),
        "preparation": {
            "method": "PDBFixer missing-heavy-atom audit plus hydrogen addition",
            "original_atom_count": original_atom_count,
            "protein_atom_count_after_repair_and_hydrogens": sum(
                1 for residue in list(modeller.topology.residues())[:protein_residue_count] for _ in residue.atoms()
            ),
            "terminal_atoms_added": terminal_atoms,
            "missing_residues_added": 0,
            "non_terminal_heavy_atoms_added": 0,
            "ph": float(config["ph"]),
        },
        "atom_count": modeller.topology.getNumAtoms(),
        "residue_count_including_solvent": modeller.topology.getNumResidues(),
        "initial_potential_energy_kj_mol": energy_kj(initial),
        "minimized_potential_energy_kj_mol": energy_kj(minimized),
        "replicates": replicates,
        "software": {
            "openmm": openmm.__version__,
            "pdbfixer": version("pdbfixer"),
            "python": sys.version,
            "platform": py_platform.platform(),
            "compute_platform": platform_name,
        },
        "provenance": {
            "commit_sha": subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=Path(__file__).resolve().parents[1],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip(),
            "worktree_modified": bool(subprocess.run(
                ["git", "status", "--porcelain", "--untracked-files=no"],
                cwd=Path(__file__).resolve().parents[1],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()),
            "config_sha256": sha256(args.config),
        },
        "evidence_quality": config["evidence_quality"],
        "claim_boundary": config["claim_boundary"],
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
