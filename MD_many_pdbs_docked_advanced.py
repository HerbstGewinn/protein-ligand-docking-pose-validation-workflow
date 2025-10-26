#!/usr/bin/env python3
"""
Molecular Dynamics Simulation Script (Batch Version)
This script performs a complete MD simulation workflow for multiple systems:
- Reads PDB IDs from a text file.
- For each system, it prepares it with PDBFixer.
- Assigns force field parameters with OpenFF, including custom ion parameters.
- Runs a multi-stage minimization and equilibration protocol.
- Runs a final production simulation.
"""

import sys
import os
from sys import stdout
import numpy as np
import mdtraj
from openff.toolkit import ForceField, Molecule, Topology
from openff.interchange import Interchange
from openff.units import unit as openff_unit
from openmm import Vec3, LangevinMiddleIntegrator, CustomExternalForce, MonteCarloBarostat
from openmm.app import PDBFile, Simulation, PDBReporter, StateDataReporter, DCDReporter
from openmm.unit import nanometer, nanometers, molar, kelvin, picosecond, picoseconds, femtosecond, angstroms, kilocalories_per_mole, atmosphere
from pdbfixer import PDBFixer
import argparse

# --- Function Definitions ---

def describe_state_of(simulation: Simulation, name: str = "State", pdb_id: str = ""):
    """Describes the current state of a simulation."""
    state = simulation.context.getState(getEnergy=True, getForces=True)
    forces = [np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in state.getForces()]
    max_force = max(forces)
    max_force_index = np.argmax(forces)
    print(
        f"  State '{name}' for {pdb_id} has energy {round(state.getPotentialEnergy()._value, 2)} kJ/mol "
        f"with maximum force {round(max_force, 2)} kJ/(mol nm) on atom {max_force_index}."
    )
    # Also, it seems you intended to use pdb_id in writePDB as well.
    # The original writePDB call inside describe_state_of was: writePDB(simulation, name)
    # You might want to change it to something like:
    writePDB(simulation, f"{pdb_id}_{name}")

def applyPDBFixer(filename):
    """Applies PDBFixer to a PDB file."""
    print(f"    Running PDBFixer on: {filename}")
    fixer = PDBFixer(filename=filename)
    fixer.findMissingResidues()
    if len(fixer.missingResidues) > 0:
        print(f"    [ERROR] PDBFixer detected missing residues in {filename}. Cannot proceed.")
        return None
    fixer.findMissingAtoms()
    if len(fixer.missingAtoms) > 0:
        print(f"    [ERROR] PDBFixer detected missing atoms in {filename}. Cannot proceed.")
        return None

    fixer.addSolvent(padding=1.2 * nanometer, ionicStrength=0.15*molar)
    return fixer

def add_harmonic_restraints(system, topology, positions, expression, strength=25):
    """Adds harmonic restraints to the system."""
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force_strength = strength * kilocalories_per_mole/angstroms**2
    force.addGlobalParameter("k", force_strength)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    mdtraj_topology = mdtraj.Topology.from_openmm(topology)
    idlist = mdtraj_topology.select(expression)
    print(f"      {len(idlist)} atoms match restraint expression '{expression}'")
    num = 0
    for i, (atom_crd, atom) in enumerate(zip(positions, topology.atoms())):
        if i in idlist:
            num += 1
            force.addParticle(i, atom_crd.value_in_unit(nanometers))

    if num > 0:
        index = system.addForce(force)
        print(f'      Restraints added to {num} atoms as force index {index}')
        return index
    return -1

def removeCustomExternalForces(system):
    """Removes all CustomExternalForce objects from a system."""
    force_indices_to_remove = []
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, CustomExternalForce):
            force_indices_to_remove.append(i)

    for index in sorted(force_indices_to_remove, reverse=True):
        print(f'      Removing CustomExternalForce at index {index}')
        system.removeForce(index)

def listForces(system):
    """Lists all forces in a system."""
    print("    Current forces in system:")
    for id_val in range(system.getNumForces()):
        force = system.getForce(id_val)
        print(f'      Force {id_val}: {force.__class__.__name__}')

def writePDB(simulation, name):
    """Writes the current simulation state to a PDB file."""
    with open(f'{name}.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

def get_unique_molecules(topology):
    uniques = {}
    for mol in topology.molecules:
        if mol.n_atoms > 500:
            if mol.hill_formula in uniques.keys():
                print("Error: two large molecules with identical hill formula")
                exit(0)
            uniques[mol.hill_formula] = mol
        else:
            if mol.to_smiles in uniques.keys():
                pass
            else:
                uniques[mol.to_smiles()] = mol
    #print(uniques.keys())
    #print(uniques.items())
    return list(uniques.values())

def molecule_already_in_list(mol, listofmols):
    """Checks if a molecule is already in a list."""
    smiles1 = mol.to_smiles()
    for mol2 in listofmols:
        if smiles1 == mol2.to_smiles():
            return True
    return False

def summarize_molecule_list(mols):
    print(len(mols), ' unique molecules:')
    for mol in mols:
        print('  ', mol.hill_formula)

def add_ion_params(forcefield, name, atomic_number, charge, rmin_half, epsilon):
    """Adds custom ion parameters to the force field."""
    library_charge_handler = forcefield.get_parameter_handler("LibraryCharges")
    smarts = f'[#{atomic_number}{charge:+}]'
    smirks = f"[{smarts}:1]"
    try:
        library_charge_handler.add_parameter({"smirks": smirks, "charge": [charge * openff_unit.elementary_charge], "id": name})
        vdw_handler = forcefield.get_parameter_handler("vdW")
        vdw_handler.add_parameter({"smirks": smirks, "rmin_half": rmin_half * openff_unit.angstrom, "epsilon": epsilon * openff_unit.kilocalories_per_mole, "id": name})
        print(f"    Custom parameters added for {name}.")
    except Exception as e:
        print(f"    [WARN] Could not add parameters for {name}. It might already exist. Error: {e}")

def process_single_pdb(pdb_id, forcefield):
    """Runs the full MD preparation and simulation workflow for a single PDB ID."""
    print(f"\n===== Starting Workflow for PDB ID: {pdb_id} =====")

    # Define dynamic filenames
    pdbfile = f"{pdb_id}_prepared_dock_advanced.pdb"
    sdf_file = f"{pdb_id}_Ligand.sdf"
    output_tag = f"{pdb_id}_MD_out_"

    # Check for required input files
    if not (os.path.exists(pdbfile) and os.path.exists(sdf_file)):
        print(f"  [ERROR] Missing input files for {pdb_id} ('{pdbfile}' or '{sdf_file}'). Skipping.")
        return False

    try:
        # Load ligand
        ligand = Molecule.from_file(sdf_file, file_format="SDF")
        uniques = [ligand]

        # Prepare system with PDBFixer
        fixer = applyPDBFixer(pdbfile)
        if fixer is None: return False

        # Build initial topology to discover all unique molecules
        complex_topology_nobox = Topology.from_pdb(pdbfile, unique_molecules=uniques)
        uniques_from_pdb = get_unique_molecules(complex_topology_nobox)
        for mol in uniques_from_pdb:
            if not molecule_already_in_list(mol, uniques):
                uniques.append(mol)

        # Add standard solvent/ions to the list for parameterization
        print('  Adding standard box components to unique molecules list...')
        for smiles in ['[H]O[H]', '[Na+]', '[Cl-]']:
            mol = Molecule.from_smiles(smiles)
            if not molecule_already_in_list(mol, uniques):
                uniques.append(mol)

        summarize_molecule_list(uniques)

        # Re-build final topology using the now solvated structure from PDBFixer
        complex_topology = Topology.from_openmm(fixer.topology, unique_molecules=uniques)
        interchange = forcefield.create_interchange(topology=complex_topology)

        # --- Equilibration and Production ---
        system = interchange.to_openmm_system()
        topology = interchange.to_openmm_topology()
        positions = fixer.positions

        # Minimize non-solute
        print("  Starting Minimization...")
        add_harmonic_restraints(system, topology, positions, 'not (water or resname NA CL)', 25)
        integrator = LangevinMiddleIntegrator(298*kelvin, 2/picosecond, 0.001*picoseconds)
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        describe_state_of(simulation, name=output_tag + 'start')
        simulation.minimizeEnergy()
        describe_state_of(simulation, name=output_tag + 'minimized')

        # Restrained NVT warmup
        print("  Starting NVT Warmup...")
        removeCustomExternalForces(system)
        add_harmonic_restraints(system, topology, positions, 'not (element H or water or resname NA CL)', 25)
        simulation.context.reinitialize(preserveState=True)
        simulation.minimizeEnergy()
        describe_state_of(simulation, output_tag + 'minimized2')
        simulation.context.setVelocitiesToTemperature(5*kelvin)
        simulation.reporters.append(StateDataReporter(file=output_tag + 'equil.csv', reportInterval=2500, step=True, time=True, potentialEnergy=True, temperature=True, density=True))

        print('  Warming up the system...')
        for temperature in range(5, 298, 20):
            integrator.setTemperature(temperature*kelvin)
            simulation.step(1000)
        integrator.setTemperature(298*kelvin)
        simulation.step(1000)
        describe_state_of(simulation, name=output_tag + 'afterNVT')

        # NPT equilibration
        print("  Starting NPT Equilibration...")
        integrator.setFriction(1.0 / picosecond)
        integrator.setStepSize(0.002 * picoseconds)
        barostat = MonteCarloBarostat(1*atmosphere, 298*kelvin, 10)
        system.addForce(barostat)
        simulation.context.reinitialize(preserveState=True)
        simulation.step(25000)
        describe_state_of(simulation, name=output_tag + 'afterNPT')

        # NPT with backbone restraints
        barostat.setFrequency(25)
        removeCustomExternalForces(system)
        add_harmonic_restraints(system, topology, positions, 'backbone or not (element H or water or protein or resname NA CL)', 25)
        simulation.context.reinitialize(preserveState=True)
        simulation.step(25000)
        describe_state_of(simulation, name=output_tag + 'after2ndNPT')

        # NPT with restraint release
        for k in range(20, 0, -5):
            print(f'    Restraint strength k = {k} kcal/mol/A^2')
            simulation.context.setParameter('k', k * kilocalories_per_mole/angstroms**2)
            simulation.step(5000)
        describe_state_of(simulation, name=output_tag + 'after3rdNPT')
        simulation.saveState(f'{pdb_id}_equilibrated.xml')
        print(f"  Equilibrated state for {pdb_id} saved to {pdb_id}_equilibrated.xml")
        checkpoint = f'{output_tag}equilibrated.chk'
        simulation.saveCheckpoint(checkpoint)
        print(f"  Equilibrated state saved to checkpoint: {checkpoint}")

        # --- Integrated Production Simulation ---
        print(f"  Starting Production Simulation...")
        removeCustomExternalForces(system)
        simulation.context.reinitialize(preserveState=True)
        listForces(system)

        production_trajectory_file = f"{output_tag}production_trajectory.dcd"
        simulation.reporters.append(DCDReporter(file=production_trajectory_file, reportInterval=50000))

        production_steps = 500000
        print(f"  Running {production_steps} steps of production dynamics...")
        simulation.step(production_steps)
        describe_state_of(simulation, name=f'after1ns', pdb_id=pdb_id)
        simulation.saveState(f'{pdb_id}_after1ns.xml')
        print(f"  Final state for {pdb_id} saved to {pdb_id}_final_state.xml")
        checkpoint = f'{output_tag}_after1ns.chk'
        simulation.saveCheckpoint(checkpoint)
        print(f"  Final state saved to checkpoint: {checkpoint}")
        print(f"  [SUCCESS] Production simulation for {pdb_id} finished.")
        return True

    except Exception as e:
        print(f"  [CRITICAL ERROR] An unexpected error occurred during workflow for {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function to loop through PDB IDs and run the MD workflow."""
    parser = argparse.ArgumentParser(description='Batch MD simulation from PDB ID list')
    parser.add_argument('-list', required=True, help='Text file containing PDB IDs (one per line)')
    args = parser.parse_args()

    # Load force field and add custom ion parameters once
    print("Loading shared force field parameters...")
    forcefield = ForceField("openff-2.2.1.offxml", "ff14sb_off_impropers_0.0.4.offxml")
    add_ion_params(forcefield, 'Ca2+LiCM', 20, 2, 1.635, 0.09788018)
    add_ion_params(forcefield, 'Mg2+LiCM', 12, 2, 1.364, 0.01055378)
    add_ion_params(forcefield, 'Zn2+LiCM', 30, 2, 1.276, 0.00354287)
    print("Force field and custom ion parameters loaded successfully.")

    if not os.path.exists(args.list):
        print(f"ERROR: PDB ID list file '{args.list}' not found!")
        sys.exit(1)

    with open(args.list, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip()]

    print(f"Found {len(pdb_ids)} PDB IDs to process.")

    successful = []
    failed = []

    for pdb_id in pdb_ids:
        if process_single_pdb(pdb_id, forcefield):
            successful.append(pdb_id)
        else:
            failed.append(pdb_id)

    print(f"\n{'='*60}")
    print("BATCH PROCESSING SUMMARY")
    print(f"Total: {len(pdb_ids)}, Successful: {len(successful)}, Failed: {len(failed)}")
    if failed:
        print("Failed PDB IDs:", ", ".join(failed))

if __name__ == "__main__":
    main()
