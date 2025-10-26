#!/usr/bin/env python
# -*- coding: utf-8 -*-

# extractLigand_ASL_batch.py
# Extracts ligand atoms based on a provided ASL from multiple prepared 
# Maestro files listed in an input file.

import argparse
import os
import sys

# --- Configuration ---
# Define the ASL to use for selecting the ligand atoms to KEEP.
# Examples:
# LIGAND_ASL = "ligand" # Use if prepwizard reliably tagged the ligand
LIGAND_ASL = "ligand" 
# LIGAND_ASL = "resname 0X5 and resnr 501" # Use for a specific ligand residue
# LIGAND_ASL = "res.ptype '0X5'" # Another way to select by residue name

# Define the suffix for input prepared complex files
INPUT_SUFFIX = "_prepared.maegz"
# Define the suffix for output extracted ligand files
OUTPUT_SUFFIX = "_Ligand_Docking_prep.maegz"
# ---------------------

# Check if SCHRODINGER environment variable is set
schrodinger_path = os.environ.get("SCHRODINGER")
if not schrodinger_path:
    print("Error: SCHRODINGER environment variable not set.")
    print("Please load the Schrodinger module environment first.")
    sys.exit(1)

# Try importing necessary Schrodinger modules
try:
    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.structutils import analyze
except ImportError as e:
    print(f"Error importing Schrodinger modules: {e}")
    print("Ensure your Schrodinger environment is correctly set up.")
    sys.exit(1)

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description=f"Extract ligand atoms (ASL: '{LIGAND_ASL}') from multiple prepared Maestro files.")
parser.add_argument("pdb_id_list_file", help="Text file containing a list of PDB IDs, one per line.")
args = parser.parse_args()

# --- Input File Existence Check ---
if not os.path.isfile(args.pdb_id_list_file):
    print(f"Error: PDB ID list file not found: {args.pdb_id_list_file}")
    sys.exit(1)

# --- Read PDB ID list ---
try:
    with open(args.pdb_id_list_file, 'r') as f:
        # Read lines, strip whitespace, and filter out empty lines
        pdb_ids = [line.strip() for line in f if line.strip()]
    if not pdb_ids:
        print(f"Error: No PDB IDs found in {args.pdb_id_list_file}")
        sys.exit(1)
    print(f"Found {len(pdb_ids)} PDB IDs to process.")
except Exception as e:
    print(f"Error reading PDB ID list file: {e}")
    sys.exit(1)

# --- Process each PDB ID ---
success_count = 0
fail_count = 0

for pdb_id in pdb_ids:
    print(f"--- Processing PDB ID: {pdb_id} ---")
    inputfile = f"{pdb_id}{INPUT_SUFFIX}"
    outputfile = f"{pdb_id}{OUTPUT_SUFFIX}"

    # --- DEBUG PRINT ADDED ---
    print(f"  Looking for input file: '{inputfile}' in current directory '{os.getcwd()}'") 
    # -------------------------

    # Check if input file exists
    if not os.path.isfile(inputfile):
        print(f"  Input file not found: {inputfile}. Skipping.")
        fail_count += 1
        continue

    # Check if output file already exists
    if os.path.isfile(outputfile):
        print(f"  Output file exists: {outputfile}. Skipping.")
        # success_count += 1 # Optionally count existing files as success
        continue

    # --- Main Logic for Extraction ---
    try:
        print(f"  Reading structure from: {inputfile}")
        st = None
        for st_i in StructureReader(inputfile):
            st = st_i
            break # Process only the first structure
        if st is None:
            print(f"  Error: No structures found in {inputfile}. Skipping.")
            fail_count += 1
            continue

        # Select atoms NOT matching the ligand ASL
        asl_to_delete = f'not ({LIGAND_ASL})'
        print(f"  Selecting atoms to DELETE using ASL: '{asl_to_delete}'")
        
        # Evaluate ASL for atoms to keep (for reporting)
        atoms_to_keep_indices = analyze.evaluate_asl(st, LIGAND_ASL)
        if not atoms_to_keep_indices:
             print(f"  Warning: ASL '{LIGAND_ASL}' matched 0 atoms in {inputfile}. Output will be empty.")
             # Decide whether to skip or write empty file - writing empty for now
        else:
             print(f"  ASL '{LIGAND_ASL}' matches {len(atoms_to_keep_indices)} atoms out of {st.atom_total} in input.")

        # Evaluate ASL for atoms to delete (returns 0-based indices)
        indices_to_delete_0based = analyze.evaluate_asl(st, asl_to_delete)

        if not indices_to_delete_0based:
            print(f"  No atoms found matching delete ASL '{asl_to_delete}'. Writing original structure (minus potential ligand).")
            # This case might happen if the input only contained the ligand
        else:
            
            print(f"  Deleting {len(indices_to_delete_0based)} non-ligand atoms...")
            st.deleteAtoms(indices_to_delete_0based) # Pass 0-based list directly
            # ---------------------------------------------------------
            # No st.repack() needed based on previous error

        # Write the remaining atoms (the ligand) to the output file
        print(f"  Writing extracted ligand to: {outputfile}")
        with StructureWriter(outputfile) as writer:
            writer.append(st)
        
        if os.path.isfile(outputfile):
             print(f"  Successfully created: {outputfile}")
             success_count += 1
        else:
             print(f"  Error: Output file was not created after writing.")
             fail_count += 1

    except Exception as e:
        print(f"  Error processing {pdb_id}: {e}")
        fail_count += 1
        # Attempt to remove potentially corrupted output file
        if os.path.isfile(outputfile):
            try:
                os.remove(outputfile)
            except OSError:
                pass # Ignore error if removal fails
        continue # Move to the next PDB ID

print("--- Batch Processing Summary ---")
print(f"Successfully processed/skipped: {success_count}")
print(f"Failed: {fail_count}")
print("--------------------------------")

if fail_count > 0:
    sys.exit(1)
else:
    sys.exit(0)

