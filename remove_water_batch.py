#!/usr/bin/env python
# -*- coding: utf-8 -*-

# remove_water_batch.py
# Removes water molecules (ASL: 'water') from multiple prepared Maestro files,
# based on a list of PDB IDs provided in a text file.

import argparse
import os
import sys

# --- Configuration ---
# Define the suffix for input prepared complex files
INPUT_SUFFIX = "_prepared.maegz"
# Define the suffix for output files without water
OUTPUT_SUFFIX = "_prepared_nowater.maegz"
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
parser = argparse.ArgumentParser(description="Remove water molecules (ASL: 'water') from multiple prepared Maestro files.")
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
skip_count = 0

for pdb_id in pdb_ids:
    print(f"--- Processing PDB ID: {pdb_id} ---")
    # Construct filenames based on PDB ID
    inputfile = f"{pdb_id}{INPUT_SUFFIX}"
    outputfile = f"{pdb_id}{OUTPUT_SUFFIX}"

    # Check if input file exists
    if not os.path.isfile(inputfile):
        print(f"  Input file not found: {inputfile}. Skipping.")
        fail_count += 1
        continue

    # Check if output file already exists
    if os.path.isfile(outputfile):
        print(f"  Output file exists: {outputfile}. Skipping.")
        skip_count += 1
        continue

    # --- Main Logic for Water Removal (same as single file script) ---
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

        # Select water atoms using ASL
        asl = "water"
        print(f"  Selecting atoms with ASL: '{asl}'")
        water_indices = analyze.evaluate_asl(st, asl)

        if not water_indices:
            print(f"  No water atoms found in {inputfile}. Writing original structure to {outputfile}.")
            # Write the original structure if no water was found
            with StructureWriter(outputfile) as writer:
                 writer.append(st)
            success_count += 1 # Count as success as file is created
        else:
            print(f"  Found {len(water_indices)} water atoms to delete.")

            st.deleteAtoms(water_indices)
            # No st.repack() needed based on previous findings

            # Write the modified structure (without water) to the output file
            print(f"  Writing structure without water to: {outputfile}")
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

print("--- Batch Water Removal Summary ---")
print(f"Successfully processed: {success_count}")
print(f"Skipped (output existed): {skip_count}")
print(f"Failed: {fail_count}")
print("-----------------------------------")

if fail_count > 0:
    sys.exit(1)
else:
    sys.exit(0)
