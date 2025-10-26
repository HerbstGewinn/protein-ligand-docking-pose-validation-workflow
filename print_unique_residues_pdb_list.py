#!/usr/bin/env python3

from openmm.app import PDBFile
import sys
import os
import csv # Added for CSV writing

# Standard residues including water
standardRes = {'HOH', 'ALA', 'SER', 'ILE', 'MET', 'TYR', 'PHE', 'LYS', 'HIS',
               'CYS', 'LEU', 'TRP', 'GLY', 'VAL', 'THR', 'ASN', 'ASP', 'GLN',
               'PRO', 'GLU', 'ARG'}

# Set to store all unique non-standard residues found across all files (for summary)
all_unique_non_standard_res_global_summary = set()

# List to store results for CSV output (list of dictionaries)
csv_output_data = []

# Check if the correct number of arguments is provided
if len(sys.argv) != 3: # Expects 3 arguments: script_name, id_list_file, pdb_dir
    print("Usage: python printUniqueResidues.py <pdb_id_list_file.txt> <pdb_files_directory>")
    print("  <pdb_id_list_file.txt>: A text file containing PDB IDs, one per line.")
    print("  <pdb_files_directory>: The directory containing the .pdb structure files.")
    sys.exit(1)

pdb_id_list_filename = sys.argv[1]
pdb_files_directory = sys.argv[2] # New argument for PDB files directory

# Check if the PDB ID list file exists
if not os.path.isfile(pdb_id_list_filename):
    print(f"Error: PDB ID list file not found: {pdb_id_list_filename}")
    sys.exit(1)

# Check if the PDB files directory exists
if not os.path.isdir(pdb_files_directory):
    print(f"Error: PDB files directory not found: {pdb_files_directory}")
    sys.exit(1)

print(f"Reading PDB IDs from: {pdb_id_list_filename}")
try:
    with open(pdb_id_list_filename, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    if not pdb_ids:
        print(f"Error: No PDB IDs found in {pdb_id_list_filename}")
        sys.exit(1)
    print(f"Found {len(pdb_ids)} PDB IDs to process.")
    print(f"Looking for PDB files in directory: {pdb_files_directory}")
except Exception as e:
    print(f"Error reading PDB ID list file: {e}")
    sys.exit(1)

# Process each PDB ID
for pdb_id in pdb_ids:
    # Construct the full PDB filename including the directory path
    # Assumes .pdb extension. Adjust if your files are e.g. PDBID_bio.pdb
    pdb_filename_base = f"{pdb_id}_prepared.pdb"
    pdb_filename_full_path = os.path.join(pdb_files_directory, pdb_filename_base)

    print(f"--- Processing: {pdb_filename_full_path} ---")

    # Data for the current PDB ID for CSV
    current_pdb_data = {"PDB_ID": pdb_id, "NonStandardResidues": ""}

    if not os.path.isfile(pdb_filename_full_path):
        print(f"  Warning: PDB file '{pdb_filename_full_path}' not found. Skipping.")
        current_pdb_data["NonStandardResidues"] = "PDB_FILE_NOT_FOUND"
        csv_output_data.append(current_pdb_data)
        continue

    try:
        # Read the PDB file
        pdb = PDBFile(pdb_filename_full_path)
        topology = pdb.getTopology()

        # Check if positions were loaded (basic check for valid PDB for OpenMM)
        if pdb.getNumFrames() == 0 or pdb.getPositions() is None:
            print(f"  Warning: No coordinate frames or positions found in '{pdb_filename_full_path}'. Skipping.")
            current_pdb_data["NonStandardResidues"] = "NO_COORDINATES_IN_PDB"
            csv_output_data.append(current_pdb_data)
            continue

        current_file_unique_res = set()
        for res in topology.residues():
            current_file_unique_res.add(res.name)

        non_standard_in_file = []
        for res_name in current_file_unique_res:
            if res_name not in standardRes:
                non_standard_in_file.append(res_name)
                all_unique_non_standard_res_global_summary.add(res_name) # Add to global set

        if non_standard_in_file:
            sorted_non_standard = sorted(non_standard_in_file)
            print(f"  Non-standard residues found in {pdb_filename_base}: {', '.join(sorted_non_standard)}")
            current_pdb_data["NonStandardResidues"] = ','.join(sorted_non_standard)
        else:
            print(f"  No non-standard residues found in {pdb_filename_base}.")
            # Keep NonStandardResidues as empty string for CSV if none found

        csv_output_data.append(current_pdb_data)

    except Exception as e:
        print(f"  Error processing PDB file {pdb_filename_base}: {e}")
        current_pdb_data["NonStandardResidues"] = f"ERROR_PROCESSING: {str(e)[:50]}" # Truncate long errors
        csv_output_data.append(current_pdb_data)
        continue

# --- Write CSV Output ---
output_csv_filename = "non_standard_residues_report.csv"
print(f"\nWriting results to CSV file: {output_csv_filename}")
try:
    with open(output_csv_filename, 'w', newline='') as csvfile:
        fieldnames = ['PDB_ID', 'NonStandardResidues']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_output_data)
    print(f"Successfully wrote {output_csv_filename}")
except Exception as e:
    print(f"Error writing CSV file: {e}")

# --- Print Global Summary (as before) ---
print("\n--- Global Summary of All Unique Non-Standard Residues Found ---")
if all_unique_non_standard_res_global_summary:
    for res_name in sorted(list(all_unique_non_standard_res_global_summary)):
        print(res_name)
else:
    print("No non-standard residues found in any of the successfully processed PDB files (or no PDB files were successfully processed).")

print("\nScript finished.")
