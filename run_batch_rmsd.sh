#!/bin/bash

# --- Configuration ---
# Define the name of the text file containing PDB IDs (one per line)
PDB_ID_LIST_FILE="pdb_ids.txt"

# Define filename patterns/suffixes based on previous steps
# Input Files for RMSD calculation:
REF_LIGAND_SUFFIX="_Ligand.maegz" # Reference ligand (extracted & prepped)
# DOCKED_POSES_SUFFIX="_dock_SP_dock_pv.maegz"   # Original suffix
DOCKED_POSES_SUFFIX="_dock_SP_dock_pv.maegz" # <-- UPDATED suffix based on user input

# Output File for RMSD results:
RMSD_CSV_SUFFIX="_rmsd.csv"

# Define rmsd.py Options
RMSD_ASL="NOT atom.ele H" # ASL for RMSD calculation (e.g., heavy atoms)
PV_OPTION="-pv second"    # Treat the second input file as a pose viewer file

# Load Schrödinger environment module (adjust if necessary)
echo "Loading Schrodinger environment..."
ml Schrodinger_Suites
ml # Optional: Display loaded modules

# Check if SCHRODINGER variable is set
if [ -z "$SCHRODINGER" ]; then
  echo "Error: SCHRODINGER environment variable is not set."
  exit 1
fi
echo "Schrodinger base directory: $SCHRODINGER"

# Check if PDB ID list file exists
if [ ! -f "$PDB_ID_LIST_FILE" ]; then
    echo "Error: Input PDB ID list file '$PDB_ID_LIST_FILE' not found."
    exit 1
fi

echo "Starting Batch RMSD Calculations from $PDB_ID_LIST_FILE..."
echo "-----------------------------------------------------"

# --- Read PDB ID list and Launch Jobs ---
# Reads the PDB ID file line by line
while IFS= read -r pdb_id || [[ -n "$pdb_id" ]]; do
    # Skip empty lines
    if [ -z "$pdb_id" ]; then
        continue
    fi

    echo "--- Processing PDB ID: $pdb_id ---"

    # --- Construct File Names ---
    REF_FILE="${pdb_id}${REF_LIGAND_SUFFIX}"
    POSE_FILE="${pdb_id}${DOCKED_POSES_SUFFIX}"
    OUTPUT_CSV="${pdb_id}${RMSD_CSV_SUFFIX}"

    echo "  Reference File: $REF_FILE"
    echo "  Pose File:      $POSE_FILE"
    echo "  Output CSV:     $OUTPUT_CSV"

    # --- Check if input files exist ---
    if [ ! -f "$REF_FILE" ]; then
        echo "  [SKIP] Error: Reference ligand file '$REF_FILE' not found."
        echo "-----------------------------------------------------"
        continue
    fi
    if [ ! -f "$POSE_FILE" ]; then
        echo "  [SKIP] Error: Docked poses file '$POSE_FILE' not found."
        echo "-----------------------------------------------------"
        continue
    fi

    # --- Check if output CSV already exists ---
    if [ -f "$OUTPUT_CSV" ]; then
        echo "  [SKIP] Output CSV file '$OUTPUT_CSV' already exists."
        echo "-----------------------------------------------------"
        continue
    fi

    # --- Run rmsd.py ---
    echo "  Running RMSD calculation..."
    # --- FIX: Removed quotes around option variables ---
    # The shell handles splitting $PV_OPTION into '-pv' and 'second'
    "$SCHRODINGER/run" rmsd.py \
        $PV_OPTION \
        -a "$RMSD_ASL" \
        -c "$OUTPUT_CSV" \
        "$REF_FILE" \
        "$POSE_FILE"
    # --------------------------------------------------

    # --- Check rmsd.py Exit Status ---
    EXIT_CODE=$? # Capture exit status immediately
    if [ $EXIT_CODE -eq 0 ] && [ -f "$OUTPUT_CSV" ]; then
        echo "  [SUCCESS] Successfully completed RMSD calculation for $pdb_id. Output: $OUTPUT_CSV"
    else
        echo "  [FAIL] Error: RMSD calculation for '$pdb_id' failed with exit code $EXIT_CODE."
        # rmsd.py often prints errors to stderr, check terminal output
        # Optionally remove potentially empty/incomplete CSV file
        # rm -f "$OUTPUT_CSV"
        # Decide if you want the script to stop on failure or continue
        # exit $EXIT_CODE # Uncomment to stop on first failure
    fi

    echo "-----------------------------------------------------"

done < "$PDB_ID_LIST_FILE" # Feed the file into the while loop

echo "All RMSD calculations submitted or checked."
exit 0
