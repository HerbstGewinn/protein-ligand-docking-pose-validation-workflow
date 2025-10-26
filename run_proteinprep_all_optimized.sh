#!/bin/bash

# Set Schrodinger path
# SCHRODINGER="${SCHRODINGER}/utilities/prepwizard" # This line is commented out, assuming SCHRODINGER is correctly set by the environment

# Load Schrödinger environment module (adjust if necessary for your system)
ml Schrodinger_Suites
ml # Display loaded modules (optional)

# Define the input file containing PDB IDs
PDB_ID_FILE="pdb_ids.txt"

# Check if the PDB ID file exists and is readable
if [ ! -f "$PDB_ID_FILE" ]; then
    echo "Error: PDB ID file '$PDB_ID_FILE' not found in the current directory."
    exit 1
fi
if [ ! -r "$PDB_ID_FILE" ]; then
    echo "Error: PDB ID file '$PDB_ID_FILE' is not readable."
    exit 1
fi

echo "Press Enter to start processing PDB IDs from $PDB_ID_FILE..."
read Y # Pause until Enter is pressed

# Loop through each PDB ID read from the file
while IFS= read -r PDB_ID || [[ -n "$PDB_ID" ]]; do
    # Skip empty lines or lines with only whitespace
    if [[ -z "${PDB_ID// }" ]]; then
        continue
    fi

    echo "Checking $PDB_ID..."

    # Define expected output filename
    OUTPUT_MAE="${PDB_ID}_prepared.maegz"

    # --- START: Added Check ---
    # Check if the output file already exists and is a regular file
    if [ -f "$OUTPUT_MAE" ]; then
        echo "${OUTPUT_MAE} already exists. Skipping $PDB_ID."
        echo "---------------------------------------------"
        continue # Skip to the next PDB_ID in the loop
    fi
    # --- END: Added Check ---

    echo "Processing $PDB_ID..."

    # Download from RCSB
    INPUT_PDB="${PDB_ID}.pdb"
    # Biological Unit should be downloaded by adding a 1 after .pdb
    # Using -q for quiet download, consider removing if you want to see wget's output
    wget -q "https://files.rcsb.org/download/${PDB_ID}.pdb1" -O "${INPUT_PDB}"

    # Check if download was successful (file exists and is not empty)
    if [ ! -s "${INPUT_PDB}" ]; then
        echo "Failed to download biological unit for ${PDB_ID} (as ${PDB_ID}.pdb1). Attempting standard PDB download."
        # Try downloading the standard PDB file if biological unit fails
        wget -q "https://files.rcsb.org/download/${PDB_ID}.pdb" -O "${INPUT_PDB}"
        if [ ! -s "${INPUT_PDB}" ]; then
            echo "Failed to download standard PDB for ${PDB_ID} as well. Skipping."
            # Optionally remove the empty file if wget created one on failure
            rm -f "${INPUT_PDB}"
            echo "---------------------------------------------"
            continue
        else
            echo "Successfully downloaded standard PDB for ${PDB_ID}."
        fi
    else
        echo "Successfully downloaded biological unit for ${PDB_ID}."
    fi


    echo "Schrodinger base directory:"
    echo "$SCHRODINGER"

    # Run PrepWizard with the same options
    # Ensure SCHRODINGER variable is correctly set by the 'ml' command
    "${SCHRODINGER}/utilities/prepwizard" \
        "${INPUT_PDB}" "${OUTPUT_MAE}" -WAIT \
        -fillsidechains -disulfides -captermini -assign_all_residues \
        -rehtreat -max_states 1 -epik_pH 7.4 -epik_pHt 1.0 \
        -antibody_cdr_scheme Kabat -samplewater -include_epik_states \
        -propka_pH 7.4 -f S-OPLS -rmsd 0.3 -keepfarwat \
        -JOBNAME "${PDB_ID}_prep"

    # Check if prepwizard was successful
    EXIT_STATUS=$?
    if [ $EXIT_STATUS -eq 0 ] && [ -s "$OUTPUT_MAE" ]; then
        echo "Successfully finished $PDB_ID"
        # Optionally remove the downloaded PDB file now that preparation is done
        # rm -f "${INPUT_PDB}"
    else
        echo "Error processing $PDB_ID with PrepWizard. Exit status: $EXIT_STATUS"
        # Optionally remove potentially corrupted output file
        # rm -f "$OUTPUT_MAE"
    fi

    echo "---------------------------------------------"
done < "$PDB_ID_FILE" # Redirect file content to the while loop

echo "All processing attempts finished!"
