# Detailed Workflows

## Docking

1.  **Create List of `pdb_ids.txt`**
    *Create a text file named `pdb_ids.txt` with one PDB ID per line.*
    *Load the necessary module:*
    ```bash
    ml Schroedinger_Suites
    ```

2.  **Protein Preparation**
    *Run protein preparation for all PDBs listed in `pdb_ids.txt`. The script automatically uses the biological unit.*
    ```bash
    ./run_proteinprep_all_optimized.sh
    ```

3.  **Remove Water Molecules**
    ```bash
    $SCHRODINGER/run python3 remove_water_batch.py pdb_ids.txt
    ```

4.  **Ligand Extraction**
    ```bash
    $SCHRODINGER/run extractLigand_ASL_batch.py pdb_ids.txt
    ```

5.  **Print Unique Residues**
    *Activate the conda environment first.*
    ```bash
    ml Miniforge3
    conda activate molsim_openmm
    ./print_unique_residues_pdb_list.py pdb_ids.txt ./
    ```

6.  **Create Grid Jobs File**
    *Manually create the `grid_jobs.csv` file. You may need to filter out unique ligands from the list of unique residues generated in the previous step.*

7.  **Run Grid Generation**
    *Automated script from Thilo.*
    ```bash
    $SCHRODINGER/run generate_glide_grids.py grid_jobs.csv -forcefield OPLS4
    ```

8.  **Run Glide Docking**
    ```bash
    ./run_docking_batch.sh
    ```

9.  **Calculate RMSDs**
    ```bash
    ./run_batch_rmsd.sh
    ```

10. **Sort Files into PDB Folders**
    ```bash
    ./organize_files.sh
    ```

11. **Organize RMSDs**
    *Collect all RMSD result files into a single folder.*
    ```bash
    ./organize_rmsds.sh
    ```

12. **Plot RMSD Results**
    ```bash
    python combined_rmsds_plot.py
    ```

---

## Decoy Docking

1.  **Prepare for Docking**
    *This workflow utilizes the `pdb_ids.txt` file again, along with the prepared protein structures that have had water molecules removed.*

2.  **Run Glide Docking for Decoys**
    ```bash
    ./run_docking_batch_advanced.sh
    ```

3.  **Calculate RMSDs**
    ```bash
    ./run_batch_rmsd.sh
    ```

4.  **Sort Files into PDB Folders**
    ```bash
    ./organize_files.sh
    ```

5.  **Organize RMSDs**
    *Collect all RMSD result files into a single folder.*
    ```bash
    ./organize_rmsds.sh
    ```

6.  **Plot RMSD Results**
    ```bash
    python plot_rmsds_first_pose_more_than_2_5_A.py
    ```

---

## MD for Correct Crystal Structures (MD_correct_complex)

1.  **Prepare for MD**
    *This workflow utilizes the `pdb_ids.txt` file again, along with the prepared protein structures that still contain water molecules.*

2.  **Activate Conda Environment**
    ```bash
    ml miniforge3
    conda activate molsim_openmm
    ```

3.  **Convert File Formats**
    *Convert extracted `Ligands.maegz` files into `.sdf` files. Convert `(PDB_ID)_prepared.maegz` files into `.pdb` files.*
    ```bash
    # Use the Schroedinger structconvert utility for this step.
    # Example: $SCHRODINGER/utilities/structconvert -imae input.maegz -osdf output.sdf
    ```

4.  **Run Batch MD**
    ```bash
    python MD_many_pdbs.py
    ```

5.  **Plot Equilibration Data**
    ```bash
    python plot_equil_data.py
    ```

6.  **Sort Files into PDB Folders**
    ```bash
    ./organize_files.sh
    ```

7.  **Run MDtraj Analysis**
    *Note: This script is currently configured for a single, hardcoded structure.*
    ```bash
    python analyze_2B7D.py
    ```

---

## MD for Correct Docked Complex (MD_docked_complex)

1.  **Prepare for MD**
    *Utilize the `pdb_ids.txt` file.*

2.  **Activate Conda Environment**
    ```bash
    ml miniforge3
    conda activate molsim_openmm
    ```

3.  **Merge Ligand and Protein**
    *Merge the docked `.maegz` files (containing the ligand) with the protein structure to create a combined PDB file: `{pdb_id}_prepared_dock.pdb`. You will also need the `Ligand.sdf` files.*
    ***Note: This step is currently performed manually in Maestro.***

4.  **Run Batch MD for Docked Structures**
    ```bash
    python MD_many_pdbs_docked.py
    ```

5.  **Plot Equilibration Data**
    ```bash
    python plot_equil_data.py
    ```

6.  **Sort Files into PDB Folders**
    ```bash
    ./organize_files.sh
    ```

7.  **Run MDtraj Analysis**
    *Note: This script is currently configured for a single, hardcoded structure.*

---

## MD for Advanced Docked Decoy (MD_docked_decoy_advanced)

1.  **Prepare for MD**
    *Utilize the `pdb_ids.txt` file.*

2.  **Activate Conda Environment**
    ```bash
    ml miniforge3
    conda activate molsim_openmm
    ```

3.  **Merge Ligand and Protein**
    *Merge the docked `.maegz` files (containing the ligand) with the protein structure to create a combined PDB file: `{pdb_id}_prepared_dock_advanced.pdb`. You will also need the `Ligand.sdf` files.*
    ***Note: This step is currently performed manually in Maestro.***

4.  **Run Batch MD for Docked Structures**
    ```bash
    python MD_many_pdbs_docked.py
    ```

5.  **Plot Equilibration Data**
    ```bash
    python plot_equil_data.py
    ```

6.  **Sort Files into PDB Folders**
    ```bash
    ./organize_files.sh
    ```

7.  **Run MDtraj Analysis**
    *Note: This script is currently configured for a single, hardcoded structure.*
    ```bash
    python analyze_2B7D.py
    ```

