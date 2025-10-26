# Protein–Ligand Docking Pose Validation Workflow
**Author:** Laurin Herbst (HerbstGewinn)  
**Project:** Master’s Thesis – *Protein–Ligand Docking Pose Validation with MD Simulation*  
**Last updated:** October 2025  

---

## 🧭 Overview
This repository contains all scripts, batch jobs, and helper files used in the master’s thesis *“Protein–Ligand Docking Pose Validation with Molecular Dynamics Simulation.”*  
The workflow couples **Schrödinger Glide** docking with **OpenMM-based molecular dynamics (MD)** simulations to assess the stability of predicted binding poses.

The goal is to determine whether short (10 ns) MD simulations can distinguish *correct* from *decoy* docking poses based on physical stability metrics such as **ligand RMSD**, **backbone RMSD**, **ligand RMSF**, and **hydrogen-bond persistence**.

---

## 🧩 Repository Structure

| File / Script | Purpose |
|----------------|----------|
| `run_proteinprep_all_optimized.sh` | Automates protein preparation in Maestro/Schrödinger using the Protein Preparation Wizard. |
| `run_docking_batch.sh` | Executes standard-precision (SP) Glide docking jobs for all systems listed in `grid_jobs.csv`. |
| `run_docking_batch_advanced.sh` | Docking with advanced parameters or decoy generation. |
| `grid_jobs.csv` | Defines PDB IDs and grid-generation settings for each complex. |
| `extractLigand_ASL_batch.py` | Extracts ligand coordinates from prepared PDB files using ASL selections. |
| `print_unique_residues_pdb_list.py` | Identifies and lists non-standard residues across all PDB inputs. |
| `remove_water_batch.py` | Removes crystallographic waters (optional pre-docking cleanup). |
| `run_batch_rmsd.sh` | Runs automated RMSD analyses on MD trajectories via `pytraj_rmsd_all_compared_original.py`. |
| `pytraj_rmsd_all_compared_original.py` | Compares RMSDs between crystal, docked, and decoy poses over time using **PyTraj**. |
| `plot_rmsds_first_pose_more_than_2_5_A.py` | Plots RMSD evolution for decoy poses with initial deviation > 2.5 Å. |
| `MD_many_pdbs.py` | Launches OpenMM MD simulations for multiple prepared complexes (original crystal poses). |
| `MD_many_pdbs_docked.py` | Runs MD simulations for top docked poses (< 2.5 Å RMSD). |
| `MD_many_pdbs_docked_advanced.py` | MD runs for intentionally perturbed decoy poses (2.5–4 Å RMSD). |
| `environment.yml` | Conda environment file listing Python dependencies (OpenMM, MDTraj, PyTraj, NumPy, Pandas, Matplotlib, Seaborn). |
| `pdb_ids.txt` | List of all protein–ligand complex IDs used in the benchmark dataset. |
| `workflow_guides.md` | Step-by-step description of the full workflow: from docking → MD → analysis. |

---

## ⚙️ Requirements

- **Python ≥ 3.9**
- **Schrödinger Suite 2023** (for Glide docking and Maestro preparation)
- **OpenMM ≥ 8.0**
- **AmberTools / PyTraj**
- **MDTraj**, **NumPy**, **Pandas**, **Matplotlib**, **Seaborn**
- **GPU (NVIDIA RTX series)** recommended for production MD runs

Install the required packages with:
```bash
conda env create -f environment.yml
conda activate docking-md
