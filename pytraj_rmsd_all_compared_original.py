import pytraj as pt
import os
import numpy as np 

# --- Ligand mapping for reference structures ---
LIGAND_MAP = {
    "3FLY": "FLY", "2D3U": "CCT", "1AJQ": "SPA", "1E66": "HUX", "2QWE": "ZMR",
    "2CET": "PGI", "1TOI": "HCI", "1N2V": "BDI", "1S39": "AQO", "2BAJ": "PQB",
    "1NDZ": "FR5", "1XGJ": "XXX", "1ZOE": "K25", "1PXO": "RPN", "1SL3": "170",
    "1NVQ": "UCN", "2AOU": "CQA", "1PBQ": "DK1", "1B8O": "IMH", "1FTM": "AMQ",
    "1FCZ": "LMU", "1F4G": "TP4", "1SQA": "UI1", "2AYR": "L4G", "1IF7": "SBR",
    "1NFY": "RTR", "1MQ6": "XLD", "2BZ6": "346", "2B7D": "C1B", "1TRD": "PGH",
    "1UTP": "PBN", "1O3F": "696", "1Y6Q": "TDI", "1ELB": "0Z4", "1K9S": "FM2",
    "1FLR": "FLU", "10GS": "VWW", "3STD": "MQ0", "1ZVX": "FIN", "2DRC": "MTX",
    "2I0D": "MUT", "1RNT": "2GP", "1H1Q": "2A6", "2GMX": "877", "2QBS": "024",
    "2ZFF": "53U", "3FMH": "533", "3FMK": "FMK", "4DJW": "0KP", "4GIH": "0X5",
    "4HW3": "19G", "1S64": "778", "1SA4": "JAN", "1SA5": "BMV", "2XBV": "XBV",
    "2XBW": "455", "2XBX": "RR8", "2XBY": "63C", "2XC0": "8NC", "2XC4": "IVK",
    "2XC5": "OYJ", "3HQW": "PF4", "3HQY": "PF6", "3HQZ": "PF8", "3HR1": "PF9",
    "3KRJ": "KRJ", "3KRL": "KRL", "3LN0": "52B", "3LN1": "HEM", "3MNR": "SD1",
    "3MQE": "416", "3NTG": "HEM", "3SFK": "ORO", "3UW4": "0DQ", "3UW5": "0DQ",
    "3ZZE": "6XP", "3ZXZ": "KRW", "4AOI": "4K0", "4AP7": "F47", "4B05": "32D",
    "4C66": "H4C", "4C67": "L5S", "4CDD": "GDI", "4FLH": "14K", "4JH0": "1MD",
    "4LKO": "1WH", "4PHW": "2W1", "5AJR": "VT1", "5F4N": "5UY", "5HLS": "62G",
    "5HM0": "62V", "5K7K": "6RJ", "5MZC": "8EQ", "5MZI": "FYK", "5MZK": "OK1",
    "5OQ4": "A3W", "5PGU": "8K4", "5PGV": "8K7", "5PGW": "8KA", "5PGX": "8KD",
    "5PGY": "8KG", "5QCK": "BUV", "5QCL": "BUY", "5QCM": "BVJ", "5QCN": "BVM",
    "5QII": "HJG", "5QQO": "NRJ", "5QQP": "NR7", "5QTY": "QLJ", "5QTU": "QEV",
    "5TQ3": "7GZ", "5TQ4": "7GY", "5TQ5": "7GX", "5TQ6": "7GV", "5TQ7": "7GT",
    "5TQ8": "7GS", "5TX5": "7MJ", "5UIR": "8BY", "5UIS": "8C1", "5UIT": "8CD",
    "5UIU": "8CG", "5UIQ": "8BV", "5VP0": "9GJ", "5VP1": "9GA", "5V8Q": "97A",
    "5WHR": "AOJ", "6BBU": "D7D", "6BBV": "D7D", "6B0Z": "FLC", "6C0S": "EEJ",
    "6E0B": "JZ8", "6E3Z": "HRV", "6ENQ": "BJB", "6H3K": "7PE", "6H6Q": "FUK",
    "6H6R": "FUE", "6HKX": "GCE", "6HKY": "GCE", "6IWI": "B0C", "6JZ0": "CKO",
    "6O3I": "LKP", "6RLN": "K8K", "6RN8": "K9T", "6RNA": "KA2", "6TCZ": "N2E",
    "6TD5": "BO2", "6USX": "M1R", "6USZ": "QH4", "6UT0": "M1X", "6W35": "SKV",
    "6W50": "SWP", "6ZGC": "H8H", "7D9O": "H0L", "7D9P": "H0R", "7D9Q": "H1R",
    "7E9B": "UNK", "7MBO": "YXG", "7M7V": "Z3J", "7N4N": "0BK", "7N5O": "PGO",
    "7N5R": "0UW", "7N5X": "PGO", "7N5Y": "0CI", "7N66": "0EW", "7N91": "1RJ",
    "7N93": "1SK", "7NQQ": "UMN", "7NQW": "UNW", "7NR3": "UO5", "7NR5": "UOH",
    "7NR8": "UOE", "7NR9": "UOW", "7QJG": "EKR", "7QJU": "EKF", "7QK4": "EJR",
    "7QPF": "EFJ", "7QPM": "EEI", "7QPQ": "EIH", "7QPV": "EHI", "7QQ4": "EIK",
    "7RFR": "4W8", "7RFS": "4WI", "7RFU": "4YG", "7RFW": "4WI", "7T78": "G2T",
    "7T79": "G1S", "7VTH": "7XB", "7VU6": "7YY", "8C9W": "U30", "8CIC": "U30",
    "8GBU": "YWE", "8HTR": "MWH", "8HTS": "N2L", "8V4U": "YDL", "9C2R": "1PE"
}

# --- Read all PDB IDs from file ---
with open("pdb_ids.txt") as f:
    pdb_ids = [line.strip() for line in f if line.strip()]

def open_trajectories(pdbname, trajname, refname):
    traj = pt.load(trajname, top=pdbname)
    ref = pt.load(refname)
    return traj, ref

def calculate_rmsd(traj, ref, traj_ligand, ref_ligand):
    # Dynamically exclude ligand from backbone mask
    bbmask = f'(@C,N,O,CA)&(!:HOH,{ref_ligand},UNK)'
    traj_ligmask = f':{traj_ligand}&!@H='
    ref_ligmask  = f':{ref_ligand}&!@H='

    # Ligand selection info
    idx_lig_traj = pt.select_atoms(traj_ligmask, traj.top)
    idx_lig_ref  = pt.select_atoms(ref_ligmask, ref.top)
    print(f"   {len(idx_lig_traj)} ligand atoms in trajectory ({traj_ligand})")
    print(f"   {len(idx_lig_ref)} ligand atoms in reference ({ref_ligand})")

    if len(idx_lig_traj) == 0 or len(idx_lig_ref) == 0:
        print("⚠️  No ligand atoms found; skipping RMSD calculation.")
        return None, None

    pt.autoimage(traj, bbmask)

    protrmsd = pt.rmsd(traj=traj, mask=bbmask, ref=ref)
    ligrmsd_raw = pt.rmsd_nofit(traj=traj, mask=traj_ligmask, ref=ref, ref_mask=ref_ligmask)

    # Handle scalar vs array RMSD return
    if np.ndim(ligrmsd_raw) == 0:
        ligrmsd = [float(ligrmsd_raw)]
    else:
        ligrmsd = list(ligrmsd_raw)

    return list(protrmsd), ligrmsd


# --- MAIN LOOP ---
ROOT = os.getcwd()
print(f"Starting in: {ROOT}")

successful = 0
failed = 0

for pid in pdb_ids:
    subdir = os.path.join(ROOT, pid)
    print(f"\n🔹 Processing {pid} ...")

    # Check folder exists
    if not os.path.isdir(subdir):
        print(f"⚠️  Folder {pid} not found — skipping.")
        failed += 1
        continue

    # Move into subfolder
    os.chdir(subdir)

    # File paths inside this folder
    pdbname = f"{pid}_after1ns.pdb"
    refname = f"_{pid}_MD_out_start_correct_test.pdb"
    trajname = f"{pid}_MD_out_production_trajectory.dcd"
    output_csv = f"{pid}_rmsd_versus_correct.csv"

    # Check files exist here
    missing = [f for f in [pdbname, refname, trajname] if not os.path.exists(f)]
    if missing:
        print(f"⚠️  Missing file(s) for {pid}: {', '.join(missing)}")
        os.chdir(ROOT)
        failed += 1
        continue

    traj_ligand = "UNK"
    ref_ligand = LIGAND_MAP.get(pid, "UNK")

    try:
        traj, ref = open_trajectories(pdbname, trajname, refname)
        protrmsd, ligrmsd = calculate_rmsd(traj, ref, traj_ligand, ref_ligand)

        if protrmsd is None or ligrmsd is None:
            print(f"❌ Skipped {pid} (no ligand atoms)")
            os.chdir(ROOT)
            failed += 1
            continue

        with open(output_csv, 'w') as f:
            f.write('Snapshot,BBrmsd,LigRMSD\n')
            for i, rmsd in enumerate(zip(protrmsd, ligrmsd)):
                f.write(f"{i:3n},{rmsd[0]:6.2f},{rmsd[1]:6.2f}\n")

        print(f"✅ Saved {output_csv} in {subdir}")
        successful += 1

    except Exception as e:
        print(f"❌ Error processing {pid}: {e}")
        failed += 1

    finally:
        # Return to main directory before next loop
        os.chdir(ROOT)

print("\n===== SUMMARY =====")
print(f"✅ Successful: {successful}")
print(f"❌ Failed: {failed}")
print("===================")
