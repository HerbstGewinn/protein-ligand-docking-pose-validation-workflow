import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def analyze_first_high_rmsd_files(directory_path):
    """
    Analyzes CSV files in a directory. For each file, it finds the first
    RMSD value chronologically that is >= 2.5Å. This value is then
    categorized, and the script plots the distribution of files per category.
    """
    # Counters for the 5 categories (representing files)
    # Note: The first 3 categories will always be 0 with this logic.
    cat1_correct_lt_1_5_count = 0
    cat2_correct_1_5_to_2_0_count = 0
    cat3_correct_2_0_to_2_5_count = 0
    cat4_partially_correct_2_5_to_4_0_count = 0
    cat5_wrong_gte_4_0_count = 0

    files_categorized = 0
    files_with_no_high_rmsd = 0
    files_with_issues = 0
    issue_details = []

    print(f"Searching for CSV files in directory: {directory_path}\n")

    try:
        all_items_in_dir = os.listdir(directory_path)
        csv_files = [f for f in all_items_in_dir if f.lower().endswith('_rmsd.csv')]

        if not csv_files:
            print(f"No '*_rmsd.csv' files found in '{directory_path}'.")
            print(f"Falling back to search for any '.csv' files...\n")
            csv_files = [f for f in all_items_in_dir if f.lower().endswith('.csv')]
            if not csv_files:
                print(f"No '.csv' files found in '{directory_path}' either. Nothing to process.")
                return
            else:
                print(f"Found {len(csv_files)} generic '.csv' file(s) to process.")
        else:
            print(f"Found {len(csv_files)} '*_rmsd.csv' file(s) to process.")

    except FileNotFoundError:
        print(f"[ERROR] Directory not found: {directory_path}")
        return
    except Exception as e:
        print(f"[ERROR] Could not list files in directory {directory_path}: {e}")
        return

    for filename in csv_files:
        filepath = os.path.join(directory_path, filename)
        print(f"Processing: {filename}...")
        try:
            df = pd.read_csv(filepath, quotechar='"', skipinitialspace=True)

            if 'RMSD' not in df.columns:
                msg = f"  Skipped: 'RMSD' column not found."
                print(msg)
                files_with_issues += 1
                issue_details.append(f"{filename}: {msg.strip()}")
                continue

            df['RMSD'] = pd.to_numeric(df['RMSD'], errors='coerce')
            df.dropna(subset=['RMSD'], inplace=True)

            if df.empty:
                msg = f"  Skipped: No valid RMSD data after cleaning."
                print(msg)
                files_with_issues += 1
                issue_details.append(f"{filename}: {msg.strip()}")
                continue

            # --- NEW RMSD SELECTION LOGIC ---
            # Filter the DataFrame to find rows where RMSD is >= 2.5
            high_rmsd_df = df[df['RMSD'] >= 2.5]

            if high_rmsd_df.empty:
                print(f"  No RMSD values >= 2.5Å found in this file. Skipping categorization.")
                files_with_no_high_rmsd += 1
                continue

            # Get the first RMSD value from the filtered DataFrame
            first_high_rmsd = high_rmsd_df['RMSD'].iloc[0]
            print(f"  First RMSD >= 2.5Å found: {first_high_rmsd:.4f} Å")

            # --- APPLY 5-CATEGORY CLASSIFICATION ---
            if first_high_rmsd < 1.5:
                cat1_correct_lt_1_5_count += 1
            elif 1.5 <= first_high_rmsd < 2.0:
                cat2_correct_1_5_to_2_0_count += 1
            elif 2.0 <= first_high_rmsd < 2.5:
                cat3_correct_2_0_to_2_5_count += 1
            elif 2.5 <= first_high_rmsd < 4.0:
                cat4_partially_correct_2_5_to_4_0_count += 1
            else: # first_high_rmsd >= 4.0
                cat5_wrong_gte_4_0_count += 1
            # --- END OF CLASSIFICATION ---

            files_categorized += 1

        except pd.errors.EmptyDataError:
            msg = f"  Skipped: File is empty."
            print(msg)
            files_with_issues += 1
            issue_details.append(f"{filename}: {msg.strip()}")
        except Exception as e:
            msg = f"  Skipped: Error processing file - {e}"
            print(msg)
            files_with_issues += 1
            issue_details.append(f"{filename}: Error - {e}")

    print(f"\n--- Analysis Summary ---")
    print(f"Total CSV files targeted: {len(csv_files)}")
    print(f"Files Categorized (had an RMSD >= 2.5Å): {files_categorized}")
    print(f"Files Skipped (no RMSD >= 2.5Å found): {files_with_no_high_rmsd}")
    print(f"Files with Other Issues (e.g., no RMSD column): {files_with_issues}")
    if issue_details:
        print("Details for files with other issues:")
        for detail in issue_details:
            print(f"  - {detail}")

    if files_categorized == 0:
        print("\nNo data to plot as no files contained an RMSD value of 2.5Å or greater.")
        return

    # --- Plotting ---
    categories = [
        'Correct\n(RMSD < 1.5Å)',
        'Correct\n(1.5Å ≤ RMSD < 2Å)',
        'Correct\n(2Å ≤ RMSD < 2.5Å)',
        'Partially Correct\n(2.5Å ≤ RMSD < 4Å)',
        'Wrong\n(RMSD ≥ 4Å)'
    ]
    counts = [
        cat1_correct_lt_1_5_count,
        cat2_correct_1_5_to_2_0_count,
        cat3_correct_2_0_to_2_5_count,
        cat4_partially_correct_2_5_to_4_0_count,
        cat5_wrong_gte_4_0_count
    ]

    sns.set_theme(style="whitegrid", palette="flare")
    plt.figure(figsize=(14, 8))

    barplot = sns.barplot(x=categories, y=counts)

    plt.title('Distribution of First Encountered "Bad" Poses (RMSD ≥ 2.5Å)', fontsize=18, fontweight='bold', pad=20)
    plt.xlabel('RMSD Category', fontsize=15, labelpad=15)
    plt.ylabel('Number of Files', fontsize=15, labelpad=15)

    max_count_val = max(counts) if counts else 0
    for i, count_val in enumerate(counts):
        text_y_offset = max_count_val * 0.015
        if count_val == 0:
            if max_count_val == 0:
                text_y_offset = 0.1

        text_display_position = count_val + text_y_offset if count_val > 0 else text_y_offset
        plt.text(i, text_display_position, str(count_val),
                 ha='center', va='bottom', fontsize=11, fontweight='bold', color='black')

    plt.ylim(0, max_count_val * 1.1 + 0.5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=11)
    sns.despine(left=True, bottom=False)
    plt.tight_layout()

    plot_filename = "first_high_rmsd_distribution_plot.png"
    try:
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"\nPlot successfully saved as {plot_filename}")
    except Exception as e:
        print(f"\n[ERROR] Could not save plot: {e}")

    plt.show()


if __name__ == '__main__':
    data_directory = '.'
    analyze_first_high_rmsd_files(data_directory)
