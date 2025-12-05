import pandas as pd
import random

# ----------------------------------------------
# User-defined settings
# ----------------------------------------------
selection_targets = {
    "fully_susceptible": 3,  # number of susceptible samples to pick
    "mono_resistant": 2,     # number of mono-resistant samples
    "poly_resistant": 5      # number of poly-resistant samples
}

# Include second-line drugs? (set True/False) - for later
include_second_line = False

# ----------------------------------------------
# Load CSV file
# ----------------------------------------------
df = pd.read_csv("./tbtamr/validation_results_tbtamr_master.csv")

# Ensure Sample_ID is consistent and set to Isolate
df["Sample_ID"] = df["Isolate"]

# ----------------------------------------------
# Keep only Current study
# ----------------------------------------------
df = df[df['Study'] == 'Current study']

# ----------------------------------------------
# First-line drugs
# ----------------------------------------------
first_line_drugs = [
    ("Rifampicin", "Rifampicin_result"),
    ("Isoniazid", "Isoniazid_result"),
    ("Pyrazinamide", "Pyrazinamide_result"),
    ("Ethambutol", "Ethambutol_result")
]

# Optional second-line drugs
second_line_drugs = [
    ("Moxifloxacin", "Moxifloxacin_result"), # to be included if include_second_line is True - for later
    ("Streptomycin", "Streptomycin_result"),
    ("Amikacin", "Amikacin_result"),
    ("Kanamycin", "Kanamycin_result"),
    ("Cycloserine", "Cycloserine_result"),
    ("Capreomycin", "Capreomycin_result"),
    ("Ethionamide", "Ethionamide_result"),
    ("Para-aminosalicylic acid", "Para-aminosalicylic acid_result")
]

# Combine if needed
drug_cols = first_line_drugs + second_line_drugs if include_second_line else first_line_drugs

# ----------------------------------------------
# Keep only True Positive or True Negative
# ----------------------------------------------
for drug, result_col in drug_cols:
    df = df[df[result_col].isin(["TP", "TN"])]

# ----------------------------------------------
# Count resistant drugs 
# ----------------------------------------------
df["AMR_first_line"] = df[[r for _, r in first_line_drugs]].apply(lambda row: (row == "TP").sum(), axis=1)

if include_second_line:
    df["AMR_second_line"] = df[[r for _, r in second_line_drugs]].apply(lambda row: (row == "TP").sum(), axis=1)
else:
    df["AMR_second_line"] = 0  # placeholder for now

# Total AMR count
df["AMR_count"] = df["AMR_first_line"] + df["AMR_second_line"]

# ----------------------------------------------
# Stratify into complexity levels
# ----------------------------------------------
fully_susceptible = df[df["AMR_count"] == 0]
mono_resistant = df[df["AMR_count"] == 1]
poly_resistant = df[df["AMR_count"] >= 2]

# ----------------------------------------------
# Random selection 
# ----------------------------------------------
def pick_random(df_group, n, label):
    if len(df_group) >= n:
        chosen = df_group.sample(n)
    else:
        chosen = df_group
    chosen = chosen.copy()
    chosen["Selection_reason"] = label
    return chosen

# ----------------------------------------------
# Pick samples
# ----------------------------------------------
samples = []
samples.append(pick_random(fully_susceptible, selection_targets["fully_susceptible"], "fully_susceptible"))
samples.append(pick_random(mono_resistant, selection_targets["mono_resistant"], "mono_resistant"))
samples.append(pick_random(poly_resistant, selection_targets["poly_resistant"], "poly_resistant"))

truth_set = pd.concat(samples).drop_duplicates()

# ----------------------------------------------
# Final random shuffle
# ----------------------------------------------
truth_set = truth_set.sample(frac=1).reset_index(drop=True)

# ----------------------------------------------
# Optional: Include mutation info or tbtamr outputs
# ----------------------------------------------
if "aaSubstitutions" in df.columns:
    truth_set["Mutations"] = truth_set["aaSubstitutions"]

# ----------------------------------------------
# Save outputs
# ----------------------------------------------
truth_set.to_csv("truth_set.csv", index=False)

# Generate a report explaining selection
report_lines = []
for _, row in truth_set.iterrows():
    reason = row["Selection_reason"]
    first_line_count = row["AMR_first_line"]
    second_line_count = row["AMR_second_line"]
    mutations = row.get("Mutations", "N/A")
    report_lines.append(
        f"Sample {row['Sample_ID']}: Selected as {reason}, "
        f"first-line resistant: {first_line_count}, "
        f"second-line resistant: {second_line_count}, mutations: {mutations}"
    )

with open("truth_set_report.txt", "w") as f:
    f.write("\n".join(report_lines))

print(f"Truth set created with {len(truth_set)} samples.")
print("Report saved to truth_set_report.txt")
