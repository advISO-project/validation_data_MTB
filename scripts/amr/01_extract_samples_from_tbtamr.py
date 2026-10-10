from pathlib import Path

import pandas as pd


# ============================================================
# Configuration
# ============================================================

INPUT_FILE = Path("./tbtamr/validation_results_tbtamr_master.csv")
OUTPUT_DIR = Path("./tbtamr/reference_dataset")

STUDY_FILTER = "Current study"

# Second-line drugs can be enabled in a later release.
INCLUDE_SECOND_LINE = False

# Fixed seed makes representative sample selection reproducible.
RANDOM_STATE = 42

SELECTION_TARGETS = {
    "fully_susceptible": 3,
    "mono_resistant": 2,
    "poly_resistant": 5,
}

# When True, stop if a category has fewer eligible samples than requested.
# When False, use all available samples and print a warning.
STRICT_SELECTION_TARGETS = False


# ============================================================
# Drug definitions
#
# dst:
#   Phenotypic drug susceptibility testing result.
#   This is retained as the reference resistance status.
#
# mechanism:
#   Resistance mechanism reported by tbtAMR.
#
# result:
#   TP/TN/FP/FN classification reported in the source dataset.
#   TP and TN are used only as the high-confidence inclusion filter.
# ============================================================

FIRST_LINE_DRUGS = {
    "rifampicin": {
        "display_name": "Rifampicin",
        "dst": "rifampicin - dst",
        "mechanism": "Rifampicin",
        "result": "Rifampicin_result",
    },
    "isoniazid": {
        "display_name": "Isoniazid",
        "dst": "isoniazid - dst",
        "mechanism": "Isoniazid",
        "result": "Isoniazid_result",
    },
    "pyrazinamide": {
        "display_name": "Pyrazinamide",
        "dst": "pyrazinamide - dst",
        "mechanism": "Pyrazinamide",
        "result": "Pyrazinamide_result",
    },
    "ethambutol": {
        "display_name": "Ethambutol",
        "dst": "ethambutol - dst",
        "mechanism": "Ethambutol",
        "result": "Ethambutol_result",
    },
}

SECOND_LINE_DRUGS = {
    "moxifloxacin": {
        "display_name": "Moxifloxacin",
        "dst": "moxifloxacin - dst",
        "mechanism": "Moxifloxacin",
        "result": "Moxifloxacin_result",
    },
    "amikacin": {
        "display_name": "Amikacin",
        "dst": "amikacin - dst",
        "mechanism": "Amikacin",
        "result": "Amikacin_result",
    },
    "kanamycin": {
        "display_name": "Kanamycin",
        "dst": "kanamycin - dst",
        "mechanism": "Kanamycin",
        "result": "Kanamycin_result",
    },
    "streptomycin": {
        "display_name": "Streptomycin",
        "dst": "streptomycin - dst",
        "mechanism": "Streptomycin",
        "result": "Streptomycin_result",
    },
    "cycloserine": {
        "display_name": "Cycloserine",
        "dst": "cycloserine - dst",
        "mechanism": "Cycloserine",
        "result": "Cycloserine_result",
    },
    "capreomycin": {
        "display_name": "Capreomycin",
        "dst": "capreomycin - dst",
        "mechanism": "Capreomycin",
        "result": "Capreomycin_result",
    },
    "ethionamide": {
        "display_name": "Ethionamide",
        "dst": "ethionamide - dst",
        "mechanism": "Ethionamide",
        "result": "Ethionamide_result",
    },
    "para_aminosalicylic_acid": {
        "display_name": "Para-aminosalicylic acid",
        "dst": "para-aminosalicylic_acid - dst",
        "mechanism": "Para-aminosalicylic acid",
        "result": "Para-aminosalicylic acid_result",
    },
}

DRUGS = dict(FIRST_LINE_DRUGS)

if INCLUDE_SECOND_LINE:
    DRUGS.update(SECOND_LINE_DRUGS)


# ============================================================
# Helper functions
# ============================================================

def validate_required_columns(
    dataframe: pd.DataFrame,
    drug_definitions: dict,
) -> None:
    """Raise a clear error if required source columns are absent."""

    required = {"Isolate", "Study"}

    for drug_info in drug_definitions.values():
        required.update(
            {
                drug_info["dst"],
                drug_info["mechanism"],
                drug_info["result"],
            }
        )

    missing = sorted(required.difference(dataframe.columns))

    if missing:
        raise ValueError(
            "The following required columns are missing from the input file:\n"
            + "\n".join(f"  - {column}" for column in missing)
        )


def clean_string_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Trim surrounding whitespace while preserving missing values."""

    cleaned = dataframe.copy()

    string_columns = cleaned.select_dtypes(
        include=["object", "string"]
    ).columns

    for column in string_columns:
        cleaned[column] = cleaned[column].apply(
            lambda value: value.strip()
            if isinstance(value, str)
            else value
        )

    return cleaned


def select_representative_samples(
    group: pd.DataFrame,
    target: int,
    label: str,
) -> pd.DataFrame:
    """Select a reproducible subset from one resistance category."""

    available = len(group)

    if available < target:
        message = (
            f"Category '{label}' requested {target} samples, "
            f"but only {available} eligible samples are available."
        )

        if STRICT_SELECTION_TARGETS:
            raise ValueError(message)

        print(f"WARNING: {message} Using all available samples.")
        selected = group.copy()
    else:
        selected = group.sample(
            n=target,
            random_state=RANDOM_STATE,
        ).copy()

    selected["Selection_reason"] = label
    selected["Selection_target"] = target
    selected["Eligible_in_category"] = available

    return selected


def build_reference_output(
    dataframe: pd.DataFrame,
    drug_definitions: dict,
) -> pd.DataFrame:
    """
    Build the released wide-format reference table.

    One row represents one isolate.
    Phenotypic DST is retained as the expected resistance status.
    Mechanisms and source TP/TN labels are retained as provenance.
    """

    output = pd.DataFrame(
        {
            "Sample_ID": dataframe["Sample_ID"],
            "Study": dataframe["Study"],
            "Resistance_class": dataframe["Resistance_class"],
            "AMR_first_line_count": dataframe["AMR_first_line_count"],
            "AMR_second_line_count": dataframe["AMR_second_line_count"],
            "AMR_total_count": dataframe["AMR_total_count"],
        }
    )

    for column in [
        "Selection_reason",
        "Selection_target",
        "Eligible_in_category",
    ]:
        if column in dataframe.columns:
            output[column] = dataframe[column]

    for drug_key, drug_info in drug_definitions.items():
        output[f"{drug_key}_expected"] = dataframe[drug_info["dst"]]
        output[f"{drug_key}_mechanism"] = dataframe[
            drug_info["mechanism"]
        ]
        output[f"{drug_key}_source_result"] = dataframe[
            drug_info["result"]
        ]

    return output


# ============================================================
# Load and validate input
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE)
df = clean_string_columns(df)

validate_required_columns(df, DRUGS)

# Use the isolate accession as the stable sample identifier.
df["Sample_ID"] = df["Isolate"]


# ============================================================
# Restrict this release to the Current study cohort
# ============================================================

source_rows = len(df)

df = df.loc[df["Study"] == STUDY_FILTER].copy()

if df.empty:
    raise ValueError(
        f"No records were found where Study == '{STUDY_FILTER}'."
    )

current_study_rows = len(df)


# ============================================================
# Apply the high-confidence filter
#
# Phenotypic DST is the reference resistance status.
# TP/TN is used only as the concordance-based confidence filter.
#
# Because this is a whole-isolate first-line AMR panel, every
# included drug must have:
#   1. a valid DST result: R or S
#   2. a concordant source result: TP or TN
# ============================================================

valid_dst_mask = pd.Series(True, index=df.index)
concordant_mask = pd.Series(True, index=df.index)

for drug_info in DRUGS.values():
    valid_dst_mask &= df[drug_info["dst"]].isin({"R", "S"})
    concordant_mask &= df[drug_info["result"]].isin({"TP", "TN"})

df["Complete_valid_DST"] = valid_dst_mask
df["Complete_concordant_profile"] = concordant_mask

high_confidence = df.loc[
    df["Complete_valid_DST"]
    & df["Complete_concordant_profile"]
].copy()

if high_confidence.empty:
    raise ValueError(
        "No samples passed the complete high-confidence filter."
    )


# ============================================================
# Derive resistance counts directly from phenotypic DST
# ============================================================

first_line_dst_columns = [
    drug_info["dst"]
    for drug_info in FIRST_LINE_DRUGS.values()
]

high_confidence["AMR_first_line_count"] = (
    high_confidence[first_line_dst_columns]
    .eq("R")
    .sum(axis=1)
)

if INCLUDE_SECOND_LINE:
    second_line_dst_columns = [
        drug_info["dst"]
        for drug_info in SECOND_LINE_DRUGS.values()
    ]

    high_confidence["AMR_second_line_count"] = (
        high_confidence[second_line_dst_columns]
        .eq("R")
        .sum(axis=1)
    )
else:
    high_confidence["AMR_second_line_count"] = 0

high_confidence["AMR_total_count"] = (
    high_confidence["AMR_first_line_count"]
    + high_confidence["AMR_second_line_count"]
)


# ============================================================
# Assign panel-selection categories
#
# These are selection categories, not formal WHO resistance classes.
# ============================================================

high_confidence["Resistance_class"] = pd.cut(
    high_confidence["AMR_total_count"],
    bins=[-1, 0, 1, float("inf")],
    labels=[
        "fully_susceptible",
        "mono_resistant",
        "poly_resistant",
    ],
).astype("string")


# ============================================================
# Create representative panel
# ============================================================

selected_groups = []

for category, target in SELECTION_TARGETS.items():
    category_group = high_confidence.loc[
        high_confidence["Resistance_class"] == category
    ]

    selected_groups.append(
        select_representative_samples(
            group=category_group,
            target=target,
            label=category,
        )
    )

reference_panel = (
    pd.concat(selected_groups, ignore_index=True)
    .drop_duplicates(subset=["Sample_ID"])
    .sample(frac=1, random_state=RANDOM_STATE)
    .reset_index(drop=True)
)


# ============================================================
# Build the two released CSV outputs
# ============================================================

complete_reference = build_reference_output(
    high_confidence,
    DRUGS,
)

reference_panel_output = build_reference_output(
    reference_panel,
    DRUGS,
)


# ============================================================
# Save outputs
# ============================================================

complete_reference_file = (
    OUTPUT_DIR
    / "tbtamr_current_study_high_confidence.csv"
)

reference_panel_file = (
    OUTPUT_DIR
    / "tbtamr_current_study_reference_panel.csv"
)

report_file = (
    OUTPUT_DIR
    / "tbtamr_current_study_reference_panel_report.txt"
)

complete_reference.to_csv(
    complete_reference_file,
    index=False,
)

reference_panel_output.to_csv(
    reference_panel_file,
    index=False,
)


# ============================================================
# Generate curation and selection report
# ============================================================

class_counts = (
    high_confidence["Resistance_class"]
    .value_counts()
    .reindex(SELECTION_TARGETS.keys(), fill_value=0)
)

selected_class_counts = (
    reference_panel["Resistance_class"]
    .value_counts()
    .reindex(SELECTION_TARGETS.keys(), fill_value=0)
)

report_lines = [
    "tbtAMR Current Study High-Confidence AMR Reference Dataset",
    "=" * 65,
    "",
    "Scope",
    "-----",
    f"Input file: {INPUT_FILE}",
    f"Study filter: {STUDY_FILTER}",
    f"Second-line drugs included: {INCLUDE_SECOND_LINE}",
    f"Random seed: {RANDOM_STATE}",
    "",
    "Curation criteria",
    "-----------------",
    "1. Only isolates labelled 'Current study' were considered.",
    "2. Phenotypic DST was retained as the reference resistance status.",
    "3. Every included drug required a phenotypic result of R or S.",
    "4. Only TP and TN source comparisons were retained.",
    "5. TP/TN was used as a concordance-based confidence filter.",
    "6. Resistance counts were calculated directly from DST == R.",
    "7. Resistance mechanisms reported by tbtAMR were retained.",
    "",
    "Cohort counts",
    "-------------",
    f"Rows in source dataset: {source_rows}",
    f"Rows from Current study: {current_study_rows}",
    (
        "Rows with complete valid DST for included drugs: "
        f"{int(valid_dst_mask.sum())}"
    ),
    (
        "Rows with complete TP/TN concordance for included drugs: "
        f"{int(concordant_mask.sum())}"
    ),
    (
        "Rows passing both high-confidence criteria: "
        f"{len(high_confidence)}"
    ),
    "",
    "Eligible high-confidence samples by category",
    "--------------------------------------------",
]

for category, count in class_counts.items():
    report_lines.append(f"{category}: {count}")

report_lines.extend(
    [
        "",
        "Representative panel",
        "--------------------",
        f"Total selected samples: {len(reference_panel)}",
    ]
)

for category, count in selected_class_counts.items():
    target = SELECTION_TARGETS[category]
    report_lines.append(
        f"{category}: selected {count}; target {target}"
    )

report_lines.extend(
    [
        "",
        "Selected samples",
        "----------------",
    ]
)

for _, row in reference_panel.iterrows():
    drug_summary = []

    for drug_key, drug_info in DRUGS.items():
        phenotype = row[drug_info["dst"]]
        mechanism = row[drug_info["mechanism"]]
        concordance = row[drug_info["result"]]

        drug_summary.append(
            f"{drug_key}={phenotype} "
            f"[{concordance}; {mechanism}]"
        )

    report_lines.append(
        f"{row['Sample_ID']} | "
        f"{row['Resistance_class']} | "
        f"{row['AMR_total_count']} resistant drug(s) | "
        + "; ".join(drug_summary)
    )

report_file.write_text(
    "\n".join(report_lines),
    encoding="utf-8",
)


# ============================================================
# Console summary
# ============================================================

print()
print("Reference dataset generation complete.")
print(f"High-confidence cohort: {len(high_confidence)} samples")
print(f"Representative panel:   {len(reference_panel)} samples")
print()
print("Outputs:")
print(f"  {complete_reference_file}")
print(f"  {reference_panel_file}")
print(f"  {report_file}")