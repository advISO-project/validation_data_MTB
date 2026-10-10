#!/usr/bin/env python3

"""
Build the contamination validation truth set from Heupink et al. (2021)
Supplementary Table 1 culture simulations.

The source sample identifiers encode the simulation truth, for example:

COV50-FASTA10-MCOV50-CONTMabs-CPER3-SEED123456

Fields:
    COV50       = dataset coverage group
    FASTA10     = simulated MTB genome
    MCOV50      = MTB coverage
    CONTMabs    = contaminant
    CPER3       = contamination percentage
    SEED...     = ART simulation seed

Only COV50 samples are retained for the initial contamination validation
dataset so that MTB coverage is held constant while contamination varies.

The source workbook is never modified.
"""

from pathlib import Path
import re
import sys

import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "metadata"
    / "heupink2021_supplementary_table1_culture_simulations.xlsx"
)

OUTPUT_DIR = PROJECT_ROOT / "contamination" / "reference_dataset"
REPORT_DIR = PROJECT_ROOT / "contamination" / "reports"

OUTPUT_FILE = OUTPUT_DIR / "contamination_truth_set.csv"
REPORT_FILE = REPORT_DIR / "contamination_truth_set_report.txt"


# ---------------------------------------------------------------------
# Study metadata
# ---------------------------------------------------------------------

SOURCE_STUDY = "Heupink et al. 2021"
BIOPROJECT = "PRJNA706121"

TARGET_COVERAGE = 50


# Contaminant codes used in the supplementary material.
CONTAMINANT_MAP = {
    "None": {
        "contaminant_species": "None",
        "contaminant_class": "none",
    },
    "Hsap": {
        "contaminant_species": "Homo sapiens",
        "contaminant_class": "human",
    },
    "Mabs": {
        "contaminant_species": "Mycobacterium abscessus",
        "contaminant_class": "NTM",
    },
    "Mavi": {
        "contaminant_species": "Mycobacterium avium",
        "contaminant_class": "NTM",
    },
    "Mint": {
        "contaminant_species": "Mycobacterium intracellulare",
        "contaminant_class": "NTM",
    },
    "Paer": {
        "contaminant_species": "Pseudomonas aeruginosa",
        "contaminant_class": "bacterial",
    },
    "Sepi": {
        "contaminant_species": "Staphylococcus epidermidis",
        "contaminant_class": "bacterial",
    },
    "Mix1": {
        "contaminant_species": "NTM mixture",
        "contaminant_class": "mixed_NTM",
    },
    "Mix2": {
        "contaminant_species": "all-contaminant mixture",
        "contaminant_class": "mixed",
    },
}


# ---------------------------------------------------------------------
# Sample-name parser
# ---------------------------------------------------------------------

SAMPLE_PATTERN = re.compile(
    r"^COV(?P<dataset_coverage>\d+)"
    r"-FASTA(?P<mtb_strain>\d+)"
    r"-MCOV(?P<mtb_coverage>\d+)"
    r"-CONT(?P<contaminant_code>[A-Za-z0-9]+)"
    r"-CPER(?P<contamination_percent>\d+)"
    r"-SEED(?P<simulation_seed>\d+)"
)


def parse_sample_id(sample_id):
    """Decode the simulation truth encoded in a sample identifier."""

    sample_id = str(sample_id).strip()

    match = SAMPLE_PATTERN.match(sample_id)

    if not match:
        return None

    parsed = match.groupdict()

    code = parsed["contaminant_code"]

    if code not in CONTAMINANT_MAP:
        raise ValueError(
            f"Unknown contaminant code '{code}' in sample: {sample_id}"
        )

    contaminant = CONTAMINANT_MAP[code]

    return {
        "source_sample_id": sample_id,
        "simulated_mtb_strain": f"FASTA{parsed['mtb_strain']}",
        "dataset_coverage": int(parsed["dataset_coverage"]),
        "mtb_coverage": int(parsed["mtb_coverage"]),
        "contaminant_code": code,
        "contaminant_species": contaminant["contaminant_species"],
        "contaminant_class": contaminant["contaminant_class"],
        "contamination_percent": int(parsed["contamination_percent"]),
        "simulation_seed": parsed["simulation_seed"],
        "source_study": SOURCE_STUDY,
        "bioproject": BIOPROJECT,
    }


# ---------------------------------------------------------------------
# Extract sample identifiers
# ---------------------------------------------------------------------

def extract_sample_ids(workbook):
    """
    Search all cells in all worksheets for sample identifiers.

    This avoids assuming a particular column name or worksheet layout in
    the supplementary workbook.
    """

    excel = pd.ExcelFile(workbook)

    sample_ids = set()

    for sheet_name in excel.sheet_names:

        df = pd.read_excel(
            workbook,
            sheet_name=sheet_name,
            header=None,
            dtype=str,
        )

        for value in df.to_numpy().ravel():

            if pd.isna(value):
                continue

            value = str(value).strip()

            match = SAMPLE_PATTERN.match(value)

            if match:
                sample_ids.add(match.group(0))

    return sorted(sample_ids)


# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------

def validate_truth_set(df):
    """Run basic integrity checks on the generated truth set."""

    problems = []

    if df.empty:
        problems.append("No COV50 samples were found.")

    duplicated = df["source_sample_id"].duplicated()

    if duplicated.any():
        problems.append(
            f"{duplicated.sum()} duplicate sample identifiers detected."
        )

    if not (df["dataset_coverage"] == TARGET_COVERAGE).all():
        problems.append("Non-COV50 samples are present.")

    if not (df["mtb_coverage"] == TARGET_COVERAGE).all():
        problems.append("Samples with MTB coverage other than 50x are present.")

    valid_percentages = {0, 1, 2, 3, 4, 5}

    observed_percentages = set(df["contamination_percent"])

    unexpected = observed_percentages - valid_percentages

    if unexpected:
        problems.append(
            f"Unexpected contamination percentages: {sorted(unexpected)}"
        )

    # A clean control should have no contaminant and 0% contamination.
    controls = df[df["contamination_percent"] == 0]

    if not controls.empty:
        invalid_controls = controls[
            controls["contaminant_code"] != "None"
        ]

        if not invalid_controls.empty:
            problems.append(
                "One or more 0% samples have a contaminant code other than None."
            )

    # Conversely, CONTNone should correspond to 0%.
    none_samples = df[df["contaminant_code"] == "None"]

    invalid_none = none_samples[
        none_samples["contamination_percent"] != 0
    ]

    if not invalid_none.empty:
        problems.append(
            "One or more CONTNone samples have non-zero contamination."
        )

    if problems:
        raise ValueError(
            "Truth-set validation failed:\n- " + "\n- ".join(problems)
        )


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def write_report(df, total_samples_found):
    """Write a human-readable provenance and dataset summary."""

    contamination_counts = (
        df["contamination_percent"]
        .value_counts()
        .sort_index()
    )

    contaminant_counts = (
        df["contaminant_code"]
        .value_counts()
        .sort_index()
    )

    cross_table = pd.crosstab(
        df["contaminant_code"],
        df["contamination_percent"],
    )

    with REPORT_FILE.open("w") as handle:

        handle.write("Contamination truth-set build report\n")
        handle.write("=" * 40 + "\n\n")

        handle.write(f"Source study: {SOURCE_STUDY}\n")
        handle.write(f"BioProject: {BIOPROJECT}\n")
        handle.write(f"Source workbook: {INPUT_FILE}\n\n")

        handle.write(
            f"Total simulation sample IDs found in workbook: "
            f"{total_samples_found}\n"
        )

        handle.write(
            f"COV{TARGET_COVERAGE} samples retained: {len(df)}\n\n"
        )

        handle.write("Samples by contamination percentage\n")
        handle.write("-" * 40 + "\n")
        handle.write(contamination_counts.to_string())
        handle.write("\n\n")

        handle.write("Samples by contaminant code\n")
        handle.write("-" * 40 + "\n")
        handle.write(contaminant_counts.to_string())
        handle.write("\n\n")

        handle.write(
            "Contaminant code x contamination percentage\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(cross_table.to_string())
        handle.write("\n")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    if not INPUT_FILE.exists():
        sys.exit(
            f"ERROR: Source workbook not found:\n{INPUT_FILE}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {INPUT_FILE}")

    sample_ids = extract_sample_ids(INPUT_FILE)

    print(f"Found {len(sample_ids)} simulation sample identifiers.")

    records = []

    for sample_id in sample_ids:

        parsed = parse_sample_id(sample_id)

        if parsed is not None:
            records.append(parsed)

    all_samples = pd.DataFrame(records)

    truth_set = all_samples[
        (all_samples["dataset_coverage"] == TARGET_COVERAGE)
        & (all_samples["mtb_coverage"] == TARGET_COVERAGE)
    ].copy()

    truth_set = truth_set.sort_values(
        [
            "contamination_percent",
            "contaminant_code",
            "simulated_mtb_strain",
            "source_sample_id",
        ]
    ).reset_index(drop=True)

    validate_truth_set(truth_set)

    truth_set.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    write_report(
        truth_set,
        total_samples_found=len(all_samples),
    )

    print()
    print("Truth set successfully created.")
    print(f"Samples retained: {len(truth_set)}")
    print(f"Truth set: {OUTPUT_FILE}")
    print(f"Report:    {REPORT_FILE}")


if __name__ == "__main__":
    main()
