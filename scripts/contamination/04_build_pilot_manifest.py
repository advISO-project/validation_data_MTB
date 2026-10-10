#!/usr/bin/env python3

"""
Build a reproducible 30-sample contamination pilot validation panel
from the 100-sample master contamination validation manifest.

Design:
    - 5 uncontaminated controls at 0%
    - 5 contaminated samples at each of 1%, 2%, 3%, 4%, and 5%
    - maximise contaminant-code diversity within each contamination level
    - balance contaminant-code representation across the complete panel
    - deterministic selection using source_sample_id as the tie-breaker

The master validation manifest is not modified.
"""

from pathlib import Path
import sys

import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

MASTER_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reference_dataset"
    / "contamination_validation_manifest.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reference_dataset"
    / "contamination_pilot_manifest.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reports"
    / "contamination_pilot_manifest_report.txt"
)


# ---------------------------------------------------------------------
# Panel design
# ---------------------------------------------------------------------

TARGET_PER_LEVEL = {
    0: 5,
    1: 5,
    2: 5,
    3: 5,
    4: 5,
    5: 5,
}

EXPECTED_PANEL_SIZE = 30


# ---------------------------------------------------------------------
# Load and validate master manifest
# ---------------------------------------------------------------------

def load_master():

    if not MASTER_FILE.exists():
        sys.exit(
            f"ERROR: Master validation manifest not found:\n{MASTER_FILE}"
        )

    df = pd.read_csv(
        MASTER_FILE,
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "source_sample_id",
        "contaminant_code",
        "contamination_percent",
        "run_accession",
        "fastq_1_ftp",
        "fastq_2_ftp",
        "fastq_1_md5",
        "fastq_2_md5",
    }

    missing = required - set(df.columns)

    if missing:
        sys.exit(
            "ERROR: Missing required columns:\n"
            + "\n".join(sorted(missing))
        )

    if len(df) != 100:
        sys.exit(
            f"ERROR: Expected 100 master samples but found {len(df)}."
        )

    df["contamination_percent"] = pd.to_numeric(
        df["contamination_percent"],
        errors="raise",
    ).astype(int)

    if df["source_sample_id"].duplicated().any():
        sys.exit(
            "ERROR: Duplicate source_sample_id values in master manifest."
        )

    if df["run_accession"].duplicated().any():
        sys.exit(
            "ERROR: Duplicate run accessions in master manifest."
        )

    return df


# ---------------------------------------------------------------------
# Select controls
# ---------------------------------------------------------------------

def select_controls(df):

    controls = df[
        df["contamination_percent"] == 0
    ].copy()

    if len(controls) < TARGET_PER_LEVEL[0]:
        sys.exit(
            "ERROR: Insufficient uncontaminated controls."
        )

    # Controls should explicitly be CONTNone.
    invalid_controls = controls[
        controls["contaminant_code"] != "None"
    ]

    if not invalid_controls.empty:
        sys.exit(
            "ERROR: 0% samples contain unexpected contaminant codes."
        )

    controls = controls.sort_values(
        "source_sample_id"
    )

    selected = controls.head(
        TARGET_PER_LEVEL[0]
    ).copy()

    selected["selection_reason"] = (
        "0% uncontaminated control"
    )

    return selected


# ---------------------------------------------------------------------
# Select contaminated samples
# ---------------------------------------------------------------------

def select_contaminated(df):

    selected_rows = []

    # Track how often each contaminant code has already been selected
    # across the pilot panel.
    global_code_counts = {}

    for percent in range(1, 6):

        candidates = df[
            df["contamination_percent"] == percent
        ].copy()

        if len(candidates) < TARGET_PER_LEVEL[percent]:
            sys.exit(
                f"ERROR: Only {len(candidates)} candidates available "
                f"at {percent}% contamination."
            )

        candidates = candidates.sort_values(
            "source_sample_id"
        )

        chosen_indices = []
        codes_used_this_level = set()

        for _ in range(TARGET_PER_LEVEL[percent]):

            remaining = candidates[
                ~candidates.index.isin(chosen_indices)
            ].copy()

            if remaining.empty:
                sys.exit(
                    f"ERROR: Unable to complete selection at "
                    f"{percent}% contamination."
                )

            # Prefer contaminant codes not yet represented at this
            # contamination percentage.
            remaining["used_this_level"] = (
                remaining["contaminant_code"]
                .isin(codes_used_this_level)
                .astype(int)
            )

            # Prefer contaminant codes with the lowest representation
            # across the panel so far.
            remaining["global_count"] = (
                remaining["contaminant_code"]
                .map(global_code_counts)
                .fillna(0)
                .astype(int)
            )

            # Deterministic priority:
            # 1. code not yet used at this contamination level
            # 2. least represented code across panel
            # 3. source_sample_id lexical ordering
            remaining = remaining.sort_values(
                [
                    "used_this_level",
                    "global_count",
                    "source_sample_id",
                ]
            )

            chosen = remaining.iloc[0]

            chosen_indices.append(chosen.name)

            code = chosen["contaminant_code"]

            codes_used_this_level.add(code)

            global_code_counts[code] = (
                global_code_counts.get(code, 0) + 1
            )

        selected = candidates.loc[
            chosen_indices
        ].copy()

        selected["selection_reason"] = (
            f"{percent}% contamination; "
            "selected for contaminant diversity and panel balance"
        )

        selected_rows.append(selected)

    return pd.concat(
        selected_rows,
        ignore_index=True,
    )


# ---------------------------------------------------------------------
# Build pilot panel
# ---------------------------------------------------------------------

def build_panel(df):

    controls = select_controls(df)
    contaminated = select_contaminated(df)

    panel = pd.concat(
        [controls, contaminated],
        ignore_index=True,
    )

    panel = panel.sort_values(
        [
            "contamination_percent",
            "contaminant_code",
            "source_sample_id",
        ]
    ).reset_index(drop=True)

    return panel


# ---------------------------------------------------------------------
# Validate selected panel
# ---------------------------------------------------------------------

def validate_panel(panel):

    problems = []

    if len(panel) != EXPECTED_PANEL_SIZE:
        problems.append(
            f"Panel contains {len(panel)} samples; "
            f"expected {EXPECTED_PANEL_SIZE}."
        )

    if panel["source_sample_id"].duplicated().any():
        problems.append(
            "Duplicate source sample IDs in pilot panel."
        )

    if panel["run_accession"].duplicated().any():
        problems.append(
            "Duplicate run accessions in pilot panel."
        )

    level_counts = (
        panel["contamination_percent"]
        .value_counts()
        .to_dict()
    )

    for percent, expected in TARGET_PER_LEVEL.items():

        observed = level_counts.get(percent, 0)

        if observed != expected:
            problems.append(
                f"{percent}% contamination: expected "
                f"{expected} samples but found {observed}."
            )

    controls = panel[
        panel["contamination_percent"] == 0
    ]

    if not (
        controls["contaminant_code"] == "None"
    ).all():
        problems.append(
            "One or more 0% controls have a contaminant code "
            "other than None."
        )

    required_complete = [
        "source_sample_id",
        "run_accession",
        "fastq_1_ftp",
        "fastq_2_ftp",
        "fastq_1_md5",
        "fastq_2_md5",
    ]

    for field in required_complete:

        if (panel[field].str.strip() == "").any():
            problems.append(
                f"Missing values detected in {field}."
            )

    if problems:
        sys.exit(
            "ERROR: Pilot panel validation failed:\n- "
            + "\n- ".join(problems)
        )


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def write_report(panel):

    percent_counts = (
        panel["contamination_percent"]
        .value_counts()
        .sort_index()
    )

    code_counts = (
        panel["contaminant_code"]
        .value_counts()
        .sort_index()
    )

    cross_table = pd.crosstab(
        panel["contaminant_code"],
        panel["contamination_percent"],
    )

    with REPORT_FILE.open("w") as handle:

        handle.write(
            "Contamination pilot manifest report\n"
        )
        handle.write("=" * 40 + "\n\n")

        handle.write(
            f"Master manifest: {MASTER_FILE}\n"
        )
        handle.write(
            f"Pilot samples selected: {len(panel)}\n\n"
        )

        handle.write(
            "Selection design\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            "5 samples at each contamination level from 0% to 5%.\n"
        )
        handle.write(
            "Contaminated samples selected to maximise contaminant-code "
            "diversity within each level and balance contaminant-code "
            "representation across the complete pilot panel.\n\n"
        )

        handle.write(
            "Samples by contamination percentage\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            percent_counts.to_string()
        )
        handle.write("\n\n")

        handle.write(
            "Samples by contaminant code\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            code_counts.to_string()
        )
        handle.write("\n\n")

        handle.write(
            "Contaminant code x contamination percentage\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            cross_table.to_string()
        )
        handle.write("\n\n")

        handle.write(
            "Selected samples\n"
        )
        handle.write("-" * 40 + "\n")

        handle.write(
            panel[
                [
                    "source_sample_id",
                    "run_accession",
                    "contamination_percent",
                    "contaminant_code",
                ]
            ].to_string(index=False)
        )

        handle.write("\n")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print("Loading master contamination validation manifest...")

    master = load_master()

    print(f"Master samples: {len(master)}")

    panel = build_panel(master)

    validate_panel(panel)

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    panel.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    write_report(panel)

    print()
    print("Pilot contamination manifest successfully created.")
    print(f"Samples selected: {len(panel)}")
    print(f"Manifest: {OUTPUT_FILE}")
    print(f"Report:   {REPORT_FILE}")


if __name__ == "__main__":
    main()
