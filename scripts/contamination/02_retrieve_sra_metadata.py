#!/usr/bin/env python3

"""
Retrieve public sequencing metadata for the Heupink et al. (2021)
contamination simulation dataset (BioProject PRJNA706121).

This script retrieves run-level metadata from the ENA Portal API.
No FASTQ files are downloaded at this stage.

Output:
    contamination/metadata/contamination_sra_metadata.tsv
    contamination/reports/contamination_sra_metadata_report.txt
"""

from pathlib import Path
import sys

import pandas as pd
import requests


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

METADATA_DIR = PROJECT_ROOT / "contamination" / "metadata"
REPORT_DIR = PROJECT_ROOT / "contamination" / "reports"

OUTPUT_FILE = METADATA_DIR / "contamination_sra_metadata.tsv"
REPORT_FILE = REPORT_DIR / "contamination_sra_metadata_report.txt"


# ---------------------------------------------------------------------
# Study / API configuration
# ---------------------------------------------------------------------

BIOPROJECT = "PRJNA706121"

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession",
    "study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "experiment_accession",
    "experiment_alias",
    "experiment_title",
    "library_layout",
    "library_strategy",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]


# ---------------------------------------------------------------------
# Retrieve metadata
# ---------------------------------------------------------------------

def retrieve_metadata():

    params = {
        "result": "read_run",
        "query": f'study_accession="{BIOPROJECT}"',
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": "0",
    }

    print(f"Retrieving ENA metadata for {BIOPROJECT}...")

    try:
        response = requests.get(
            ENA_API,
            params=params,
            timeout=120,
        )

        response.raise_for_status()

    except requests.RequestException as exc:
        sys.exit(
            "ERROR: Failed to retrieve metadata from ENA.\n"
            f"{exc}"
        )

    if not response.text.strip():
        sys.exit(
            f"ERROR: ENA returned no metadata for {BIOPROJECT}."
        )

    OUTPUT_FILE.write_text(response.text)

    print(f"Metadata saved: {OUTPUT_FILE}")


# ---------------------------------------------------------------------
# Validate and summarise
# ---------------------------------------------------------------------

def validate_metadata():

    df = pd.read_csv(
        OUTPUT_FILE,
        sep="\t",
        dtype=str,
    )

    if df.empty:
        sys.exit("ERROR: Retrieved metadata table is empty.")

    if "run_accession" not in df.columns:
        sys.exit(
            "ERROR: run_accession column is missing from ENA output."
        )

    duplicate_runs = df["run_accession"].duplicated().sum()

    missing_runs = df["run_accession"].isna().sum()

    return df, duplicate_runs, missing_runs


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def write_report(df, duplicate_runs, missing_runs):

    with REPORT_FILE.open("w") as handle:

        handle.write("Contamination SRA/ENA metadata retrieval report\n")
        handle.write("=" * 48 + "\n\n")

        handle.write(f"BioProject: {BIOPROJECT}\n")
        handle.write("Source: ENA Portal API\n\n")

        handle.write(f"Runs retrieved: {len(df)}\n")
        handle.write(
            f"Unique run accessions: "
            f"{df['run_accession'].nunique()}\n"
        )
        handle.write(
            f"Duplicate run accessions: {duplicate_runs}\n"
        )
        handle.write(
            f"Missing run accessions: {missing_runs}\n\n"
        )

        for field in [
            "sample_alias",
            "sample_title",
            "sample_description",
            "experiment_alias",
            "experiment_title",
        ]:
            if field in df.columns:

                populated = df[field].notna().sum()

                handle.write(
                    f"{field} populated: "
                    f"{populated}/{len(df)}\n"
                )

        handle.write("\n")

        if "library_layout" in df.columns:
            handle.write("Library layout\n")
            handle.write("-" * 30 + "\n")
            handle.write(
                df["library_layout"]
                .fillna("MISSING")
                .value_counts()
                .to_string()
            )
            handle.write("\n\n")

        if "instrument_platform" in df.columns:
            handle.write("Instrument platform\n")
            handle.write("-" * 30 + "\n")
            handle.write(
                df["instrument_platform"]
                .fillna("MISSING")
                .value_counts()
                .to_string()
            )
            handle.write("\n\n")

        handle.write("Example records\n")
        handle.write("-" * 30 + "\n")

        preview_fields = [
            field for field in [
                "run_accession",
                "sample_accession",
                "sample_alias",
                "experiment_alias",
                "experiment_title",
            ]
            if field in df.columns
        ]

        handle.write(
            df[preview_fields]
            .head(10)
            .to_string(index=False)
        )

        handle.write("\n")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    retrieve_metadata()

    df, duplicate_runs, missing_runs = validate_metadata()

    write_report(
        df,
        duplicate_runs,
        missing_runs,
    )

    print()
    print("Metadata retrieval successfully completed.")
    print(f"Runs retrieved: {len(df)}")
    print(f"Metadata: {OUTPUT_FILE}")
    print(f"Report:   {REPORT_FILE}")


if __name__ == "__main__":
    main()
