#!/usr/bin/env python3

"""
Build the contamination validation manifest by linking the published
Heupink et al. (2021) simulation truth set to ENA sequencing runs.

Join:
    contamination_truth_set.csv: source_sample_id
    contamination_sra_metadata.tsv: sample_alias

Expected result:
    100 truth samples
    100 unique ENA run matches
    paired-end FASTQs for every sample

No FASTQ files are downloaded by this script.
"""

from pathlib import Path
import sys

import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

TRUTH_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reference_dataset"
    / "contamination_truth_set.csv"
)

ENA_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "metadata"
    / "contamination_sra_metadata.tsv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reference_dataset"
    / "contamination_validation_manifest.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reports"
    / "contamination_validation_manifest_report.txt"
)

EXPECTED_TRUTH_SAMPLES = 100


# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------

def load_inputs():

    if not TRUTH_FILE.exists():
        sys.exit(
            f"ERROR: Truth set not found:\n{TRUTH_FILE}"
        )

    if not ENA_FILE.exists():
        sys.exit(
            f"ERROR: ENA metadata not found:\n{ENA_FILE}"
        )

    truth = pd.read_csv(
        TRUTH_FILE,
        dtype=str,
        keep_default_na=False,
    )

    ena = pd.read_csv(
        ENA_FILE,
        sep="\t",
        dtype=str,
    )

    return truth, ena


# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------

def validate_inputs(truth, ena):

    required_truth = {
        "source_sample_id",
        "simulated_mtb_strain",
        "dataset_coverage",
        "mtb_coverage",
        "contaminant_code",
        "contaminant_species",
        "contaminant_class",
        "contamination_percent",
        "simulation_seed",
        "source_study",
        "bioproject",
    }

    required_ena = {
        "run_accession",
        "sample_accession",
        "sample_alias",
        "library_layout",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    }

    missing_truth = required_truth - set(truth.columns)
    missing_ena = required_ena - set(ena.columns)

    if missing_truth:
        sys.exit(
            "ERROR: Missing truth-set columns:\n"
            + "\n".join(sorted(missing_truth))
        )

    if missing_ena:
        sys.exit(
            "ERROR: Missing ENA metadata columns:\n"
            + "\n".join(sorted(missing_ena))
        )

    if len(truth) != EXPECTED_TRUTH_SAMPLES:
        sys.exit(
            f"ERROR: Expected {EXPECTED_TRUTH_SAMPLES} truth samples "
            f"but found {len(truth)}."
        )

    if truth["source_sample_id"].duplicated().any():
        sys.exit(
            "ERROR: Duplicate source_sample_id values in truth set."
        )


# ---------------------------------------------------------------------
# Restrict ENA metadata to truth-set samples
# ---------------------------------------------------------------------

def select_matching_ena_records(truth, ena):

    truth_ids = set(truth["source_sample_id"])

    matching = ena[
        ena["sample_alias"].isin(truth_ids)
    ].copy()

    return matching


# ---------------------------------------------------------------------
# Validate mapping
# ---------------------------------------------------------------------

def validate_mapping(truth, matching):

    problems = []

    truth_ids = set(truth["source_sample_id"])
    matched_ids = set(matching["sample_alias"])

    unmatched = sorted(truth_ids - matched_ids)

    if unmatched:
        problems.append(
            f"{len(unmatched)} truth samples have no ENA match."
        )

    # Each truth sample should correspond to exactly one run.
    run_counts = (
        matching.groupby("sample_alias")["run_accession"]
        .nunique()
    )

    multiple_runs = run_counts[run_counts != 1]

    if not multiple_runs.empty:
        problems.append(
            f"{len(multiple_runs)} samples do not map to exactly one "
            "unique run accession."
        )

    if matching["run_accession"].duplicated().any():
        problems.append(
            "Duplicate run accessions detected among matched samples."
        )

    non_paired = matching[
        matching["library_layout"] != "PAIRED"
    ]

    if not non_paired.empty:
        problems.append(
            f"{len(non_paired)} matched runs are not PAIRED."
        )

    # Validate FASTQ information.
    missing_fastq = matching[
        matching["fastq_ftp"].isna()
        | (matching["fastq_ftp"].str.strip() == "")
    ]

    if not missing_fastq.empty:
        problems.append(
            f"{len(missing_fastq)} matched runs have no FASTQ URLs."
        )

    invalid_fastq_count = []

    for _, row in matching.iterrows():

        if pd.isna(row["fastq_ftp"]):
            continue

        urls = [
            x for x in row["fastq_ftp"].split(";")
            if x.strip()
        ]

        if len(urls) != 2:
            invalid_fastq_count.append(
                row["run_accession"]
            )

    if invalid_fastq_count:
        problems.append(
            f"{len(invalid_fastq_count)} runs do not have exactly "
            "two FASTQ URLs."
        )

    if len(matching) != len(truth):
        problems.append(
            f"Expected {len(truth)} matched records but found "
            f"{len(matching)}."
        )

    if problems:
        sys.exit(
            "ERROR: Validation manifest mapping failed:\n- "
            + "\n- ".join(problems)
        )


# ---------------------------------------------------------------------
# Build manifest
# ---------------------------------------------------------------------

def build_manifest(truth, matching):

    ena_fields = [
        "sample_alias",
        "run_accession",
        "sample_accession",
        "secondary_sample_accession",
        "experiment_accession",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]

    ena_fields = [
        field for field in ena_fields
        if field in matching.columns
    ]

    ena_subset = matching[ena_fields].copy()

    manifest = truth.merge(
        ena_subset,
        left_on="source_sample_id",
        right_on="sample_alias",
        how="left",
        validate="one_to_one",
    )

    # sample_alias duplicates source_sample_id after a successful join.
    manifest = manifest.drop(
        columns=["sample_alias"]
    )

    # Split paired FASTQ metadata into explicit R1/R2 fields.
    manifest[
        ["fastq_1_ftp", "fastq_2_ftp"]
    ] = manifest["fastq_ftp"].str.split(
        ";",
        n=1,
        expand=True,
    )

    manifest[
        ["fastq_1_md5", "fastq_2_md5"]
    ] = manifest["fastq_md5"].str.split(
        ";",
        n=1,
        expand=True,
    )

    manifest[
        ["fastq_1_bytes", "fastq_2_bytes"]
    ] = manifest["fastq_bytes"].str.split(
        ";",
        n=1,
        expand=True,
    )

    # Add protocol for FTP URLs so they can be used directly later.
    for field in ["fastq_1_ftp", "fastq_2_ftp"]:

        manifest[field] = manifest[field].apply(
            lambda x: (
                f"ftp://{x}"
                if pd.notna(x)
                and not str(x).startswith(
                    ("ftp://", "http://", "https://")
                )
                else x
            )
        )

    # Raw combined ENA fields are no longer needed.
    manifest = manifest.drop(
        columns=[
            "fastq_ftp",
            "fastq_md5",
            "fastq_bytes",
        ]
    )

    return manifest


# ---------------------------------------------------------------------
# Final manifest validation
# ---------------------------------------------------------------------

def validate_manifest(manifest):

    problems = []

    if len(manifest) != EXPECTED_TRUTH_SAMPLES:
        problems.append(
            f"Manifest contains {len(manifest)} rows; "
            f"expected {EXPECTED_TRUTH_SAMPLES}."
        )

    if manifest["source_sample_id"].duplicated().any():
        problems.append(
            "Duplicate source sample IDs in final manifest."
        )

    if manifest["run_accession"].duplicated().any():
        problems.append(
            "Duplicate run accessions in final manifest."
        )

    required_complete = [
        "source_sample_id",
        "run_accession",
        "sample_accession",
        "fastq_1_ftp",
        "fastq_2_ftp",
        "fastq_1_md5",
        "fastq_2_md5",
    ]

    for field in required_complete:

        if manifest[field].isna().any():
            problems.append(
                f"Missing values detected in {field}."
            )

    if problems:
        sys.exit(
            "ERROR: Final manifest validation failed:\n- "
            + "\n- ".join(problems)
        )


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def write_report(truth, matching, manifest):

    contamination_counts = (
        manifest["contamination_percent"]
        .value_counts()
        .sort_index()
    )

    contaminant_counts = (
        manifest["contaminant_code"]
        .value_counts()
        .sort_index()
    )

    with REPORT_FILE.open("w") as handle:

        handle.write(
            "Contamination validation manifest report\n"
        )
        handle.write("=" * 42 + "\n\n")

        handle.write(
            f"Truth-set samples: {len(truth)}\n"
        )

        handle.write(
            f"ENA records matched: {len(matching)}\n"
        )

        handle.write(
            f"Final manifest rows: {len(manifest)}\n"
        )

        handle.write(
            f"Unique run accessions: "
            f"{manifest['run_accession'].nunique()}\n"
        )

        handle.write(
            f"Paired FASTQ sets: {len(manifest)}\n\n"
        )

        handle.write(
            "Mapping key:\n"
            "source_sample_id = ENA sample_alias\n\n"
        )

        handle.write(
            "Samples by contamination percentage\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            contamination_counts.to_string()
        )
        handle.write("\n\n")

        handle.write(
            "Samples by contaminant code\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            contaminant_counts.to_string()
        )
        handle.write("\n\n")

        handle.write("Example mappings\n")
        handle.write("-" * 40 + "\n")

        handle.write(
            manifest[
                [
                    "source_sample_id",
                    "run_accession",
                    "sample_accession",
                ]
            ]
            .head(10)
            .to_string(index=False)
        )

        handle.write("\n")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print("Loading contamination truth set and ENA metadata...")

    truth, ena = load_inputs()

    validate_inputs(truth, ena)

    print(f"Truth samples: {len(truth)}")
    print(f"ENA runs available: {len(ena)}")

    matching = select_matching_ena_records(
        truth,
        ena,
    )

    print(
        f"ENA records matching truth samples: "
        f"{len(matching)}"
    )

    validate_mapping(
        truth,
        matching,
    )

    manifest = build_manifest(
        truth,
        matching,
    )

    validate_manifest(manifest)

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    manifest.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    write_report(
        truth,
        matching,
        manifest,
    )

    print()
    print("Validation manifest successfully created.")
    print(f"Samples:  {len(manifest)}")
    print(
        f"Unique runs: "
        f"{manifest['run_accession'].nunique()}"
    )
    print(f"Manifest: {OUTPUT_FILE}")
    print(f"Report:   {REPORT_FILE}")


if __name__ == "__main__":
    main()
