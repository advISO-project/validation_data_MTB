#!/usr/bin/env python3

"""
Build the mixed-infection validation manifest.

This script joins:

1. The curated truth set derived from known sample composition
   reported in Sobkowiak et al. (2018), Table 1.

2. ENA sequencing metadata retrieved for those exact run accessions.

The resulting manifest provides the traceable link between the
expected validation result and the corresponding public sequence
files.

No truth labels are assigned or modified by this script.
No external resources are queried.
"""

from pathlib import Path
import csv
import hashlib
import sys
from collections import Counter


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

TRUTH_FILE = Path(
    "mixed_infection/reference_dataset/mixed_infection_truth_set.csv"
)

ENA_FILE = Path(
    "mixed_infection/metadata/mixed_infection_ena_metadata.tsv"
)

OUTPUT_FILE = Path(
    "mixed_infection/reference_dataset/"
    "mixed_infection_validation_manifest.csv"
)

REPORT_FILE = Path(
    "mixed_infection/reports/"
    "mixed_infection_validation_manifest_report.txt"
)

CHECKSUM_FILE = Path(
    "mixed_infection/metadata/"
    "mixed_infection_validation_manifest.sha256"
)


# ---------------------------------------------------------------------
# Expected dataset properties
# ---------------------------------------------------------------------

EXPECTED_TOTAL = 48

EXPECTED_STATUS_COUNTS = {
    "mixed": 36,
    "pure": 12,
}

EXPECTED_PROPORTION_COUNTS = {
    "0.70": 12,
    "0.90": 12,
    "0.95": 12,
    "1.00": 12,
}

EXPECTED_STUDY = "PRJEB2794"
EXPECTED_LAYOUT = "PAIRED"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_table(path, delimiter):
    if not path.exists():
        fail(f"Input file not found: {path}")

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter=delimiter,
            )
        )


def split_ena_pair(value, field_name, accession):
    """
    Split a semicolon-separated ENA paired-file field.

    Exactly two non-empty values are required.
    """

    values = [
        item.strip()
        for item in value.split(";")
        if item.strip()
    ]

    if len(values) != 2:
        fail(
            f"{accession}: expected exactly 2 values in "
            f"{field_name}, found {len(values)}."
        )

    return values[0], values[1]


def basename_from_ftp(path):
    return path.rstrip("/").split("/")[-1]


# ---------------------------------------------------------------------
# Read inputs
# ---------------------------------------------------------------------

truth_records = read_table(
    TRUTH_FILE,
    delimiter=",",
)

ena_records = read_table(
    ENA_FILE,
    delimiter="\t",
)


if len(truth_records) != EXPECTED_TOTAL:
    fail(
        f"Expected {EXPECTED_TOTAL} truth records, "
        f"found {len(truth_records)}."
    )

if len(ena_records) != EXPECTED_TOTAL:
    fail(
        f"Expected {EXPECTED_TOTAL} ENA records, "
        f"found {len(ena_records)}."
    )


# ---------------------------------------------------------------------
# Validate accession uniqueness
# ---------------------------------------------------------------------

truth_accessions = [
    row["run_accession"].strip()
    for row in truth_records
]

ena_accessions = [
    row["run_accession"].strip()
    for row in ena_records
]


if len(set(truth_accessions)) != EXPECTED_TOTAL:
    fail("Truth set does not contain 48 unique run accessions.")

if len(set(ena_accessions)) != EXPECTED_TOTAL:
    fail("ENA metadata does not contain 48 unique run accessions.")


missing_from_ena = sorted(
    set(truth_accessions) - set(ena_accessions)
)

unexpected_in_ena = sorted(
    set(ena_accessions) - set(truth_accessions)
)


if missing_from_ena:
    fail(
        "Truth-set accessions missing from ENA metadata: "
        + ", ".join(missing_from_ena)
    )

if unexpected_in_ena:
    fail(
        "Unexpected accessions in ENA metadata: "
        + ", ".join(unexpected_in_ena)
    )


# ---------------------------------------------------------------------
# Index ENA records
# ---------------------------------------------------------------------

ena_by_run = {
    row["run_accession"].strip(): row
    for row in ena_records
}


# ---------------------------------------------------------------------
# Join truth + ENA metadata
# ---------------------------------------------------------------------

manifest = []

for truth in truth_records:

    accession = truth["run_accession"].strip()
    ena = ena_by_run[accession]

    # -------------------------------------------------------------
    # Validate expected ENA properties
    # -------------------------------------------------------------

    if ena["study_accession"].strip() != EXPECTED_STUDY:
        fail(
            f"{accession}: unexpected study accession "
            f"{ena['study_accession']!r}; "
            f"expected {EXPECTED_STUDY}."
        )

    if ena["library_layout"].strip() != EXPECTED_LAYOUT:
        fail(
            f"{accession}: unexpected library layout "
            f"{ena['library_layout']!r}; "
            f"expected {EXPECTED_LAYOUT}."
        )

    # -------------------------------------------------------------
    # Split paired FASTQ metadata
    # -------------------------------------------------------------

    fastq_1, fastq_2 = split_ena_pair(
        ena["fastq_ftp"],
        "fastq_ftp",
        accession,
    )

    md5_1, md5_2 = split_ena_pair(
        ena["fastq_md5"],
        "fastq_md5",
        accession,
    )

    bytes_1, bytes_2 = split_ena_pair(
        ena["fastq_bytes"],
        "fastq_bytes",
        accession,
    )

    # -------------------------------------------------------------
    # Check expected filenames
    # -------------------------------------------------------------

    expected_name_1 = f"{accession}_1.fastq.gz"
    expected_name_2 = f"{accession}_2.fastq.gz"

    observed_name_1 = basename_from_ftp(fastq_1)
    observed_name_2 = basename_from_ftp(fastq_2)

    if observed_name_1 != expected_name_1:
        fail(
            f"{accession}: unexpected read-1 filename "
            f"{observed_name_1!r}; "
            f"expected {expected_name_1!r}."
        )

    if observed_name_2 != expected_name_2:
        fail(
            f"{accession}: unexpected read-2 filename "
            f"{observed_name_2!r}; "
            f"expected {expected_name_2!r}."
        )

    # -------------------------------------------------------------
    # Validate byte fields
    # -------------------------------------------------------------

    try:
        bytes_1_int = int(bytes_1)
        bytes_2_int = int(bytes_2)
    except ValueError:
        fail(
            f"{accession}: non-integer FASTQ byte value."
        )

    if bytes_1_int <= 0 or bytes_2_int <= 0:
        fail(
            f"{accession}: FASTQ byte values must be > 0."
        )

    # -------------------------------------------------------------
    # Construct joined record
    # -------------------------------------------------------------

    manifest.append(
        {
            # Identity
            "run_accession": accession,

            # Expected result / truth
            "expected_status":
                truth["expected_status"],

            "published_major_proportion":
                truth["published_major_proportion"],

            "expected_minor_proportion":
                truth["expected_minor_proportion"],

            "major_lineage":
                truth["major_lineage"],

            "minor_lineage":
                truth["minor_lineage"],

            "major_spoligotype":
                truth["major_spoligotype"],

            "minor_spoligotype":
                truth["minor_spoligotype"],

            "truth_basis":
                truth["truth_basis"],

            "primary_truth_source":
                truth["primary_truth_source"],

            # ENA identity
            "study_accession":
                ena["study_accession"],

            "sample_accession":
                ena["sample_accession"],

            "experiment_accession":
                ena["experiment_accession"],

            "scientific_name":
                ena["scientific_name"],

            "library_layout":
                ena["library_layout"],

            "instrument_platform":
                ena["instrument_platform"],

            "instrument_model":
                ena["instrument_model"],

            # FASTQ provenance
            "fastq_1":
                fastq_1,

            "fastq_2":
                fastq_2,

            "md5_1":
                md5_1,

            "md5_2":
                md5_2,

            "bytes_1":
                bytes_1,

            "bytes_2":
                bytes_2,
        }
    )


# ---------------------------------------------------------------------
# Dataset-level integrity checks
# ---------------------------------------------------------------------

if len(manifest) != EXPECTED_TOTAL:
    fail(
        f"Joined manifest contains {len(manifest)} records; "
        f"expected {EXPECTED_TOTAL}."
    )


status_counts = Counter(
    row["expected_status"]
    for row in manifest
)

if dict(status_counts) != EXPECTED_STATUS_COUNTS:
    fail(
        "Unexpected mixed/pure distribution.\n"
        f"Observed: {dict(status_counts)}\n"
        f"Expected: {EXPECTED_STATUS_COUNTS}"
    )


proportion_counts = Counter(
    row["published_major_proportion"]
    for row in manifest
)

if dict(proportion_counts) != EXPECTED_PROPORTION_COUNTS:
    fail(
        "Unexpected major-proportion distribution.\n"
        f"Observed: {dict(proportion_counts)}\n"
        f"Expected: {EXPECTED_PROPORTION_COUNTS}"
    )


# ---------------------------------------------------------------------
# Check ENA identifiers are populated
# ---------------------------------------------------------------------

for row in manifest:

    accession = row["run_accession"]

    for field in [
        "sample_accession",
        "experiment_accession",
        "scientific_name",
        "instrument_platform",
        "instrument_model",
    ]:
        if not row[field].strip():
            fail(
                f"{accession}: missing required ENA field "
                f"{field}."
            )


# ---------------------------------------------------------------------
# Check sample and experiment mappings
# ---------------------------------------------------------------------

sample_accessions = [
    row["sample_accession"]
    for row in manifest
]

experiment_accessions = [
    row["experiment_accession"]
    for row in manifest
]


duplicate_samples = sorted(
    accession
    for accession, count
    in Counter(sample_accessions).items()
    if count > 1
)

duplicate_experiments = sorted(
    accession
    for accession, count
    in Counter(experiment_accessions).items()
    if count > 1
)


if duplicate_samples:
    fail(
        "Duplicate ENA sample accessions detected: "
        + ", ".join(duplicate_samples)
    )

if duplicate_experiments:
    fail(
        "Duplicate ENA experiment accessions detected: "
        + ", ".join(duplicate_experiments)
    )


# ---------------------------------------------------------------------
# Write manifest
# ---------------------------------------------------------------------

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

fieldnames = [
    "run_accession",

    "expected_status",
    "published_major_proportion",
    "expected_minor_proportion",
    "major_lineage",
    "minor_lineage",
    "major_spoligotype",
    "minor_spoligotype",
    "truth_basis",
    "primary_truth_source",

    "study_accession",
    "sample_accession",
    "experiment_accession",
    "scientific_name",
    "library_layout",
    "instrument_platform",
    "instrument_model",

    "fastq_1",
    "fastq_2",
    "md5_1",
    "md5_2",
    "bytes_1",
    "bytes_2",
]


with OUTPUT_FILE.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(manifest)


# ---------------------------------------------------------------------
# SHA256
# ---------------------------------------------------------------------

output_bytes = OUTPUT_FILE.read_bytes()

sha256 = hashlib.sha256(
    output_bytes
).hexdigest()

CHECKSUM_FILE.write_text(
    f"{sha256}  {OUTPUT_FILE}\n",
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

total_fastq_bytes = sum(
    int(row["bytes_1"]) + int(row["bytes_2"])
    for row in manifest
)

total_fastq_gb = (
    total_fastq_bytes / (1024 ** 3)
)


report = [
    "Mixed-infection validation manifest report",
    "==========================================",
    "",
    f"Truth input: {TRUTH_FILE}",
    f"ENA input: {ENA_FILE}",
    "",
    f"Truth records: {len(truth_records)}",
    f"ENA records: {len(ena_records)}",
    f"Joined records: {len(manifest)}",
    "",
    f"Missing from ENA: {len(missing_from_ena)}",
    f"Unexpected in ENA: {len(unexpected_in_ena)}",
    "",
    "Expected status:",
    f"  mixed: {status_counts['mixed']}",
    f"  pure: {status_counts['pure']}",
    "",
    "Published major-strain proportions:",
]

for proportion in [
    "0.70",
    "0.90",
    "0.95",
    "1.00",
]:
    report.append(
        f"  {proportion}: "
        f"{proportion_counts[proportion]}"
    )


report.extend(
    [
        "",
        f"Study accession: {EXPECTED_STUDY}",
        f"Library layout: {EXPECTED_LAYOUT}",
        "",
        f"Unique sample accessions: "
        f"{len(set(sample_accessions))}",
        f"Unique experiment accessions: "
        f"{len(set(experiment_accessions))}",
        "",
        "FASTQ metadata:",
        f"  paired FASTQ sets: {len(manifest)}",
        f"  FASTQ files represented: {len(manifest) * 2}",
        f"  expected download size: {total_fastq_gb:.2f} GiB",
        "",
        f"Output: {OUTPUT_FILE}",
        f"SHA256: {sha256}",
        "",
        "STATUS: PASS",
    ]
)


REPORT_FILE.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

print("Mixed-infection validation manifest constructed")
print("-----------------------------------------------")
print(f"Truth records:       {len(truth_records)}")
print(f"ENA records:         {len(ena_records)}")
print(f"Joined records:      {len(manifest)}")
print(f"Mixed:               {status_counts['mixed']}")
print(f"Pure:                {status_counts['pure']}")
print(f"FASTQ files:         {len(manifest) * 2}")
print(f"Expected size:       {total_fastq_gb:.2f} GiB")
print()
print(f"Manifest: {OUTPUT_FILE}")
print(f"Report:   {REPORT_FILE}")
print(f"SHA256:   {sha256}")
print("STATUS: PASS")
