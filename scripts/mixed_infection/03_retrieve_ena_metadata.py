#!/usr/bin/env python3

"""
Retrieve ENA sequencing metadata for the 48 accessions in the
mixed-infection truth set.

Input:
    mixed_infection/reference_dataset/mixed_infection_truth_set.csv

Output:
    mixed_infection/metadata/mixed_infection_ena_metadata.tsv

Purpose:
    Link the published truth-set accessions to their public ENA
    sequencing records.

This script does not assign or modify truth labels.
"""

from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import sys
import urllib.parse
import urllib.request


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

TRUTH_FILE = Path(
    "mixed_infection/reference_dataset/mixed_infection_truth_set.csv"
)

OUTPUT_FILE = Path(
    "mixed_infection/metadata/mixed_infection_ena_metadata.tsv"
)

CHECKSUM_FILE = Path(
    "mixed_infection/metadata/mixed_infection_ena_metadata.sha256"
)

REPORT_FILE = Path(
    "mixed_infection/reports/mixed_infection_ena_retrieval_report.txt"
)


# ---------------------------------------------------------------------
# ENA configuration
# ---------------------------------------------------------------------

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "scientific_name",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def retrieve_run(accession):
    """
    Retrieve one ENA read-run record by run accession.
    """

    params = {
        "result": "read_run",
        "query": f'run_accession="{accession}"',
        "fields": ",".join(FIELDS),
        "format": "tsv",
    }

    url = ENA_API + "?" + urllib.parse.urlencode(params)

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "ISO-in-a-Box validation dataset metadata retrieval"
        },
    )

    with urllib.request.urlopen(request) as response:
        text = response.read().decode("utf-8")

    lines = [
        line
        for line in text.splitlines()
        if line.strip()
    ]

    if len(lines) != 2:
        fail(
            f"{accession}: expected exactly one ENA record, "
            f"but response contained {max(len(lines) - 1, 0)}."
        )

    reader = csv.DictReader(lines, delimiter="\t")
    records = list(reader)

    if len(records) != 1:
        fail(
            f"{accession}: could not resolve exactly one ENA record."
        )

    record = records[0]

    if record["run_accession"] != accession:
        fail(
            f"{accession}: ENA returned unexpected run accession "
            f"{record['run_accession']}."
        )

    return record


# ---------------------------------------------------------------------
# Read truth-set accessions
# ---------------------------------------------------------------------

if not TRUTH_FILE.exists():
    fail(f"Truth set not found: {TRUTH_FILE}")

with TRUTH_FILE.open(
    newline="",
    encoding="utf-8",
) as handle:

    truth_records = list(csv.DictReader(handle))


if len(truth_records) != 48:
    fail(
        f"Expected 48 truth-set records, "
        f"found {len(truth_records)}."
    )


accessions = [
    row["run_accession"].strip()
    for row in truth_records
]


if len(set(accessions)) != 48:
    fail("Truth set does not contain 48 unique run accessions.")


print("ENA metadata retrieval")
print("----------------------")
print(f"Truth-set accessions: {len(accessions)}")
print()


# ---------------------------------------------------------------------
# Retrieve exact accessions
# ---------------------------------------------------------------------

ena_records = []

for number, accession in enumerate(accessions, start=1):

    print(
        f"[{number:02d}/{len(accessions)}] "
        f"Retrieving {accession}"
    )

    record = retrieve_run(accession)

    ena_records.append(record)


# ---------------------------------------------------------------------
# Final accession integrity check
# ---------------------------------------------------------------------

returned_accessions = [
    record["run_accession"]
    for record in ena_records
]

missing = sorted(
    set(accessions) - set(returned_accessions)
)

unexpected = sorted(
    set(returned_accessions) - set(accessions)
)

duplicates = sorted(
    {
        accession
        for accession in returned_accessions
        if returned_accessions.count(accession) > 1
    }
)


if missing:
    fail(
        "Truth-set accessions missing from ENA response: "
        + ", ".join(missing)
    )

if unexpected:
    fail(
        "Unexpected ENA accessions returned: "
        + ", ".join(unexpected)
    )

if duplicates:
    fail(
        "Duplicate ENA accessions returned: "
        + ", ".join(duplicates)
    )


# ---------------------------------------------------------------------
# Write metadata
# ---------------------------------------------------------------------

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

with OUTPUT_FILE.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(ena_records)


# ---------------------------------------------------------------------
# Checksum
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
# Summaries
# ---------------------------------------------------------------------

layout_counts = {}
platform_counts = {}
model_counts = {}
study_counts = {}

missing_fastq = []

for record in ena_records:

    layout = record["library_layout"] or "<missing>"
    platform = record["instrument_platform"] or "<missing>"
    model = record["instrument_model"] or "<missing>"
    study = record["study_accession"] or "<missing>"

    layout_counts[layout] = layout_counts.get(layout, 0) + 1
    platform_counts[platform] = platform_counts.get(platform, 0) + 1
    model_counts[model] = model_counts.get(model, 0) + 1
    study_counts[study] = study_counts.get(study, 0) + 1

    if not record["fastq_ftp"]:
        missing_fastq.append(record["run_accession"])


retrieved_utc = datetime.now(timezone.utc).isoformat()


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

report = [
    "Mixed-infection ENA metadata retrieval report",
    "=============================================",
    "",
    f"Input truth set: {TRUTH_FILE}",
    f"Retrieved UTC: {retrieved_utc}",
    "",
    f"Requested accessions: {len(accessions)}",
    f"ENA records retrieved: {len(ena_records)}",
    f"Unique ENA run accessions: {len(set(returned_accessions))}",
    f"Missing requested accessions: {len(missing)}",
    f"Unexpected accessions: {len(unexpected)}",
    "",
    "Study accessions:",
]

for value in sorted(study_counts):
    report.append(
        f"  {value}: {study_counts[value]}"
    )

report.extend(
    [
        "",
        "Library layout:",
    ]
)

for value in sorted(layout_counts):
    report.append(
        f"  {value}: {layout_counts[value]}"
    )

report.extend(
    [
        "",
        "Instrument platform:",
    ]
)

for value in sorted(platform_counts):
    report.append(
        f"  {value}: {platform_counts[value]}"
    )

report.extend(
    [
        "",
        "Instrument model:",
    ]
)

for value in sorted(model_counts):
    report.append(
        f"  {value}: {model_counts[value]}"
    )

report.extend(
    [
        "",
        f"Records without FASTQ links: {len(missing_fastq)}",
    ]
)

if missing_fastq:
    for accession in missing_fastq:
        report.append(
            f"  {accession}"
        )

report.extend(
    [
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
# Console output
# ---------------------------------------------------------------------

print()
print("Retrieval complete")
print("------------------")
print(f"Requested:         {len(accessions)}")
print(f"Retrieved:         {len(ena_records)}")
print(f"Unique accessions: {len(set(returned_accessions))}")
print(f"Missing:           {len(missing)}")
print(f"Unexpected:        {len(unexpected)}")
print(f"No FASTQ links:    {len(missing_fastq)}")
print()
print(f"Metadata: {OUTPUT_FILE}")
print(f"Report:   {REPORT_FILE}")
print(f"SHA256:   {sha256}")
print("STATUS: PASS")
