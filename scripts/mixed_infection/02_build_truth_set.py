#!/usr/bin/env python3

"""
Build the mixed-infection truth set from the frozen Table 1 of
Sobkowiak et al. (2018).

Truth basis
-----------
Expected mixed-infection status is derived from the published known
sample composition, not from the computational results reported in
the same table.

Published source fields retained:
    - Sample identifier
    - Known major strain proportion
    - Major/minor strain lineage
    - Major/minor strain spoligotype

Derived fields:
    - Expected minor strain proportion = 1 - major strain proportion
    - Expected status:
          major proportion < 1.0 -> mixed
          major proportion = 1.0 -> pure

Computational-result columns from the 2018 analysis are deliberately
excluded from truth assignment.

Source:
Sobkowiak B et al.
Identifying mixed Mycobacterium tuberculosis infections from whole
genome sequence data.
BMC Genomics. 2018;19:613.
DOI: 10.1186/s12864-018-4988-z
PMCID: PMC6092779
Table 1
"""

from pathlib import Path
from html.parser import HTMLParser
import csv
import hashlib
import re
import sys


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

SOURCE_FILE = Path(
    "mixed_infection/raw/sobkowiak2018_table1.html"
)

OUTPUT_FILE = Path(
    "mixed_infection/reference_dataset/mixed_infection_truth_set.csv"
)

REPORT_FILE = Path(
    "mixed_infection/reports/mixed_infection_truth_set_report.txt"
)

CHECKSUM_FILE = Path(
    "mixed_infection/metadata/mixed_infection_truth_set.sha256"
)


# ---------------------------------------------------------------------
# Expected dataset structure
# ---------------------------------------------------------------------

EXPECTED_TOTAL = 48

EXPECTED_PROPORTION_COUNTS = {
    0.70: 12,
    0.90: 12,
    0.95: 12,
    1.00: 12,
}

EXPECTED_STATUS_COUNTS = {
    "mixed": 36,
    "pure": 12,
}


# ---------------------------------------------------------------------
# Minimal HTML table parser
# ---------------------------------------------------------------------

class TableParser(HTMLParser):
    """
    Extract text from HTML table rows and cells.
    """

    def __init__(self):
        super().__init__()
        self.in_row = False
        self.in_cell = False
        self.current_cell = []
        self.current_row = []
        self.rows = []

    def handle_starttag(self, tag, attrs):
        tag = tag.lower()

        if tag == "tr":
            self.in_row = True
            self.current_row = []

        elif tag in {"td", "th"} and self.in_row:
            self.in_cell = True
            self.current_cell = []

    def handle_data(self, data):
        if self.in_cell:
            self.current_cell.append(data)

    def handle_endtag(self, tag):
        tag = tag.lower()

        if tag in {"td", "th"} and self.in_cell:
            value = " ".join(
                "".join(self.current_cell).split()
            )
            self.current_row.append(value)
            self.in_cell = False

        elif tag == "tr" and self.in_row:
            if self.current_row:
                self.rows.append(self.current_row)

            self.in_row = False


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def split_pair(value):
    """
    Split a published major/minor value.

    Mixed samples contain values such as:
        1/3
        EAI6-BGD1/CAS1-Kili

    Pure samples contain only one value:
        4
        T

    Returns:
        major, minor

    Minor is blank for pure samples.
    """

    parts = value.split("/", 1)

    major = parts[0].strip()

    if len(parts) == 2:
        minor = parts[1].strip()
    else:
        minor = ""

    return major, minor


# ---------------------------------------------------------------------
# Validate source
# ---------------------------------------------------------------------

if not SOURCE_FILE.exists():
    fail(f"Source file not found: {SOURCE_FILE}")

html = SOURCE_FILE.read_text(
    encoding="utf-8",
    errors="replace",
)

parser = TableParser()
parser.feed(html)


# ---------------------------------------------------------------------
# Identify the 48 sample rows
# ---------------------------------------------------------------------

sample_pattern = re.compile(r"^ERR\d+$")

sample_rows = [
    row
    for row in parser.rows
    if row and sample_pattern.fullmatch(row[0].strip())
]

if len(sample_rows) != EXPECTED_TOTAL:
    fail(
        f"Expected {EXPECTED_TOTAL} sample rows, "
        f"found {len(sample_rows)}."
    )


# ---------------------------------------------------------------------
# Parse truth fields
# ---------------------------------------------------------------------

records = []

for row in sample_rows:

    if len(row) < 4:
        fail(
            f"Sample row has fewer than four columns: {row}"
        )

    run_accession = row[0].strip()

    try:
        major_proportion = float(row[1])
    except ValueError:
        fail(
            f"Invalid major proportion for "
            f"{run_accession}: {row[1]}"
        )

    published_lineages = row[2].strip()
    published_spoligotypes = row[3].strip()

    major_lineage, minor_lineage = split_pair(
        published_lineages
    )

    major_spoligotype, minor_spoligotype = split_pair(
        published_spoligotypes
    )

    # -------------------------------------------------------------
    # Derive expected status from known composition
    # -------------------------------------------------------------

    if major_proportion == 1.0:
        expected_status = "pure"
        minor_proportion = 0.0
    elif major_proportion < 1.0:
        expected_status = "mixed"
        minor_proportion = round(
            1.0 - major_proportion,
            2,
        )
    else:
        fail(
            f"Unexpected major proportion > 1.0 "
            f"for {run_accession}: {major_proportion}"
        )

    records.append(
        {
            "run_accession": run_accession,
            "expected_status": expected_status,
            "published_major_proportion": (
                f"{major_proportion:.2f}"
            ),
            "expected_minor_proportion": (
                f"{minor_proportion:.2f}"
            ),
            "major_lineage": major_lineage,
            "minor_lineage": minor_lineage,
            "major_spoligotype": major_spoligotype,
            "minor_spoligotype": minor_spoligotype,
            "truth_basis": "known_sample_composition",
            "primary_truth_source": (
                "Sobkowiak_et_al_2018_Table_1"
            ),
        }
    )


# ---------------------------------------------------------------------
# Integrity checks
# ---------------------------------------------------------------------

accessions = [
    record["run_accession"]
    for record in records
]

if len(set(accessions)) != EXPECTED_TOTAL:
    duplicates = sorted(
        {
            accession
            for accession in accessions
            if accessions.count(accession) > 1
        }
    )

    fail(
        "Duplicate run accessions detected: "
        + ", ".join(duplicates)
    )


# Check expected proportion groups

proportion_counts = {}

for record in records:
    proportion = float(
        record["published_major_proportion"]
    )

    proportion_counts[proportion] = (
        proportion_counts.get(proportion, 0) + 1
    )

if proportion_counts != EXPECTED_PROPORTION_COUNTS:
    fail(
        "Unexpected major-proportion distribution.\n"
        f"Observed: {proportion_counts}\n"
        f"Expected: {EXPECTED_PROPORTION_COUNTS}"
    )


# Check mixed/pure counts

status_counts = {}

for record in records:
    status = record["expected_status"]

    status_counts[status] = (
        status_counts.get(status, 0) + 1
    )

if status_counts != EXPECTED_STATUS_COUNTS:
    fail(
        "Unexpected status distribution.\n"
        f"Observed: {status_counts}\n"
        f"Expected: {EXPECTED_STATUS_COUNTS}"
    )


# Check pure samples have no minor constituent

for record in records:

    if record["expected_status"] == "pure":

        if record["minor_lineage"]:
            fail(
                f"Pure sample {record['run_accession']} "
                "unexpectedly contains a minor lineage."
            )

        if record["minor_spoligotype"]:
            fail(
                f"Pure sample {record['run_accession']} "
                "unexpectedly contains a minor spoligotype."
            )


# Check mixed samples contain both constituents

for record in records:

    if record["expected_status"] == "mixed":

        if not record["minor_lineage"]:
            fail(
                f"Mixed sample {record['run_accession']} "
                "has no minor lineage."
            )

        if not record["minor_spoligotype"]:
            fail(
                f"Mixed sample {record['run_accession']} "
                "has no minor spoligotype."
            )


# ---------------------------------------------------------------------
# Write canonical truth set
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
    writer.writerows(records)


# ---------------------------------------------------------------------
# SHA256 of generated truth set
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
# Write report
# ---------------------------------------------------------------------

report_lines = [
    "Mixed-infection truth set construction report",
    "=============================================",
    "",
    f"Source: {SOURCE_FILE}",
    "Primary evidence: Sobkowiak et al. 2018 Table 1",
    "Truth basis: known sample composition",
    "",
    f"Total samples: {len(records)}",
    f"Unique run accessions: {len(set(accessions))}",
    "",
    "Expected status:",
    f"  mixed: {status_counts.get('mixed', 0)}",
    f"  pure: {status_counts.get('pure', 0)}",
    "",
    "Published major-strain proportions:",
]

for proportion in sorted(proportion_counts):
    report_lines.append(
        f"  {proportion:.2f}: "
        f"{proportion_counts[proportion]}"
    )

report_lines.extend(
    [
        "",
        "Truth assignment policy:",
        "  0.70 -> mixed (expected minor proportion 0.30)",
        "  0.90 -> mixed (expected minor proportion 0.10)",
        "  0.95 -> mixed (expected minor proportion 0.05)",
        "  1.00 -> pure  (expected minor proportion 0.00)",
        "",
        "The computational-result columns reported in the",
        "2018 Table 1 were not used to assign expected status.",
        "",
        f"Output: {OUTPUT_FILE}",
        f"SHA256: {sha256}",
        "",
        "STATUS: PASS",
    ]
)

REPORT_FILE.write_text(
    "\n".join(report_lines) + "\n",
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

print("Mixed-infection truth set constructed")
print("--------------------------------------")
print(f"Source rows:        {len(sample_rows)}")
print(f"Truth-set samples:  {len(records)}")
print(f"Unique accessions:  {len(set(accessions))}")
print(f"Mixed:              {status_counts['mixed']}")
print(f"Pure:               {status_counts['pure']}")
print()
print("Major-proportion groups:")

for proportion in sorted(proportion_counts):
    print(
        f"  {proportion:.2f}: "
        f"{proportion_counts[proportion]}"
    )

print()
print(f"Truth set: {OUTPUT_FILE}")
print(f"Report:    {REPORT_FILE}")
print(f"SHA256:    {sha256}")
print("STATUS: PASS")
