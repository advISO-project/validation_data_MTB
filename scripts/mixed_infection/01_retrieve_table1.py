#!/usr/bin/env python3

"""
Retrieve and freeze Table 1 from Sobkowiak et al. (2018).

Source:
Sobkowiak B et al. Identifying mixed Mycobacterium tuberculosis
infections from whole genome sequence data.
BMC Genomics. 2018;19:613.

PMCID: PMC6092779
PMID: 30107785
DOI: 10.1186/s12864-018-4988-z

This script retrieves the published PMC Table 1 HTML and stores it
unchanged as a raw source resource.

No truth labels or derived values are generated at this stage.
"""

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import urllib.request


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

TABLE_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/"
    "PMC6092779/table/Tab1/"
)

RAW_DIR = Path("mixed_infection/raw")
METADATA_DIR = Path("mixed_infection/metadata")

RAW_FILE = RAW_DIR / "sobkowiak2018_table1.html"
CHECKSUM_FILE = METADATA_DIR / "sobkowiak2018_table1.sha256"
SOURCE_FILE = METADATA_DIR / "sobkowiak2018_table1_source.txt"


# ---------------------------------------------------------------------
# Prepare directories
# ---------------------------------------------------------------------

RAW_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Retrieve published table
# ---------------------------------------------------------------------

print("Retrieving Sobkowiak et al. 2018 Table 1")
print(f"Source: {TABLE_URL}")

request = urllib.request.Request(
    TABLE_URL,
    headers={
        "User-Agent": (
            "ISO-in-a-Box validation dataset resource retrieval"
        )
    },
)

with urllib.request.urlopen(request) as response:
    data = response.read()

if not data:
    raise RuntimeError("Retrieved Table 1 file is empty.")

RAW_FILE.write_bytes(data)


# ---------------------------------------------------------------------
# SHA256 checksum
# ---------------------------------------------------------------------

sha256 = hashlib.sha256(data).hexdigest()

CHECKSUM_FILE.write_text(
    f"{sha256}  {RAW_FILE}\n",
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Source provenance
# ---------------------------------------------------------------------

retrieved_utc = datetime.now(timezone.utc).isoformat()

source_metadata = f"""Resource: Sobkowiak et al. 2018 Table 1
Article title: Identifying mixed Mycobacterium tuberculosis infections from whole genome sequence data
Journal: BMC Genomics
Year: 2018
Volume: 19
Article: 613
DOI: 10.1186/s12864-018-4988-z
PMCID: PMC6092779
PMID: 30107785
Table: Table 1
Source URL: {TABLE_URL}
Retrieved UTC: {retrieved_utc}
Raw file: {RAW_FILE}
SHA256: {sha256}
Purpose: Published source resource for construction of the mixed-infection validation truth dataset.
"""


SOURCE_FILE.write_text(
    source_metadata,
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

print()
print("Retrieval complete")
print("------------------")
print(f"Raw file:       {RAW_FILE}")
print(f"Bytes:          {len(data):,}")
print(f"SHA256:         {sha256}")
print(f"Checksum file:  {CHECKSUM_FILE}")
print(f"Source record:  {SOURCE_FILE}")
