#!/usr/bin/env python3

"""
Download and verify FASTQ files for the mixed-infection validation set.

Input:
    mixed_infection/reference_dataset/
        mixed_infection_validation_manifest.csv

Output:
    mixed_infection/fastqs/

Reports:
    mixed_infection/reports/
        mixed_infection_fastq_download_manifest.csv
        mixed_infection_fastq_download_report.txt

Behaviour
---------
- Reads FASTQ locations and ENA MD5 checksums from the frozen
  validation manifest.
- Checks available disk space before downloading.
- Downloads paired FASTQs for all 48 runs.
- Verifies every downloaded file against the ENA MD5.
- Correct existing files are retained on reruns.
- Incorrect/incomplete existing files are replaced.
- Temporary downloads use a .part suffix.
- A failed download is not treated as complete.
"""

from pathlib import Path
import csv
import hashlib
import shutil
import sys
import urllib.request
from datetime import datetime, timezone


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

MANIFEST_FILE = Path(
    "mixed_infection/reference_dataset/"
    "mixed_infection_validation_manifest.csv"
)

FASTQ_DIR = Path(
    "mixed_infection/fastqs"
)

DOWNLOAD_MANIFEST = Path(
    "mixed_infection/reports/"
    "mixed_infection_fastq_download_manifest.csv"
)

REPORT_FILE = Path(
    "mixed_infection/reports/"
    "mixed_infection_fastq_download_report.txt"
)


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

EXPECTED_SAMPLES = 48
EXPECTED_FILES = 96

# Additional free-space buffer beyond the amount still requiring
# download.
SPACE_BUFFER_GIB = 2.0

CHUNK_SIZE = 1024 * 1024  # 1 MiB


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def md5sum(path):
    digest = hashlib.md5()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def ftp_to_https(path):
    """
    ENA returns paths such as:

        ftp.sra.ebi.ac.uk/vol1/fastq/...

    Use HTTPS for retrieval.
    """

    path = path.strip()

    if path.startswith("https://"):
        return path

    if path.startswith("ftp://"):
        path = path[len("ftp://"):]

    return "https://" + path


def download_file(url, destination):
    """
    Download to a temporary .part file and rename only after
    successful transfer.
    """

    temporary = Path(str(destination) + ".part")

    if temporary.exists():
        temporary.unlink()

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "ISO-in-a-Box validation dataset FASTQ retrieval"
        },
    )

    with urllib.request.urlopen(request) as response:
        with temporary.open("wb") as output:

            while True:
                chunk = response.read(CHUNK_SIZE)

                if not chunk:
                    break

                output.write(chunk)

    temporary.replace(destination)


def gib(number_bytes):
    return number_bytes / (1024 ** 3)


# ---------------------------------------------------------------------
# Read validation manifest
# ---------------------------------------------------------------------

if not MANIFEST_FILE.exists():
    fail(
        f"Validation manifest not found: {MANIFEST_FILE}"
    )


with MANIFEST_FILE.open(
    newline="",
    encoding="utf-8",
) as handle:

    records = list(csv.DictReader(handle))


if len(records) != EXPECTED_SAMPLES:
    fail(
        f"Expected {EXPECTED_SAMPLES} samples, "
        f"found {len(records)}."
    )


FASTQ_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

DOWNLOAD_MANIFEST.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ---------------------------------------------------------------------
# Expand sample manifest into one record per FASTQ
# ---------------------------------------------------------------------

files = []

for row in records:

    accession = row["run_accession"].strip()

    for read_number in ("1", "2"):

        remote_path = row[f"fastq_{read_number}"].strip()
        expected_md5 = row[f"md5_{read_number}"].strip().lower()

        try:
            expected_bytes = int(
                row[f"bytes_{read_number}"]
            )
        except ValueError:
            fail(
                f"{accession}: invalid bytes_{read_number}."
            )

        filename = Path(remote_path).name

        expected_filename = (
            f"{accession}_{read_number}.fastq.gz"
        )

        if filename != expected_filename:
            fail(
                f"{accession}: unexpected FASTQ filename "
                f"{filename!r}; expected "
                f"{expected_filename!r}."
            )

        if len(expected_md5) != 32:
            fail(
                f"{accession}: invalid MD5 for read "
                f"{read_number}: {expected_md5!r}"
            )

        files.append(
            {
                "run_accession": accession,
                "read": read_number,
                "remote_path": remote_path,
                "url": ftp_to_https(remote_path),
                "filename": filename,
                "local_path": FASTQ_DIR / filename,
                "expected_md5": expected_md5,
                "expected_bytes": expected_bytes,
            }
        )


if len(files) != EXPECTED_FILES:
    fail(
        f"Expected {EXPECTED_FILES} FASTQ files, "
        f"found {len(files)} in manifest."
    )


# ---------------------------------------------------------------------
# Pre-scan existing files
# ---------------------------------------------------------------------

print("Checking existing FASTQ files...")
print()

bytes_still_needed = 0

for item in files:

    path = item["local_path"]

    if not path.exists():
        item["initial_state"] = "missing"
        bytes_still_needed += item["expected_bytes"]
        continue

    observed_size = path.stat().st_size

    if observed_size != item["expected_bytes"]:
        item["initial_state"] = "invalid_size"
        bytes_still_needed += item["expected_bytes"]
        continue

    observed_md5 = md5sum(path)

    if observed_md5 == item["expected_md5"]:
        item["initial_state"] = "verified_existing"
    else:
        item["initial_state"] = "invalid_md5"
        bytes_still_needed += item["expected_bytes"]


# ---------------------------------------------------------------------
# Disk-space check
# ---------------------------------------------------------------------

disk = shutil.disk_usage(FASTQ_DIR)

free_bytes = disk.free

buffer_bytes = int(
    SPACE_BUFFER_GIB * (1024 ** 3)
)

required_bytes = (
    bytes_still_needed + buffer_bytes
)


print("Disk-space check")
print("----------------")
print(
    f"FASTQ data still required: "
    f"{gib(bytes_still_needed):.2f} GiB"
)
print(
    f"Safety buffer:             "
    f"{SPACE_BUFFER_GIB:.2f} GiB"
)
print(
    f"Required free space:       "
    f"{gib(required_bytes):.2f} GiB"
)
print(
    f"Available free space:      "
    f"{gib(free_bytes):.2f} GiB"
)
print()


if free_bytes < required_bytes:
    fail(
        "Insufficient free disk space for remaining "
        "FASTQ downloads plus safety buffer."
    )


# ---------------------------------------------------------------------
# Download / verify
# ---------------------------------------------------------------------

download_results = []

for number, item in enumerate(files, start=1):

    accession = item["run_accession"]
    read_number = item["read"]
    path = item["local_path"]

    print(
        f"[{number:02d}/{EXPECTED_FILES}] "
        f"{accession} read {read_number}"
    )

    initial_state = item["initial_state"]

    # -------------------------------------------------------------
    # Already verified
    # -------------------------------------------------------------

    if initial_state == "verified_existing":

        observed_bytes = path.stat().st_size
        observed_md5 = item["expected_md5"]

        status = "verified_existing"

        print("  existing file: MD5 OK")

    else:

        # ---------------------------------------------------------
        # Remove invalid existing file before replacement
        # ---------------------------------------------------------

        if path.exists():

            print(
                f"  existing file state: {initial_state}"
            )
            print("  removing invalid existing file")

            path.unlink()

        # ---------------------------------------------------------
        # Download
        # ---------------------------------------------------------

        print(
            f"  downloading "
            f"({gib(item['expected_bytes']):.2f} GiB)"
        )

        try:
            download_file(
                item["url"],
                path,
            )

        except Exception as error:

            part_file = Path(str(path) + ".part")

            if part_file.exists():
                part_file.unlink()

            print(
                f"  FAILED: {error}",
                file=sys.stderr,
            )

            download_results.append(
                {
                    "run_accession": accession,
                    "read": read_number,
                    "filename": item["filename"],
                    "expected_bytes":
                        item["expected_bytes"],
                    "observed_bytes": "",
                    "expected_md5":
                        item["expected_md5"],
                    "observed_md5": "",
                    "status": "download_failed",
                }
            )

            continue

        # ---------------------------------------------------------
        # Size verification
        # ---------------------------------------------------------

        observed_bytes = path.stat().st_size

        if observed_bytes != item["expected_bytes"]:

            print(
                "  FAILED: downloaded size does not match ENA"
            )

            download_results.append(
                {
                    "run_accession": accession,
                    "read": read_number,
                    "filename": item["filename"],
                    "expected_bytes":
                        item["expected_bytes"],
                    "observed_bytes":
                        observed_bytes,
                    "expected_md5":
                        item["expected_md5"],
                    "observed_md5": "",
                    "status": "size_failed",
                }
            )

            continue

        # ---------------------------------------------------------
        # MD5 verification
        # ---------------------------------------------------------

        print("  verifying MD5...")

        observed_md5 = md5sum(path)

        if observed_md5 != item["expected_md5"]:

            print("  FAILED: MD5 mismatch")

            status = "md5_failed"

        else:

            print("  MD5 OK")
            status = "downloaded_verified"


    download_results.append(
        {
            "run_accession": accession,
            "read": read_number,
            "filename": item["filename"],
            "expected_bytes":
                item["expected_bytes"],
            "observed_bytes":
                observed_bytes,
            "expected_md5":
                item["expected_md5"],
            "observed_md5":
                observed_md5,
            "status":
                status,
        }
    )


# ---------------------------------------------------------------------
# Write per-file download manifest
# ---------------------------------------------------------------------

download_fields = [
    "run_accession",
    "read",
    "filename",
    "expected_bytes",
    "observed_bytes",
    "expected_md5",
    "observed_md5",
    "status",
]


with DOWNLOAD_MANIFEST.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=download_fields,
    )

    writer.writeheader()
    writer.writerows(download_results)


# ---------------------------------------------------------------------
# Final verification
#
# Re-check every file from disk rather than trusting only the
# download loop.
# ---------------------------------------------------------------------

print()
print("Performing final full-set verification...")
print()


verified_files = 0
failed_files = []

for item in files:

    path = item["local_path"]

    if not path.exists():
        failed_files.append(
            (item["filename"], "missing")
        )
        continue

    if path.stat().st_size != item["expected_bytes"]:
        failed_files.append(
            (item["filename"], "size_mismatch")
        )
        continue

    observed_md5 = md5sum(path)

    if observed_md5 != item["expected_md5"]:
        failed_files.append(
            (item["filename"], "md5_mismatch")
        )
        continue

    verified_files += 1


# ---------------------------------------------------------------------
# Determine complete paired samples
# ---------------------------------------------------------------------

verified_by_sample = {}

for item in files:

    accession = item["run_accession"]
    path = item["local_path"]

    valid = False

    if path.exists():

        if path.stat().st_size == item["expected_bytes"]:

            if md5sum(path) == item["expected_md5"]:
                valid = True

    verified_by_sample.setdefault(
        accession,
        {},
    )[item["read"]] = valid


complete_samples = sum(
    1
    for reads in verified_by_sample.values()
    if reads.get("1") and reads.get("2")
)


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

retrieved_utc = datetime.now(
    timezone.utc
).isoformat()


report = [
    "Mixed-infection FASTQ download report",
    "=====================================",
    "",
    f"Manifest: {MANIFEST_FILE}",
    f"Completed UTC: {retrieved_utc}",
    "",
    f"Expected samples: {EXPECTED_SAMPLES}",
    f"Expected FASTQ files: {EXPECTED_FILES}",
    f"Verified FASTQ files: {verified_files}",
    f"Complete paired samples: {complete_samples}",
    f"Failed/missing FASTQ files: {len(failed_files)}",
    "",
    f"Download directory: {FASTQ_DIR}",
    f"Per-file manifest: {DOWNLOAD_MANIFEST}",
]


if failed_files:

    report.extend(
        [
            "",
            "Failed/missing files:",
        ]
    )

    for filename, reason in failed_files:
        report.append(
            f"  {filename}: {reason}"
        )


if (
    verified_files == EXPECTED_FILES
    and complete_samples == EXPECTED_SAMPLES
):

    final_status = "PASS"

else:

    final_status = "INCOMPLETE"


report.extend(
    [
        "",
        f"STATUS: {final_status}",
    ]
)


REPORT_FILE.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8",
)


# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

print("FASTQ retrieval summary")
print("-----------------------")
print(f"Expected files:          {EXPECTED_FILES}")
print(f"Verified files:          {verified_files}")
print(f"Complete paired samples: {complete_samples}")
print(f"Failed/missing files:    {len(failed_files)}")
print()
print(f"Report:   {REPORT_FILE}")
print(f"Manifest: {DOWNLOAD_MANIFEST}")
print(f"STATUS:   {final_status}")


if final_status != "PASS":
    sys.exit(1)
