#!/usr/bin/env python3

"""
Download FASTQ files for the contamination pilot validation panel.

Input:
    contamination/reference_dataset/contamination_pilot_manifest.csv

Output:
    contamination/fastqs/
    contamination/reports/contamination_fastq_download_report.txt

For each FASTQ:
    - use the ENA URL stored in the pilot manifest
    - verify the ENA MD5 checksum
    - skip an existing file if its checksum is correct
    - replace an existing file if its checksum is incorrect
    - download to a temporary .part file first
    - only rename to the final filename after checksum verification

Expected pilot panel:
    30 samples
    60 FASTQ files
"""

from pathlib import Path
from urllib.parse import urlparse
import hashlib
import shutil
import sys
import time

import pandas as pd
import requests


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

MANIFEST_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reference_dataset"
    / "contamination_pilot_manifest.csv"
)

FASTQ_DIR = (
    PROJECT_ROOT
    / "contamination"
    / "fastqs"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "contamination"
    / "reports"
    / "contamination_fastq_download_report.txt"
)


# ---------------------------------------------------------------------
# Expected panel size
# ---------------------------------------------------------------------

EXPECTED_SAMPLES = 30
EXPECTED_FASTQS = 60


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def calculate_md5(path, chunk_size=1024 * 1024):

    md5 = hashlib.md5()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            md5.update(chunk)

    return md5.hexdigest()


def normalise_url(url):

    """
    ENA metadata may contain ftp:// URLs.

    Convert ENA FTP URLs to HTTPS for downloading with requests.
    """

    if url.startswith("ftp://ftp.sra.ebi.ac.uk/"):
        return url.replace(
            "ftp://ftp.sra.ebi.ac.uk/",
            "https://ftp.sra.ebi.ac.uk/",
            1,
        )

    return url


def filename_from_url(url):

    parsed = urlparse(url)

    filename = Path(parsed.path).name

    if not filename:
        sys.exit(
            f"ERROR: Could not determine filename from URL:\n{url}"
        )

    return filename


# ---------------------------------------------------------------------
# Load manifest
# ---------------------------------------------------------------------

def load_manifest():

    if not MANIFEST_FILE.exists():
        sys.exit(
            f"ERROR: Pilot manifest not found:\n{MANIFEST_FILE}"
        )

    df = pd.read_csv(
        MANIFEST_FILE,
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "source_sample_id",
        "run_accession",
        "fastq_1_ftp",
        "fastq_2_ftp",
        "fastq_1_md5",
        "fastq_2_md5",
        "fastq_1_bytes",
        "fastq_2_bytes",
    }

    missing = required - set(df.columns)

    if missing:
        sys.exit(
            "ERROR: Missing required manifest columns:\n"
            + "\n".join(sorted(missing))
        )

    if len(df) != EXPECTED_SAMPLES:
        sys.exit(
            f"ERROR: Expected {EXPECTED_SAMPLES} pilot samples "
            f"but found {len(df)}."
        )

    if df["run_accession"].duplicated().any():
        sys.exit(
            "ERROR: Duplicate run accessions detected in pilot manifest."
        )

    return df


# ---------------------------------------------------------------------
# Build download table
# ---------------------------------------------------------------------

def build_download_table(manifest):

    records = []

    for _, row in manifest.iterrows():

        for read_number in (1, 2):

            url = row[f"fastq_{read_number}_ftp"]
            expected_md5 = row[f"fastq_{read_number}_md5"]
            expected_bytes = row[f"fastq_{read_number}_bytes"]

            if not url:
                sys.exit(
                    f"ERROR: Missing FASTQ URL for "
                    f"{row['run_accession']} R{read_number}."
                )

            if not expected_md5:
                sys.exit(
                    f"ERROR: Missing MD5 for "
                    f"{row['run_accession']} R{read_number}."
                )

            download_url = normalise_url(url)
            filename = filename_from_url(download_url)

            records.append(
                {
                    "source_sample_id": row["source_sample_id"],
                    "run_accession": row["run_accession"],
                    "read": f"R{read_number}",
                    "url": download_url,
                    "filename": filename,
                    "expected_md5": expected_md5.lower(),
                    "expected_bytes": expected_bytes,
                }
            )

    downloads = pd.DataFrame(records)

    if len(downloads) != EXPECTED_FASTQS:
        sys.exit(
            f"ERROR: Expected {EXPECTED_FASTQS} FASTQ records "
            f"but constructed {len(downloads)}."
        )

    if downloads["filename"].duplicated().any():
        duplicates = downloads[
            downloads["filename"].duplicated(keep=False)
        ]["filename"].tolist()

        sys.exit(
            "ERROR: Duplicate FASTQ filenames detected:\n"
            + "\n".join(duplicates)
        )

    return downloads


# ---------------------------------------------------------------------
# Download one FASTQ
# ---------------------------------------------------------------------

def download_file(record):

    final_path = FASTQ_DIR / record["filename"]
    temp_path = FASTQ_DIR / f"{record['filename']}.part"

    expected_md5 = record["expected_md5"]

    # -------------------------------------------------------------
    # Existing final file
    # -------------------------------------------------------------

    if final_path.exists():

        print(f"Checking existing file: {final_path.name}")

        observed_md5 = calculate_md5(final_path)

        if observed_md5 == expected_md5:

            print("  checksum OK — skipping")

            return {
                **record,
                "status": "existing_verified",
                "observed_md5": observed_md5,
                "final_bytes": final_path.stat().st_size,
            }

        print("  checksum FAILED — removing existing file")

        final_path.unlink()

    # -------------------------------------------------------------
    # Remove stale partial download
    # -------------------------------------------------------------

    if temp_path.exists():

        print(
            f"  removing incomplete partial file: "
            f"{temp_path.name}"
        )

        temp_path.unlink()

    # -------------------------------------------------------------
    # Download
    # -------------------------------------------------------------

    print(
        f"Downloading {record['run_accession']} "
        f"{record['read']}: {record['filename']}"
    )

    try:

        with requests.get(
            record["url"],
            stream=True,
            timeout=(30, 300),
        ) as response:

            response.raise_for_status()

            with temp_path.open("wb") as handle:

                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):

                    if chunk:
                        handle.write(chunk)

    except requests.RequestException as exc:

        if temp_path.exists():
            temp_path.unlink()

        raise RuntimeError(
            f"Download failed for {record['filename']}: {exc}"
        )

    # -------------------------------------------------------------
    # Verify checksum
    # -------------------------------------------------------------

    observed_md5 = calculate_md5(temp_path)

    if observed_md5 != expected_md5:

        temp_path.unlink()

        raise RuntimeError(
            f"MD5 mismatch for {record['filename']}\n"
            f"Expected: {expected_md5}\n"
            f"Observed: {observed_md5}"
        )

    # -------------------------------------------------------------
    # Optional byte-size check
    # -------------------------------------------------------------

    final_bytes = temp_path.stat().st_size

    if record["expected_bytes"]:

        try:
            expected_bytes = int(record["expected_bytes"])

            if final_bytes != expected_bytes:

                temp_path.unlink()

                raise RuntimeError(
                    f"File-size mismatch for {record['filename']}\n"
                    f"Expected: {expected_bytes}\n"
                    f"Observed: {final_bytes}"
                )

        except ValueError:
            raise RuntimeError(
                f"Invalid expected byte count for "
                f"{record['filename']}: "
                f"{record['expected_bytes']}"
            )

    # -------------------------------------------------------------
    # Finalise
    # -------------------------------------------------------------

    shutil.move(
        str(temp_path),
        str(final_path),
    )

    print("  checksum OK")

    return {
        **record,
        "status": "downloaded_verified",
        "observed_md5": observed_md5,
        "final_bytes": final_bytes,
    }


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def write_report(results, elapsed_seconds):

    result_df = pd.DataFrame(results)

    status_counts = (
        result_df["status"]
        .value_counts()
        .sort_index()
    )

    total_bytes = result_df["final_bytes"].astype(int).sum()

    with REPORT_FILE.open("w") as handle:

        handle.write(
            "Contamination FASTQ download report\n"
        )
        handle.write("=" * 40 + "\n\n")

        handle.write(
            f"Pilot manifest: {MANIFEST_FILE}\n"
        )
        handle.write(
            f"Samples: {EXPECTED_SAMPLES}\n"
        )
        handle.write(
            f"FASTQ files: {len(result_df)}\n"
        )
        handle.write(
            f"Total verified bytes: {total_bytes}\n"
        )
        handle.write(
            f"Elapsed seconds: {elapsed_seconds:.1f}\n\n"
        )

        handle.write(
            "Download status\n"
        )
        handle.write("-" * 40 + "\n")
        handle.write(
            status_counts.to_string()
        )
        handle.write("\n\n")

        handle.write(
            "Verified files\n"
        )
        handle.write("-" * 40 + "\n")

        handle.write(
            result_df[
                [
                    "run_accession",
                    "read",
                    "filename",
                    "status",
                    "expected_md5",
                    "observed_md5",
                    "final_bytes",
                ]
            ].to_string(index=False)
        )

        handle.write("\n")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    FASTQ_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("Loading contamination pilot manifest...")

    manifest = load_manifest()

    downloads = build_download_table(manifest)

    print(f"Samples: {len(manifest)}")
    print(f"FASTQ files to verify/download: {len(downloads)}")
    print(f"Destination: {FASTQ_DIR}")
    print()

    results = []

    start_time = time.time()

    try:

        for _, record in downloads.iterrows():

            result = download_file(
                record.to_dict()
            )

            results.append(result)

    except Exception as exc:

        print()
        print("ERROR: FASTQ download process stopped.")
        print(str(exc))
        print()
        print(
            "Files already successfully downloaded and verified "
            "have been retained."
        )
        print(
            "Rerun this script to continue; verified files will "
            "be skipped."
        )

        sys.exit(1)

    elapsed = time.time() - start_time

    write_report(
        results,
        elapsed,
    )

    print()
    print("FASTQ download successfully completed.")
    print(f"Verified FASTQ files: {len(results)}")
    print(f"FASTQ directory: {FASTQ_DIR}")
    print(f"Report: {REPORT_FILE}")


if __name__ == "__main__":
    main()
