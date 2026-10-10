#!/usr/bin/env python3
"""
Download and verify paired FASTQ files for the selected lineage pilot panel.

Key behavior
------------
- Processes only manifest rows with manifest_status == "ready".
- Downloads over HTTPS.
- Does NOT attempt ranged resume requests.
- Removes partial files before each retry.
- Captures curl/wget error output.
- Verifies every completed file against the ENA MD5 checksum.
- Skips existing files only when their checksum is already valid.
- Writes an updated manifest and download report.

This script does not inspect barcode markers or confirm lineage labels.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "broad_lineage",
    "fastq_1_url",
    "fastq_2_url",
    "fastq_1_md5",
    "fastq_2_md5",
    "manifest_status",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download and MD5-verify paired FASTQ files for the lineage pilot."
        )
    )

    parser.add_argument(
        "manifest",
        type=Path,
        help="Path to pilot_single_lineage_fastq_manifest.csv or TSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lineage/pilot/fastqs"),
        help="Directory for downloaded FASTQs.",
    )

    parser.add_argument(
        "--report-dir",
        type=Path,
        default=Path("lineage/pilot/reports"),
        help="Directory for updated manifest and report.",
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Maximum fresh download attempts per file. Default: 3.",
    )

    parser.add_argument(
        "--retry-delay",
        type=float,
        default=5.0,
        help="Seconds between retries. Default: 5.",
    )

    parser.add_argument(
        "--connect-timeout",
        type=int,
        default=30,
        help="Connection timeout in seconds. Default: 30.",
    )

    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path, dtype=str).fillna("")

    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")

    raise ValueError(
        f"Unsupported manifest format {suffix!r}. Use CSV or TSV."
    )


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def md5_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)

    return digest.hexdigest()


def choose_downloader() -> str:
    if shutil.which("curl"):
        return "curl"

    if shutil.which("wget"):
        return "wget"

    raise RuntimeError("Neither curl nor wget is available in PATH.")


def filename_from_url(url: str) -> str:
    name = Path(url.split("?", 1)[0]).name

    if not name:
        raise ValueError(f"Could not derive filename from URL: {url}")

    return name


def verify_existing_file(
    path: Path,
    expected_md5: str,
) -> tuple[bool, str]:
    if not path.exists():
        return False, "missing"

    if not path.is_file():
        return False, "path_is_not_a_file"

    if path.stat().st_size == 0:
        return False, "empty_file"

    if not expected_md5:
        return False, "missing_expected_md5"

    observed = md5_file(path)

    if observed.lower() == expected_md5.lower():
        return True, observed

    return False, observed


def remove_partial(path: Path) -> None:
    if path.exists():
        path.unlink()


def run_download(
    downloader: str,
    url: str,
    destination: Path,
    connect_timeout: int,
) -> subprocess.CompletedProcess[str]:
    destination.parent.mkdir(parents=True, exist_ok=True)

    if downloader == "curl":
        command = [
            "curl",
            "--fail",
            "--location",
            "--show-error",
            "--connect-timeout",
            str(connect_timeout),
            "--output",
            str(destination),
            url,
        ]
    elif downloader == "wget":
        command = [
            "wget",
            "--timeout",
            str(connect_timeout),
            "--output-document",
            str(destination),
            url,
        ]
    else:
        raise ValueError(f"Unsupported downloader: {downloader}")

    print("  $ " + " ".join(command))

    return subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def classify_downloader_failure(
    returncode: int,
    stderr: str,
) -> str:
    message = clean_text(stderr)

    if returncode == 18:
        category = "transfer_interrupted"
    elif returncode == 22 and "403" in message:
        category = "http_403"
    elif returncode == 22 and "404" in message:
        category = "http_404"
    elif returncode == 28:
        category = "timeout"
    elif returncode == 6:
        category = "dns_resolution_failure"
    elif returncode == 7:
        category = "connection_failure"
    else:
        category = f"downloader_exit_{returncode}"

    if message:
        return f"{category}: {message}"

    return category


def download_and_verify(
    *,
    downloader: str,
    url: str,
    destination: Path,
    expected_md5: str,
    retries: int,
    retry_delay: float,
    connect_timeout: int,
) -> dict[str, object]:
    valid, observed = verify_existing_file(
        destination,
        expected_md5,
    )

    if valid:
        return {
            "status": "already_present_valid",
            "observed_md5": observed,
            "attempts": 0,
            "message": "Existing file passed MD5 verification.",
        }

    remove_partial(destination)

    last_message = ""

    for attempt in range(1, retries + 1):
        # Always start from a clean file. No Range/resume request is used.
        remove_partial(destination)

        result = run_download(
            downloader=downloader,
            url=url,
            destination=destination,
            connect_timeout=connect_timeout,
        )

        if result.returncode != 0:
            last_message = classify_downloader_failure(
                returncode=result.returncode,
                stderr=result.stderr,
            )
        elif not destination.exists():
            last_message = "downloader_finished_but_file_missing"
        else:
            valid, observed = verify_existing_file(
                destination,
                expected_md5,
            )

            if valid:
                return {
                    "status": "downloaded_valid",
                    "observed_md5": observed,
                    "attempts": attempt,
                    "message": (
                        "Download completed and MD5 matched."
                    ),
                }

            last_message = (
                "md5_mismatch: "
                f"observed={observed}; expected={expected_md5}"
            )

        remove_partial(destination)

        if attempt < retries:
            print(
                f"  Attempt {attempt}/{retries} failed: "
                f"{last_message}"
            )
            print(
                f"  Retrying from scratch in {retry_delay:g} seconds."
            )
            time.sleep(retry_delay)

    return {
        "status": "failed",
        "observed_md5": "",
        "attempts": retries,
        "message": last_message or "unknown_download_failure",
    }


def prepare_status_columns(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "local_fastq_1",
        "local_fastq_2",
        "fastq_1_download_status",
        "fastq_2_download_status",
        "fastq_1_observed_md5",
        "fastq_2_observed_md5",
        "fastq_1_download_attempts",
        "fastq_2_download_attempts",
        "fastq_1_download_message",
        "fastq_2_download_message",
        "sample_download_status",
    ]

    for column in columns:
        if column not in df.columns:
            df[column] = ""

    return df


def summarize_sample_status(row: pd.Series) -> str:
    valid_statuses = {
        "already_present_valid",
        "downloaded_valid",
    }

    r1 = clean_text(row["fastq_1_download_status"])
    r2 = clean_text(row["fastq_2_download_status"])

    if r1 in valid_statuses and r2 in valid_statuses:
        return "complete"

    if "failed" in {r1, r2}:
        return "failed"

    return "incomplete"


def write_report(
    path: Path,
    manifest_path: Path,
    output_dir: Path,
    results: pd.DataFrame,
    downloader: str,
) -> None:
    sample_counts = (
        results["sample_download_status"]
        .value_counts(dropna=False)
    )

    file_statuses = pd.concat(
        [
            results["fastq_1_download_status"],
            results["fastq_2_download_status"],
        ],
        ignore_index=True,
    ).value_counts(dropna=False)

    total_size = 0

    for column in ["local_fastq_1", "local_fastq_2"]:
        for value in results[column]:
            local_path = Path(clean_text(value))

            if local_path.exists() and local_path.is_file():
                total_size += local_path.stat().st_size

    lines = [
        "LINEAGE PILOT FASTQ DOWNLOAD REPORT",
        "=" * 37,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Input manifest: {manifest_path}",
        f"Output directory: {output_dir}",
        f"Downloader: {downloader}",
        "",
        "SAMPLE STATUS",
        f"complete: {int(sample_counts.get('complete', 0))}",
        f"failed: {int(sample_counts.get('failed', 0))}",
        f"incomplete: {int(sample_counts.get('incomplete', 0))}",
        "",
        "FILE STATUS",
        (
            "already_present_valid: "
            f"{int(file_statuses.get('already_present_valid', 0))}"
        ),
        (
            "downloaded_valid: "
            f"{int(file_statuses.get('downloaded_valid', 0))}"
        ),
        f"failed: {int(file_statuses.get('failed', 0))}",
        "",
        "LOCAL STORAGE",
        f"Verified bytes: {total_size}",
        f"Approximate GiB: {total_size / (1024 ** 3):.3f}",
        "",
        "FAILED OR INCOMPLETE SAMPLES",
    ]

    incomplete = results[
        results["sample_download_status"] != "complete"
    ]

    if incomplete.empty:
        lines.append("None")
    else:
        for row in incomplete.itertuples(index=False):
            lines.append(
                f"{row.sample_id}\t"
                f"{row.sample_download_status}\t"
                f"R1={row.fastq_1_download_status}\t"
                f"{row.fastq_1_download_message}\t"
                f"R2={row.fastq_2_download_status}\t"
                f"{row.fastq_2_download_message}"
            )

    lines.extend(
        [
            "",
            "NOTE",
            (
                "Every retry started from a clean local file. "
                "No ranged resume request was used. Files were accepted "
                "only after matching the ENA-provided MD5 checksum."
            ),
            (
                "This step does not confirm candidate lineage labels."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    manifest_path = args.manifest.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    report_dir = args.report_dir.expanduser().resolve()

    if not manifest_path.exists():
        print(
            f"ERROR: manifest not found: {manifest_path}",
            file=sys.stderr,
        )
        return 1

    if args.retries < 1:
        print(
            "ERROR: --retries must be at least 1.",
            file=sys.stderr,
        )
        return 1

    try:
        manifest = read_table(manifest_path)

        missing = REQUIRED_COLUMNS.difference(manifest.columns)

        if missing:
            raise ValueError(
                f"Manifest is missing required columns: {sorted(missing)}"
            )

        ready = manifest[
            manifest["manifest_status"].str.strip().eq("ready")
        ].copy()

        if ready.empty:
            raise ValueError(
                "Manifest contains no rows with manifest_status=ready."
            )

        if ready["sample_id"].duplicated().any():
            duplicates = sorted(
                ready.loc[
                    ready["sample_id"].duplicated(keep=False),
                    "sample_id",
                ].unique()
            )

            raise ValueError(
                f"Duplicate ready sample IDs: {duplicates}"
            )

        downloader = choose_downloader()

        output_dir.mkdir(parents=True, exist_ok=True)
        report_dir.mkdir(parents=True, exist_ok=True)

        ready = prepare_status_columns(ready)

        total = len(ready)

        for position, row_index in enumerate(
            ready.index,
            start=1,
        ):
            sample_id = clean_text(
                ready.at[row_index, "sample_id"]
            )

            print()
            print(
                f"=== [{position}/{total}] Processing {sample_id} ==="
            )

            for read_number in (1, 2):
                url_column = f"fastq_{read_number}_url"
                md5_column = f"fastq_{read_number}_md5"

                url = clean_text(
                    ready.at[row_index, url_column]
                )
                expected_md5 = clean_text(
                    ready.at[row_index, md5_column]
                )

                destination = Path("")

                if not url:
                    result = {
                        "status": "failed",
                        "observed_md5": "",
                        "attempts": 0,
                        "message": f"missing_{url_column}",
                    }
                elif not expected_md5:
                    result = {
                        "status": "failed",
                        "observed_md5": "",
                        "attempts": 0,
                        "message": f"missing_{md5_column}",
                    }
                else:
                    filename = filename_from_url(url)
                    destination = output_dir / filename

                    print(
                        f"  R{read_number}: {filename}"
                    )

                    result = download_and_verify(
                        downloader=downloader,
                        url=url,
                        destination=destination,
                        expected_md5=expected_md5,
                        retries=args.retries,
                        retry_delay=args.retry_delay,
                        connect_timeout=args.connect_timeout,
                    )

                ready.at[
                    row_index,
                    f"local_fastq_{read_number}",
                ] = str(destination) if str(destination) else ""

                ready.at[
                    row_index,
                    f"fastq_{read_number}_download_status",
                ] = result["status"]

                ready.at[
                    row_index,
                    f"fastq_{read_number}_observed_md5",
                ] = result["observed_md5"]

                ready.at[
                    row_index,
                    f"fastq_{read_number}_download_attempts",
                ] = result["attempts"]

                ready.at[
                    row_index,
                    f"fastq_{read_number}_download_message",
                ] = result["message"]

                print(
                    f"    {result['status']}: {result['message']}"
                )

            ready.at[
                row_index,
                "sample_download_status",
            ] = summarize_sample_status(
                ready.loc[row_index]
            )

            print(
                "  Sample status: "
                f"{ready.at[row_index, 'sample_download_status']}"
            )

        updated_manifest_path = (
            report_dir
            / "pilot_single_lineage_fastq_download_manifest.csv"
        )

        report_path = (
            report_dir
            / "pilot_single_lineage_fastq_download_report.txt"
        )

        ready = ready.sort_values(
            ["broad_lineage", "candidate_lineage", "sample_id"],
            kind="stable",
        ).reset_index(drop=True)

        ready.to_csv(
            updated_manifest_path,
            index=False,
        )

        write_report(
            path=report_path,
            manifest_path=manifest_path,
            output_dir=output_dir,
            results=ready,
            downloader=downloader,
        )

    except (
        ValueError,
        RuntimeError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    complete = int(
        ready["sample_download_status"].eq("complete").sum()
    )
    failed = int(
        ready["sample_download_status"].eq("failed").sum()
    )

    print()
    print("FASTQ processing complete.")
    print(f"Complete samples: {complete}")
    print(f"Failed samples:   {failed}")
    print(f"Download manifest: {updated_manifest_path}")
    print(f"Report:            {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
