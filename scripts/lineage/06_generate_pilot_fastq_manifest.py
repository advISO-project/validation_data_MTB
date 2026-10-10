#!/usr/bin/env python3
"""
Generate an ENA FASTQ manifest for the selected single-lineage pilot panel.

Input
-----
pilot_single_lineage_candidates.csv

For each selected run accession, the script queries the ENA file-report API
and retrieves read-run metadata, including FASTQ locations and checksums.

Outputs
-------
pilot_single_lineage_fastq_manifest.csv
pilot_single_lineage_fastq_manifest_report.txt

This script does not download FASTQs and does not confirm lineage labels.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


ENA_FILE_REPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
)

ENA_FIELDS = [
    "run_accession",
    "sample_accession",
    "study_accession",
    "experiment_accession",
    "scientific_name",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]

REQUIRED_INPUT_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "broad_lineage",
    "pilot_selection_reason",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve pilot run accessions to ENA FASTQ locations and create "
            "a reproducible manifest."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Path to pilot_single_lineage_candidates.csv or TSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for outputs. Default: current directory.",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="HTTP request timeout in seconds. Default: 60.",
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of ENA request attempts per sample. Default: 3.",
    )

    parser.add_argument(
        "--request-delay",
        type=float,
        default=0.2,
        help="Delay between ENA requests in seconds. Default: 0.2.",
    )

    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path, dtype=str)

    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str)

    raise ValueError(
        f"Unsupported input format {suffix!r}. Use CSV or TSV."
    )


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""

    return str(value).strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def split_semicolon_field(value: object) -> list[str]:
    text = clean_text(value)

    if not text:
        return []

    return [
        item.strip()
        for item in text.split(";")
        if item.strip()
    ]


def to_https(url: str) -> str:
    url = clean_text(url)

    if not url:
        return ""

    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]

    if url.startswith("http://") or url.startswith("https://"):
        return url

    return "https://" + url


def fetch_ena_record(
    run_accession: str,
    session: requests.Session,
    timeout: int,
    retries: int,
) -> tuple[pd.DataFrame | None, str]:
    params = {
        "accession": run_accession,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
    }

    last_error = ""

    for attempt in range(1, retries + 1):
        try:
            response = session.get(
                ENA_FILE_REPORT_URL,
                params=params,
                timeout=timeout,
                headers={
                    "Accept": "text/tab-separated-values",
                    "User-Agent": (
                        "validation-data-MTB-pilot-manifest/1.0"
                    ),
                },
            )
            response.raise_for_status()

            if not response.text.strip():
                last_error = "empty ENA response"
                continue

            dataframe = pd.read_csv(
                io.StringIO(response.text),
                sep="\t",
                dtype=str,
            ).fillna("")

            if dataframe.empty:
                last_error = "ENA returned no read_run rows"
                continue

            return dataframe, ""

        except (
            requests.RequestException,
            pd.errors.ParserError,
        ) as error:
            last_error = (
                f"attempt {attempt}/{retries}: {error}"
            )

            if attempt < retries:
                time.sleep(min(attempt, 3))

    return None, last_error


def classify_fastq_record(
    library_layout: str,
    urls: list[str],
) -> tuple[str, str]:
    """
    Return manifest_status and status_reason.
    """

    layout = clean_text(library_layout).upper()

    if not urls:
        return "unresolved", "no_fastq_urls"

    if layout == "PAIRED":
        if len(urls) == 2:
            return "ready", "paired_fastq_urls_resolved"

        if len(urls) == 1:
            return (
                "review",
                "paired_library_has_one_fastq_url",
            )

        return (
            "review",
            f"paired_library_has_{len(urls)}_fastq_urls",
        )

    if layout == "SINGLE":
        return (
            "excluded_nonpaired",
            "single_end_library",
        )

    return (
        "review",
        f"unrecognised_library_layout:{layout or 'missing'}",
    )


def prepare_manifest_row(
    pilot_row: pd.Series,
    ena_row: pd.Series,
) -> dict[str, object]:
    urls = [
        to_https(url)
        for url in split_semicolon_field(
            ena_row.get("fastq_ftp", "")
        )
    ]

    md5_values = split_semicolon_field(
        ena_row.get("fastq_md5", "")
    )

    byte_values = split_semicolon_field(
        ena_row.get("fastq_bytes", "")
    )

    status, reason = classify_fastq_record(
        library_layout=ena_row.get("library_layout", ""),
        urls=urls,
    )

    fastq_1 = urls[0] if len(urls) >= 1 else ""
    fastq_2 = urls[1] if len(urls) >= 2 else ""

    md5_1 = md5_values[0] if len(md5_values) >= 1 else ""
    md5_2 = md5_values[1] if len(md5_values) >= 2 else ""

    bytes_1 = byte_values[0] if len(byte_values) >= 1 else ""
    bytes_2 = byte_values[1] if len(byte_values) >= 2 else ""

    return {
        "sample_id": clean_text(pilot_row["sample_id"]),
        "candidate_lineage": clean_text(
            pilot_row["candidate_lineage"]
        ),
        "broad_lineage": clean_text(
            pilot_row["broad_lineage"]
        ),
        "stage1_selection_policy": clean_text(
            pilot_row.get("stage1_selection_policy", "")
        ),
        "pilot_selection_reason": clean_text(
            pilot_row["pilot_selection_reason"]
        ),
        "run_accession": clean_text(
            ena_row.get("run_accession", "")
        ),
        "sample_accession": clean_text(
            ena_row.get("sample_accession", "")
        ),
        "study_accession": clean_text(
            ena_row.get("study_accession", "")
        ),
        "experiment_accession": clean_text(
            ena_row.get("experiment_accession", "")
        ),
        "scientific_name": clean_text(
            ena_row.get("scientific_name", "")
        ),
        "library_layout": clean_text(
            ena_row.get("library_layout", "")
        ),
        "instrument_platform": clean_text(
            ena_row.get("instrument_platform", "")
        ),
        "instrument_model": clean_text(
            ena_row.get("instrument_model", "")
        ),
        "fastq_1_url": fastq_1,
        "fastq_2_url": fastq_2,
        "fastq_1_md5": md5_1,
        "fastq_2_md5": md5_2,
        "fastq_1_bytes": bytes_1,
        "fastq_2_bytes": bytes_2,
        "fastq_url_count": len(urls),
        "manifest_status": status,
        "manifest_status_reason": reason,
        "candidate_status": (
            "candidate_not_yet_barcode_confirmed"
        ),
    }


def unresolved_manifest_row(
    pilot_row: pd.Series,
    reason: str,
) -> dict[str, object]:
    return {
        "sample_id": clean_text(pilot_row["sample_id"]),
        "candidate_lineage": clean_text(
            pilot_row["candidate_lineage"]
        ),
        "broad_lineage": clean_text(
            pilot_row["broad_lineage"]
        ),
        "stage1_selection_policy": clean_text(
            pilot_row.get("stage1_selection_policy", "")
        ),
        "pilot_selection_reason": clean_text(
            pilot_row["pilot_selection_reason"]
        ),
        "run_accession": clean_text(
            pilot_row["sample_id"]
        ),
        "sample_accession": "",
        "study_accession": "",
        "experiment_accession": "",
        "scientific_name": "",
        "library_layout": "",
        "instrument_platform": "",
        "instrument_model": "",
        "fastq_1_url": "",
        "fastq_2_url": "",
        "fastq_1_md5": "",
        "fastq_2_md5": "",
        "fastq_1_bytes": "",
        "fastq_2_bytes": "",
        "fastq_url_count": 0,
        "manifest_status": "unresolved",
        "manifest_status_reason": reason,
        "candidate_status": (
            "candidate_not_yet_barcode_confirmed"
        ),
    }


def write_report(
    path: Path,
    input_path: Path,
    input_sha256: str,
    manifest: pd.DataFrame,
) -> None:
    status_counts = (
        manifest["manifest_status"]
        .value_counts(dropna=False)
    )

    broad_ready = (
        manifest[
            manifest["manifest_status"] == "ready"
        ]
        .groupby("broad_lineage")["sample_id"]
        .nunique()
        .sort_index()
    )

    total_bytes = pd.to_numeric(
        manifest["fastq_1_bytes"],
        errors="coerce",
    ).fillna(0).sum()

    total_bytes += pd.to_numeric(
        manifest["fastq_2_bytes"],
        errors="coerce",
    ).fillna(0).sum()

    lines = [
        "PILOT SINGLE-LINEAGE FASTQ MANIFEST REPORT",
        "=" * 45,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Input pilot: {input_path}",
        f"Input SHA-256: {input_sha256}",
        f"Pilot rows: {len(manifest)}",
        "",
        "MANIFEST STATUS",
    ]

    for status in [
        "ready",
        "review",
        "excluded_nonpaired",
        "unresolved",
    ]:
        lines.append(
            f"{status}: {int(status_counts.get(status, 0))}"
        )

    lines.extend(
        [
            "",
            "READY-SAMPLE BROAD-LINEAGE COUNTS",
        ]
    )

    if broad_ready.empty:
        lines.append("None")
    else:
        for broad_lineage, count in broad_ready.items():
            lines.append(
                f"{broad_lineage}: {int(count)}"
            )

    lines.extend(
        [
            "",
            "ESTIMATED DOWNLOAD SIZE",
            f"Total bytes from ENA metadata: {int(total_bytes)}",
            f"Approximate GiB: {total_bytes / (1024 ** 3):.3f}",
            "",
            "SAMPLES REQUIRING REVIEW",
        ]
    )

    review_rows = manifest[
        manifest["manifest_status"] != "ready"
    ]

    if review_rows.empty:
        lines.append("None")
    else:
        for row in review_rows.itertuples(index=False):
            lines.append(
                f"{row.sample_id}\t"
                f"{row.manifest_status}\t"
                f"{row.manifest_status_reason}"
            )

    lines.extend(
        [
            "",
            "NOTE",
            (
                "This manifest resolves candidate pilot accessions to ENA "
                "file locations. It does not download files or confirm "
                "candidate lineage labels."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_path.exists():
        print(
            f"ERROR: input file not found: {input_path}",
            file=sys.stderr,
        )
        return 1

    try:
        pilot = read_table(input_path)

        missing = REQUIRED_INPUT_COLUMNS.difference(
            pilot.columns
        )

        if missing:
            raise ValueError(
                "Input pilot is missing required columns: "
                f"{sorted(missing)}"
            )

        if pilot["sample_id"].duplicated().any():
            duplicates = sorted(
                pilot.loc[
                    pilot["sample_id"].duplicated(keep=False),
                    "sample_id",
                ].unique()
            )

            raise ValueError(
                f"Duplicate pilot sample IDs: {duplicates}"
            )

        manifest_rows: list[dict[str, object]] = []

        session = requests.Session()

        for index, pilot_row in pilot.iterrows():
            run_accession = clean_text(
                pilot_row["sample_id"]
            )

            print(
                f"[{index + 1}/{len(pilot)}] "
                f"Resolving {run_accession}"
            )

            ena_df, error = fetch_ena_record(
                run_accession=run_accession,
                session=session,
                timeout=args.timeout,
                retries=args.retries,
            )

            if ena_df is None:
                manifest_rows.append(
                    unresolved_manifest_row(
                        pilot_row,
                        reason=error or "ENA lookup failed",
                    )
                )
                time.sleep(args.request_delay)
                continue

            exact_rows = ena_df[
                ena_df["run_accession"]
                .astype(str)
                .str.strip()
                .str.upper()
                .eq(run_accession.upper())
            ]

            if len(exact_rows) != 1:
                manifest_rows.append(
                    unresolved_manifest_row(
                        pilot_row,
                        reason=(
                            "ENA lookup returned "
                            f"{len(exact_rows)} exact run rows"
                        ),
                    )
                )
                time.sleep(args.request_delay)
                continue

            manifest_rows.append(
                prepare_manifest_row(
                    pilot_row=pilot_row,
                    ena_row=exact_rows.iloc[0],
                )
            )

            time.sleep(args.request_delay)

        manifest = pd.DataFrame(manifest_rows)

        input_sha256 = sha256_file(input_path)

    except (
        ValueError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = (
        output_dir
        / "pilot_single_lineage_fastq_manifest.csv"
    )

    report_path = (
        output_dir
        / "pilot_single_lineage_fastq_manifest_report.txt"
    )

    manifest.to_csv(
        manifest_path,
        index=False,
    )

    write_report(
        path=report_path,
        input_path=input_path,
        input_sha256=input_sha256,
        manifest=manifest,
    )

    status_counts = manifest[
        "manifest_status"
    ].value_counts()

    print()
    print(f"Manifest rows: {len(manifest)}")
    print(
        f"Ready:         "
        f"{int(status_counts.get('ready', 0))}"
    )
    print(
        f"Review:        "
        f"{int(status_counts.get('review', 0))}"
    )
    print(
        f"Non-paired:    "
        f"{int(status_counts.get('excluded_nonpaired', 0))}"
    )
    print(
        f"Unresolved:    "
        f"{int(status_counts.get('unresolved', 0))}"
    )
    print()
    print(f"Manifest: {manifest_path}")
    print(f"Report:   {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
