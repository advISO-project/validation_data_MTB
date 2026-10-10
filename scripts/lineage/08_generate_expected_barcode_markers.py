#!/usr/bin/env python3
"""
Generate expected terminal-lineage barcode markers for the selected pilot.

Inputs
------
1. pilot_single_lineage_candidates.csv
2. lineage_barcode_reference.tsv

For each pilot sample, the script retrieves only markers explicitly assigned
to the sample's exact candidate lineage.

Example
-------
Candidate lineage: 4.2.1.1

Markers included:
    lineage == 4.2.1.1

Markers excluded from the primary expectation set:
    4
    4.2
    4.2.1

Outputs
-------
pilot_expected_terminal_barcode_markers.tsv
pilot_expected_terminal_barcode_summary.csv
pilot_expected_terminal_barcode_report.txt

This script does not inspect FASTQs, BAMs, or VCFs and does not confirm
candidate lineages.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_PILOT_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "broad_lineage",
}

REQUIRED_BARCODE_COLUMNS = {
    "marker_id",
    "position_1based",
    "reference_accession",
    "reference_allele",
    "marker_allele",
    "lineage",
    "broad_lineage",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate expected barcode markers assigned only to each "
            "pilot sample's exact candidate lineage."
        )
    )

    parser.add_argument(
        "--pilot",
        type=Path,
        default=Path(
            "lineage/metadata/fastlin/"
            "pilot_single_lineage_candidates.csv"
        ),
        help="Pilot candidate CSV or TSV.",
    )

    parser.add_argument(
        "--barcode",
        type=Path,
        default=Path(
            "references/tbdb/processed/"
            "lineage_barcode_reference.tsv"
        ),
        help="Canonical lineage barcode TSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lineage/pilot/barcode_expectations"),
        help="Output directory.",
    )

    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path, dtype=str).fillna("")

    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")

    raise ValueError(
        f"Unsupported file format {suffix!r}: {path}"
    )


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""

    return str(value).strip()


def normalise_lineage(value: object) -> str:
    lineage = clean_text(value)

    if not lineage:
        return ""

    if lineage.lower().startswith("lineage"):
        lineage = lineage[len("lineage"):].strip()

    aliases = {
        "L1": "1",
        "L2": "2",
        "L3": "3",
        "L4": "4",
        "L5": "5",
        "L6": "6",
        "L7": "7",
        "L8": "8",
        "L9": "9",
    }

    return aliases.get(lineage.upper(), lineage)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def validate_inputs(
    pilot: pd.DataFrame,
    barcode: pd.DataFrame,
) -> None:
    missing_pilot = REQUIRED_PILOT_COLUMNS.difference(
        pilot.columns
    )

    if missing_pilot:
        raise ValueError(
            "Pilot input is missing required columns: "
            f"{sorted(missing_pilot)}"
        )

    missing_barcode = REQUIRED_BARCODE_COLUMNS.difference(
        barcode.columns
    )

    if missing_barcode:
        raise ValueError(
            "Barcode input is missing required columns: "
            f"{sorted(missing_barcode)}"
        )

    if pilot["sample_id"].duplicated().any():
        duplicates = sorted(
            pilot.loc[
                pilot["sample_id"].duplicated(keep=False),
                "sample_id",
            ].unique()
        )

        raise ValueError(
            f"Duplicate pilot sample IDs found: {duplicates}"
        )

    if barcode["marker_id"].duplicated().any():
        duplicates = sorted(
            barcode.loc[
                barcode["marker_id"].duplicated(keep=False),
                "marker_id",
            ].unique()
        )

        raise ValueError(
            f"Duplicate barcode marker IDs found: {duplicates[:10]}"
        )


def build_terminal_expectations(
    pilot: pd.DataFrame,
    barcode: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    barcode = barcode.copy()
    barcode["lineage"] = barcode["lineage"].map(
        normalise_lineage
    )

    expected_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for pilot_row in pilot.itertuples(index=False):
        sample_id = clean_text(pilot_row.sample_id)
        candidate_lineage = normalise_lineage(
            pilot_row.candidate_lineage
        )
        broad_lineage = clean_text(
            pilot_row.broad_lineage
        )

        if not candidate_lineage:
            raise ValueError(
                f"Sample {sample_id} has an empty candidate lineage."
            )

        matching = barcode[
            barcode["lineage"] == candidate_lineage
        ].copy()

        marker_count = len(matching)

        if marker_count > 0:
            for marker in matching.itertuples(index=False):
                expected_rows.append(
                    {
                        "sample_id": sample_id,
                        "candidate_lineage": candidate_lineage,
                        "candidate_broad_lineage": broad_lineage,
                        "supporting_lineage_level": candidate_lineage,
                        "marker_scope": "terminal_lineage_only",
                        "marker_id": marker.marker_id,
                        "position_1based": marker.position_1based,
                        "reference_accession": (
                            marker.reference_accession
                        ),
                        "reference_allele": (
                            marker.reference_allele
                        ),
                        "expected_marker_allele": (
                            marker.marker_allele
                        ),
                        "marker_broad_lineage": (
                            marker.broad_lineage
                        ),
                    }
                )

        summary_rows.append(
            {
                "sample_id": sample_id,
                "candidate_lineage": candidate_lineage,
                "broad_lineage": broad_lineage,
                "expected_terminal_marker_count": marker_count,
                "terminal_lineage_has_markers": marker_count > 0,
                "expectation_status": (
                    "ready"
                    if marker_count > 0
                    else "no_terminal_markers"
                ),
            }
        )

    expected = pd.DataFrame(expected_rows)
    summary = pd.DataFrame(summary_rows)

    if not expected.empty:
        duplicate_pairs = expected.duplicated(
            subset=["sample_id", "marker_id"],
            keep=False,
        )

        if duplicate_pairs.any():
            duplicates = expected.loc[
                duplicate_pairs,
                ["sample_id", "marker_id"],
            ]

            raise ValueError(
                "Duplicate sample-marker expectations were generated.\n"
                f"{duplicates.head(20).to_string(index=False)}"
            )

        expected = expected.sort_values(
            [
                "sample_id",
                "position_1based",
                "marker_id",
            ],
            kind="stable",
        ).reset_index(drop=True)

    summary = summary.sort_values(
        [
            "broad_lineage",
            "candidate_lineage",
            "sample_id",
        ],
        kind="stable",
    ).reset_index(drop=True)

    return expected, summary


def write_report(
    path: Path,
    *,
    pilot_path: Path,
    barcode_path: Path,
    pilot_sha256: str,
    barcode_sha256: str,
    expected: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    total_samples = len(summary)

    ready_samples = int(
        summary["expectation_status"].eq("ready").sum()
    )

    no_marker_samples = int(
        summary["expectation_status"]
        .eq("no_terminal_markers")
        .sum()
    )

    lines = [
        "PILOT TERMINAL-LINEAGE BARCODE MARKER REPORT",
        "=" * 45,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Pilot input: {pilot_path}",
        f"Pilot SHA-256: {pilot_sha256}",
        f"Barcode input: {barcode_path}",
        f"Barcode SHA-256: {barcode_sha256}",
        "",
        "SUMMARY",
        f"Pilot samples: {total_samples}",
        f"Samples with terminal-lineage markers: {ready_samples}",
        f"Samples with no terminal-lineage markers: {no_marker_samples}",
        f"Expected sample-marker rows: {len(expected)}",
        "",
        "PER-SAMPLE TERMINAL EXPECTATIONS",
    ]

    for row in summary.itertuples(index=False):
        lines.append(
            f"{row.sample_id}\t"
            f"{row.candidate_lineage}\t"
            f"terminal_markers="
            f"{row.expected_terminal_marker_count}\t"
            f"status={row.expectation_status}"
        )

    review = summary[
        summary["expectation_status"] != "ready"
    ]

    lines.extend(
        [
            "",
            "SAMPLES REQUIRING EXPECTATION REVIEW",
        ]
    )

    if review.empty:
        lines.append("None")
    else:
        for row in review.itertuples(index=False):
            lines.append(
                f"{row.sample_id}\t"
                f"{row.candidate_lineage}\t"
                f"status={row.expectation_status}"
            )

    lines.extend(
        [
            "",
            "INTERPRETATION",
            (
                "Expected markers were taken only from the frozen canonical "
                "barcode reference."
            ),
            (
                "For each pilot sample, only markers explicitly assigned to "
                "the exact candidate lineage were included."
            ),
            (
                "Markers assigned to ancestral lineage levels were excluded "
                "from the primary expectation set."
            ),
            (
                "This step defines expected evidence only. It does not inspect "
                "sequencing data or confirm any lineage label."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    pilot_path = args.pilot.expanduser().resolve()
    barcode_path = args.barcode.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not pilot_path.exists():
        print(
            f"ERROR: pilot input not found: {pilot_path}",
            file=sys.stderr,
        )
        return 1

    if not barcode_path.exists():
        print(
            f"ERROR: barcode input not found: {barcode_path}",
            file=sys.stderr,
        )
        return 1

    try:
        pilot = read_table(pilot_path)
        barcode = read_table(barcode_path)

        validate_inputs(pilot, barcode)

        expected, summary = build_terminal_expectations(
            pilot=pilot,
            barcode=barcode,
        )

        pilot_sha256 = sha256_file(pilot_path)
        barcode_sha256 = sha256_file(barcode_path)

    except (
        ValueError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    expected_path = (
        output_dir
        / "pilot_expected_terminal_barcode_markers.tsv"
    )

    summary_path = (
        output_dir
        / "pilot_expected_terminal_barcode_summary.csv"
    )

    report_path = (
        output_dir
        / "pilot_expected_terminal_barcode_report.txt"
    )

    expected.to_csv(
        expected_path,
        sep="\t",
        index=False,
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    write_report(
        path=report_path,
        pilot_path=pilot_path,
        barcode_path=barcode_path,
        pilot_sha256=pilot_sha256,
        barcode_sha256=barcode_sha256,
        expected=expected,
        summary=summary,
    )

    print(f"Pilot samples:        {len(summary)}")
    print(f"Expected marker rows: {len(expected)}")
    print(
        "Samples ready:       "
        f"{int(summary['expectation_status'].eq('ready').sum())}"
    )
    print()
    print(f"Expected markers: {expected_path}")
    print(f"Summary:          {summary_path}")
    print(f"Report:           {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
