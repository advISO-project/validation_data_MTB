#!/usr/bin/env python3
"""
Merge TB-Profiler and fastlin sample-level lineage results from
fastlin Supplementary Material 1.

This script performs data integration only. It does not:
- select candidate samples,
- exclude mixtures,
- apply quality thresholds,
- confirm barcode SNPs, or
- assign truth labels.

Default input sheets:
    TB-profiler all
    fastlin default - 3139 SRA

Outputs:
    candidate_dataset_master.csv
    candidate_dataset_merge_report.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


DEFAULT_TBPROFILER_SHEET = "TB-profiler all"
DEFAULT_FASTLIN_SHEET = "fastlin default - 3139 SRA"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge TB-Profiler and fastlin lineage results by SRA accession. "
            "No samples are filtered."
        )
    )
    parser.add_argument(
        "workbook",
        type=Path,
        help="Path to the fastlin supplementary Excel workbook.",
    )
    parser.add_argument(
        "--tbprofiler-sheet",
        default=DEFAULT_TBPROFILER_SHEET,
        help=f"TB-Profiler sheet name (default: {DEFAULT_TBPROFILER_SHEET!r}).",
    )
    parser.add_argument(
        "--fastlin-sheet",
        default=DEFAULT_FASTLIN_SHEET,
        help=f"fastlin sheet name (default: {DEFAULT_FASTLIN_SHEET!r}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for output files (default: current directory).",
    )
    parser.add_argument(
        "--output-prefix",
        default="candidate_dataset_master",
        help="Output filename prefix (default: candidate_dataset_master).",
    )
    return parser.parse_args()


def clean_text(value: object) -> str | pd.NA:
    """Return stripped text, or pd.NA for blank/missing values."""
    if pd.isna(value):
        return pd.NA

    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "n/a"}:
        return pd.NA

    return text


def normalize_sample_id(value: object) -> str | pd.NA:
    """Normalize an SRA accession for matching."""
    text = clean_text(value)
    if pd.isna(text):
        return pd.NA
    return str(text).upper()


def remove_fastlin_abundance(lineage_value: object) -> str | pd.NA:
    """
    Remove fastlin abundance values from lineage output.

    Examples:
        '4.6.5 (83)' -> '4.6.5'
        '4.1 (6), 1.1.2 (59)' -> '4.1;1.1.2'

    Multiple calls are retained and separated with semicolons.
    """
    text = clean_text(lineage_value)
    if pd.isna(text):
        return pd.NA

    calls = []
    for item in str(text).split(","):
        item = re.sub(r"\s*\([^()]*\)\s*$", "", item).strip()
        if item:
            calls.append(item)

    return ";".join(calls) if calls else pd.NA


def normalize_lineage_list(value: object) -> str | pd.NA:
    """
    Normalize a lineage string without changing lineage names.

    Both TB-Profiler and fastlin multiple-lineage results are represented
    as semicolon-separated calls. Original order is retained.
    """
    text = clean_text(value)
    if pd.isna(text):
        return pd.NA

    parts = [
        part.strip()
        for part in re.split(r"[;,]", str(text))
        if part.strip()
    ]
    return ";".join(parts) if parts else pd.NA


def count_lineage_calls(value: object) -> int:
    """Count semicolon-separated lineage calls."""
    if pd.isna(value):
        return 0
    return len([part for part in str(value).split(";") if part.strip()])


def compare_lineages(row: pd.Series) -> str:
    """
    Compare normalized TB-Profiler and fastlin lineage results.

    This is descriptive only and is not a candidate-selection decision.
    """
    tb = row["tbprofiler_lineage"]
    fastlin = row["fastlin_lineage"]

    if row["_merge"] == "left_only":
        return "tbprofiler_only"
    if row["_merge"] == "right_only":
        return "fastlin_only"
    if pd.isna(tb) and pd.isna(fastlin):
        return "both_no_call"
    if pd.isna(tb):
        return "tbprofiler_no_call"
    if pd.isna(fastlin):
        return "fastlin_no_call"
    if tb == fastlin:
        return "exact_agreement"
    return "disagreement"


def read_and_prepare_tbprofiler(
    workbook: Path,
    sheet_name: str,
) -> pd.DataFrame:
    df = pd.read_excel(workbook, sheet_name=sheet_name)

    required = {"SRA ID", "lineages", "strain proportions"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"TB-Profiler sheet {sheet_name!r} is missing columns: "
            f"{sorted(missing)}"
        )

    df = df.loc[:, ["SRA ID", "lineages", "strain proportions"]].copy()
    df = df.rename(
        columns={
            "SRA ID": "sample_id",
            "lineages": "tbprofiler_lineage_raw",
            "strain proportions": "tbprofiler_strain_proportions",
        }
    )

    df["sample_id"] = df["sample_id"].map(normalize_sample_id)
    df["tbprofiler_lineage_raw"] = df["tbprofiler_lineage_raw"].map(clean_text)
    df["tbprofiler_lineage"] = (
        df["tbprofiler_lineage_raw"].map(normalize_lineage_list)
    )
    df["tbprofiler_lineage_count"] = (
        df["tbprofiler_lineage"].map(count_lineage_calls)
    )

    if df["sample_id"].isna().any():
        raise ValueError("TB-Profiler sheet contains missing sample IDs.")
    if df["sample_id"].duplicated().any():
        duplicates = sorted(
            df.loc[df["sample_id"].duplicated(keep=False), "sample_id"].unique()
        )
        raise ValueError(
            "TB-Profiler sheet contains duplicate sample IDs: "
            f"{duplicates[:10]}"
        )

    return df


def read_and_prepare_fastlin(
    workbook: Path,
    sheet_name: str,
) -> pd.DataFrame:
    df = pd.read_excel(workbook, sheet_name=sheet_name)

    required = {
        "#sample",
        "nb_files",
        "k_cov",
        "mixture",
        "lineages",
        "log_barcodes",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"fastlin sheet {sheet_name!r} is missing columns: "
            f"{sorted(missing)}"
        )

    df = df.loc[
        :,
        [
            "#sample",
            "nb_files",
            "k_cov",
            "mixture",
            "lineages",
            "log_barcodes",
        ],
    ].copy()

    df = df.rename(
        columns={
            "#sample": "sample_id",
            "nb_files": "fastlin_nb_files",
            "k_cov": "fastlin_k_cov",
            "mixture": "fastlin_mixture",
            "lineages": "fastlin_lineage_raw",
            "log_barcodes": "fastlin_log_barcodes",
        }
    )

    df["sample_id"] = df["sample_id"].map(normalize_sample_id)
    df["fastlin_mixture"] = (
        df["fastlin_mixture"].map(clean_text).astype("string").str.lower()
    )
    df["fastlin_lineage_raw"] = df["fastlin_lineage_raw"].map(clean_text)
    df["fastlin_lineage"] = (
        df["fastlin_lineage_raw"]
        .map(remove_fastlin_abundance)
        .map(normalize_lineage_list)
    )
    df["fastlin_lineage_count"] = (
        df["fastlin_lineage"].map(count_lineage_calls)
    )

    if df["sample_id"].isna().any():
        raise ValueError("fastlin sheet contains missing sample IDs.")
    if df["sample_id"].duplicated().any():
        duplicates = sorted(
            df.loc[df["sample_id"].duplicated(keep=False), "sample_id"].unique()
        )
        raise ValueError(
            "fastlin sheet contains duplicate sample IDs: "
            f"{duplicates[:10]}"
        )

    return df


def build_master_table(
    tbprofiler: pd.DataFrame,
    fastlin: pd.DataFrame,
) -> pd.DataFrame:
    master = tbprofiler.merge(
        fastlin,
        on="sample_id",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )

    master["dataset_presence"] = master["_merge"].map(
        {
            "both": "both",
            "left_only": "tbprofiler_only",
            "right_only": "fastlin_only",
        }
    )
    master["lineage_comparison"] = master.apply(compare_lineages, axis=1)

    # Descriptive flags only; no rows are filtered.
    master["tbprofiler_multiple_lineages"] = (
        master["tbprofiler_lineage_count"].fillna(0).astype(int) > 1
    )
    master["fastlin_multiple_lineages"] = (
        master["fastlin_lineage_count"].fillna(0).astype(int) > 1
    )
    master["either_method_multiple_lineages"] = (
        master["tbprofiler_multiple_lineages"]
        | master["fastlin_multiple_lineages"]
    )

    column_order = [
        "sample_id",
        "dataset_presence",
        "lineage_comparison",
        "tbprofiler_lineage_raw",
        "tbprofiler_lineage",
        "tbprofiler_lineage_count",
        "tbprofiler_strain_proportions",
        "tbprofiler_multiple_lineages",
        "fastlin_lineage_raw",
        "fastlin_lineage",
        "fastlin_lineage_count",
        "fastlin_mixture",
        "fastlin_multiple_lineages",
        "either_method_multiple_lineages",
        "fastlin_nb_files",
        "fastlin_k_cov",
        "fastlin_log_barcodes",
    ]

    master = master.loc[:, column_order]
    return master.sort_values("sample_id", kind="stable").reset_index(drop=True)


def write_report(
    report_path: Path,
    workbook: Path,
    tbprofiler_sheet: str,
    fastlin_sheet: str,
    tbprofiler: pd.DataFrame,
    fastlin: pd.DataFrame,
    master: pd.DataFrame,
) -> None:
    presence_counts = master["dataset_presence"].value_counts(dropna=False)
    comparison_counts = master["lineage_comparison"].value_counts(dropna=False)

    lines = [
        "FASTLIN / TB-PROFILER MERGE REPORT",
        "=" * 39,
        "",
        f"Input workbook: {workbook}",
        f"TB-Profiler sheet: {tbprofiler_sheet}",
        f"fastlin sheet: {fastlin_sheet}",
        "",
        "INPUT COUNTS",
        f"TB-Profiler rows: {len(tbprofiler)}",
        f"fastlin rows: {len(fastlin)}",
        f"Master rows after outer merge: {len(master)}",
        "",
        "DATASET PRESENCE",
    ]

    for category in ["both", "tbprofiler_only", "fastlin_only"]:
        lines.append(f"{category}: {int(presence_counts.get(category, 0))}")

    lines.extend(["", "LINEAGE COMPARISON"])
    comparison_order = [
        "exact_agreement",
        "disagreement",
        "tbprofiler_no_call",
        "fastlin_no_call",
        "both_no_call",
        "tbprofiler_only",
        "fastlin_only",
    ]
    for category in comparison_order:
        lines.append(f"{category}: {int(comparison_counts.get(category, 0))}")

    lines.extend(
        [
            "",
            "DESCRIPTIVE COUNTS",
            (
                "fastlin mixture=yes: "
                f"{int(master['fastlin_mixture'].eq('yes').sum())}"
            ),
            (
                "TB-Profiler multiple-lineage rows: "
                f"{int(master['tbprofiler_multiple_lineages'].sum())}"
            ),
            (
                "fastlin multiple-lineage rows: "
                f"{int(master['fastlin_multiple_lineages'].sum())}"
            ),
            "",
            "NOTE",
            (
                "No samples were filtered and no truth labels were assigned. "
                "The comparison fields are descriptive only."
            ),
        ]
    )

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_arguments()

    workbook = args.workbook.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not workbook.exists():
        print(f"ERROR: workbook not found: {workbook}", file=sys.stderr)
        return 1
    if workbook.suffix.lower() not in {".xlsx", ".xls"}:
        print(
            f"ERROR: expected an Excel workbook, received: {workbook.name}",
            file=sys.stderr,
        )
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        tbprofiler = read_and_prepare_tbprofiler(
            workbook,
            args.tbprofiler_sheet,
        )
        fastlin = read_and_prepare_fastlin(
            workbook,
            args.fastlin_sheet,
        )
        master = build_master_table(tbprofiler, fastlin)
    except (ValueError, KeyError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    csv_path = output_dir / f"{args.output_prefix}.csv"
    report_path = output_dir / f"{args.output_prefix}_merge_report.txt"

    master.to_csv(csv_path, index=False)
    write_report(
        report_path=report_path,
        workbook=workbook,
        tbprofiler_sheet=args.tbprofiler_sheet,
        fastlin_sheet=args.fastlin_sheet,
        tbprofiler=tbprofiler,
        fastlin=fastlin,
        master=master,
    )

    print(f"Wrote {len(master):,} rows.")
    print(f"CSV:    {csv_path}")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
