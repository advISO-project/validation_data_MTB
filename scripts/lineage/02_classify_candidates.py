#!/usr/bin/env python3
"""
Classify merged fastlin and TB-Profiler lineage results into:

1. Single-lineage candidates
2. Mixed-lineage candidates
3. Excluded samples

This script does not confirm lineage truth and does not select a final panel.
It only identifies candidate samples for later independent barcode confirmation.

Expected input:
    candidate_dataset_master.csv

Default outputs:
    single_lineage_candidates.csv
    mixed_lineage_candidates.csv
    excluded_candidates.csv
    candidate_classification_master.csv
    candidate_classification_report.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "sample_id",
    "dataset_presence",
    "lineage_comparison",
    "tbprofiler_lineage",
    "tbprofiler_lineage_count",
    "tbprofiler_multiple_lineages",
    "fastlin_lineage",
    "fastlin_lineage_count",
    "fastlin_mixture",
    "fastlin_multiple_lineages",
    "fastlin_k_cov",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify merged fastlin and TB-Profiler results into "
            "single-lineage candidates, mixed-lineage candidates, "
            "and excluded samples."
        )
    )
    parser.add_argument(
        "merged_dataset",
        type=Path,
        help="Path to candidate_dataset_master.csv or .tsv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for output files (default: current directory).",
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t")

    raise ValueError(
        f"Unsupported input format {suffix!r}. Use a CSV or TSV file."
    )


def lineage_set(value: object) -> frozenset[str]:
    """Convert a semicolon-separated lineage string into a set."""
    if pd.isna(value):
        return frozenset()

    return frozenset(
        lineage.strip()
        for lineage in str(value).split(";")
        if lineage.strip()
    )


def lineages_match_as_sets(row: pd.Series) -> bool:
    """
    Compare lineage components without depending on their reported order.

    Example:
        TB-Profiler: 2.2;4.1
        fastlin:     4.1;2.2

    These are treated as the same lineage combination.
    """
    tb_set = lineage_set(row["tbprofiler_lineage"])
    fastlin_set = lineage_set(row["fastlin_lineage"])

    return bool(tb_set) and tb_set == fastlin_set


def classify_row(row: pd.Series) -> tuple[str, str]:
    """
    Return:
        classification
        classification_reason
    """

    if row["dataset_presence"] != "both":
        return (
            "excluded",
            f"not_present_in_both_datasets:{row['dataset_presence']}",
        )

    tb_count = int(row["tbprofiler_lineage_count"])
    fastlin_count = int(row["fastlin_lineage_count"])
    mixture = str(row["fastlin_mixture"]).strip().lower()

    if tb_count == 0 and fastlin_count == 0:
        return "excluded", "both_methods_no_call"

    if tb_count == 0:
        return "excluded", "tbprofiler_no_call"

    if fastlin_count == 0:
        return "excluded", "fastlin_no_call"

    # Conservative single-lineage candidate rule:
    # - present in both datasets
    # - exactly one lineage from each method
    # - exact lineage agreement
    # - fastlin does not flag a mixture
    if (
        tb_count == 1
        and fastlin_count == 1
        and row["lineage_comparison"] == "exact_agreement"
        and mixture == "no"
    ):
        return (
            "single_lineage_candidate",
            "single_lineage_concordant_non_mixed",
        )

    # Conservative mixed-lineage candidate rule:
    # - present in both datasets
    # - fastlin explicitly flags a mixture
    # - both methods report more than one lineage
    # - both methods report the same lineage components
    if (
        mixture == "yes"
        and tb_count > 1
        and fastlin_count > 1
        and bool(row["lineage_components_match"])
    ):
        return (
            "mixed_lineage_candidate",
            "mixed_lineage_components_concordant",
        )

    # Remaining rows are excluded from candidate panels but retained
    # with a clear reason.
    if mixture == "yes" and tb_count <= 1:
        return (
            "excluded",
            "fastlin_mixture_but_tbprofiler_single_or_no_call",
        )

    if mixture == "yes" and fastlin_count > 1 and tb_count > 1:
        return (
            "excluded",
            "mixed_lineage_components_disagree",
        )

    if mixture == "yes":
        return "excluded", "fastlin_mixture_not_concordant"

    if tb_count > 1 or fastlin_count > 1:
        return "excluded", "multiple_lineages_without_concordant_mixture"

    if row["lineage_comparison"] == "disagreement":
        return "excluded", "single_lineage_disagreement"

    if row["lineage_comparison"] == "tbprofiler_no_call":
        return "excluded", "tbprofiler_no_call"

    if row["lineage_comparison"] == "fastlin_no_call":
        return "excluded", "fastlin_no_call"

    if row["lineage_comparison"] == "both_no_call":
        return "excluded", "both_methods_no_call"

    return "excluded", "does_not_meet_candidate_rules"


def write_report(
    path: Path,
    input_path: Path,
    classified: pd.DataFrame,
) -> None:
    class_counts = classified["candidate_classification"].value_counts()
    reason_counts = classified["classification_reason"].value_counts()

    lines = [
        "FASTLIN / TB-PROFILER CANDIDATE CLASSIFICATION REPORT",
        "=" * 56,
        "",
        f"Input dataset: {input_path}",
        f"Total rows classified: {len(classified)}",
        "",
        "CANDIDATE CLASSIFICATION",
        (
            "single_lineage_candidate: "
            f"{int(class_counts.get('single_lineage_candidate', 0))}"
        ),
        (
            "mixed_lineage_candidate: "
            f"{int(class_counts.get('mixed_lineage_candidate', 0))}"
        ),
        f"excluded: {int(class_counts.get('excluded', 0))}",
        "",
        "CLASSIFICATION REASONS",
    ]

    for reason, count in reason_counts.items():
        lines.append(f"{reason}: {int(count)}")

    lines.extend(
        [
            "",
            "RULES",
            (
                "Single-lineage candidate: present in both datasets; one "
                "lineage reported by each method; exact agreement; "
                "fastlin mixture=no."
            ),
            (
                "Mixed-lineage candidate: present in both datasets; "
                "fastlin mixture=yes; both methods report multiple lineages; "
                "the same lineage components are reported by both methods."
            ),
            (
                "Excluded: any sample that does not meet either candidate rule."
            ),
            "",
            "NOTE",
            (
                "Candidate classification is not truth assignment. "
                "All candidate samples still require independent barcode "
                "evidence and QC confirmation."
            ),
        ]
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_arguments()

    input_path = args.merged_dataset.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        return 1

    try:
        df = read_table(input_path)
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        print(
            "ERROR: input dataset is missing required columns: "
            f"{sorted(missing)}",
            file=sys.stderr,
        )
        return 1

    # Ensure count columns are numeric and safely handle missing values.
    for column in [
        "tbprofiler_lineage_count",
        "fastlin_lineage_count",
    ]:
        df[column] = (
            pd.to_numeric(df[column], errors="coerce")
            .fillna(0)
            .astype(int)
        )

    df["lineage_components_match"] = df.apply(
        lineages_match_as_sets,
        axis=1,
    )

    classifications = df.apply(classify_row, axis=1)
    df["candidate_classification"] = [
        result[0] for result in classifications
    ]
    df["classification_reason"] = [
        result[1] for result in classifications
    ]

    output_dir.mkdir(parents=True, exist_ok=True)

    master_path = output_dir / "candidate_classification_master.csv"
    single_path = output_dir / "single_lineage_candidates.csv"
    mixed_path = output_dir / "mixed_lineage_candidates.csv"
    excluded_path = output_dir / "excluded_candidates.csv"
    report_path = output_dir / "candidate_classification_report.txt"

    df.to_csv(master_path, index=False)

    df.loc[
        df["candidate_classification"] == "single_lineage_candidate"
    ].to_csv(single_path, index=False)

    df.loc[
        df["candidate_classification"] == "mixed_lineage_candidate"
    ].to_csv(mixed_path, index=False)

    df.loc[
        df["candidate_classification"] == "excluded"
    ].to_csv(excluded_path, index=False)

    write_report(
        report_path,
        input_path,
        df,
    )

    print(f"Classified {len(df):,} rows.")
    print(f"Master:   {master_path}")
    print(f"Single:   {single_path}")
    print(f"Mixed:    {mixed_path}")
    print(f"Excluded: {excluded_path}")
    print(f"Report:   {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
