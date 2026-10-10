#!/usr/bin/env python3
"""
Summarise the distribution of single-lineage reference candidates.

Input
-----
single_lineage_candidates.csv

Outputs
-------
single_lineage_candidate_summary.csv
broad_lineage_candidate_summary.csv
single_lineage_candidate_summary_report.txt

This script does not:
- select samples;
- download sequencing data;
- examine barcode markers;
- confirm candidate lineages; or
- assign final truth labels.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "sample_id",
    "tbprofiler_lineage",
    "fastlin_lineage",
    "candidate_classification",
}


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Summarise exact-lineage and broad-lineage representation "
            "among single-lineage reference candidates."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Path to single_lineage_candidates.csv or TSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for output files. Default: current directory.",
    )

    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV or tab-separated input file as strings."""

    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path, dtype=str)

    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str)

    raise ValueError(
        f"Unsupported input format {suffix!r}. "
        "Use a CSV or tab-separated file."
    )


def clean_text(value: object) -> str:
    """Convert a value to stripped text, returning an empty string for NA."""

    if pd.isna(value):
        return ""

    return str(value).strip()


def normalise_lineage(value: object) -> str:
    """
    Normalise a single lineage label without changing its biological meaning.

    Examples
    --------
    lineage4.3.4    -> 4.3.4
    L4              -> 4
    4.3.4           -> 4.3.4
    La1.1           -> La1.1
    La1.7.X-unk4    -> La1.7.X-unk4
    """

    lineage = clean_text(value)

    if not lineage:
        return ""

    if lineage.lower().startswith("lineage"):
        lineage = lineage[len("lineage"):].strip()

    broad_aliases = {
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

    return broad_aliases.get(lineage.upper(), lineage)


def extract_broad_lineage(lineage: str) -> str:
    """
    Derive the broad lineage group from a normalised lineage label.

    Numeric lineages are grouped by their first numerical component.

    Examples
    --------
    4.3.4.2.1     -> 4
    4.2           -> 4
    1.2.2.2       -> 1
    8             -> 8

    La-prefixed lineages are grouped by the first La component.

    Examples
    --------
    La1.8.1       -> La1
    La1.7.X-unk4  -> La1
    La1.2.BCG     -> La1
    La2           -> La2
    La3           -> La3

    Unrecognised labels are retained with an UNRESOLVED prefix so that they
    are not silently assigned to an incorrect group.
    """

    lineage = clean_text(lineage)

    if not lineage:
        return "UNRESOLVED:EMPTY"

    numeric_match = re.match(
        r"^(\d+)(?:\.|$)",
        lineage,
    )

    if numeric_match:
        return numeric_match.group(1)

    la_match = re.match(
        r"^(La\d+)(?:\.|$)",
        lineage,
        flags=re.IGNORECASE,
    )

    if la_match:
        broad = la_match.group(1)

        # Preserve standard capitalisation: La1, La2, La3.
        return "La" + broad[2:]

    return f"UNRESOLVED:{lineage}"


def lineage_depth(lineage: str) -> int:
    """
    Return the number of lineage hierarchy components.

    Examples
    --------
    4             -> 1
    4.3.4         -> 3
    4.3.4.2.1     -> 5
    La1.1         -> 2
    La1.7.X-unk4  -> 3
    La2           -> 1
    """

    lineage = clean_text(lineage)

    if not lineage:
        return 0

    return len(lineage.split("."))


def validate_single_lineage_candidates(df: pd.DataFrame) -> None:
    """Validate that the input is the expected candidate dataset."""

    missing = REQUIRED_COLUMNS.difference(df.columns)

    if missing:
        raise ValueError(
            "Input is missing required columns: "
            f"{sorted(missing)}"
        )

    invalid_classification = df[
        df["candidate_classification"]
        .fillna("")
        .ne("single_lineage_candidate")
    ]

    if not invalid_classification.empty:
        raise ValueError(
            "Input contains rows that are not classified as "
            "'single_lineage_candidate'. "
            f"Unexpected rows: {len(invalid_classification)}"
        )


def build_exact_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Create one summary row per exact lineage label."""

    summary = (
        df.groupby(
            [
                "candidate_lineage",
                "broad_lineage",
                "lineage_depth",
            ],
            dropna=False,
        )
        .agg(
            candidate_count=("sample_id", "nunique"),
        )
        .reset_index()
    )

    total_candidates = df["sample_id"].nunique()

    summary["candidate_percentage"] = (
        summary["candidate_count"]
        / total_candidates
        * 100
    ).round(3)

    summary = summary.sort_values(
        by=[
            "broad_lineage",
            "lineage_depth",
            "candidate_lineage",
        ],
        kind="stable",
    ).reset_index(drop=True)

    return summary


def build_broad_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Create one summary row per broad lineage group."""

    summary = (
        df.groupby(
            "broad_lineage",
            dropna=False,
        )
        .agg(
            candidate_count=("sample_id", "nunique"),
            exact_lineage_count=("candidate_lineage", "nunique"),
            minimum_lineage_depth=("lineage_depth", "min"),
            maximum_lineage_depth=("lineage_depth", "max"),
        )
        .reset_index()
    )

    total_candidates = df["sample_id"].nunique()

    summary["candidate_percentage"] = (
        summary["candidate_count"]
        / total_candidates
        * 100
    ).round(3)

    summary = summary.sort_values(
        by=[
            "candidate_count",
            "broad_lineage",
        ],
        ascending=[
            False,
            True,
        ],
        kind="stable",
    ).reset_index(drop=True)

    return summary


def write_report(
    report_path: Path,
    input_path: Path,
    df: pd.DataFrame,
    exact_summary: pd.DataFrame,
    broad_summary: pd.DataFrame,
) -> None:
    """Write a human-readable candidate-distribution report."""

    total_candidates = df["sample_id"].nunique()
    exact_lineage_count = df["candidate_lineage"].nunique()
    broad_lineage_count = df["broad_lineage"].nunique()

    exact_counts = exact_summary["candidate_count"]

    unresolved = broad_summary[
        broad_summary["broad_lineage"].str.startswith(
            "UNRESOLVED:",
            na=False,
        )
    ]

    lines = [
        "SINGLE-LINEAGE CANDIDATE SUMMARY REPORT",
        "=" * 43,
        "",
        f"Input dataset: {input_path}",
        f"Total candidate rows: {len(df)}",
        f"Unique candidate samples: {total_candidates}",
        (
            "Exact lineage labels represented: "
            f"{exact_lineage_count}"
        ),
        (
            "Broad lineage groups represented: "
            f"{broad_lineage_count}"
        ),
        "",
        "EXACT-LINEAGE REPRESENTATION",
        (
            "Minimum candidates in an exact lineage: "
            f"{int(exact_counts.min())}"
        ),
        (
            "Maximum candidates in an exact lineage: "
            f"{int(exact_counts.max())}"
        ),
        (
            "Median candidates per exact lineage: "
            f"{exact_counts.median():.1f}"
        ),
        "",
        "BROAD-LINEAGE COUNTS",
    ]

    for row in broad_summary.itertuples(index=False):
        lines.append(
            f"{row.broad_lineage}: "
            f"{int(row.candidate_count)} candidates; "
            f"{int(row.exact_lineage_count)} exact lineage labels; "
            f"{float(row.candidate_percentage):.3f}% of candidates"
        )

    lines.extend(
        [
            "",
            "RAREST EXACT LINEAGES",
        ]
    )

    rarest = exact_summary.sort_values(
        by=[
            "candidate_count",
            "candidate_lineage",
        ],
        ascending=[
            True,
            True,
        ],
        kind="stable",
    ).head(20)

    for row in rarest.itertuples(index=False):
        lines.append(
            f"{row.candidate_lineage}: "
            f"{int(row.candidate_count)}"
        )

    lines.extend(
        [
            "",
            "MOST COMMON EXACT LINEAGES",
        ]
    )

    most_common = exact_summary.sort_values(
        by=[
            "candidate_count",
            "candidate_lineage",
        ],
        ascending=[
            False,
            True,
        ],
        kind="stable",
    ).head(20)

    for row in most_common.itertuples(index=False):
        lines.append(
            f"{row.candidate_lineage}: "
            f"{int(row.candidate_count)}"
        )

    lines.extend(
        [
            "",
            "GROUPING RULES",
            (
                "Numeric lineage labels are grouped using the first "
                "numerical component."
            ),
            "Example: 4.3.4.2.1 -> broad lineage 4.",
            (
                "La-prefixed lineage labels are grouped using the first "
                "La component."
            ),
            "Example: La1.7.X-unk4 -> broad lineage La1.",
            (
                "Labels that do not match either rule are retained as "
                "UNRESOLVED labels for manual review."
            ),
        ]
    )

    if unresolved.empty:
        lines.extend(
            [
                "",
                "UNRESOLVED GROUPING LABELS",
                "None",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "UNRESOLVED GROUPING LABELS",
            ]
        )

        for row in unresolved.itertuples(index=False):
            lines.append(
                f"{row.broad_lineage}: "
                f"{int(row.candidate_count)} candidates"
            )

    lines.extend(
        [
            "",
            "NOTE",
            (
                "This report describes candidate representation only. "
                "No samples were selected, downloaded, confirmed, or "
                "assigned final truth labels."
            ),
        ]
    )

    report_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    """Run the candidate-summary workflow."""

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
        df = read_table(input_path)

        validate_single_lineage_candidates(df)

        df = df.copy()

        df["sample_id"] = df["sample_id"].map(clean_text)

        df["candidate_lineage"] = (
            df["tbprofiler_lineage"]
            .map(normalise_lineage)
        )

        df["fastlin_lineage_normalised"] = (
            df["fastlin_lineage"]
            .map(normalise_lineage)
        )

        # The previous classification script required exact agreement.
        # Verify that this remains true after normalisation.
        lineage_mismatches = df[
            df["candidate_lineage"]
            != df["fastlin_lineage_normalised"]
        ]

        if not lineage_mismatches.empty:
            example_columns = [
                "sample_id",
                "tbprofiler_lineage",
                "fastlin_lineage",
                "candidate_lineage",
                "fastlin_lineage_normalised",
            ]

            examples = lineage_mismatches[
                example_columns
            ].head(10)

            raise ValueError(
                "Normalised fastlin and TB-Profiler lineage labels differ "
                f"for {len(lineage_mismatches)} supposedly concordant "
                "candidates.\n"
                f"Examples:\n{examples.to_string(index=False)}"
            )

        if (df["sample_id"] == "").any():
            raise ValueError(
                "At least one candidate has an empty sample ID."
            )

        if (df["candidate_lineage"] == "").any():
            raise ValueError(
                "At least one candidate has an empty lineage label."
            )

        if df["sample_id"].duplicated().any():
            duplicate_ids = sorted(
                df.loc[
                    df["sample_id"].duplicated(keep=False),
                    "sample_id",
                ].unique()
            )

            raise ValueError(
                "Duplicate candidate sample IDs found: "
                f"{duplicate_ids[:10]}"
            )

        df["broad_lineage"] = (
            df["candidate_lineage"]
            .map(extract_broad_lineage)
        )

        df["lineage_depth"] = (
            df["candidate_lineage"]
            .map(lineage_depth)
        )

        exact_summary = build_exact_summary(df)
        broad_summary = build_broad_summary(df)

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

    exact_summary_path = (
        output_dir
        / "single_lineage_candidate_summary.csv"
    )

    broad_summary_path = (
        output_dir
        / "broad_lineage_candidate_summary.csv"
    )

    report_path = (
        output_dir
        / "single_lineage_candidate_summary_report.txt"
    )

    exact_summary.to_csv(
        exact_summary_path,
        index=False,
    )

    broad_summary.to_csv(
        broad_summary_path,
        index=False,
    )

    write_report(
        report_path=report_path,
        input_path=input_path,
        df=df,
        exact_summary=exact_summary,
        broad_summary=broad_summary,
    )

    unresolved_count = int(
        broad_summary["broad_lineage"]
        .str.startswith(
            "UNRESOLVED:",
            na=False,
        )
        .sum()
    )

    print(
        f"Summarised {len(df):,} "
        "single-lineage candidates."
    )

    print(
        "Exact-lineage summary: "
        f"{exact_summary_path}"
    )

    print(
        "Broad-lineage summary: "
        f"{broad_summary_path}"
    )

    print(
        f"Report:                {report_path}"
    )

    print(
        "Broad lineage groups:  "
        f"{len(broad_summary)}"
    )

    print(
        "Unresolved groups:      "
        f"{unresolved_count}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())