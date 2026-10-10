#!/usr/bin/env python3
"""
Apply the terminal-lineage confirmation policy to pilot marker observations.

Inputs
------
1. pilot_terminal_marker_observations.tsv
2. pilot_terminal_marker_observation_summary.csv

Outputs
-------
pilot_terminal_truth_candidates.tsv
pilot_terminal_truth_summary.csv
pilot_terminal_truth_report.txt
terminal_confirmation_policy.json

Decision policy
---------------
rejected:
    reference_allele_observed > 0 OR other_allele_observed > 0

review:
    mixed_alleles > 0 OR ambiguous_allele_observation > 0
    OR callable_marker_count != expected_allele_observed

confirmed_low_depth:
    all callable markers support the expected terminal allele,
    no conflicts are present, and at least one marker is low depth
    or has no coverage

confirmed:
    every expected terminal marker supports the expected allele and
    no non-callable or conflicting marker is present

This script interprets existing evidence only. It does not inspect BAMs,
FASTQs, or ancestral-lineage markers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_OBSERVATION_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "candidate_broad_lineage",
    "marker_id",
    "observation_status",
}

REQUIRED_SUMMARY_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "candidate_broad_lineage",
    "expected_marker_count",
    "expected_allele_observed",
    "reference_allele_observed",
    "mixed_alleles",
    "other_allele_observed",
    "ambiguous_allele_observation",
    "low_depth",
    "no_coverage",
    "callable_marker_count",
    "callable_marker_fraction",
    "expected_marker_support_fraction",
    "median_depth",
    "minimum_depth",
    "maximum_depth",
    "mean_expected_allele_fraction",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply the predefined terminal-lineage confirmation policy "
            "to pilot marker observations."
        )
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path(
            "lineage/pilot/barcode_observations/"
            "pilot_terminal_marker_observations.tsv"
        ),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(
            "lineage/pilot/barcode_observations/"
            "pilot_terminal_marker_observation_summary.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lineage/pilot/truth_candidates"),
    )
    parser.add_argument(
        "--policy-version",
        default="1.0",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_observations(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    missing = REQUIRED_OBSERVATION_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            f"Observation table missing columns: {sorted(missing)}"
        )

    if df.duplicated(["sample_id", "marker_id"], keep=False).any():
        raise ValueError(
            "Duplicate sample-marker rows found in observation table."
        )

    return df


def read_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    missing = REQUIRED_SUMMARY_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            f"Summary table missing columns: {sorted(missing)}"
        )

    if df["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample rows found in summary table.")

    integer_columns = [
        "expected_marker_count",
        "expected_allele_observed",
        "reference_allele_observed",
        "mixed_alleles",
        "other_allele_observed",
        "ambiguous_allele_observation",
        "low_depth",
        "no_coverage",
        "callable_marker_count",
    ]

    for column in integer_columns:
        df[column] = pd.to_numeric(
            df[column], errors="raise"
        ).astype(int)

    return df


def validate_consistency(
    observations: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    if set(observations["sample_id"]) != set(summary["sample_id"]):
        raise ValueError(
            "Observation and summary sample sets do not match."
        )

    counts = (
        observations.groupby("sample_id")["marker_id"]
        .count()
        .rename("observation_rows")
    )

    check = summary.set_index("sample_id").join(counts)

    bad = check[
        check["expected_marker_count"] != check["observation_rows"]
    ]

    if not bad.empty:
        raise ValueError(
            "Expected marker counts do not match observation rows for: "
            f"{bad.index.tolist()}"
        )


def determine_status(row: pd.Series) -> tuple[str, str]:
    reference_count = int(row["reference_allele_observed"])
    other_count = int(row["other_allele_observed"])
    mixed_count = int(row["mixed_alleles"])
    ambiguous_count = int(row["ambiguous_allele_observation"])
    low_depth_count = int(row["low_depth"])
    no_coverage_count = int(row["no_coverage"])
    supported = int(row["expected_allele_observed"])
    callable_count = int(row["callable_marker_count"])
    expected_count = int(row["expected_marker_count"])

    if reference_count > 0 or other_count > 0:
        reasons = []
        if reference_count > 0:
            reasons.append(f"{reference_count}_reference_dominant_marker")
        if other_count > 0:
            reasons.append(f"{other_count}_unexpected_allele_marker")
        return "rejected", ";".join(reasons)

    if mixed_count > 0 or ambiguous_count > 0:
        reasons = []
        if mixed_count > 0:
            reasons.append(f"{mixed_count}_mixed_marker")
        if ambiguous_count > 0:
            reasons.append(f"{ambiguous_count}_ambiguous_marker")
        return "review", ";".join(reasons)

    if supported != callable_count:
        return (
            "review",
            (
                "callable_markers_not_all_expected:"
                f"supported={supported};callable={callable_count}"
            ),
        )

    if low_depth_count > 0 or no_coverage_count > 0:
        reasons = []
        if low_depth_count > 0:
            reasons.append(f"{low_depth_count}_low_depth_marker")
        if no_coverage_count > 0:
            reasons.append(f"{no_coverage_count}_no_coverage_marker")
        return "confirmed_low_depth", ";".join(reasons)

    if supported == expected_count:
        return "confirmed", "all_terminal_markers_supported"

    return "review", "unclassified_terminal_marker_pattern"


def build_truth_candidates(
    observations: pd.DataFrame,
    summary: pd.DataFrame,
    policy_version: str,
) -> pd.DataFrame:
    output = summary.copy()

    decisions = output.apply(
        determine_status,
        axis=1,
        result_type="expand",
    )
    decisions.columns = [
        "truth_candidate_status",
        "decision_reason",
    ]

    output = pd.concat([output, decisions], axis=1)

    output["conflicting_marker_count"] = (
        output["reference_allele_observed"]
        + output["other_allele_observed"]
    )

    output["non_callable_marker_count"] = (
        output["low_depth"] + output["no_coverage"]
    )

    non_callable = (
        observations[
            observations["observation_status"].isin(
                ["low_depth", "no_coverage"]
            )
        ]
        .groupby("sample_id")["marker_id"]
        .apply(lambda values: ";".join(sorted(values)))
        .rename("non_callable_marker_ids")
    )

    output = output.merge(
        non_callable,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    output["non_callable_marker_ids"] = (
        output["non_callable_marker_ids"].fillna("")
    )
    output["terminal_confirmation_scope"] = (
        "exact_candidate_lineage_only"
    )
    output["terminal_confirmation_policy_version"] = policy_version

    columns = [
        "sample_id",
        "candidate_lineage",
        "candidate_broad_lineage",
        "expected_marker_count",
        "expected_allele_observed",
        "callable_marker_count",
        "callable_marker_fraction",
        "expected_marker_support_fraction",
        "low_depth",
        "no_coverage",
        "non_callable_marker_count",
        "non_callable_marker_ids",
        "reference_allele_observed",
        "other_allele_observed",
        "mixed_alleles",
        "ambiguous_allele_observation",
        "conflicting_marker_count",
        "median_depth",
        "minimum_depth",
        "maximum_depth",
        "mean_expected_allele_fraction",
        "truth_candidate_status",
        "decision_reason",
        "terminal_confirmation_scope",
        "terminal_confirmation_policy_version",
    ]

    return output[columns].sort_values(
        [
            "candidate_broad_lineage",
            "candidate_lineage",
            "sample_id",
        ],
        kind="stable",
    ).reset_index(drop=True)


def build_status_summary(
    truth_candidates: pd.DataFrame,
) -> pd.DataFrame:
    summary = (
        truth_candidates.groupby("truth_candidate_status")
        .agg(
            sample_count=("sample_id", "nunique"),
            broad_lineage_count=(
                "candidate_broad_lineage", "nunique"
            ),
            exact_lineage_count=("candidate_lineage", "nunique"),
        )
        .reset_index()
    )

    total = truth_candidates["sample_id"].nunique()
    summary["sample_percentage"] = (
        summary["sample_count"] / total * 100
    ).round(3)

    order = {
        "confirmed": 0,
        "confirmed_low_depth": 1,
        "review": 2,
        "rejected": 3,
    }

    summary["_order"] = (
        summary["truth_candidate_status"].map(order).fillna(99)
    )

    return summary.sort_values(
        ["_order", "truth_candidate_status"]
    ).drop(columns="_order").reset_index(drop=True)


def write_policy(
    path: Path,
    observations_path: Path,
    summary_path: Path,
    policy_version: str,
) -> None:
    policy = {
        "policy_name": "terminal_lineage_confirmation",
        "policy_version": policy_version,
        "confirmation_scope": "exact_candidate_lineage_only",
        "input_observations": str(observations_path),
        "input_observations_sha256": sha256_file(observations_path),
        "input_summary": str(summary_path),
        "input_summary_sha256": sha256_file(summary_path),
        "status_precedence": [
            "rejected",
            "review",
            "confirmed_low_depth",
            "confirmed",
        ],
        "rules": {
            "rejected": (
                "reference_allele_observed > 0 OR "
                "other_allele_observed > 0"
            ),
            "review": (
                "mixed_alleles > 0 OR "
                "ambiguous_allele_observation > 0 OR "
                "callable_marker_count != expected_allele_observed"
            ),
            "confirmed_low_depth": (
                "all callable markers support expected allele; "
                "no conflicts; low_depth + no_coverage > 0"
            ),
            "confirmed": (
                "all expected terminal markers support expected allele; "
                "no non-callable or conflicting markers"
            ),
        },
        "notes": [
            (
                "The status field records the decision; evidence counts and "
                "marker identifiers are stored in separate columns."
            ),
            (
                "Low-depth and uncovered markers are non-callable evidence, "
                "not contradictory evidence."
            ),
            "Ancestral-path markers are not evaluated.",
        ],
    }

    path.write_text(
        json.dumps(policy, indent=2) + "\n",
        encoding="utf-8",
    )


def write_report(
    path: Path,
    truth_candidates: pd.DataFrame,
    observations_path: Path,
    summary_path: Path,
    policy_path: Path,
    policy_version: str,
) -> None:
    counts = truth_candidates[
        "truth_candidate_status"
    ].value_counts()

    lines = [
        "PILOT TERMINAL-LINEAGE TRUTH-CANDIDATE REPORT",
        "=" * 48,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Observation input: {observations_path}",
        f"Observation SHA-256: {sha256_file(observations_path)}",
        f"Summary input: {summary_path}",
        f"Summary SHA-256: {sha256_file(summary_path)}",
        f"Policy file: {policy_path}",
        f"Policy version: {policy_version}",
        "",
        "STATUS COUNTS",
        f"confirmed: {int(counts.get('confirmed', 0))}",
        (
            "confirmed_low_depth: "
            f"{int(counts.get('confirmed_low_depth', 0))}"
        ),
        f"review: {int(counts.get('review', 0))}",
        f"rejected: {int(counts.get('rejected', 0))}",
        "",
        "PER-SAMPLE DECISIONS",
    ]

    for row in truth_candidates.itertuples(index=False):
        lines.append(
            f"{row.sample_id}\t"
            f"{row.candidate_lineage}\t"
            f"expected={row.expected_marker_count}\t"
            f"supported={row.expected_allele_observed}\t"
            f"low_depth={row.low_depth}\t"
            f"no_coverage={row.no_coverage}\t"
            f"conflicts={row.conflicting_marker_count}\t"
            f"status={row.truth_candidate_status}\t"
            f"reason={row.decision_reason}"
        )

    attention = truth_candidates[
        truth_candidates["truth_candidate_status"] != "confirmed"
    ]

    lines.extend(
        [
            "",
            "SAMPLES REQUIRING QUALIFIED INTERPRETATION",
        ]
    )

    if attention.empty:
        lines.append("None")
    else:
        for row in attention.itertuples(index=False):
            lines.append(
                f"{row.sample_id}\t"
                f"{row.candidate_lineage}\t"
                f"{row.truth_candidate_status}\t"
                f"{row.decision_reason}\t"
                f"non_callable_marker_ids="
                f"{row.non_callable_marker_ids or 'none'}"
            )

    lines.extend(
        [
            "",
            "INTERPRETATION",
            (
                "Confirmation used only markers explicitly assigned to each "
                "exact candidate lineage in the frozen barcode."
            ),
            (
                "Low-depth and uncovered positions were retained as "
                "non-callable evidence and did not count as contradictions."
            ),
            (
                "These remain truth candidates until the panel is reviewed "
                "and formally frozen."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    observations_path = args.observations.expanduser().resolve()
    summary_path = args.summary.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    for path, label in [
        (observations_path, "observation table"),
        (summary_path, "observation summary"),
    ]:
        if not path.is_file():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 1

    try:
        observations = read_observations(observations_path)
        summary = read_summary(summary_path)
        validate_consistency(observations, summary)

        truth_candidates = build_truth_candidates(
            observations,
            summary,
            args.policy_version,
        )
        status_summary = build_status_summary(truth_candidates)

        output_dir.mkdir(parents=True, exist_ok=True)

        truth_path = (
            output_dir / "pilot_terminal_truth_candidates.tsv"
        )
        status_path = (
            output_dir / "pilot_terminal_truth_summary.csv"
        )
        report_path = (
            output_dir / "pilot_terminal_truth_report.txt"
        )
        policy_path = (
            output_dir / "terminal_confirmation_policy.json"
        )

        truth_candidates.to_csv(
            truth_path,
            sep="\t",
            index=False,
        )
        status_summary.to_csv(status_path, index=False)

        write_policy(
            policy_path,
            observations_path,
            summary_path,
            args.policy_version,
        )
        write_report(
            report_path,
            truth_candidates,
            observations_path,
            summary_path,
            policy_path,
            args.policy_version,
        )

    except (
        ValueError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    counts = truth_candidates[
        "truth_candidate_status"
    ].value_counts()

    print(f"Truth candidates:       {len(truth_candidates)}")
    print(f"Confirmed:              {int(counts.get('confirmed', 0))}")
    print(
        "Confirmed low depth:    "
        f"{int(counts.get('confirmed_low_depth', 0))}"
    )
    print(f"Review:                 {int(counts.get('review', 0))}")
    print(f"Rejected:               {int(counts.get('rejected', 0))}")
    print()
    print(f"Truth candidates: {truth_path}")
    print(f"Status summary:   {status_path}")
    print(f"Policy:           {policy_path}")
    print(f"Report:           {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
