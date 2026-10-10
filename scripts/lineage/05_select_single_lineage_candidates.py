#!/usr/bin/env python3
"""
Select a balanced single-lineage candidate pool and a reproducible pilot panel.

Stage 1: Balanced candidate pool
--------------------------------
Selection is applied independently within each exact lineage:

    1 candidate      -> keep 1
    2-5 candidates   -> keep all
    6-20 candidates  -> randomly select 5
    >20 candidates   -> randomly select 10

Stage 2: Pilot panel
--------------------
From the balanced pool:

    1. Select one sample from every broad lineage group.
    2. Fill remaining slots from exact lineages not yet represented.
    3. Only allow duplicate exact lineages if the pilot target cannot
       otherwise be reached.
    4. Use a fixed random seed for reproducibility.

This script does not:
- download FASTQs;
- inspect barcode markers;
- confirm candidate labels; or
- assign final truth labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "sample_id",
    "tbprofiler_lineage",
    "fastlin_lineage",
    "candidate_classification",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a balanced exact-lineage candidate pool and a "
            "broad-lineage-stratified pilot panel."
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
        help="Directory for outputs. Default: current directory.",
    )

    parser.add_argument(
        "--pilot-size",
        type=int,
        default=20,
        help="Target number of pilot samples. Default: 20.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed used for all selections. Default: 42.",
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


def extract_broad_lineage(lineage: str) -> str:
    lineage = clean_text(lineage)

    if not lineage:
        return "UNRESOLVED:EMPTY"

    if lineage[0].isdigit():
        return lineage.split(".")[0]

    lower = lineage.lower()

    if lower.startswith("la"):
        first_component = lineage.split(".")[0]

        if (
            len(first_component) >= 3
            and first_component[2:].isdigit()
        ):
            return "La" + first_component[2:]

    return f"UNRESOLVED:{lineage}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def validate_input(df: pd.DataFrame) -> None:
    missing = REQUIRED_COLUMNS.difference(df.columns)

    if missing:
        raise ValueError(
            f"Input is missing required columns: {sorted(missing)}"
        )

    invalid = df[
        df["candidate_classification"]
        .fillna("")
        .ne("single_lineage_candidate")
    ]

    if not invalid.empty:
        raise ValueError(
            "Input contains rows not classified as "
            "'single_lineage_candidate'. "
            f"Unexpected rows: {len(invalid)}"
        )


def stage1_target(available: int) -> int:
    if available <= 0:
        return 0

    if available == 1:
        return 1

    if 2 <= available <= 5:
        return available

    if 6 <= available <= 20:
        return 5

    return 10


def deterministic_sample(
    sample_ids: list[str],
    n: int,
    seed: int,
    grouping_key: str,
) -> list[str]:
    if n >= len(sample_ids):
        return sorted(sample_ids)

    stable_material = f"{seed}|{grouping_key}".encode("utf-8")
    stable_seed = int(
        hashlib.sha256(stable_material).hexdigest()[:16],
        16,
    )

    rng = random.Random(stable_seed)

    return sorted(rng.sample(sorted(sample_ids), n))


def build_balanced_pool(
    df: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_ids: list[str] = []
    summary_rows: list[dict[str, object]] = []

    grouped = df.groupby(
        "candidate_lineage",
        sort=True,
        dropna=False,
    )

    for lineage, group in grouped:
        sample_ids = sorted(group["sample_id"].tolist())
        available = len(sample_ids)
        target = stage1_target(available)

        chosen = deterministic_sample(
            sample_ids=sample_ids,
            n=target,
            seed=seed,
            grouping_key=f"stage1|{lineage}",
        )

        selected_ids.extend(chosen)

        if available == 1:
            policy = "keep_1"
        elif available <= 5:
            policy = "keep_all_2_to_5"
        elif available <= 20:
            policy = "select_5_from_6_to_20"
        else:
            policy = "select_10_from_over_20"

        summary_rows.append(
            {
                "candidate_lineage": lineage,
                "broad_lineage": group["broad_lineage"].iloc[0],
                "available_candidates": available,
                "selection_target": target,
                "selected_candidates": len(chosen),
                "selection_policy": policy,
                "selected_sample_ids": ";".join(chosen),
            }
        )

    balanced = df[df["sample_id"].isin(selected_ids)].copy()

    policy_lookup = {
        row["candidate_lineage"]: row["selection_policy"]
        for row in summary_rows
    }

    balanced["balanced_pool_selected"] = True
    balanced["balanced_pool_seed"] = seed
    balanced["stage1_selection_policy"] = (
        balanced["candidate_lineage"].map(policy_lookup)
    )

    balanced = balanced.sort_values(
        ["broad_lineage", "candidate_lineage", "sample_id"],
        kind="stable",
    ).reset_index(drop=True)

    summary = pd.DataFrame(summary_rows).sort_values(
        ["broad_lineage", "candidate_lineage"],
        kind="stable",
    ).reset_index(drop=True)

    return balanced, summary


def choose_one_per_broad_lineage(
    balanced: pd.DataFrame,
    seed: int,
) -> tuple[list[str], set[str]]:
    selected_ids: list[str] = []
    broad_coverage_ids: set[str] = set()

    for broad_lineage in sorted(
        balanced["broad_lineage"].dropna().unique()
    ):
        group = balanced[
            balanced["broad_lineage"] == broad_lineage
        ]

        chosen = deterministic_sample(
            sample_ids=group["sample_id"].tolist(),
            n=1,
            seed=seed,
            grouping_key=f"pilot_broad|{broad_lineage}",
        )[0]

        selected_ids.append(chosen)
        broad_coverage_ids.add(chosen)

    return selected_ids, broad_coverage_ids


def choose_unique_exact_lineage_fill(
    balanced: pd.DataFrame,
    selected_ids: list[str],
    remaining_slots: int,
    seed: int,
) -> list[str]:
    """
    Fill pilot slots using exact lineages not already represented.

    One candidate is selected per previously unrepresented exact lineage.
    """

    if remaining_slots <= 0:
        return []

    selected = balanced[
        balanced["sample_id"].isin(selected_ids)
    ]

    represented_lineages = set(
        selected["candidate_lineage"]
    )

    eligible = balanced[
        ~balanced["candidate_lineage"].isin(represented_lineages)
        & ~balanced["sample_id"].isin(selected_ids)
    ]

    exact_lineages = sorted(
        eligible["candidate_lineage"].unique()
    )

    if not exact_lineages:
        return []

    lineages_to_select = deterministic_sample(
        sample_ids=exact_lineages,
        n=min(remaining_slots, len(exact_lineages)),
        seed=seed,
        grouping_key="pilot_unique_exact_lineages",
    )

    additional_ids: list[str] = []

    for lineage in lineages_to_select:
        lineage_group = eligible[
            eligible["candidate_lineage"] == lineage
        ]

        chosen = deterministic_sample(
            sample_ids=lineage_group["sample_id"].tolist(),
            n=1,
            seed=seed,
            grouping_key=f"pilot_exact|{lineage}",
        )[0]

        additional_ids.append(chosen)

    return additional_ids


def choose_duplicate_fill(
    balanced: pd.DataFrame,
    selected_ids: list[str],
    remaining_slots: int,
    seed: int,
) -> list[str]:
    """
    Fill any final slots when unique exact-lineage coverage is exhausted.
    """

    if remaining_slots <= 0:
        return []

    remaining_ids = sorted(
        set(balanced["sample_id"]) - set(selected_ids)
    )

    if remaining_slots > len(remaining_ids):
        raise ValueError(
            f"Pilot requires {remaining_slots} additional samples, "
            f"but only {len(remaining_ids)} remain."
        )

    return deterministic_sample(
        sample_ids=remaining_ids,
        n=remaining_slots,
        seed=seed,
        grouping_key="pilot_duplicate_exact_fill",
    )


def build_pilot(
    balanced: pd.DataFrame,
    pilot_size: int,
    seed: int,
) -> pd.DataFrame:
    broad_groups = sorted(
        balanced["broad_lineage"].dropna().unique()
    )

    unresolved = [
        group
        for group in broad_groups
        if str(group).startswith("UNRESOLVED:")
    ]

    if unresolved:
        raise ValueError(
            "Pilot selection cannot proceed with unresolved broad groups: "
            f"{unresolved}"
        )

    if pilot_size < len(broad_groups):
        raise ValueError(
            f"Pilot size {pilot_size} is smaller than the number of broad "
            f"lineage groups ({len(broad_groups)})."
        )

    selected_ids, broad_coverage_ids = choose_one_per_broad_lineage(
        balanced=balanced,
        seed=seed,
    )

    unique_exact_fill_ids = choose_unique_exact_lineage_fill(
        balanced=balanced,
        selected_ids=selected_ids,
        remaining_slots=pilot_size - len(selected_ids),
        seed=seed,
    )

    selected_ids.extend(unique_exact_fill_ids)

    duplicate_fill_ids = choose_duplicate_fill(
        balanced=balanced,
        selected_ids=selected_ids,
        remaining_slots=pilot_size - len(selected_ids),
        seed=seed,
    )

    selected_ids.extend(duplicate_fill_ids)

    unique_exact_fill_set = set(unique_exact_fill_ids)
    duplicate_fill_set = set(duplicate_fill_ids)

    pilot = balanced[
        balanced["sample_id"].isin(selected_ids)
    ].copy()

    def reason(sample_id: str) -> str:
        if sample_id in broad_coverage_ids:
            return "broad_lineage_coverage"
        if sample_id in unique_exact_fill_set:
            return "unique_exact_lineage_fill"
        if sample_id in duplicate_fill_set:
            return "duplicate_exact_lineage_fill"
        raise ValueError(
            f"Selected sample has no pilot selection reason: {sample_id}"
        )

    pilot["pilot_selected"] = True
    pilot["pilot_seed"] = seed
    pilot["pilot_target_size"] = pilot_size
    pilot["pilot_selection_reason"] = (
        pilot["sample_id"].map(reason)
    )

    reason_order = {
        "broad_lineage_coverage": 0,
        "unique_exact_lineage_fill": 1,
        "duplicate_exact_lineage_fill": 2,
    }

    pilot["_reason_order"] = (
        pilot["pilot_selection_reason"].map(reason_order)
    )

    pilot = pilot.sort_values(
        [
            "_reason_order",
            "broad_lineage",
            "candidate_lineage",
            "sample_id",
        ],
        kind="stable",
    ).drop(columns="_reason_order").reset_index(drop=True)

    return pilot


def write_report(
    report_path: Path,
    *,
    input_path: Path,
    input_sha256: str,
    input_df: pd.DataFrame,
    balanced: pd.DataFrame,
    stage1_summary: pd.DataFrame,
    pilot: pd.DataFrame,
    pilot_size: int,
    seed: int,
) -> None:
    broad_counts = (
        pilot.groupby("broad_lineage")["sample_id"]
        .nunique()
        .sort_index()
    )

    pilot_reason_counts = (
        pilot["pilot_selection_reason"]
        .value_counts()
    )

    stage1_policy_counts = Counter(
        stage1_summary["selection_policy"]
    )

    duplicate_exact_lineage_count = int(
        pilot["candidate_lineage"].duplicated(keep=False).sum()
    )

    lines = [
        "SINGLE-LINEAGE CANDIDATE SELECTION REPORT",
        "=" * 45,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Input dataset: {input_path}",
        f"Input SHA-256: {input_sha256}",
        f"Random seed: {seed}",
        "",
        "INPUT",
        f"Eligible single-lineage candidates: {len(input_df)}",
        (
            "Exact lineage labels represented: "
            f"{input_df['candidate_lineage'].nunique()}"
        ),
        (
            "Broad lineage groups represented: "
            f"{input_df['broad_lineage'].nunique()}"
        ),
        "",
        "STAGE 1: BALANCED CANDIDATE POOL",
        (
            "Policy: 1 available -> keep 1; "
            "2-5 -> keep all; 6-20 -> select 5; "
            ">20 -> select 10."
        ),
        f"Balanced pool size: {len(balanced)}",
        (
            "Exact lineage labels retained: "
            f"{balanced['candidate_lineage'].nunique()}"
        ),
        (
            "Broad lineage groups retained: "
            f"{balanced['broad_lineage'].nunique()}"
        ),
        "",
        "Stage 1 policy applications:",
    ]

    for policy in [
        "keep_1",
        "keep_all_2_to_5",
        "select_5_from_6_to_20",
        "select_10_from_over_20",
    ]:
        lines.append(
            f"{policy}: {int(stage1_policy_counts.get(policy, 0))} "
            "exact lineage labels"
        )

    lines.extend(
        [
            "",
            "STAGE 2: PILOT PANEL",
            f"Pilot target size: {pilot_size}",
            f"Pilot selected size: {len(pilot)}",
            (
                "Broad lineage groups represented: "
                f"{pilot['broad_lineage'].nunique()}"
            ),
            (
                "Exact lineage labels represented: "
                f"{pilot['candidate_lineage'].nunique()}"
            ),
            (
                "Broad-lineage coverage selections: "
                f"{int(pilot_reason_counts.get('broad_lineage_coverage', 0))}"
            ),
            (
                "Unique exact-lineage fill selections: "
                f"{int(pilot_reason_counts.get('unique_exact_lineage_fill', 0))}"
            ),
            (
                "Duplicate exact-lineage fill selections: "
                f"{int(pilot_reason_counts.get('duplicate_exact_lineage_fill', 0))}"
            ),
            (
                "Pilot rows belonging to duplicated exact lineages: "
                f"{duplicate_exact_lineage_count}"
            ),
            "",
            "PILOT BROAD-LINEAGE COUNTS",
        ]
    )

    for broad_lineage, count in broad_counts.items():
        lines.append(
            f"{broad_lineage}: {int(count)}"
        )

    lines.extend(
        [
            "",
            "PILOT SAMPLES",
        ]
    )

    for row in pilot.itertuples(index=False):
        lines.append(
            f"{row.sample_id}\t"
            f"{row.candidate_lineage}\t"
            f"{row.broad_lineage}\t"
            f"{row.stage1_selection_policy}\t"
            f"{row.pilot_selection_reason}"
        )

    lines.extend(
        [
            "",
            "NOTE",
            (
                "The balanced pool and pilot are candidate selections only. "
                "No selected sample has yet been confirmed using the canonical "
                "barcode reference, and no final truth labels were assigned."
            ),
        ]
    )

    report_path.write_text(
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

    if args.pilot_size <= 0:
        print(
            "ERROR: --pilot-size must be greater than zero.",
            file=sys.stderr,
        )
        return 1

    try:
        df = read_table(input_path)
        validate_input(df)

        df = df.copy()

        df["sample_id"] = df["sample_id"].map(clean_text)
        df["candidate_lineage"] = (
            df["tbprofiler_lineage"].map(normalise_lineage)
        )
        df["fastlin_lineage_normalised"] = (
            df["fastlin_lineage"].map(normalise_lineage)
        )

        mismatch = df[
            df["candidate_lineage"]
            != df["fastlin_lineage_normalised"]
        ]

        if not mismatch.empty:
            raise ValueError(
                "Concordant input contains normalized lineage mismatches: "
                f"{len(mismatch)} rows."
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
            duplicates = sorted(
                df.loc[
                    df["sample_id"].duplicated(keep=False),
                    "sample_id",
                ].unique()
            )

            raise ValueError(
                "Duplicate sample IDs found: "
                f"{duplicates[:10]}"
            )

        df["broad_lineage"] = (
            df["candidate_lineage"].map(extract_broad_lineage)
        )

        balanced, stage1_summary = build_balanced_pool(
            df=df,
            seed=args.seed,
        )

        pilot = build_pilot(
            balanced=balanced,
            pilot_size=args.pilot_size,
            seed=args.seed,
        )

        input_sha256 = sha256_file(input_path)

    except (
        ValueError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    balanced_path = (
        output_dir / "balanced_single_lineage_candidates.csv"
    )
    stage1_summary_path = (
        output_dir / "balanced_single_lineage_selection_summary.csv"
    )
    pilot_path = (
        output_dir / "pilot_single_lineage_candidates.csv"
    )
    report_path = (
        output_dir / "single_lineage_selection_report.txt"
    )
    policy_path = (
        output_dir / "single_lineage_selection_policy.json"
    )

    balanced.to_csv(balanced_path, index=False)
    stage1_summary.to_csv(stage1_summary_path, index=False)
    pilot.to_csv(pilot_path, index=False)

    policy = {
        "input": str(input_path),
        "input_sha256": input_sha256,
        "random_seed": args.seed,
        "stage1_policy": {
            "1": 1,
            "2-5": "all",
            "6-20": 5,
            ">20": 10,
        },
        "pilot_size": args.pilot_size,
        "pilot_policy": [
            "select one candidate per broad lineage",
            "fill remaining slots using previously unrepresented exact lineages",
            (
                "allow duplicate exact lineages only when unique exact "
                "lineages are exhausted"
            ),
        ],
    }

    policy_path.write_text(
        json.dumps(policy, indent=2) + "\n",
        encoding="utf-8",
    )

    write_report(
        report_path=report_path,
        input_path=input_path,
        input_sha256=input_sha256,
        input_df=df,
        balanced=balanced,
        stage1_summary=stage1_summary,
        pilot=pilot,
        pilot_size=args.pilot_size,
        seed=args.seed,
    )

    print(f"Input candidates:       {len(df):,}")
    print(f"Balanced pool:          {len(balanced):,}")
    print(f"Pilot panel:            {len(pilot):,}")
    print(f"Pilot broad groups:     {pilot['broad_lineage'].nunique()}")
    print(f"Pilot exact lineages:   {pilot['candidate_lineage'].nunique()}")
    print()
    print(f"Balanced pool:          {balanced_path}")
    print(f"Stage 1 summary:        {stage1_summary_path}")
    print(f"Pilot panel:            {pilot_path}")
    print(f"Policy record:          {policy_path}")
    print(f"Selection report:       {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
