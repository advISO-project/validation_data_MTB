#!/usr/bin/env python3
"""
Extract allele counts at expected terminal-lineage barcode positions.

Inputs
------
1. pilot_alignment_manifest.csv
2. pilot_expected_terminal_barcode_markers.tsv
3. Frozen reference FASTA

For each expected sample-marker pair, the script:
- locates the sample BAM;
- runs samtools mpileup at the expected marker positions;
- parses the pileup base string;
- counts A, C, G, T, N, deletion placeholders, and other symbols;
- calculates expected-marker and reference-allele fractions;
- assigns a descriptive observation status.

Outputs
-------
pilot_terminal_marker_observations.tsv
pilot_terminal_marker_observation_summary.csv
pilot_terminal_marker_observation_report.txt

This script does not confirm or reject lineage labels.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_ALIGNMENT_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "broad_lineage",
    "bam_path",
    "bai_path",
    "alignment_status",
}

REQUIRED_EXPECTATION_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "candidate_broad_lineage",
    "marker_id",
    "position_1based",
    "reference_accession",
    "reference_allele",
    "expected_marker_allele",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract per-base allele counts from pilot BAMs at expected "
            "terminal-lineage barcode positions."
        )
    )

    parser.add_argument(
        "--alignment-manifest",
        type=Path,
        default=Path(
            "lineage/pilot/alignment/reports/"
            "pilot_alignment_manifest.csv"
        ),
        help="Pilot alignment manifest.",
    )

    parser.add_argument(
        "--expectations",
        type=Path,
        default=Path(
            "lineage/pilot/barcode_expectations/"
            "pilot_expected_terminal_barcode_markers.tsv"
        ),
        help="Expected terminal-marker table.",
    )

    parser.add_argument(
        "--reference",
        type=Path,
        default=Path(
            "references/reference/raw/NC_000962.3.fasta"
        ),
        help="Frozen reference FASTA.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lineage/pilot/barcode_observations"),
        help="Output directory.",
    )

    parser.add_argument(
        "--minimum-base-quality",
        type=int,
        default=20,
        help="Minimum base quality used by samtools mpileup. Default: 20.",
    )

    parser.add_argument(
        "--minimum-mapping-quality",
        type=int,
        default=20,
        help="Minimum mapping quality used by samtools mpileup. Default: 20.",
    )

    parser.add_argument(
        "--low-depth-threshold",
        type=int,
        default=10,
        help=(
            "Depth below which a callable position is labelled low_depth. "
            "Default: 10."
        ),
    )

    parser.add_argument(
        "--dominant-allele-threshold",
        type=float,
        default=0.90,
        help=(
            "Allele fraction required for a dominant-allele observation. "
            "Default: 0.90."
        ),
    )

    parser.add_argument(
        "--mixed-minor-threshold",
        type=float,
        default=0.10,
        help=(
            "Minor allele fraction used to label mixed_alleles. "
            "Default: 0.10."
        ),
    )

    return parser.parse_args()


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


def require_program(program: str) -> str:
    result = subprocess.run(
        ["which", program],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    path = clean_text(result.stdout)

    if result.returncode != 0 or not path:
        raise RuntimeError(
            f"Required program was not found in PATH: {program}"
        )

    return path


def read_alignment_manifest(path: Path) -> pd.DataFrame:
    dataframe = pd.read_csv(path, dtype=str).fillna("")

    missing = REQUIRED_ALIGNMENT_COLUMNS.difference(
        dataframe.columns
    )

    if missing:
        raise ValueError(
            "Alignment manifest is missing required columns: "
            f"{sorted(missing)}"
        )

    valid_statuses = {
        "complete",
        "already_present_valid",
    }

    dataframe = dataframe[
        dataframe["alignment_status"].isin(valid_statuses)
    ].copy()

    if dataframe.empty:
        raise ValueError(
            "No valid aligned samples were found in the alignment manifest."
        )

    if dataframe["sample_id"].duplicated().any():
        raise ValueError(
            "Duplicate sample IDs were found in the alignment manifest."
        )

    return dataframe


def read_expectations(path: Path) -> pd.DataFrame:
    dataframe = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    ).fillna("")

    missing = REQUIRED_EXPECTATION_COLUMNS.difference(
        dataframe.columns
    )

    if missing:
        raise ValueError(
            "Expectation table is missing required columns: "
            f"{sorted(missing)}"
        )

    dataframe["position_1based"] = pd.to_numeric(
        dataframe["position_1based"],
        errors="raise",
    ).astype(int)

    duplicate_pairs = dataframe.duplicated(
        subset=["sample_id", "marker_id"],
        keep=False,
    )

    if duplicate_pairs.any():
        raise ValueError(
            "Duplicate sample-marker expectations were found."
        )

    return dataframe


def parse_pileup_bases(
    base_string: str,
    reference_allele: str,
) -> Counter:
    """
    Parse the samtools mpileup read-bases field.

    Pileup conventions handled:
    - '.' and ',' represent the reference allele;
    - '^' plus the following character marks read start and mapping quality;
    - '$' marks read end;
    - '+' or '-' followed by an integer and sequence represents an indel;
    - '*' and '#' are deletion placeholders;
    - '<' and '>' are reference skips;
    - A/C/G/T/N in either case are explicit base observations.
    """

    counts: Counter = Counter(
        {
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0,
            "N": 0,
            "deletion_placeholder": 0,
            "reference_skip": 0,
            "other": 0,
        }
    )

    ref = reference_allele.upper()
    index = 0

    while index < len(base_string):
        char = base_string[index]

        if char in ".,":
            if ref in {"A", "C", "G", "T", "N"}:
                counts[ref] += 1
            else:
                counts["other"] += 1
            index += 1
            continue

        if char == "^":
            index += 2
            continue

        if char == "$":
            index += 1
            continue

        if char in "+-":
            index += 1
            match = re.match(r"\d+", base_string[index:])

            if not match:
                counts["other"] += 1
                continue

            length_text = match.group(0)
            indel_length = int(length_text)
            index += len(length_text)
            index += indel_length
            continue

        if char in "*#":
            counts["deletion_placeholder"] += 1
            index += 1
            continue

        if char in "<>":
            counts["reference_skip"] += 1
            index += 1
            continue

        upper = char.upper()

        if upper in {"A", "C", "G", "T", "N"}:
            counts[upper] += 1
        else:
            counts["other"] += 1

        index += 1

    return counts


def run_mpileup(
    *,
    samtools: str,
    bam_path: Path,
    reference_path: Path,
    contig: str,
    positions: list[int],
    minimum_base_quality: int,
    minimum_mapping_quality: int,
) -> dict[int, dict[str, object]]:
    """
    Run one mpileup command for all expected positions in one sample.
    """

    observations: dict[int, dict[str, object]] = {}

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".bed",
        prefix="barcode_positions_",
        delete=False,
        encoding="utf-8",
    ) as bed_handle:
        bed_path = Path(bed_handle.name)

        for position in sorted(set(positions)):
            start_0based = position - 1
            bed_handle.write(
                f"{contig}\t{start_0based}\t{position}\n"
            )

    command = [
        samtools,
        "mpileup",
        "-f",
        str(reference_path),
        "-l",
        str(bed_path),
        "-Q",
        str(minimum_base_quality),
        "-q",
        str(minimum_mapping_quality),
        "-B",
        "-a",
        str(bam_path),
    ]

    try:
        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    finally:
        bed_path.unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(
            "samtools mpileup failed for "
            f"{bam_path.name}: {clean_text(result.stderr)}"
        )

    for line in result.stdout.splitlines():
        if not line.strip():
            continue

        parts = line.split("\t")

        if len(parts) < 6:
            raise ValueError(
                f"Unexpected mpileup output line: {line}"
            )

        line_contig = parts[0]
        position = int(parts[1])
        reference_allele = parts[2].upper()
        reported_depth = int(parts[3])
        base_string = parts[4]
        quality_string = parts[5]

        counts = parse_pileup_bases(
            base_string=base_string,
            reference_allele=reference_allele,
        )

        counted_depth = sum(
            counts[base]
            for base in ["A", "C", "G", "T", "N"]
        )

        observations[position] = {
            "mpileup_contig": line_contig,
            "mpileup_reference_allele": reference_allele,
            "reported_depth": reported_depth,
            "counted_base_depth": counted_depth,
            "quality_character_count": len(quality_string),
            **counts,
        }

    return observations


def classify_observation(
    *,
    total_depth: int,
    expected_fraction: float,
    reference_fraction: float,
    counts: Counter,
    expected_allele: str,
    reference_allele: str,
    low_depth_threshold: int,
    dominant_allele_threshold: float,
    mixed_minor_threshold: float,
) -> str:
    if total_depth == 0:
        return "no_coverage"

    if total_depth < low_depth_threshold:
        return "low_depth"

    nucleotide_counts = {
        base: counts.get(base, 0)
        for base in ["A", "C", "G", "T"]
    }

    sorted_counts = sorted(
        nucleotide_counts.items(),
        key=lambda item: item[1],
        reverse=True,
    )

    dominant_base, dominant_count = sorted_counts[0]
    second_count = sorted_counts[1][1]
    dominant_fraction = dominant_count / total_depth
    second_fraction = second_count / total_depth

    if (
        expected_fraction >= dominant_allele_threshold
        and dominant_base == expected_allele
    ):
        return "expected_allele_observed"

    if (
        reference_fraction >= dominant_allele_threshold
        and dominant_base == reference_allele
    ):
        return "reference_allele_observed"

    if (
        dominant_fraction < dominant_allele_threshold
        and second_fraction >= mixed_minor_threshold
    ):
        return "mixed_alleles"

    if (
        dominant_fraction >= dominant_allele_threshold
        and dominant_base not in {
            expected_allele,
            reference_allele,
        }
    ):
        return "other_allele_observed"

    return "ambiguous_allele_observation"


def build_observations(
    *,
    alignments: pd.DataFrame,
    expectations: pd.DataFrame,
    samtools: str,
    reference_path: Path,
    minimum_base_quality: int,
    minimum_mapping_quality: int,
    low_depth_threshold: int,
    dominant_allele_threshold: float,
    mixed_minor_threshold: float,
) -> pd.DataFrame:
    alignment_lookup = alignments.set_index(
        "sample_id",
        drop=False,
    )

    missing_samples = sorted(
        set(expectations["sample_id"])
        - set(alignment_lookup.index)
    )

    if missing_samples:
        raise ValueError(
            "Expected samples are missing from the valid alignment manifest: "
            f"{missing_samples}"
        )

    output_rows: list[dict[str, object]] = []

    sample_groups = expectations.groupby(
        "sample_id",
        sort=True,
    )

    for sample_number, (sample_id, sample_expectations) in enumerate(
        sample_groups,
        start=1,
    ):
        alignment_row = alignment_lookup.loc[sample_id]
        bam_path = Path(clean_text(alignment_row["bam_path"]))
        bai_path = Path(clean_text(alignment_row["bai_path"]))

        if not bam_path.is_file():
            raise FileNotFoundError(
                f"BAM not found for {sample_id}: {bam_path}"
            )

        if not bai_path.is_file():
            raise FileNotFoundError(
                f"BAM index not found for {sample_id}: {bai_path}"
            )

        contigs = sample_expectations[
            "reference_accession"
        ].unique()

        if len(contigs) != 1:
            raise ValueError(
                f"Sample {sample_id} has multiple expected reference "
                f"accessions: {contigs.tolist()}"
            )

        contig = clean_text(contigs[0])
        positions = sample_expectations[
            "position_1based"
        ].astype(int).tolist()

        print(
            f"[{sample_number}/{len(sample_groups)}] "
            f"Extracting {len(positions)} positions for {sample_id}"
        )

        pileup = run_mpileup(
            samtools=samtools,
            bam_path=bam_path,
            reference_path=reference_path,
            contig=contig,
            positions=positions,
            minimum_base_quality=minimum_base_quality,
            minimum_mapping_quality=minimum_mapping_quality,
        )

        for expectation in sample_expectations.itertuples(
            index=False
        ):
            position = int(expectation.position_1based)
            expected_allele = clean_text(
                expectation.expected_marker_allele
            ).upper()
            reference_allele = clean_text(
                expectation.reference_allele
            ).upper()

            pileup_row = pileup.get(
                position,
                {
                    "mpileup_contig": contig,
                    "mpileup_reference_allele": reference_allele,
                    "reported_depth": 0,
                    "counted_base_depth": 0,
                    "quality_character_count": 0,
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "N": 0,
                    "deletion_placeholder": 0,
                    "reference_skip": 0,
                    "other": 0,
                },
            )

            counts = Counter(
                {
                    base: int(pileup_row.get(base, 0))
                    for base in ["A", "C", "G", "T", "N"]
                }
            )

            total_depth = sum(
                counts[base]
                for base in ["A", "C", "G", "T", "N"]
            )

            expected_count = counts.get(expected_allele, 0)
            reference_count = counts.get(reference_allele, 0)

            expected_fraction = (
                expected_count / total_depth
                if total_depth > 0
                else 0.0
            )

            reference_fraction = (
                reference_count / total_depth
                if total_depth > 0
                else 0.0
            )

            observation_status = classify_observation(
                total_depth=total_depth,
                expected_fraction=expected_fraction,
                reference_fraction=reference_fraction,
                counts=counts,
                expected_allele=expected_allele,
                reference_allele=reference_allele,
                low_depth_threshold=low_depth_threshold,
                dominant_allele_threshold=dominant_allele_threshold,
                mixed_minor_threshold=mixed_minor_threshold,
            )

            output_rows.append(
                {
                    "sample_id": sample_id,
                    "candidate_lineage": (
                        expectation.candidate_lineage
                    ),
                    "candidate_broad_lineage": (
                        expectation.candidate_broad_lineage
                    ),
                    "marker_id": expectation.marker_id,
                    "position_1based": position,
                    "reference_accession": (
                        expectation.reference_accession
                    ),
                    "reference_allele": reference_allele,
                    "expected_marker_allele": expected_allele,
                    "bam_path": str(bam_path),
                    "minimum_base_quality": minimum_base_quality,
                    "minimum_mapping_quality": (
                        minimum_mapping_quality
                    ),
                    "reported_depth": int(
                        pileup_row["reported_depth"]
                    ),
                    "counted_base_depth": total_depth,
                    "A_count": counts["A"],
                    "C_count": counts["C"],
                    "G_count": counts["G"],
                    "T_count": counts["T"],
                    "N_count": counts["N"],
                    "deletion_placeholder_count": int(
                        pileup_row["deletion_placeholder"]
                    ),
                    "reference_skip_count": int(
                        pileup_row["reference_skip"]
                    ),
                    "other_symbol_count": int(
                        pileup_row["other"]
                    ),
                    "expected_allele_count": expected_count,
                    "expected_allele_fraction": round(
                        expected_fraction,
                        6,
                    ),
                    "reference_allele_count": reference_count,
                    "reference_allele_fraction": round(
                        reference_fraction,
                        6,
                    ),
                    "observation_status": observation_status,
                }
            )

    observations = pd.DataFrame(output_rows)

    observations = observations.sort_values(
        [
            "sample_id",
            "position_1based",
            "marker_id",
        ],
        kind="stable",
    ).reset_index(drop=True)

    return observations


def build_summary(
    observations: pd.DataFrame,
) -> pd.DataFrame:
    status_counts = (
        observations.groupby(
            [
                "sample_id",
                "candidate_lineage",
                "candidate_broad_lineage",
            ]
        )["observation_status"]
        .value_counts()
        .unstack(fill_value=0)
        .reset_index()
    )

    for status in [
        "expected_allele_observed",
        "reference_allele_observed",
        "mixed_alleles",
        "other_allele_observed",
        "ambiguous_allele_observation",
        "low_depth",
        "no_coverage",
    ]:
        if status not in status_counts.columns:
            status_counts[status] = 0

    totals = (
        observations.groupby("sample_id")
        .agg(
            expected_marker_count=("marker_id", "count"),
            callable_marker_count=(
                "observation_status",
                lambda values: sum(
                    value not in {
                        "low_depth",
                        "no_coverage",
                    }
                    for value in values
                ),
            ),
            median_depth=("counted_base_depth", "median"),
            minimum_depth=("counted_base_depth", "min"),
            maximum_depth=("counted_base_depth", "max"),
            mean_expected_allele_fraction=(
                "expected_allele_fraction",
                "mean",
            ),
        )
        .reset_index()
    )

    summary = status_counts.merge(
        totals,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    summary["callable_marker_fraction"] = (
        summary["callable_marker_count"]
        / summary["expected_marker_count"]
    ).round(6)

    summary["expected_marker_support_fraction"] = (
        summary["expected_allele_observed"]
        / summary["expected_marker_count"]
    ).round(6)

    summary = summary.sort_values(
        [
            "candidate_broad_lineage",
            "candidate_lineage",
            "sample_id",
        ],
        kind="stable",
    ).reset_index(drop=True)

    return summary


def write_report(
    path: Path,
    *,
    alignment_manifest_path: Path,
    expectation_path: Path,
    reference_path: Path,
    observations: pd.DataFrame,
    summary: pd.DataFrame,
    minimum_base_quality: int,
    minimum_mapping_quality: int,
    low_depth_threshold: int,
    dominant_allele_threshold: float,
    mixed_minor_threshold: float,
) -> None:
    overall_counts = observations[
        "observation_status"
    ].value_counts()

    lines = [
        "PILOT TERMINAL-MARKER OBSERVATION REPORT",
        "=" * 42,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Alignment manifest: {alignment_manifest_path}",
        (
            "Alignment manifest SHA-256: "
            f"{sha256_file(alignment_manifest_path)}"
        ),
        f"Expectation table: {expectation_path}",
        (
            "Expectation table SHA-256: "
            f"{sha256_file(expectation_path)}"
        ),
        f"Reference: {reference_path}",
        f"Reference SHA-256: {sha256_file(reference_path)}",
        "",
        "PILEUP FILTERS",
        f"Minimum mapping quality: {minimum_mapping_quality}",
        f"Minimum base quality: {minimum_base_quality}",
        f"Low-depth threshold: {low_depth_threshold}",
        f"Dominant-allele threshold: {dominant_allele_threshold}",
        f"Mixed minor-allele threshold: {mixed_minor_threshold}",
        "",
        "OVERALL OBSERVATION COUNTS",
    ]

    for status in [
        "expected_allele_observed",
        "reference_allele_observed",
        "mixed_alleles",
        "other_allele_observed",
        "ambiguous_allele_observation",
        "low_depth",
        "no_coverage",
    ]:
        lines.append(
            f"{status}: {int(overall_counts.get(status, 0))}"
        )

    lines.extend(
        [
            "",
            "PER-SAMPLE SUMMARY",
        ]
    )

    for row in summary.itertuples(index=False):
        lines.append(
            f"{row.sample_id}\t"
            f"{row.candidate_lineage}\t"
            f"markers={row.expected_marker_count}\t"
            f"expected={row.expected_allele_observed}\t"
            f"reference={row.reference_allele_observed}\t"
            f"mixed={row.mixed_alleles}\t"
            f"other={row.other_allele_observed}\t"
            f"low_depth={row.low_depth}\t"
            f"no_coverage={row.no_coverage}\t"
            f"median_depth={row.median_depth}"
        )

    lines.extend(
        [
            "",
            "NOTE",
            (
                "Observation statuses are descriptive only. "
                "No confirmation, review, or rejection rule was applied."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    alignment_manifest_path = (
        args.alignment_manifest.expanduser().resolve()
    )
    expectation_path = args.expectations.expanduser().resolve()
    reference_path = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if args.minimum_base_quality < 0:
        print(
            "ERROR: --minimum-base-quality cannot be negative.",
            file=sys.stderr,
        )
        return 1

    if args.minimum_mapping_quality < 0:
        print(
            "ERROR: --minimum-mapping-quality cannot be negative.",
            file=sys.stderr,
        )
        return 1

    if args.low_depth_threshold < 1:
        print(
            "ERROR: --low-depth-threshold must be at least 1.",
            file=sys.stderr,
        )
        return 1

    if not 0 < args.dominant_allele_threshold <= 1:
        print(
            "ERROR: --dominant-allele-threshold must be in (0, 1].",
            file=sys.stderr,
        )
        return 1

    if not 0 <= args.mixed_minor_threshold < 1:
        print(
            "ERROR: --mixed-minor-threshold must be in [0, 1).",
            file=sys.stderr,
        )
        return 1

    for path, description in [
        (alignment_manifest_path, "alignment manifest"),
        (expectation_path, "expectation table"),
        (reference_path, "reference FASTA"),
    ]:
        if not path.is_file():
            print(
                f"ERROR: {description} not found: {path}",
                file=sys.stderr,
            )
            return 1

    try:
        samtools = require_program("samtools")

        alignments = read_alignment_manifest(
            alignment_manifest_path
        )
        expectations = read_expectations(expectation_path)

        observations = build_observations(
            alignments=alignments,
            expectations=expectations,
            samtools=samtools,
            reference_path=reference_path,
            minimum_base_quality=args.minimum_base_quality,
            minimum_mapping_quality=args.minimum_mapping_quality,
            low_depth_threshold=args.low_depth_threshold,
            dominant_allele_threshold=(
                args.dominant_allele_threshold
            ),
            mixed_minor_threshold=args.mixed_minor_threshold,
        )

        summary = build_summary(observations)

        output_dir.mkdir(parents=True, exist_ok=True)

        observations_path = (
            output_dir
            / "pilot_terminal_marker_observations.tsv"
        )

        summary_path = (
            output_dir
            / "pilot_terminal_marker_observation_summary.csv"
        )

        report_path = (
            output_dir
            / "pilot_terminal_marker_observation_report.txt"
        )

        observations.to_csv(
            observations_path,
            sep="\t",
            index=False,
        )

        summary.to_csv(
            summary_path,
            index=False,
        )

        write_report(
            report_path,
            alignment_manifest_path=alignment_manifest_path,
            expectation_path=expectation_path,
            reference_path=reference_path,
            observations=observations,
            summary=summary,
            minimum_base_quality=args.minimum_base_quality,
            minimum_mapping_quality=(
                args.minimum_mapping_quality
            ),
            low_depth_threshold=args.low_depth_threshold,
            dominant_allele_threshold=(
                args.dominant_allele_threshold
            ),
            mixed_minor_threshold=args.mixed_minor_threshold,
        )

    except (
        FileNotFoundError,
        RuntimeError,
        ValueError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print()
    print(f"Observation rows: {len(observations)}")
    print(f"Samples:          {summary['sample_id'].nunique()}")
    print(f"Observations:     {observations_path}")
    print(f"Summary:          {summary_path}")
    print(f"Report:           {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
