#!/usr/bin/env python3

"""
Build the canonical lineage barcode reference used by the validation project.

Inputs
------
1. Frozen TBDB barcode BED file.
2. Frozen single-record reference FASTA.
3. Metadata files documenting both frozen resources.

Outputs
-------
1. lineage_barcode_reference.tsv
2. lineage_barcode_reference_report.txt

This script does not assign lineages to samples. It only creates a verified,
project-specific reference table from frozen external inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


EXPECTED_TBDB_COLUMNS = 8
VALID_BASES = set("ACGTN")


@dataclass(frozen=True)
class FastaRecord:
    accession: str
    header: str
    sequence: str


@dataclass(frozen=True)
class BarcodeMarker:
    row_number: int
    source_chrom: str
    start_0based: int
    end_0based: int
    lineage_raw: str
    marker_allele: str
    broad_lineage: str
    family: str
    rd: str


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a canonical lineage barcode reference from a frozen "
            "TBDB BED file and frozen reference FASTA."
        )
    )

    parser.add_argument(
        "--barcode",
        type=Path,
        default=Path("references/tbdb/raw/barcode.bed"),
        help="Frozen TBDB barcode BED file.",
    )

    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("references/reference/raw/NC_000962.3.fasta"),
        help="Frozen reference FASTA.",
    )

    parser.add_argument(
        "--tbdb-metadata",
        type=Path,
        default=Path("references/tbdb/metadata/source.yaml"),
        help="TBDB provenance metadata.",
    )

    parser.add_argument(
        "--reference-metadata",
        type=Path,
        default=Path("references/reference/metadata/source.yaml"),
        help="Reference provenance metadata.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "references/tbdb/processed/lineage_barcode_reference.tsv"
        ),
        help="Canonical processed barcode output.",
    )

    parser.add_argument(
        "--report",
        type=Path,
        default=Path(
            "references/tbdb/processed/"
            "lineage_barcode_reference_report.txt"
        ),
        help="Build report output.",
    )

    return parser.parse_args()


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of a file."""

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def parse_simple_yaml(path: Path) -> dict[str, str]:
    """
    Read the flat key-value YAML files used in this project.

    This intentionally supports only the simple metadata format currently
    used by the project. It avoids introducing a PyYAML dependency.
    """

    metadata: dict[str, str] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            if ":" not in line:
                raise ValueError(
                    f"{path}: line {line_number} is not a key-value pair: "
                    f"{raw_line.rstrip()}"
                )

            key, value = line.split(":", maxsplit=1)
            key = key.strip()
            value = value.strip()

            if (
                len(value) >= 2
                and value[0] == value[-1]
                and value[0] in {'"', "'"}
            ):
                value = value[1:-1]

            metadata[key] = value

    return metadata


def verify_metadata_checksum(
    resource_path: Path,
    metadata: dict[str, str],
    metadata_path: Path,
) -> str:
    """Verify that a file matches the checksum recorded in metadata."""

    expected = metadata.get("sha256")

    if not expected:
        raise ValueError(
            f"No sha256 value found in metadata file: {metadata_path}"
        )

    observed = sha256_file(resource_path)

    if observed.lower() != expected.lower():
        raise ValueError(
            "Checksum verification failed.\n"
            f"Resource: {resource_path}\n"
            f"Expected: {expected}\n"
            f"Observed: {observed}"
        )

    return observed


def read_single_fasta(path: Path) -> FastaRecord:
    """Read a FASTA containing exactly one sequence record."""

    header: str | None = None
    sequence_parts: list[str] = []
    sequence_count = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            # Blank lines are ignored, but the original file remains unchanged.
            if not line:
                continue

            if line.startswith(">"):
                sequence_count += 1

                if sequence_count > 1:
                    raise ValueError(
                        f"{path} contains more than one FASTA record."
                    )

                header = line[1:].strip()
                continue

            if header is None:
                raise ValueError(
                    f"{path}: sequence content occurs before the FASTA header "
                    f"at line {line_number}."
                )

            sequence_line = line.upper()
            unexpected = set(sequence_line) - VALID_BASES

            if unexpected:
                raise ValueError(
                    f"{path}: unexpected sequence characters on line "
                    f"{line_number}: {sorted(unexpected)}"
                )

            sequence_parts.append(sequence_line)

    if sequence_count != 1 or header is None:
        raise ValueError(
            f"{path} must contain exactly one FASTA record."
        )

    sequence = "".join(sequence_parts)

    if not sequence:
        raise ValueError(f"{path} contains an empty sequence.")

    accession = header.split()[0]

    return FastaRecord(
        accession=accession,
        header=header,
        sequence=sequence,
    )


def read_barcode(path: Path) -> list[BarcodeMarker]:
    """Read and validate the eight-column TBDB barcode BED file."""

    markers: list[BarcodeMarker] = []
    errors: list[str] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row_number, row in enumerate(reader, start=1):
            if not row or all(not value.strip() for value in row):
                continue

            if len(row) != EXPECTED_TBDB_COLUMNS:
                errors.append(
                    f"Row {row_number}: expected "
                    f"{EXPECTED_TBDB_COLUMNS} columns, found {len(row)}."
                )
                continue

            (
                source_chrom,
                start_text,
                end_text,
                lineage_raw,
                marker_allele,
                broad_lineage,
                family,
                rd,
            ) = (value.strip() for value in row)

            try:
                start_0based = int(start_text)
                end_0based = int(end_text)
            except ValueError:
                errors.append(
                    f"Row {row_number}: non-integer BED coordinate "
                    f"({start_text}, {end_text})."
                )
                continue

            marker_allele = marker_allele.upper()

            if start_0based < 0:
                errors.append(
                    f"Row {row_number}: start coordinate is negative."
                )

            if end_0based <= start_0based:
                errors.append(
                    f"Row {row_number}: end must be greater than start."
                )

            if end_0based - start_0based != 1:
                errors.append(
                    f"Row {row_number}: expected a one-base interval, found "
                    f"{start_0based}-{end_0based}."
                )

            if marker_allele not in VALID_BASES:
                errors.append(
                    f"Row {row_number}: invalid marker allele "
                    f"{marker_allele!r}."
                )

            if not lineage_raw:
                errors.append(
                    f"Row {row_number}: lineage label is empty."
                )

            markers.append(
                BarcodeMarker(
                    row_number=row_number,
                    source_chrom=source_chrom,
                    start_0based=start_0based,
                    end_0based=end_0based,
                    lineage_raw=lineage_raw,
                    marker_allele=marker_allele,
                    broad_lineage=broad_lineage,
                    family=family,
                    rd=rd,
                )
            )

    if errors:
        preview = "\n".join(errors[:20])

        if len(errors) > 20:
            preview += (
                f"\n...plus {len(errors) - 20} additional errors."
            )

        raise ValueError(
            f"Barcode validation failed with {len(errors)} error(s):\n"
            f"{preview}"
        )

    if not markers:
        raise ValueError(f"No barcode markers were read from {path}.")

    return markers


def normalise_lineage(lineage_raw: str) -> str:
    """
    Normalise the common TBDB 'lineage' prefix.

    Examples
    --------
    lineage4.2.2.1 -> 4.2.2.1
    lineage1       -> 1
    La1.8          -> La1.8
    """

    lineage = lineage_raw.strip()

    if lineage.lower().startswith("lineage"):
        lineage = lineage[len("lineage"):]

    return lineage


def clean_annotation(value: str) -> str:
    """Convert placeholder annotation values to an empty string."""

    cleaned = value.strip()

    if cleaned.lower() in {"none", "na", "n/a", "null", "."}:
        return ""

    return cleaned


def build_rows(
    markers: Iterable[BarcodeMarker],
    reference: FastaRecord,
    tbdb_commit: str,
    tbdb_sha256: str,
    reference_sha256: str,
) -> tuple[list[dict[str, object]], Counter]:
    """Create canonical output records and summary statistics."""

    rows: list[dict[str, object]] = []
    statistics: Counter = Counter()

    reference_length = len(reference.sequence)

    for marker_number, marker in enumerate(markers, start=1):
        if marker.end_0based > reference_length:
            raise ValueError(
                f"Barcode row {marker.row_number} is outside the reference "
                f"sequence: end={marker.end_0based}, "
                f"reference_length={reference_length}."
            )

        reference_allele = reference.sequence[
            marker.start_0based:marker.end_0based
        ]

        if len(reference_allele) != 1:
            raise ValueError(
                f"Barcode row {marker.row_number} did not resolve to exactly "
                f"one reference base."
            )

        marker_matches_reference = (
            marker.marker_allele == reference_allele
        )

        if marker_matches_reference:
            statistics["marker_matches_reference"] += 1
            marker_relation = "matches_reference"
        else:
            statistics["marker_differs_from_reference"] += 1
            marker_relation = "differs_from_reference"

        lineage = normalise_lineage(marker.lineage_raw)

        rows.append(
            {
                "marker_id": f"TBDB{marker_number:04d}",
                "source_row": marker.row_number,
                "source_chrom": marker.source_chrom,
                "start_0based": marker.start_0based,
                "end_0based": marker.end_0based,
                "position_1based": marker.start_0based + 1,
                "reference_accession": reference.accession,
                "reference_allele": reference_allele,
                "marker_allele": marker.marker_allele,
                "marker_relation_to_reference": marker_relation,
                "lineage_raw": marker.lineage_raw,
                "lineage": lineage,
                "broad_lineage": clean_annotation(
                    marker.broad_lineage
                ),
                "family": clean_annotation(marker.family),
                "rd": clean_annotation(marker.rd),
                "tbdb_commit": tbdb_commit,
                "tbdb_sha256": tbdb_sha256,
                "reference_sha256": reference_sha256,
            }
        )

    statistics["total_markers"] = len(rows)
    statistics["unique_positions"] = len(
        {row["position_1based"] for row in rows}
    )
    statistics["unique_lineage_labels"] = len(
        {row["lineage"] for row in rows}
    )

    position_counts = Counter(
        row["position_1based"] for row in rows
    )
    statistics["positions_with_multiple_markers"] = sum(
        1 for count in position_counts.values() if count > 1
    )

    exact_record_counts = Counter(
        (
            row["start_0based"],
            row["end_0based"],
            row["lineage"],
            row["marker_allele"],
        )
        for row in rows
    )
    statistics["duplicate_marker_records"] = sum(
        count - 1
        for count in exact_record_counts.values()
        if count > 1
    )

    return rows, statistics


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write canonical barcode records."""

    if not rows:
        raise ValueError("Cannot write an empty barcode reference.")

    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys())

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    *,
    barcode_path: Path,
    reference_path: Path,
    output_path: Path,
    reference: FastaRecord,
    tbdb_commit: str,
    tbdb_sha256: str,
    reference_sha256: str,
    statistics: Counter,
) -> None:
    """Write a human-readable audit report."""

    path.parent.mkdir(parents=True, exist_ok=True)

    generated_at = datetime.now(timezone.utc).isoformat()

    report_lines = [
        "Lineage barcode reference build report",
        "======================================",
        "",
        f"Generated at UTC: {generated_at}",
        "",
        "Inputs",
        "------",
        f"TBDB barcode: {barcode_path}",
        f"TBDB commit: {tbdb_commit}",
        f"TBDB SHA-256: {tbdb_sha256}",
        f"Reference FASTA: {reference_path}",
        f"Reference accession: {reference.accession}",
        f"Reference length: {len(reference.sequence)}",
        f"Reference SHA-256: {reference_sha256}",
        "",
        "Output",
        "------",
        f"Canonical barcode table: {output_path}",
        "",
        "Build summary",
        "-------------",
        f"Total markers: {statistics['total_markers']}",
        f"Unique positions: {statistics['unique_positions']}",
        (
            "Positions represented by multiple markers: "
            f"{statistics['positions_with_multiple_markers']}"
        ),
        (
            "Unique normalised lineage labels: "
            f"{statistics['unique_lineage_labels']}"
        ),
        (
            "Marker allele matches reference: "
            f"{statistics['marker_matches_reference']}"
        ),
        (
            "Marker allele differs from reference: "
            f"{statistics['marker_differs_from_reference']}"
        ),
        (
            "Duplicate marker records: "
            f"{statistics['duplicate_marker_records']}"
        ),
        "",
        "Validation status",
        "-----------------",
        "PASS",
        "",
        "Notes",
        "-----",
        (
            "BED intervals were interpreted as zero-based, half-open "
            "coordinates."
        ),
        (
            "The one-based position was calculated as start_0based + 1."
        ),
        (
            "The TBDB allele is retained as marker_allele and is not "
            "automatically assumed to be an alternate allele."
        ),
        (
            "No sample lineage classifications were made during this build."
        ),
    ]

    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(report_lines))
        handle.write("\n")


def require_file(path: Path, description: str) -> None:
    """Fail with a clear message if an input file is missing."""

    if not path.is_file():
        raise FileNotFoundError(
            f"{description} not found: {path}"
        )


def main() -> int:
    args = parse_arguments()

    try:
        require_file(args.barcode, "Frozen TBDB barcode")
        require_file(args.reference, "Frozen reference FASTA")
        require_file(args.tbdb_metadata, "TBDB metadata")
        require_file(
            args.reference_metadata,
            "Reference metadata",
        )

        tbdb_metadata = parse_simple_yaml(args.tbdb_metadata)
        reference_metadata = parse_simple_yaml(
            args.reference_metadata
        )

        tbdb_sha256 = verify_metadata_checksum(
            args.barcode,
            tbdb_metadata,
            args.tbdb_metadata,
        )
        reference_sha256 = verify_metadata_checksum(
            args.reference,
            reference_metadata,
            args.reference_metadata,
        )

        tbdb_commit = tbdb_metadata.get("git_commit", "")

        if not tbdb_commit:
            raise ValueError(
                f"No git_commit found in {args.tbdb_metadata}."
            )

        reference = read_single_fasta(args.reference)

        expected_accession = reference_metadata.get(
            "source_accession"
        )

        if (
            expected_accession
            and reference.accession != expected_accession
        ):
            raise ValueError(
                "Reference accession mismatch.\n"
                f"Metadata: {expected_accession}\n"
                f"FASTA: {reference.accession}"
            )

        expected_length = reference_metadata.get(
            "sequence_length"
        )

        if expected_length:
            try:
                expected_length_int = int(expected_length)
            except ValueError as exc:
                raise ValueError(
                    "sequence_length in reference metadata is not "
                    f"an integer: {expected_length}"
                ) from exc

            if len(reference.sequence) != expected_length_int:
                raise ValueError(
                    "Reference length mismatch.\n"
                    f"Metadata: {expected_length_int}\n"
                    f"FASTA: {len(reference.sequence)}"
                )

        markers = read_barcode(args.barcode)

        rows, statistics = build_rows(
            markers=markers,
            reference=reference,
            tbdb_commit=tbdb_commit,
            tbdb_sha256=tbdb_sha256,
            reference_sha256=reference_sha256,
        )

        write_tsv(args.output, rows)

        write_report(
            args.report,
            barcode_path=args.barcode,
            reference_path=args.reference,
            output_path=args.output,
            reference=reference,
            tbdb_commit=tbdb_commit,
            tbdb_sha256=tbdb_sha256,
            reference_sha256=reference_sha256,
            statistics=statistics,
        )

    except (
        FileNotFoundError,
        PermissionError,
        ValueError,
        OSError,
        csv.Error,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("Lineage barcode reference built successfully.")
    print(f"Canonical table: {args.output}")
    print(f"Build report:    {args.report}")
    print(f"Markers written: {statistics['total_markers']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
