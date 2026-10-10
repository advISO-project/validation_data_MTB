#!/usr/bin/env python3

"""
Generate the final AMR validation-set manifests.

Inputs
------
1. tbtAMR current-study reference panel
2. FASTQ retrieval manifest

Outputs
-------
1. Canonical CSV validation manifest
2. Galaxy list:pair collection YAML

The script verifies that the selected reference panel and FASTQ manifest
contain exactly the same samples and that all selected samples have valid
paired FASTQ files.
"""

from pathlib import Path
import csv
import sys

try:
    import yaml
except ImportError:
    sys.exit(
        "ERROR: PyYAML is required. Install with: pip install pyyaml"
    )


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]

REFERENCE_PANEL = (
    REPO_ROOT
    / "tbtamr"
    / "reference_dataset"
    / "tbtamr_current_study_reference_panel.csv"
)

FASTQ_MANIFEST = (
    REPO_ROOT
    / "tbtamr"
    / "reference_dataset"
    / "tbtamr_current_study_fastq_manifest.csv"
)

OUTPUT_DIR = REPO_ROOT / "tbtamr" / "reference_dataset"

VALIDATION_MANIFEST = (
    OUTPUT_DIR / "tbtamr_current_study_validation_manifest.csv"
)

GALAXY_MANIFEST = (
    OUTPUT_DIR / "tbtamr_current_study_galaxy_collection.yaml"
)


# ---------------------------------------------------------------------
# Required fields
# ---------------------------------------------------------------------

EXPECTED_FIELDS = [
    "rifampicin_expected",
    "isoniazid_expected",
    "pyrazinamide_expected",
    "ethambutol_expected",
]

OUTPUT_FIELDS = [
    "Sample_ID",
    "Resistance_class",
    "Selection_reason",
    "rifampicin_expected",
    "isoniazid_expected",
    "pyrazinamide_expected",
    "ethambutol_expected",
    "R1",
    "R2",
]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def read_csv(path):
    if not path.exists():
        sys.exit(f"ERROR: Input file does not exist: {path}")

    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def index_by_sample(rows, source_name):
    indexed = {}

    for row in rows:
        sample = row.get("Sample_ID", "").strip()

        if not sample:
            sys.exit(
                f"ERROR: Missing Sample_ID in {source_name}"
            )

        if sample in indexed:
            sys.exit(
                f"ERROR: Duplicate Sample_ID {sample} in {source_name}"
            )

        indexed[sample] = row

    return indexed


def repository_relative_path(path_string):
    path = Path(path_string).expanduser()

    if not path.is_absolute():
        path = (REPO_ROOT / path).resolve()
    else:
        path = path.resolve()

    try:
        relative = path.relative_to(REPO_ROOT.resolve())
    except ValueError:
        sys.exit(
            f"ERROR: FASTQ is outside repository: {path}"
        )

    return relative.as_posix(), path


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    panel_rows = read_csv(REFERENCE_PANEL)
    fastq_rows = read_csv(FASTQ_MANIFEST)

    panel = index_by_sample(panel_rows, "reference panel")
    fastqs = index_by_sample(fastq_rows, "FASTQ manifest")

    panel_samples = set(panel)
    fastq_samples = set(fastqs)

    if panel_samples != fastq_samples:

        missing_fastqs = sorted(panel_samples - fastq_samples)
        unexpected_fastqs = sorted(fastq_samples - panel_samples)

        if missing_fastqs:
            print(
                "ERROR: Selected samples missing from FASTQ manifest:",
                ", ".join(missing_fastqs),
                file=sys.stderr,
            )

        if unexpected_fastqs:
            print(
                "ERROR: FASTQ manifest contains samples not in "
                "reference panel:",
                ", ".join(unexpected_fastqs),
                file=sys.stderr,
            )

        sys.exit(1)

    output_rows = []
    galaxy_elements = []

    # Preserve reference-panel ordering
    for panel_row in panel_rows:

        sample = panel_row["Sample_ID"].strip()
        fastq_row = fastqs[sample]

        # -------------------------------------------------------------
        # Validate expected phenotypes
        # -------------------------------------------------------------

        for field in EXPECTED_FIELDS:
            value = panel_row.get(field, "").strip()

            if value not in {"R", "S"}:
                sys.exit(
                    f"ERROR: {sample}: {field} must be R or S; "
                    f"found '{value}'"
                )

        # -------------------------------------------------------------
        # Validate FASTQ retrieval status
        # -------------------------------------------------------------

        status = fastq_row.get("Status", "").strip()

        if status not in {"downloaded", "already_present"}:
            sys.exit(
                f"ERROR: {sample}: invalid FASTQ status '{status}'"
            )

        r1_value = fastq_row.get("R1", "").strip()
        r2_value = fastq_row.get("R2", "").strip()

        if not r1_value or not r2_value:
            sys.exit(
                f"ERROR: {sample}: R1 or R2 path is missing"
            )

        r1_relative, r1_absolute = repository_relative_path(r1_value)
        r2_relative, r2_absolute = repository_relative_path(r2_value)

        for label, path in [
            ("R1", r1_absolute),
            ("R2", r2_absolute),
        ]:
            if not path.exists():
                sys.exit(
                    f"ERROR: {sample}: {label} does not exist: {path}"
                )

            if path.stat().st_size == 0:
                sys.exit(
                    f"ERROR: {sample}: {label} is empty: {path}"
                )

        # -------------------------------------------------------------
        # Canonical validation manifest
        # -------------------------------------------------------------

        output_rows.append(
            {
                "Sample_ID": sample,
                "Resistance_class": panel_row[
                    "Resistance_class"
                ],
                "Selection_reason": panel_row[
                    "Selection_reason"
                ],
                "rifampicin_expected": panel_row[
                    "rifampicin_expected"
                ],
                "isoniazid_expected": panel_row[
                    "isoniazid_expected"
                ],
                "pyrazinamide_expected": panel_row[
                    "pyrazinamide_expected"
                ],
                "ethambutol_expected": panel_row[
                    "ethambutol_expected"
                ],
                "R1": r1_relative,
                "R2": r2_relative,
            }
        )

        # -------------------------------------------------------------
        # Galaxy list:pair collection
        # -------------------------------------------------------------

        galaxy_elements.append(
            {
                "class": "Collection",
                "collection_type": "pair",
                "identifier": sample,
                "elements": [
                    {
                        "class": "File",
                        "identifier": "forward",
                        "path": r1_relative,
                    },
                    {
                        "class": "File",
                        "identifier": "reverse",
                        "path": r2_relative,
                    },
                ],
            }
        )

    # -----------------------------------------------------------------
    # Write CSV
    # -----------------------------------------------------------------

    with VALIDATION_MANIFEST.open("w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
        )

        writer.writeheader()
        writer.writerows(output_rows)

    # -----------------------------------------------------------------
    # Write Galaxy YAML
    # -----------------------------------------------------------------

    galaxy_manifest = {
        "name": "tbtamr_current_study_amr_validation_set",
        "description": (
            "Paired-end FASTQ collection for the tbtAMR "
            "current-study AMR validation reference panel"
        ),
        "version": "1.0.0",
        "samples": {
            "class": "Collection",
            "collection_type": "list:pair",
            "elements": galaxy_elements,
        },
    }

    with GALAXY_MANIFEST.open("w") as handle:
        yaml.safe_dump(
            galaxy_manifest,
            handle,
            sort_keys=False,
        )

    print("AMR validation manifest generation complete")
    print(f"Samples: {len(output_rows)}")
    print(f"CSV:    {VALIDATION_MANIFEST}")
    print(f"Galaxy: {GALAXY_MANIFEST}")


if __name__ == "__main__":
    main()