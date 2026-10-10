#!/usr/bin/env python3

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


# ============================================================
# Configuration
# ============================================================

REFERENCE_PANEL_FILE = Path(
    "./tbtamr/reference_dataset/"
    "tbtamr_current_study_reference_panel.csv"
)

RUN_COLUMN = "Sample_ID"

FASTQ_DIR = Path("./tbtamr/fastqs")

DOWNLOAD_MANIFEST = Path(
    "./tbtamr/reference_dataset/"
    "tbtamr_current_study_fastq_manifest.csv"
)

# Number of temporary threads used by fasterq-dump.
FASTERQ_THREADS = 4

# Number of threads used by pigz, if pigz is installed.
PIGZ_THREADS = 4

# Temporary directory used by fasterq-dump.
TEMP_DIR = Path("./tbtamr/tmp/fasterq")

# If True, overwrite incomplete or invalid outputs.
OVERWRITE_INCOMPLETE = True


# ============================================================
# Utility functions
# ============================================================

def require_command(command: str) -> str:
    """Return the executable path or stop with a clear error."""

    executable = shutil.which(command)

    if executable is None:
        raise RuntimeError(
            f"Required command '{command}' was not found in PATH."
        )

    return executable


def run_command(command: list[str]) -> tuple[bool, str]:
    """
    Run an external command.

    Output is shown in the terminal while the command runs.
    A short error message is returned if the command fails.
    """

    print("  $ " + " ".join(command))

    try:
        subprocess.run(
            command,
            check=True,
        )
        return True, ""

    except subprocess.CalledProcessError as error:
        message = (
            f"Command failed with exit code {error.returncode}: "
            + " ".join(command)
        )
        return False, message


def gzip_with_python(input_file: Path) -> Path:
    """Compress one FASTQ using Python gzip."""

    output_file = Path(f"{input_file}.gz")

    with input_file.open("rb") as source:
        with gzip.open(output_file, "wb") as destination:
            shutil.copyfileobj(source, destination)

    input_file.unlink()

    return output_file


def compress_fastq(input_file: Path, pigz_path: str | None) -> tuple[bool, str]:
    """
    Compress a FASTQ.

    pigz is used when available; otherwise Python gzip is used.
    """

    if not input_file.exists():
        return False, f"FASTQ not found before compression: {input_file}"

    if pigz_path:
        success, error = run_command(
            [
                pigz_path,
                "-f",
                "-p",
                str(PIGZ_THREADS),
                str(input_file),
            ]
        )

        if not success:
            return False, error

        output_file = Path(f"{input_file}.gz")

    else:
        print(f"  Compressing with Python gzip: {input_file.name}")

        try:
            output_file = gzip_with_python(input_file)

        except OSError as error:
            return False, f"Compression failed for {input_file}: {error}"

    if not output_file.exists() or output_file.stat().st_size == 0:
        return False, f"Compressed output is missing or empty: {output_file}"

    return True, ""


def remove_if_exists(path: Path) -> None:
    """Remove a file if it exists."""

    if path.exists():
        path.unlink()


def outputs_are_complete(r1: Path, r2: Path) -> bool:
    """Return True when both paired FASTQ files exist and are non-empty."""

    return (
        r1.exists()
        and r2.exists()
        and r1.stat().st_size > 0
        and r2.stat().st_size > 0
    )


def file_size_bytes(path: Path) -> int | None:
    """Return file size when the file exists."""

    if not path.exists():
        return None

    return path.stat().st_size


# ============================================================
# Main
# ============================================================

def main() -> int:
    """Download and verify FASTQs for the selected reference panel."""

    if not REFERENCE_PANEL_FILE.exists():
        print(
            f"ERROR: Reference panel not found: {REFERENCE_PANEL_FILE}",
            file=sys.stderr,
        )
        return 1

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    DOWNLOAD_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    try:
        prefetch_path = require_command("prefetch")
        fasterq_path = require_command("fasterq-dump")

    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    pigz_path = shutil.which("pigz")

    if pigz_path:
        print(f"Using pigz for compression: {pigz_path}")
    else:
        print("pigz not found; using Python gzip compression.")

    reference_panel = pd.read_csv(REFERENCE_PANEL_FILE)

    if RUN_COLUMN not in reference_panel.columns:
        print(
            f"ERROR: Column '{RUN_COLUMN}' was not found in "
            f"{REFERENCE_PANEL_FILE}.",
            file=sys.stderr,
        )
        return 1

    runs = (
        reference_panel[RUN_COLUMN]
        .dropna()
        .astype(str)
        .str.strip()
    )

    runs = runs[runs.ne("")].drop_duplicates().tolist()

    if not runs:
        print("ERROR: No run accessions were found.", file=sys.stderr)
        return 1

    print(f"Found {len(runs)} unique runs to process.")

    manifest_records: list[dict] = []

    for position, run_id in enumerate(runs, start=1):
        print()
        print(f"=== [{position}/{len(runs)}] Processing {run_id} ===")

        r1_gz = FASTQ_DIR / f"{run_id}_1.fastq.gz"
        r2_gz = FASTQ_DIR / f"{run_id}_2.fastq.gz"

        r1_fastq = FASTQ_DIR / f"{run_id}_1.fastq"
        r2_fastq = FASTQ_DIR / f"{run_id}_2.fastq"

        status = "failed"
        error_message = ""

        if outputs_are_complete(r1_gz, r2_gz):
            print("  Complete paired FASTQs already exist; skipping download.")
            status = "already_present"

        else:
            if OVERWRITE_INCOMPLETE:
                print("  Removing incomplete or partial FASTQ outputs.")
                for path in [r1_gz, r2_gz, r1_fastq, r2_fastq]:
                    remove_if_exists(path)

            print("  Running prefetch...")

            success, error_message = run_command(
                [
                    prefetch_path,
                    run_id,
                ]
            )

            if success:
                print("  Running fasterq-dump...")

                success, error_message = run_command(
                    [
                        fasterq_path,
                        run_id,
                        "--split-files",
                        "--threads",
                        str(FASTERQ_THREADS),
                        "--temp",
                        str(TEMP_DIR),
                        "--outdir",
                        str(FASTQ_DIR),
                    ]
                )

            if success:
                if not r1_fastq.exists() or not r2_fastq.exists():
                    success = False
                    error_message = (
                        "fasterq-dump did not produce both expected paired files: "
                        f"{r1_fastq.name}, {r2_fastq.name}"
                    )

            if success:
                print("  Compressing R1...")
                success, error_message = compress_fastq(
                    r1_fastq,
                    pigz_path,
                )

            if success:
                print("  Compressing R2...")
                success, error_message = compress_fastq(
                    r2_fastq,
                    pigz_path,
                )

            if success and outputs_are_complete(r1_gz, r2_gz):
                status = "downloaded"
                print(f"  Completed: {run_id}")

            elif success:
                status = "failed"
                error_message = (
                    "Final paired FASTQ verification failed."
                )

        manifest_records.append(
            {
                "Sample_ID": run_id,
                "R1": str(r1_gz.resolve()) if r1_gz.exists() else "",
                "R2": str(r2_gz.resolve()) if r2_gz.exists() else "",
                "R1_size_bytes": file_size_bytes(r1_gz),
                "R2_size_bytes": file_size_bytes(r2_gz),
                "Status": status,
                "Error": error_message,
            }
        )

        # Write after every sample so progress is preserved if interrupted.
        pd.DataFrame(manifest_records).to_csv(
            DOWNLOAD_MANIFEST,
            index=False,
        )

        if status == "failed":
            print(f"  ERROR: {error_message}")

    manifest = pd.DataFrame(manifest_records)

    downloaded_count = manifest["Status"].isin(
        ["downloaded", "already_present"]
    ).sum()

    failed_count = manifest["Status"].eq("failed").sum()

    print()
    print("FASTQ processing complete.")
    print(f"Successful or already present: {downloaded_count}")
    print(f"Failed:                        {failed_count}")
    print(f"Manifest: {DOWNLOAD_MANIFEST}")

    if failed_count:
        print()
        print("Failed samples:")
        for sample_id in manifest.loc[
            manifest["Status"] == "failed",
            "Sample_ID",
        ]:
            print(f"  - {sample_id}")

        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())