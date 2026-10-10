#!/usr/bin/env python3
"""
Align verified lineage-pilot FASTQs to the frozen project reference.

Input
-----
pilot_single_lineage_fastq_download_manifest.csv

Requirements
------------
- bwa
- samtools

Outputs
-------
lineage/pilot/alignment/
├── reference/
├── bam/
├── logs/
├── qc/
└── reports/

For each complete sample, the script:
- checks that both local FASTQ files exist;
- aligns paired reads with bwa mem;
- sorts the alignment with samtools sort;
- indexes the BAM;
- validates the BAM with samtools quickcheck;
- writes samtools flagstat output;
- writes an updated alignment manifest and summary report.

This script does not inspect barcode loci or confirm lineage labels.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "sample_id",
    "candidate_lineage",
    "broad_lineage",
    "local_fastq_1",
    "local_fastq_2",
    "sample_download_status",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align verified lineage-pilot FASTQs and produce sorted, "
            "indexed BAM files."
        )
    )

    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "lineage/pilot/reports/"
            "pilot_single_lineage_fastq_download_manifest.csv"
        ),
        help="Verified FASTQ download manifest.",
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
        default=Path("lineage/pilot/alignment"),
        help="Alignment output directory.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads used by bwa mem and samtools sort. Default: 4.",
    )

    parser.add_argument(
        "--sort-memory",
        default="1G",
        help="Memory per samtools sort thread. Default: 1G.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild BAMs even when valid outputs already exist.",
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
    path = shutil.which(program)

    if not path:
        raise RuntimeError(
            f"Required program was not found in PATH: {program}"
        )

    return path


def run_command(
    command: list[str],
    *,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> subprocess.CompletedProcess:
    stdout_handle = None
    stderr_handle = None

    try:
        if stdout_path is not None:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_handle = stdout_path.open("w", encoding="utf-8")

        if stderr_path is not None:
            stderr_path.parent.mkdir(parents=True, exist_ok=True)
            stderr_handle = stderr_path.open("w", encoding="utf-8")

        return subprocess.run(
            command,
            check=False,
            stdout=stdout_handle or subprocess.PIPE,
            stderr=stderr_handle or subprocess.PIPE,
            text=True,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()


def tool_version(program: str, args: list[str]) -> str:
    """
    Return the first line of a tool's version output.

    The command output is captured as bytes and decoded with replacement
    characters so unexpected build metadata cannot stop report generation.
    """

    result = subprocess.run(
        [program, *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    text = result.stdout.decode(
        "utf-8",
        errors="replace",
    ).strip()

    if not text:
        return "unknown"

    return text.splitlines()[0]


def ensure_bwa_index(
    bwa: str,
    reference: Path,
    index_prefix: Path,
    log_path: Path,
) -> None:
    expected = [
        Path(f"{index_prefix}.amb"),
        Path(f"{index_prefix}.ann"),
        Path(f"{index_prefix}.bwt"),
        Path(f"{index_prefix}.pac"),
        Path(f"{index_prefix}.sa"),
    ]

    if all(path.is_file() for path in expected):
        print("BWA index already present.")
        return

    index_prefix.parent.mkdir(parents=True, exist_ok=True)

    command = [
        bwa,
        "index",
        "-p",
        str(index_prefix),
        str(reference),
    ]

    print("Building BWA index.")
    print("  $ " + " ".join(command))

    result = run_command(
        command,
        stderr_path=log_path,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "bwa index failed. See log: "
            f"{log_path}"
        )

    missing = [
        str(path)
        for path in expected
        if not path.is_file()
    ]

    if missing:
        raise RuntimeError(
            "bwa index completed but expected files are missing: "
            f"{missing}"
        )


def bam_is_valid(
    samtools: str,
    bam_path: Path,
    bai_path: Path,
) -> bool:
    if not bam_path.is_file() or not bai_path.is_file():
        return False

    result = subprocess.run(
        [samtools, "quickcheck", "-v", str(bam_path)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    return result.returncode == 0


def align_sample(
    *,
    sample_id: str,
    read1: Path,
    read2: Path,
    index_prefix: Path,
    bam_path: Path,
    bai_path: Path,
    align_log: Path,
    sort_log: Path,
    flagstat_path: Path,
    bwa: str,
    samtools: str,
    threads: int,
    sort_memory: str,
) -> tuple[str, str]:
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    flagstat_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_bam = bam_path.with_suffix(".bam.tmp")

    if temporary_bam.exists():
        temporary_bam.unlink()

    bwa_command = [
        bwa,
        "mem",
        "-t",
        str(threads),
        "-R",
        (
            f"@RG\\tID:{sample_id}\\tSM:{sample_id}"
            f"\\tPL:ILLUMINA"
        ),
        str(index_prefix),
        str(read1),
        str(read2),
    ]

    sort_command = [
        samtools,
        "sort",
        "-@",
        str(threads),
        "-m",
        sort_memory,
        "-o",
        str(temporary_bam),
        "-",
    ]

    print("  Aligning and sorting.")

    with align_log.open("w", encoding="utf-8") as align_handle, \
         sort_log.open("w", encoding="utf-8") as sort_handle:

        bwa_process = subprocess.Popen(
            bwa_command,
            stdout=subprocess.PIPE,
            stderr=align_handle,
        )

        assert bwa_process.stdout is not None

        sort_process = subprocess.Popen(
            sort_command,
            stdin=bwa_process.stdout,
            stdout=subprocess.DEVNULL,
            stderr=sort_handle,
        )

        bwa_process.stdout.close()

        sort_returncode = sort_process.wait()
        bwa_returncode = bwa_process.wait()

    if bwa_returncode != 0:
        temporary_bam.unlink(missing_ok=True)
        return (
            "failed",
            f"bwa_mem_failed_exit_{bwa_returncode}",
        )

    if sort_returncode != 0:
        temporary_bam.unlink(missing_ok=True)
        return (
            "failed",
            f"samtools_sort_failed_exit_{sort_returncode}",
        )

    temporary_bam.replace(bam_path)

    index_result = subprocess.run(
        [
            samtools,
            "index",
            "-@",
            str(threads),
            str(bam_path),
            str(bai_path),
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if index_result.returncode != 0:
        return (
            "failed",
            (
                "samtools_index_failed: "
                f"{clean_text(index_result.stderr)}"
            ),
        )

    quickcheck_result = subprocess.run(
        [samtools, "quickcheck", "-v", str(bam_path)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if quickcheck_result.returncode != 0:
        return (
            "failed",
            (
                "samtools_quickcheck_failed: "
                f"{clean_text(quickcheck_result.stderr)}"
            ),
        )

    flagstat_result = run_command(
        [samtools, "flagstat", "-@", str(threads), str(bam_path)],
        stdout_path=flagstat_path,
    )

    if flagstat_result.returncode != 0:
        return (
            "failed",
            "samtools_flagstat_failed",
        )

    return "complete", "alignment_sort_index_and_quickcheck_passed"


def write_report(
    path: Path,
    *,
    manifest_path: Path,
    manifest_sha256: str,
    reference_path: Path,
    reference_sha256: str,
    bwa_version: str,
    samtools_version: str,
    threads: int,
    results: pd.DataFrame,
) -> None:
    counts = results["alignment_status"].value_counts()

    lines = [
        "LINEAGE PILOT ALIGNMENT REPORT",
        "=" * 31,
        "",
        f"Generated at UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Input manifest: {manifest_path}",
        f"Input manifest SHA-256: {manifest_sha256}",
        f"Reference: {reference_path}",
        f"Reference SHA-256: {reference_sha256}",
        f"BWA version: {bwa_version}",
        f"samtools version: {samtools_version}",
        f"Threads: {threads}",
        "",
        "ALIGNMENT STATUS",
        f"complete: {int(counts.get('complete', 0))}",
        f"already_present_valid: "
        f"{int(counts.get('already_present_valid', 0))}",
        f"failed: {int(counts.get('failed', 0))}",
        "",
        "FAILED SAMPLES",
    ]

    failed = results[results["alignment_status"] == "failed"]

    if failed.empty:
        lines.append("None")
    else:
        for row in failed.itertuples(index=False):
            lines.append(
                f"{row.sample_id}\t{row.alignment_message}"
            )

    lines.extend(
        [
            "",
            "OUTPUTS",
            "One coordinate-sorted BAM, BAM index, flagstat file, "
            "and alignment logs were generated per successful sample.",
            "",
            "NOTE",
            (
                "This step prepares aligned evidence only. "
                "No barcode position was interrogated and no lineage label "
                "was confirmed or rejected."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_arguments()

    manifest_path = args.manifest.expanduser().resolve()
    reference_path = args.reference.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if args.threads < 1:
        print(
            "ERROR: --threads must be at least 1.",
            file=sys.stderr,
        )
        return 1

    if not manifest_path.is_file():
        print(
            f"ERROR: manifest not found: {manifest_path}",
            file=sys.stderr,
        )
        return 1

    if not reference_path.is_file():
        print(
            f"ERROR: reference not found: {reference_path}",
            file=sys.stderr,
        )
        return 1

    try:
        bwa = require_program("bwa")
        samtools = require_program("samtools")

        manifest = pd.read_csv(
            manifest_path,
            dtype=str,
        ).fillna("")

        missing = REQUIRED_COLUMNS.difference(manifest.columns)

        if missing:
            raise ValueError(
                f"Manifest is missing required columns: {sorted(missing)}"
            )

        complete_manifest = manifest[
            manifest["sample_download_status"]
            .str.strip()
            .eq("complete")
        ].copy()

        if complete_manifest.empty:
            raise ValueError(
                "No complete samples were found in the download manifest."
            )

        if complete_manifest["sample_id"].duplicated().any():
            raise ValueError(
                "Duplicate sample IDs were found in the alignment input."
            )

        reference_dir = output_dir / "reference"
        bam_dir = output_dir / "bam"
        logs_dir = output_dir / "logs"
        qc_dir = output_dir / "qc"
        reports_dir = output_dir / "reports"

        for directory in [
            reference_dir,
            bam_dir,
            logs_dir,
            qc_dir,
            reports_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)

        index_prefix = reference_dir / "NC_000962.3"

        ensure_bwa_index(
            bwa=bwa,
            reference=reference_path,
            index_prefix=index_prefix,
            log_path=logs_dir / "bwa_index.log",
        )

        results = complete_manifest.copy()

        for column in [
            "bam_path",
            "bai_path",
            "flagstat_path",
            "alignment_status",
            "alignment_message",
        ]:
            if column not in results.columns:
                results[column] = ""

        total = len(results)

        for position, row_index in enumerate(
            results.index,
            start=1,
        ):
            sample_id = clean_text(
                results.at[row_index, "sample_id"]
            )

            read1 = Path(
                clean_text(results.at[row_index, "local_fastq_1"])
            )
            read2 = Path(
                clean_text(results.at[row_index, "local_fastq_2"])
            )

            bam_path = bam_dir / f"{sample_id}.sorted.bam"
            bai_path = bam_dir / f"{sample_id}.sorted.bam.bai"
            flagstat_path = qc_dir / f"{sample_id}.flagstat.txt"

            print()
            print(f"=== [{position}/{total}] {sample_id} ===")

            if not read1.is_file() or not read2.is_file():
                status = "failed"
                message = (
                    "missing_local_fastq: "
                    f"R1_exists={read1.is_file()}, "
                    f"R2_exists={read2.is_file()}"
                )
            elif (
                not args.force
                and bam_is_valid(
                    samtools=samtools,
                    bam_path=bam_path,
                    bai_path=bai_path,
                )
            ):
                status = "already_present_valid"
                message = "Existing BAM and BAI passed quickcheck."
                print("  Existing valid BAM found; skipping.")
            else:
                status, message = align_sample(
                    sample_id=sample_id,
                    read1=read1,
                    read2=read2,
                    index_prefix=index_prefix,
                    bam_path=bam_path,
                    bai_path=bai_path,
                    align_log=logs_dir / f"{sample_id}.bwa.log",
                    sort_log=logs_dir / f"{sample_id}.sort.log",
                    flagstat_path=flagstat_path,
                    bwa=bwa,
                    samtools=samtools,
                    threads=args.threads,
                    sort_memory=args.sort_memory,
                )

                print(f"  {status}: {message}")

            results.at[row_index, "bam_path"] = str(bam_path)
            results.at[row_index, "bai_path"] = str(bai_path)
            results.at[row_index, "flagstat_path"] = str(
                flagstat_path
            )
            results.at[row_index, "alignment_status"] = status
            results.at[row_index, "alignment_message"] = message

        results = results.sort_values(
            ["broad_lineage", "candidate_lineage", "sample_id"],
            kind="stable",
        ).reset_index(drop=True)

        alignment_manifest_path = (
            reports_dir / "pilot_alignment_manifest.csv"
        )
        report_path = (
            reports_dir / "pilot_alignment_report.txt"
        )

        results.to_csv(
            alignment_manifest_path,
            index=False,
        )

        write_report(
            report_path,
            manifest_path=manifest_path,
            manifest_sha256=sha256_file(manifest_path),
            reference_path=reference_path,
            reference_sha256=sha256_file(reference_path),
            bwa_version=tool_version(bwa, []),
            samtools_version=tool_version(
                samtools,
                ["--version"],
            ),
            threads=args.threads,
            results=results,
        )

    except (
        ValueError,
        RuntimeError,
        OSError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    status_counts = results["alignment_status"].value_counts()

    print()
    print("Alignment processing complete.")
    print(
        "Complete:               "
        f"{int(status_counts.get('complete', 0))}"
    )
    print(
        "Already present valid:  "
        f"{int(status_counts.get('already_present_valid', 0))}"
    )
    print(
        "Failed:                 "
        f"{int(status_counts.get('failed', 0))}"
    )
    print(f"Manifest: {alignment_manifest_path}")
    print(f"Report:   {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())