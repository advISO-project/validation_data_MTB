#!/usr/bin/env python3
import pandas as pd
import subprocess
import gzip
import shutil
from pathlib import Path

# -----------------------------
# CONFIG
# -----------------------------
truth_set_file = "./tbtamr/truth_set.csv"   # CSV with the Isolate column (SRR IDs)
run_column = "Isolate"             # column containing SRR IDs
out_dir = Path("fastqs")           # output directory
out_dir.mkdir(exist_ok=True)

# -----------------------------
# LOAD TRUTH SET
# -----------------------------
df = pd.read_csv(truth_set_file)
runs = df[run_column].dropna().unique()

print(f"Found {len(runs)} runs to download.")

# -----------------------------
# FUNCTIONS
# -----------------------------
def run_cmd(cmd):
    """Run a shell command safely."""
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR running command: {' '.join(cmd)}")
        print(f"  {e.stderr.decode()}")
        return False


def gzip_file(filepath):
    """Compress FASTQ to .gz and remove original file."""
    gz_path = str(filepath) + ".gz"
    with open(filepath, "rb") as f_in:
        with gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    filepath.unlink()  # delete original
    return gz_path


# -----------------------------
# MAIN LOOP
# -----------------------------
for run_id in runs:
    print(f"\n=== Processing {run_id} ===")

    # Check if both paired files already exist
    r1 = out_dir / f"{run_id}_1.fastq.gz"
    r2 = out_dir / f"{run_id}_2.fastq.gz"

    if r1.exists() and r2.exists():
        print("  FASTQs already exist — skipping.")
        continue

    # -----------------------------
    # 1️⃣ Prefetch SRR
    # -----------------------------
    print("  Running prefetch...")
    if not run_cmd(["prefetch", run_id]):
        print("  Skipping due to prefetch error.")
        continue

    # -----------------------------
    # 2️⃣ fasterq-dump
    # -----------------------------
    print("  Running fasterq-dump...")
    if not run_cmd(["fasterq-dump", run_id, "--split-files", "-O", str(out_dir)]):
        print("  Skipping due to fasterq-dump error.")
        continue

    # Output from fasterq-dump:
    #   SRRxxx_1.fastq
    #   SRRxxx_2.fastq
    fq1 = out_dir / f"{run_id}_1.fastq"
    fq2 = out_dir / f"{run_id}_2.fastq"

    # -----------------------------
    # 3️⃣ Gzip compression
    # -----------------------------
    if fq1.exists():
        print("  Compressing R1...")
        gzip_file(fq1)

    if fq2.exists():
        print("  Compressing R2...")
        gzip_file(fq2)

    print(f"  Done: {run_id}")

print("\nAll downloads complete.")
