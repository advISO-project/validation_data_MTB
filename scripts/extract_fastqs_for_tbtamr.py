#!/usr/bin/env python3
import pandas as pd
import requests
from pathlib import Path

# -----------------------------
# CONFIG
# -----------------------------
truth_set_file = "truth_set.csv"  
run_column = "Isolate"           
out_dir = Path("fastqs")         
out_dir.mkdir(exist_ok=True)

ena_api_url = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?"
    "accession={run}&result=read_run&fields=fastq_ftp&format=tsv"
)

# -----------------------------
# LOAD TRUTH SET
# -----------------------------
df = pd.read_csv(truth_set_file)
runs = df[run_column].dropna().unique()

# -----------------------------
# DOWNLOAD FASTQs
# -----------------------------
for run_id in runs:
    print(f"\nProcessing {run_id} ...")

    # Query ENA API for fastq_ftp
    try:
        r = requests.get(ena_api_url.format(run=run_id), timeout=10)
        r.raise_for_status()
    except Exception as e:
        print(f"  Failed to query ENA for {run_id}: {e}")
        continue

    lines = r.text.strip().split("\n")
    if len(lines) < 2:
        print(f"  No FASTQ URL found for {run_id}")
        continue

    # Get fastq_ftp column (second line)
    fastq_field = lines[1].split("\t")[0].strip()
    if not fastq_field:
        print(f"  No FASTQ URL found for {run_id}")
        continue

    urls = fastq_field.split(";")  # split paired-end if present
    for url in urls:
        url = url.strip()
        if not url.startswith("http"):
            url = "https://" + url  # prepend https for requests

        filename = out_dir / url.split("/")[-1]
        if filename.exists():
            print(f"  {filename.name} already exists, skipping.")
            continue

        print(f"  Downloading {filename.name} ...")
        try:
            with requests.get(url, stream=True, timeout=60) as resp:
                resp.raise_for_status()
                with open(filename, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=1024*1024):
                        f.write(chunk)
            print(f"  {filename.name} downloaded successfully.")
        except Exception as e:
            print(f"  Failed to download {filename.name}: {e}")
