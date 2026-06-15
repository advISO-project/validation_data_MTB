#!/usr/bin/env python3

import pandas as pd
import subprocess
from pathlib import Path

INPUT = Path("lineage/pilot_fastq_candidates.csv")
OUTDIR = Path("lineage/fastqs")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

for _, row in df.iterrows():
    run = row["run_accession"]
    urls = str(row["fastq_ftp"]).split(";")

    print(f"\nDownloading {run}")

    for url in urls:
        if url.startswith("ftp.sra.ebi.ac.uk"):
            url = "https://" + url
        elif not url.startswith("http"):
            url = "https://" + url

        filename = OUTDIR / Path(url).name

        if filename.exists():
            print(f"Skipping existing {filename}")
            continue

        subprocess.run(["wget", "-c", "-O", str(filename), url], check=True)

print("\nDone.")
