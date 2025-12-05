import pandas as pd
import requests
from pathlib import Path
import yaml


# -----------------------------
# Config
# -----------------------------
truth_set_file = "tbtamr/truth_set.csv"         
fastq_dir = Path("tbtamr/fastqs")        # folder with downloaded FASTQs
manifest_file = "tbtamr/truth_set_manifest.csv"
run_column = "Isolate"                    # column with SRR IDs

# ENA API URL template
ena_api_url = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?"
    "accession={run}&result=read_run&fields=fastq_ftp&format=tsv"
)

# -----------------------------
# Load truth set
# -----------------------------
df = pd.read_csv(truth_set_file)

# -----------------------------
# Build manifest rows
# -----------------------------
fastq_rows = []

for _, row in df.iterrows():
    isolate = row[run_column]
    sample_id = row.get("Sample_ID", isolate)  # fallback to isolate if Sample_ID not in CSV

    # Find local FASTQs
    files = sorted(fastq_dir.glob(f"{isolate}*.fastq.gz"))  # sort so _1 and _2 are consistent
    if not files:
        print(f"Warning: No local FASTQs found for {isolate}")

    # Query ENA for FTP URLs
    try:
        r = requests.get(ena_api_url.format(run=isolate), timeout=10)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        if len(lines) < 2:
            raise ValueError("No FTP URLs found")

        # Parse header and values
        header = lines[0].split("\t")
        values = lines[1].split("\t")
        if "fastq_ftp" in header:
            idx = header.index("fastq_ftp")
            ftp_field = values[idx]
            ftp_urls = ftp_field.split(";")  # paired-end URLs
        else:
            ftp_urls = [""] * len(files)

    except Exception as e:
        print(f"Warning: Could not retrieve FTP URLs for {isolate}: {e}")
        ftp_urls = [""] * len(files)

    # Match local files with FTP URLs
    for i, f in enumerate(files):
        ftp_url = f"https://{ftp_urls[i]}" if i < len(ftp_urls) else ""
        fastq_rows.append({
            "Sample_ID": sample_id,
            "Isolate": isolate,
            "Study": row["Study"],
            "AMR_first_line": row["AMR_first_line"],
            "AMR_second_line": row["AMR_second_line"],
            "AMR_count": row["AMR_count"],
            "Selection_reason": row["Selection_reason"],
            "Mutations": row.get("Mutations", ""),
            "fastq_path": str(f),
            "fastq_ftp": ftp_url
        })

# -----------------------------
# Save manifest CSV
# -----------------------------
manifest_df = pd.DataFrame(fastq_rows)
manifest_df.to_csv(manifest_file, index=False)

print(f"Manifest created: {manifest_file} with {len(manifest_df)} entries")


# -----------------------------
# Save YAML collection descriptor
# -----------------------------
samples = []
for isolate, group in manifest_df.groupby("Isolate"):
    pair_elements = []
    for _, row in group.iterrows():
        pair_elements.append({
            "class": "File",
            "path": str(Path(row["fastq_path"]).relative_to(fastq_dir.parent)),  # include tbtamr/fastqs/ prefix
            "location": row["fastq_ftp"]
        })
    samples.append({
        "class": "Collection",
        "collection_type": "pair",
        "identifier": isolate,
        "elements": pair_elements
    })

descriptor = {
    "name": "validation_set_1",
    "description": "A set of validation files derived from tbtamr",
    "version": "1.0.0",
    "samples": {
        "class": "Collection",
        "collection_type": "list:pair",
        "elements": samples
    }
}

yaml_file = "tbtamr/validation_set_manifest.yaml"
with open(yaml_file, "w") as f:
    yaml.dump(descriptor, f, sort_keys=False)

print(f"YAML collection descriptor created: {yaml_file}")
