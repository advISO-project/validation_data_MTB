#!/usr/bin/env python3

import pandas as pd
import requests
from pathlib import Path
from io import StringIO

INPUT = Path("lineage/candidate_samples.csv")
OUTPUT = Path("lineage/candidate_samples_with_reads.csv")

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession",
    "sample_accession",
    "study_accession",
    "experiment_accession",
    "scientific_name",
    "fastq_ftp",
    "library_layout",
    "instrument_platform",
    "instrument_model",
]


def query_ena(accession):
    params = {
        "result": "read_run",
        "query": f'assembly_accession="{accession}" OR submitted_ftp LIKE "%{accession}%" OR sample_alias="{accession}"',
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": 100,
    }

    r = requests.get(ENA_URL, params=params, timeout=30)

    if r.status_code != 200 or not r.text.strip():
        return pd.DataFrame()

    try:
        return pd.read_csv(StringIO(r.text), sep="\t")
    except Exception:
        return pd.DataFrame()


def main():
    candidates = pd.read_csv(INPUT)

    rows = []

    for _, row in candidates.iterrows():
        accession = row["accession"]
        print(f"Checking {accession}")

        ena_hits = query_ena(accession)

        if ena_hits.empty:
            out = row.to_dict()
            out.update({
                "has_fastq": "no",
                "run_accession": "",
                "fastq_ftp": "",
                "read_lookup_notes": "No ENA read_run match found"
            })
            rows.append(out)
        else:
            for _, hit in ena_hits.iterrows():
                out = row.to_dict()
                out.update({
                    "has_fastq": "yes",
                    "run_accession": hit.get("run_accession", ""),
                    "fastq_ftp": hit.get("fastq_ftp", ""),
                    "read_lookup_notes": "ENA read_run match found"
                })
                rows.append(out)

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT, index=False)

    print(f"\nWrote: {OUTPUT}")
    print(result["has_fastq"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
