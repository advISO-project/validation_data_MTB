#!/usr/bin/env python3

import pandas as pd
import requests
from pathlib import Path

PROJECT_FILE = "lineage/metadata/coll2014/coll2014_projects.txt"
OUTPUT = "lineage/metadata/coll2014/coll2014_project_metadata.tsv"

FIELDS = [
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "run_accession",
    "scientific_name",
    "country",
    "collection_date",
    "instrument_platform",
    "instrument_model",
    "library_layout",
    "fastq_ftp"
]

BASE_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
)


def fetch_project(project):
    params = {
        "accession": project,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv"
    }

    r = requests.get(BASE_URL, params=params)

    r.raise_for_status()

    return pd.read_csv(
        pd.io.common.StringIO(r.text),
        sep="\t"
    )


def main():

    with open(PROJECT_FILE) as f:
        projects = [x.strip() for x in f if x.strip()]

    dfs = []

    for project in projects:
        print(f"Fetching {project}")

        try:
            df = fetch_project(project)
            df["source_project"] = project
            dfs.append(df)

        except Exception as e:
            print(f"ERROR: {project}: {e}")

    combined = pd.concat(dfs, ignore_index=True)

    Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)

    combined.to_csv(OUTPUT, sep="\t", index=False)

    print(f"\nWrote {len(combined)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
