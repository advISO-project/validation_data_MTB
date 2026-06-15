#!/usr/bin/env python3

import pandas as pd

barcode = "references/tbdb/parsed_barcode_markers.tsv"
df = pd.read_csv(barcode, sep="\t")

df["major_lineage"] = (
    df["lineage_marker"]
    .astype(str)
    .str.replace("lineage", "", regex=False)
    .str.extract(r"^([0-9]+|La[0-9]+|M)", expand=False)
)

print("\nMarker counts by major lineage/group:")
print(df["major_lineage"].value_counts(dropna=False).sort_index())

print("\nMarker counts by broad lineage:")
print(df["broad_lineage"].value_counts(dropna=False))

print("\nTotal markers:")
print(len(df))
