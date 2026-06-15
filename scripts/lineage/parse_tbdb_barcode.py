#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def parse_barcode(barcode_file):
    rows = []

    with open(barcode_file) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            rows.append({
                "chrom": parts[0] if len(parts) > 0 else "",
                "start": parts[1] if len(parts) > 1 else "",
                "end": parts[2] if len(parts) > 2 else "",
                "lineage_marker": parts[3] if len(parts) > 3 else "",
                "marker_allele": parts[4] if len(parts) > 4 else "",
                "broad_lineage": parts[5] if len(parts) > 5 else "",
                "family": parts[6] if len(parts) > 6 else "",
                "rd": parts[7] if len(parts) > 7 else "",
                "raw_record": line
            })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Parse TBDB barcode BED file into a readable TSV."
    )
    parser.add_argument(
        "--barcode",
        default="references/tbdb/barcode.bed",
        help="Path to TBDB barcode BED file"
    )
    parser.add_argument(
        "--output",
        default="references/tbdb/parsed_barcode_markers.tsv",
        help="Output TSV file"
    )

    args = parser.parse_args()

    barcode_file = Path(args.barcode)
    output_file = Path(args.output)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    df = parse_barcode(barcode_file)
    df.to_csv(output_file, sep="\t", index=False)

    print(f"Parsed {len(df)} barcode markers")
    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
