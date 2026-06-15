#!/bin/bash

mkdir -p lineage/tbprofiler_results

for r1 in lineage/fastqs/*_1.fastq.gz
do
    sample=$(basename "$r1" _1.fastq.gz)

    r2="lineage/fastqs/${sample}_2.fastq.gz"

    echo ""
    echo "===================================="
    echo "Running $sample"
    echo "===================================="

    tb-profiler profile \
        --read1 "$r1" \
        --read2 "$r2" \
        --prefix "$sample" \
        --dir lineage/tbprofiler_results

done

echo ""
echo "TB-Profiler batch complete"
