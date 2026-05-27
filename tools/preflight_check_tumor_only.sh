#!/bin/bash

# Preflight check for tumor only MoCaSeq
# This script analyzes your runtime environment prior to pipeline launch.
# Make sure to export sampletable and fastq_directory!

echo "Using sample table at: $sampletable"

while IFS=, read -r name tfile; do
    tf_fastq="$fastq_directory/$tfile.R1.fastq.gz"
    tr_fastq="$fastq_directory/$tfile.R2.fastq.gz"
    [[ ! -e "$tf_fastq" ]] && echo "$tf_fastq NOT FOUND!"
    [[ ! -e "$tr_fastq" ]] && echo "$tr_fastq NOT FOUND!"
    [[ -r "$tf_fastq" && -r "$tr_fastq" ]] && echo "$name ok" || echo "$name: FILE(s) NOT READABLE!"
done < $sampletable
