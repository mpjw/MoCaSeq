#!/bin/bash

# Script to check if complete samples are really completed

samples=$(find . -maxdepth 1 -type d -a \( -name "ASHPC_*" -o -name "PCSI_*" -o -name "RAMP_*" \) | sort)

# this parameter is intended to be used with a specific sample id or the wildcard
if [ -n "$1" ]; then
    samples=$@
    echo "processing samples: $samples"
fi

echo ""
for s in $samples ; do
    sample=$(basename $s)
    printf ${sample}

    # determine if sample ran on nextflow or shell
    if [[ -d "${sample}/results/QC/" ]]; then
        echo ": detected bash pipeline version"
        # check Mutect output (post processed files and more than header contet)
        wc -l "${sample}/results/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
        # check CNVKit and HMMCopy segments
        wc -l "${sample}/results/CNVKit/single/${sample}.Tumor.cns"
        wc -l "${sample}/results/Copywriter/${sample}.Copywriter.segments.Mode.txt"
        # QC information when sample finished
        tail -n 2 "${sample}/results/QC/${sample}.report.txt"

    else # [[ -n $(ls "${sample}/results/Mutect2/" | grep '.Tumor.') ]]; then
        echo ": detected nextflow version tumor"
        # check Mutect output (post processed files and more than header contet)
        wc -l "${sample}/results/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
        # check CNVKit and HMMCopy segments
        wc -l "${sample}/results/CNVKit/${sample}.Tumor.CNVKit.cns"
        wc -l "${sample}/results/HMMCopy/${sample}.Tumor.HMMCopy.1000.segments.txt"

    fi

    # check bam file size (expected ~100GB)
    # match any file that does not end on i to exclude .bai but accept .bam or .bam.c4gh
    du -h "${sample}/results/bam/${sample}"*[^i]
    
    echo ""
done

