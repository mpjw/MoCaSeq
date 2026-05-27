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
    if [[ -n $(ls "${sample}/results/bam/" | grep Normal) ]] && [[ -d "${sample}/results/CNVKit/matched/" ]] && [[ -d "${sample}/results/QC/" ]]; then
        echo ": detected bash pipeline version"
        # check Mutect output (post processed files and more than header contet)
        wc -l ${sample}/results/Mutect2/${sample}.Mutect2.NoCommonSNPs.OnlyImpact.txt*
        # check CNVKit and HMMCopy segments
        wc -l "${sample}/results/CNVKit/matched/${sample}.cns"
        wc -l ${sample}/results/HMMCopy/${sample}.HMMCopy.1000.segments.txt
        # check number of LOH variants
        du -h ${sample}/results/LOH/${sample}.VariantsForLOH.txt*
        # QC information when sample finished
        tail -n 2 "${sample}/results/QC/${sample}.report.txt"

    else
        if [[ -n $(echo $sample | grep '_R') ]] ; then
            echo ": detected nextflow version normal"
            # check Mutect output (post processed files and more than header contet)
            wc -l ${sample}/results/Mutect2/${sample}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt*
            # check CNVKit and HMMCopy segments
            wc -l "${sample}/results/CNVKit/${sample}.Normal.CNVKit.cns"
            wc -l ${sample}/results/HMMCopy/${sample}.Normal.HMMCopy.1000.segments.txt

        elif  [[ -n $(ls "${sample}/results/Mutect2/" | grep '.matched.') ]]; then
            echo ": detected nextflow version tumor"
            # check Mutect output (post processed files and more than header contet)
            wc -l ${sample}/results/Mutect2/${sample}.matched.Mutect2.NoCommonSNPs.OnlyImpact.txt*
            # check CNVKit and HMMCopy segments
            wc -l "${sample}/results/CNVKit/${sample}.matched.CNVKit.cns"
            wc -l ${sample}/results/HMMCopy/${sample}.HMMCopy.1000.segments.txt
            # check number of LOH variants
            du -h ${sample}/results/LOH/${sample}.VariantsForLOH.txt*

        else
            echo ": skipping this folder !!! SOMETHING WENT WRONG !!!"
            echo ""
            echo ""
            echo ""

            continue
        fi
    fi

    # check bam file size (expected ~100GB)
    # match any file that does not end on i to exclude .bai but accept .bam or .bam.c4gh
    du -h ${sample}/results/bam/${sample}*[^i]
    
    echo ""
    echo ""
    echo ""
done

# WS3 code from Niklas
# echo ""
# for d in */ ; do
#     if [ "$d" == "ref/" ] || [ "$d" == "raw/" ] || [ "$d" == "temp/" ] || [ "$d" == "done/" ] || [ "$d" == "BatchSheets/" ] ; then
#         continue
#     fi
#         sample=$(basename $d)
#     echo ${sample}
#     wc -l "${sample}/results/Mutect2/${sample}.Mutect2.NoCommonSNPs.OnlyImpact.txt"
#     wc -l "${sample}/results/CNVKit/matched/${sample}.cns"
#     tail -n 2 "${sample}/results/QC/${sample}.report.txt"
#     echo ""
#     echo ""
#     echo ""
# done