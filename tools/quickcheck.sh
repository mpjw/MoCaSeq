#!/bin/bash

# Quickcheck script for output from Mouse Cancer Sequencing (MoCaSeq) pipeline

set -e

VERBOSITY=0
SEQUENCING_TYPE="WGS"
EXPERIMENT_TYPE="tumor_only"
PATH_TO_RESULTS="$(PWD)"
SAMPLES=()

usage() {
	echo "
Usage: $(basename "$0") [options] [sample1/ sample2/ ...]
Options:
  -v, --verbosity			Verbose output level 0 (default) to 3
  -s, --sequencing_type 	Type of sequencing data. One of "WGS", "lcWGS", or "WES"
  -e, --experiment_type 	Type of cancer experiment. One of "tumor_only" or "matched"
  -h, --help				Show this help message

Performs a fast check on results from MoCaSeq pipeline. The results are assumed
to be located in the current working directory. To check specific samples just
list the folder names after the options. This will check '$PWD/sampleX'.
"
}

# argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbosity)
            VERBOSITY="$2"
            shift 2
            ;;
        -s|--sequencing_type)
            SEQUENCING_TYPE="$2"
            shift 2
            ;;
        -e|--experiment_type)
            EXPERIMENT_TYPE="$2"
            shift 2
        -h|--help)
            usage
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)
            SAMPLES+=("$1")
            shift
            ;;
    esac
done

# check sample dirs exists and are readable
VALID_SAMPLES=()
for s in "${SAMPLES[@]}"; do
    if [[ -d "$s" ]]; then
		if [[ -r "$s" ]]; then
			VALID_SAMPLES+=("$s")
		else
			echo "Sample '$PATH_TO_RESULTS/$s' not readable. Skipping." >&2
    else
        echo "Sample '$s' not found in '$PATH_TO_RESULTS'. Skipping." >&2
    fi
done

# exit if there are no valid sample dirs
if [[ ${#VALID_SAMPLES[@]} -eq 0 ]]; then
    echo "No valid samples found to process!"
    exit 0
fi

if [[ $VERBOSITY -gt 1 ]]; then
	echo "Assuming sequencing type: $SEQUENCING_TYPE"
	echo "Assuming experiment type: $EXPERIMENT_TYPE"
	echo "Processing samples: $VALID_SAMPLES"
fi

# main quick check loop
for sample in "${VALID_SAMPLES[@]}"; do
	sample_status=

	# detect pipeline version (i.e. bash or nextflow)
	pipeline_version=
	if [[ -e "${sample}/results/QC/${sample}.report.txt" ]]; then
		pipeline_version="bash"
	else
		pipeline_version="nextflow"
	fi

	if [[ "$pipeline_version" == "bash"]]; then
		# bash checks
		if [[ "$EXPERIMENT_TYPE" == "matched" ]]; then

			# matched data
			wc -l "${sample}/results/Mutect2/${sample}.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
			# Copy Number callers (i.e. CNVKit, HMMCopy or Copywriter) segments
			wc -l "${sample}/results/CNVKit/matched/${sample}.cns"
			wc -l "${sample}/results/HMMCopy/${sample}.HMMCopy.1000.segments.txt"
			wc -l "${sample}/results/Copywriter/${sample}.Copywriter.segments.Mode.txt"

			if [[ $VERBOSITY -gt 1 ]]; then
				# check normal only if provided and verbosity

				# Mutect output (post processed files and more than header content)
				wc -l "${sample}/results/Mutect2/${sample}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
				# Copy Number callers (i.e. CNVKit, HMMCopy or Copywriter) segments
				wc -l "${sample}/results/CNVKit/single/${sample}.Normal.cns"
				wc -l "${sample}/results/HMMCopy/${sample}.Normal.1000.wig"
			fi

			# check normal bam regardless of verbosity
			du -h "${sample}/results/bam/${sample}.Normal."*[^i]
		elif [[ $VERBOSITY -gt 1 ]]; then
			# Mutect output (post processed files and more than header content)
			wc -l "${sample}/results/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
			# Copy Number callers (i.e. CNVKit, HMMCopy or Copywriter) segments
			wc -l "${sample}/results/CNVKit/single/${sample}.Tumor.cns"
			wc -l "${sample}/results/HMMCopy/${sample}.Tumor.1000.wig"
			wc -l "${sample}/results/Copywriter/${sample}.Copywriter.segments.Mode.txt"
			
			# QC information when sample finished
			tail -n 2 "${sample}/results/QC/${sample}.report.txt"
		fi

		# alignment file size
		du -h "${sample}/results/bam/${sample}.Tumor."*[^i]
	else
		# nextflow checks
		if [[ -n $(echo $sample | grep '_R') ]] ; then
            echo ": detected nextflow version normal"
            # check Mutect output (post processed files and more than header contet)
            wc -l ${sample}/results/Mutect2/${sample}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt*
            # check CNVKit and HMMCopy segments
            wc -l "${sample}/results/CNVKit/${sample}.Normal.CNVKit.cns"
            wc -l ${sample}/results/HMMCopy/${sample}.Normal.HMMCopy.1000.segments.txt
		fi
        # check Mutect output (post processed files and more than header contet)
        wc -l "${sample}/results/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt"*
        # check CNVKit and HMMCopy segments
        wc -l "${sample}/results/CNVKit/${sample}.Tumor.CNVKit.cns"
        wc -l "${sample}/results/HMMCopy/${sample}.Tumor.HMMCopy.1000.segments.txt"
		du -h "${sample}/results/bam/${sample}.Tumor."*[^i]
		if [[ "$EXPERIMENT_TYPE" == "matched" ]]; then
			du -h "${sample}/results/bam/${sample}.Normal."*[^i]
		fi
	fi

	# build output string
	output_details="$sample: OK"
	if [[ -z "$sample_status" ]]; then
		(( VERBOSITY > 0 )) output_details="$output_details ($pipeline_version)"
		# (( $VERBOSITY > 1 )) output_details="$output_details\n"
	else
		(( VERBOSITY > 0 )) output_details="$output_details ($pipeline_version)\n$sample_status"
	fi
	echo "$output_details"

done
