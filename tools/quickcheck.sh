#!/bin/bash

# Quick check script for output from Mouse Cancer Sequencing (MoCaSeq) pipeline

VERBOSITY=0
SEQUENCING_TYPE="WGS"
EXPERIMENT_TYPE="tumor_only"
IN_DIR="$PWD"
SAMPLES=()

usage() {
	echo "
Usage: $(basename "$0") [options] [sample1/ sample2/ ...]
Options:
  -v, --verbosity			Verbose output level 0 (default) to 3
  -s, --sequencing_type 	Type of sequencing data. One of "WGS", "lcWGS", or "WES"
  -e, --experiment_type 	Type of cancer experiment. One of "tumor_only" or "matched"
  -i, --input_dir		 	Path to input directory.
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
			;;
		-i|--input_dir)
            IN_DIR="$2"
            shift 2
			;;
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

# check input dir exists and is readable
if [[ ! -r "$IN_DIR" ]]; then
	echo "Cannot read files from $IN_DIR! Please provide valid input_dir"
	exit 1
fi

# check sample dirs exists and are readable
VALID_SAMPLES=()
for s in "${SAMPLES[@]}"; do
    if [[ -d "$IN_DIR/$s" ]]; then
		if [[ -r "$IN_DIR/$s" ]]; then
			VALID_SAMPLES+=("$(basename "$s")")
		else
			echo "Sample '$IN_DIR/$s' not readable. Skipping." >&2
		fi
    else
        echo "Sample '$s' not found in '$IN_DIR'. Skipping." >&2
    fi
done

# exit if there are no valid sample dirs
if [[ ${#VALID_SAMPLES[@]} -eq 0 ]]; then
    echo "No valid samples found to process!"
    exit 0
fi

if [[ $VERBOSITY -gt 2 ]]; then
	echo "Assuming sequencing type: $SEQUENCING_TYPE"
	echo "Assuming experiment type: $EXPERIMENT_TYPE"
	echo "Processing samples: ${SAMPLES[*]}"
fi

check_result_files() {
	check_tool="$1"
	results_name="$2"
	file_glob="$3"
	fail_if_result_file_missing="$4"
	if [[ $VERBOSITY -gt 3 ]]; then
		echo "checking $results_name glob: $file_glob with $check_tool"
	fi

	if compgen -G "$file_glob" > /dev/null; then
		file_details="$file_details\n$($check_tool $file_glob)"
	else
		if [[ $fail_if_result_file_missing == "fail_if_missing" ]]; then
			sample_status="FAILED"
			results_name="$results_name (fail)"
		fi
		missing_files="$missing_files $results_name,"
		file_globs_tried="$file_globs_tried\n$check_tool $file_glob"
	fi
}

cd "$IN_DIR"
# main quick check loop
for sample in "${VALID_SAMPLES[@]}"; do
	sample_path="$sample/results"
	file_details="Results found in $IN_DIR/$sample/:"
	file_globs_tried="File globs tried in $IN_DIR/$sample/:"
	sample_status="OK"
	missing_files=
	QC_log_tail=

	# detect pipeline version (i.e. bash or nextflow)
	pipeline_version=
	if [[ -e "${sample_path}/QC/${sample}.report.txt" ]]; then
		pipeline_version="bash"
	else
		pipeline_version="nextflow"
	fi

	if [[ $VERBOSITY -gt 2 ]]; then
		echo "Checking $sample (pipeline version: $pipeline_version)"
	fi

	if [[ "$pipeline_version" == "bash" ]]; then
		# bash checks
		if [[ "$EXPERIMENT_TYPE" == "matched" ]]; then
			# matched data
			if [[ $VERBOSITY -gt 2 ]]; then
				echo "* checking $sample matched"
			fi
			# Mutect output (post processed files and more than header content)
			check_result_files 'wc -l' matched.Mutect2 "${sample_path}/Mutect2/${sample}.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
			# Copy Number callers (i.e. CNVKit, HMMCopy or Copywriter) segments
			check_result_files 'wc -l' matched.CNVKit "${sample_path}/CNVKit/matched/${sample}.cns"
			check_result_files 'wc -l' matched.HMMCopy "${sample_path}/HMMCopy/${sample}.HMMCopy.1000.segments.txt"
			check_result_files 'wc -l' matched.Copywriter "${sample_path}/Copywriter/${sample}.Copywriter.segments.Mode.txt"

			# check normal data as well
			if [[ $VERBOSITY -gt 2 ]]; then
				echo "* checking $sample Normal"
			fi
			check_result_files 'wc -l' Normal.Mutect2 "${sample_path}/Mutect2/${sample}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
			check_result_files 'wc -l' Normal.CNVKit "${sample_path}/CNVKit/single/${sample}.Normal.cns"
			check_result_files 'wc -l' Normal.HMMCopy "${sample_path}/HMMCopy/${sample}.Normal.1000.wig"
			# normal alignment file size
			check_result_files 'du -h' Normal.bam "${sample_path}/bam/${sample}.Normal."*[^i] fail_if_missing
		fi
		# check tumor only
		if [[ $VERBOSITY -gt 2 ]]; then
			echo "* checking $sample Tumor"
		fi
		check_result_files 'wc -l' Tumor.Mutect2 "${sample_path}/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
		check_result_files 'wc -l' Tumor.CNVKit "${sample_path}/CNVKit/single/${sample}.Tumor.cns"
		check_result_files 'wc -l' Tumor.HMMCopy "${sample_path}/HMMCopy/${sample}.Tumor.1000.wig"
		# tumor alignment file size
		check_result_files 'du -h' Tumor.bam "${sample_path}/bam/${sample}.Tumor."*[^i] fail_if_missing

		# QC information when sample finished
		QC_log_tail="$(tail -n 2 ${sample_path}/QC/${sample}.report.txt)"
	else
		# nextflow checks
		if echo "$sample" | grep -q '_R'; then
            # check Mutect output (post processed files and more than header content)
            check_result_files 'wc -l' Normal.Mutect2 "${sample_path}/Mutect2/${sample}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
            # check CNVKit and HMMCopy segments
            check_result_files 'wc -l' Normal.CNVKit "${sample_path}/CNVKit/${sample}.Normal.CNVKit.cns"
            check_result_files 'wc -l' Normal.HMMCopy "${sample_path}/HMMCopy/${sample}.Normal.HMMCopy.1000.segments.txt"
			check_result_files 'du -h' Normal.bam "${sample_path}/bam/${sample}.Normal."*[^i] fail_if_missing
		else
			# check Mutect output (post processed files and more than header content)
			check_result_files 'wc -l' Tumor.Mutect2 "${sample_path}/Mutect2/${sample}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
			# check CNVKit and HMMCopy segments
			check_result_files 'wc -l' Tumor.CNVKit "${sample_path}/CNVKit/${sample}.Tumor.CNVKit.cns"
			check_result_files 'wc -l' Tumor.HMMCopy "${sample_path}/HMMCopy/${sample}.Tumor.HMMCopy.1000.segments.txt"
			check_result_files 'du -h' Tumor.bam "${sample_path}/bam/${sample}.Tumor."*[^i] fail_if_missing
			if [[ "$EXPERIMENT_TYPE" == "matched" ]]; then
				check_result_files 'wc -l' 'Mutect2 (Tumor)' "${sample_path}/Mutect2/${sample}.matched.Mutect2.NoCommonSNPs.OnlyImpact.txt*" fail_if_missing
				check_result_files 'wc -l' 'CNVKit (Tumor)' "${sample_path}/CNVKit/${sample}.matched.CNVKit.cns"
				check_result_files 'wc -l' 'HMMCopy (Tumor)' "${sample_path}/HMMCopy/${sample}.matched.HMMCopy.1000.segments.txt"
			fi
		fi
	fi

	# build output string
	output_string="$sample ($pipeline_version)"
	if [[ "$sample_status" == "OK" ]]; then
		output_string="$output_string: $sample_status"
		if [[ $VERBOSITY -gt 0 ]]; then
			output_string="$output_string\nnot found:$missing_files"
		fi
		if [[ $VERBOSITY -gt 2 ]]; then
			output_string="$output_string\n$file_globs_tried"
		fi
		if [[ $VERBOSITY -gt 1 ]]; then
			output_string="$output_string\n$file_details"
		fi
	else
		output_string="$output_string: FAILED!"
		output_string="$output_string\nmissing:$missing_files"
		if [[ $VERBOSITY -gt 0 ]]; then
			output_string="$output_string\n$file_globs_tried"
		fi
		if [[ $VERBOSITY -gt 1 ]]; then
			output_string="$output_string\n$file_details"
		fi
	fi

	if [[ "$pipeline_version" == "bash" && $VERBOSITY -gt 1 ]]; then
		output_string="$output_string\n$QC_log_tail"
	fi

	printf "$output_string\n"

done
