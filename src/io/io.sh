#!/bin/bash

# Code for input, output and paths in MoCaSeq 

# Builds path to results based on
# * tool requested
# * pipeline version (i.e. bash or nextflow)
function get_result_prefix() {
	pipeline_version="$1"
	sample_name="$2"
	tool="$3"
	matched="$4" # default: true
	tumor="$5" # default: true

	# TODO: sanity check version and tool

	# set sample type to Tumor or Normal
	(( tumor > 0 )) && type="Tumor" || type="Normal"

	if [[ "$pipeline_version" == "bash" ]]; then
		# bash version paths here
		echo "$sample_name/results/$tool/$sample_name.$tool.$type"
	elif [[ "$pipeline_version" == "bash" ]]; then
		#  nextflow paths here
		echo "$sample_name/results/$tool/$sample_name.$tool.$type"
	else
		echo "Unknown pipeline version: ${pipeline_version} expected 'bash' or 'nextflow'"
		exit 1
	fi

}

