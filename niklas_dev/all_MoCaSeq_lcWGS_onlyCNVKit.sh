#!/bin/bash

##########################################################################################
##
## all_MoCaSeq_lcWGS.sh
##
## Main workflow for lcWGS
##
##########################################################################################

usage()
{
	echo "  Usage: $0 "
	echo "	-n, --name               Name of the sample."
	echo "	-nf, --fastq_normal_fw   Path to first normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-nr, --fastq_normal_rev  Path to second normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-tf, --fastq_tumor_fw    Path to first tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-tr, --fastq_tumor_rev   Path to second tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-s, --species   	 Determine the species. Choose from Mouse or Human. Defaults to Mouse"
	echo "	-e, --ends               Determine sequencing mode. Choose from PE or SE. Defaults to SE."
	echo "	-t, --threads            Number of CPU threads. Optional. Defaults to 8."
	echo "	-r, --RAM                Amount of Gb RAM. Optional. Defaults to 32."
	echo "	-temp, --temp_dir        Path to temporary directory. Optional. Defaults to current working directory."
	echo "	-p, --phred              If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional."
	echo "	--help                   Show this help."
  exit 1
}

# default parameters
fastq_normal_1=
fastq_normal_2=
fastq_tumor_1=
fastq_tumor_2=
sequencing_type=lcWGS
species=Mouse
quality_control=yes
threads=8
RAM=32
temp_dir=/var/pipeline/temp
phred=phred33
runmode=MS
types="Tumor Normal"
config_file=
GATK=4.1.3.0
chromosomes=19
ends=SE


# HARDCODED!!!!!
runmode="MS"
types="Tumor Normal"



# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-n|--name) shift;name="$1";;
	-nf|--fastq_normal_fw) shift;fastq_normal_1="$1";;
	-nr|--fastq_normal_rev) shift;fastq_normal_2="$1";;
	-tf|--fastq_tumor_fw) shift;fastq_tumor_1="$1";;
	-tr|--fastq_tumor_rev) shift;fastq_tumor_2="$1";;
	-s|--species) shift;species="$1";;
	-e|--ends) shift;ends="$1";;
	-t|--threads) shift;threads="$1";;
	-r|--RAM) shift;RAM="$1";;
	-temp|--temp_dir) shift;temp_dir="$1";;
    -p|--phred) shift;phred="$1";;
    --test) shift;test="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

config_file=/opt/MoCaSeq/config.sh

#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo '---- Starting Mouse Cancer Genome Analysis ----'
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/results/CNVKit


if [ ! -d $temp_dir ]; then
  mkdir -p $temp_dir/
fi

if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

# this can not be more than 30GB
tmpRAM=30
MAX_RECORDS_IN_RAM=$(expr $tmpRAM \* 250000)
HASH_TABLE_SIZE=$((tmpRAM*1000000000/500))

echo '---- Run CNVKit ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Run CNVKit ----' | tee -a $name/results/QC/$name.report.txt
# types is a space separated string, so it needs the ""
sh $repository_dir/CNV_RunCNVKit.sh $name $runmode WGS $config_file $species $threads "$types"
Rscript $repository_dir/CNV_PlotCNVKit.R $name $species $repository_dir "$types"
rm GRCm38.p6.bed # remove tmp files

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0
