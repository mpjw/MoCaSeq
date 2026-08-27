#!/bin/bash

##########################################################################################
##
## MoCaSeq.sh
##
## Main workflow
##
##########################################################################################

usage()
{
	echo "  Usage: $0 "
	echo "	-n, --name               Name of the sample."
	echo "	-s, --species            Set to 'Mouse' or 'Human'. Defaults to Mouse."
	echo "	-nf, --fastq_normal_fw   Path to first normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-nr, --fastq_normal_rev  Path to second normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-tf, --fastq_tumor_fw    Path to first tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-tr, --fastq_tumor_rev   Path to second tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-nb, --bam_normal        Path to normal BAM. Do NOT use in combination with -nf or -nr. When used, -rm MUST be specified."
	echo "	-tb, --bam_tumor         Path to tumor BAM. Do NOT use in combination with -tf or -tr. When used, -rm MUST be specified."
	echo "	-rm, --repeat_mapping    If -nb or -tb are specified, determines whether mapping is re-done ('yes') or whether the complete mapping procedure is skipped ('no')."
	echo "	-st, --sequencing_type   Set to 'WES' or 'WGS'. Defaults to WES."
	echo "	-c, --config             Path to configuration file. Optional."
	echo "	-qc, --quality_control   Determines whether QC is done ('yes') or skipped ('no'). Optional."
	echo "	-t, --threads            Number of CPU threads. Optional. Defaults to 8."
	echo "	-r, --RAM                Amount of Gb RAM. Optional. Defaults to 32."
	echo "	-temp, --temp_dir        Path to temporary directory. Optional. Defaults to current working directory."
	echo "	-art, --artefact         Set to 'no' for no filter or 'yes' to automatically filter read-orientation bias artifacts using GATK (ob-priors). This includes OxoG oxidation artefacts and FFPE artefacts. Optional. Defaults to yes."
	echo "	-filt, --filtering       Set to 'soft' (AF >= 0.05, , Variant in Tumor >= 2, Variant in Normal <= 1, Coverage >= 5, dbSNP common for human), 'hard' (AF >= 0.1, Variant in Tumor >= 3, Variant in Normal = 0, Coverage >= 10, dbSNP all for human) or 'none' (no filters). Optional. Defaults to 'soft'."
	echo "	-p, --phred              If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional."
	echo "	-mu, --Mutect2           Set to 'yes' or 'no'. Needed for LOH analysis and Titan. Greatly increases runtime for WGS. Optional. Defaults to 'yes'."
	echo "	-ck, --CNVKit            Set to 'yes' or 'no'. Needed for LOH analysis. Optional. Defaults to 'yes'."
	echo "	-de, --Delly             Set to 'yes' or 'no'. Needed for chromothripsis inference. Do not use for WES. Optional. Defaults to 'no'. Only use in matched-sample mode."
	echo "	-ti, --Titan             Set to 'yes' or 'no'. Runs TITAN to model subclonal copy number alterations, predict LOH and estimate tumor purity. Greatly increases runtime for WGS. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-abs, --Absolute         Set to 'yes' or 'no'. Runs ABSOLUTE to estimate purity/ploidy and compute copy-numbers. Optional. Can also include information from somatic mutation data, for this set Mutect2 to 'yes'. Defaults to 'yes' for WES and WGS."
	echo "	-fac, --Facets           Set to 'yes' or 'no'. Runs the allele-specific copy-number caller FACETS with sample purity estimations. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-bt, --BubbleTree             Set to 'yes' or 'no'. Runs the analysis of tumoral aneuploidy and clonality. Optional. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and WGS. Only use in matched-sample mode."
	echo "	-gatk, --GATKVersion     Set to '4.1.0.0', '4.1.3.0', '4.1.4.1', '4.1.7.0' or '4.2.0.0', determining which GATK version is used. Optional. Defaults to 4.2.0.0 (high recommended to identify all mutations)"
	echo "	--test                   If set to 'yes': Will download reference files (if needed) and start a test run. All other parameters will be ignored"
	echo "	--memstats               If integer > 0 specified, will write timestamped memory usage and cumulative CPU time usage of the docker container to ./results/memstats.txt every <integer> seconds. Defaults to '0'."
	echo "  --para                   Run Mutect2 in parallel"
	echo "	--help                   Show this help."
  exit 1
}


# default parameters
fastq_normal_1=
fastq_normal_2=
fastq_tumor_1=
fastq_tumor_2=
bam_normal=
bam_tumor=
repeat_mapping=yes
sequencing_type=WES
quality_control=yes
threads=8
RAM=32
temp_dir=/var/pipeline/temp
artefact_type=yes
filtering=soft
phred=
Mutect2=yes
CNVKit=yes
# Titan and Facets will be set below
BubbleTree=yes
Absolute=yes
Delly=no
runmode=MS
GATK=4.2.0.0
test=no
memstats=0
config_file=
species=Mouse
para=

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-n|--name) shift;name="$1";;
	-s|--species) shift;species="$1";;
	-nf|--fastq_normal_fw) shift;fastq_normal_1="$1";;
	-nr|--fastq_normal_rev) shift;fastq_normal_2="$1";;
	-tf|--fastq_tumor_fw) shift;fastq_tumor_1="$1";;
	-tr|--fastq_tumor_rev) shift;fastq_tumor_2="$1";;
	-nb|--bam_normal) shift;bam_normal="$1";;
	-tb|--bam_tumor) shift;bam_tumor="$1";;
	-rm|--repeat_mapping) shift;repeat_mapping="$1";;
	-rq|--quality_control) shift;quality_control="$1";;
	-st|--sequencing_type) shift;sequencing_type="$1";;
	-c|--config) shift;config_file="$1";;
	-t|--threads) shift;threads="$1";;
	-r|--RAM) shift;RAM="$1";;
	-temp|--temp_dir) shift;temp_dir="$1";;
	-art|--artefact) shift;artefact_type="$1";;
	-filt|--filtering) shift;filtering="$1";;
    -p|--phred) shift;phred="$1";;
    -mu|--Mutect2) shift;Mutect2="$1";;
    -ck|--CNVKit) shift;CNVKit="$1";;
    -de|--Delly) shift;Delly="$1";;
    -ti|--Titan) shift;Titan="$1";;
		-abs|--Absolute) shift;Absolute="$1";;
		-fac|--Facets) shift;Facets="$1";;
		-bt|--BubbleTree) shift;BubbleTree="$1";;
    -gatk|--GATKVersion) shift;GATK="$1";;
    --memstats) shift;memstats="$1";;
    --test) shift;test="$1";;
		--para) shift;para="yes";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

if [ -z $config_file ]; then
	config_file=/opt/MoCaSeq/config.sh
fi

# init variables of test run
test_dir=${config_file%/*}/test

# function to check if a file exists and is not empty
function CheckFile {
  file=$1

  if [ ! -s "$file" ]
  then
   echo "ERROR! File not created successfully: $file"
   exit 1
   fi
}


#test=yes
if [ $test = 'yes' ]; then
	name=MoCaSeq_Test
	species=Mouse
	fastq_normal_1=$test_dir/Mouse.Normal.R1.fastq.gz
	fastq_normal_2=$test_dir/Mouse.Normal.R2.fastq.gz
	fastq_tumor_1=$test_dir/Mouse.Tumor.R1.fastq.gz
	fastq_tumor_2=$test_dir/Mouse.Tumor.R2.fastq.gz
	sequencing_type=WES
	bam_normal=
	bam_tumor=
	repeat_mapping=yes
	quality_control=yes
	threads=4
	RAM=8
	temp_dir=/var/pipeline/temp
	artefact_type=no
	filtering=hard
	phred=
	Mutect2=yes
	CNVKit=no
	Titan=no
	Absolute=no
	Facets=no
	BubbleTree=no
	Delly=no
	runmode=MS
fi

# set some species specific arguments
if [ $species = 'Mouse' ]; then
	echo 'Species set to Mouse'
	chromosomes=19
	echo "Species $species_lowercase"
elif [ $species = 'Human' ]; then
	echo 'Species set to Human'
	chromosomes=22
else echo "Invalid species input (${species}). Choose Mouse or Human"; exit 1
fi

# TODO: add parameter check (in general, for all values)
if [ $filtering = 'no' ]; then
	filtering="none"
fi

species_lowercase=${species,,} # set the species to lowercase to match some scripts input format (e.g. Chromothripsis)

if [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Tumor"
elif [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Normal"
elif [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="MS"
	types="Tumor Normal"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	runmode="SS"
	types="Tumor"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Normal"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	runmode="MS"
	types="Tumor Normal"
else echo 'Invalid combination of input files. Either use -tf/-tr/-nf/-nr OR -tb/-nb'; #exit 1
fi

# SET PARAMETERS FOR PURITY ANALYSIS

# TITAN ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $Titan ]; then
	Titan=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Titan ]; then
	Titan=no
fi

# If TITAN set to 'yes', forces Mutect2 to 'yes'. set to 'no' if runmode SS (needs tumor+normal WIG files)
if [ $Titan = 'yes' ] && [ $runmode = 'MS' ]; then
	Mutect2=yes
	echo 'PARAMETER CHANGE: TITAN needs LOH data from Mutect2. Mutect2 was set to "yes".'
	Titan=yes
else
	if [ $Titan = 'yes' ]; then
	echo 'PARAMETER CHANGE: TITAN needs matched tumor/normal and can not be used on single sample runs. TITAN was set to "no".'
	fi
	Titan=no
fi


# Absolute ist default 'yes' for WES and default 'no' for WGS
# Absolute needs segmented copy ratios data data ((/Copywriter/CNAprofiles/segment.Rdata or HMMCopy/HMMCopy.20000.segments.txt"))
# Absolute can optionally use somatic mutation data (Mutect2.vep.maf.fn), this will be determined by the value of Mutect2
if [ $sequencing_type = 'WES' ] && [ -z $Absolute ]; then
	Absolute=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Absolute ]; then
	Absolute=no
fi


# Facets ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $Facets ]; then
	Facets=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Facets ]; then
	Facets=no
fi

# Facets needs tumor and normal BAM files, so it can only be used in MS mode
if [ $Facets = 'yes' ] && [ $runmode = 'SS' ]; then
	Facets=no
	echo 'PARAMETER CHANGE: FACETS needs matched tumor/normal and can not be used on single sample runs. Facets was set to "no".'
fi

# BubbleTree ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $BubbleTree ]; then
	BubbleTree=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $BubbleTree ]; then
	BubbleTree=no
fi

# BubbleTree needs LOH and paired samples (tumor+normal WIG files)
# If BubbleTree set to 'yes', forces Mutect2 to 'yes'. set to 'no' if runmode SS
if [ $BubbleTree = 'yes' ] && [ $runmode = 'MS' ]; then

	if [ $Mutect2 = 'no' ]; then
	Mutect2=yes
	echo 'PARAMETER CHANGE: BubbleTree needs LOH data from Mutect2. Mutect2 was set to "yes".'
	fi

	BubbleTree=yes
else

	if [ $BubbleTree = 'yes' ]; then
	echo 'PARAMETER CHANGE: BubbleTree needs matched tumor/normal and can not be used on single sample runs. BubbleTree was set to "no".'
	fi

	BubbleTree=no
fi

echo $Absolute
echo "xx"
echo $BubbleTree
echo "xx"
echo $Facets

exit 1
