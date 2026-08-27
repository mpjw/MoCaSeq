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
	echo "	-abs, --Absolute         Set to 'yes' or 'no'. Runs ABSOLUTE to estimate purity/ploidy and compute copy-numbers. Optional. Can also include information from somatic mutation data, for this set Mutect2 to 'yes'."
	echo "	-fac, --Facets           Set to 'yes' or 'no'. Runs the allele-specific copy-number caller FACETS with sample purity estimations. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-bt, --BubbleTree             Set to 'yes' or 'no'. Runs the analysis of tumoral aneuploidy and clonality. Optional. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
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
Titan=no
Absolute=no
Facets=no
BubbleTree=no
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

# CHECK IF ALL NEEDED FILE EXIST
if [ "$fastq_normal_1" != "" ] && [ ! -f "$fastq_normal_1" ]; then
echo "ERROR, File not found: $fastq_normal_1"
exit 1
fi

if [ "$fastq_normal_2" != "" ] && [ ! -f "$fastq_normal_2" ]; then
echo "ERROR, File not found: $fastq_normal_2"
exit 1
fi

if [ "$fastq_tumor_1" != "" ] && [ ! -f "$fastq_tumor_1" ]; then
echo "ERROR, File not found: $fastq_tumor_1"
exit 1
fi

if [ "$fastq_tumor_2" != "" ] && [ ! -f "$fastq_tumor_2" ]; then
echo "ERROR, File not found: $fastq_tumor_2"
exit 1
fi

if [ "$bam_normal" != "" ] && [ ! -f "$bam_normal" ]; then
echo "ERROR, File not found: $bam_normal"
exit 1
fi

if [ "$bam_tumor" != "" ] && [ ! -f "$bam_tumor" ]; then
echo "ERROR, File not found: $bam_tumor"
exit 1
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


#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo "---- Starting ${species} Cancer Genome Analysis ----"
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/
mkdir -p $name/pipeline/
mkdir -p $name/fastq/
mkdir -p $name/results/QC
mkdir -p $name/results/bam
mkdir -p $name/results/Manta
mkdir -p $name/results/Strelka
mkdir -p $name/results/msisensor

# log memory and cpu usage
logstats(){
	echo -e "date \t timestamp \t memory_usage_bytes \t cumulative_cpu_nanoseconds \t cores" > $name/results/memstats.txt
	while sleep $memstats; do ($repository_dir/Meta_logstats.sh >> $name/results/memstats.txt &) ; done
}

if [ $memstats -gt 0 ]; then
	logstats &
fi


echo '---- Checking for available reference files ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

# Niklas: removed the "-z ...txt" since it should be "! -f" to work and therefore was not used anyways
if grep --quiet DONE $genome_dir/GetReferenceData.txt; then
	echo '---- Reference files not found - Files will be downloaded ----' | tee -a $name/results/QC/${name}.report.txt
	# echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/${name}.report.txt
  # rm -rf $genome_dir
	# if [ $species = 'Mouse' ]; then
	# sh $repository_dir/Preparation_GetReferenceDataMouse.sh $config_file $temp_dir
	# elif [ $species = 'Human' ]; then
	# sh $repository_dir/Preparation_GetReferenceDataHuman.sh $config_file $temp_dir
	# fi
else
	echo '---- Reference files found! ----' | tee -a ${name}/results/QC/${name}.report.txt
fi

# check if all needed files are given in the reference folder
# || exit 1 will exit if the secondary script fails and itself calls "exit 1"
$repository_dir/CheckReferenceFiles.sh $FileList  || exit 1


# CREATE SUBDIRS
if [ ! -d $temp_dir ]; then
  mkdir -p $temp_dir/
fi

if [ $sequencing_type = 'WES' ]; then
	mkdir -p $name/results/Copywriter
fi

if [ $runmode = 'MS' ]; then
	mkdir -p $name/results/HMMCopy
fi

if [ $Mutect2 = 'yes' ]; then
	mkdir -p $name/results/Mutect2
fi

if [ $CNVKit = 'yes' ]; then
	CNVKit_folder=$name/results/CNVKit
	mkdir -p ${CNVKit_folder}
fi

if [ $Mutect2 = 'yes' ] && [ $runmode = 'MS' ]; then
	mkdir -p $name/results/LOH
fi

if [ $Titan = 'yes' ] && [ $Mutect2 = 'yes' ] && [ $runmode = 'MS' ]; then
	mkdir -p $name/results/Titan
fi

if [ $Absolute = 'yes' ]; then
	mkdir -p $name/results/ABSOLUTE
fi

if [ $Facets = 'yes' ]; then
	mkdir -p $name/results/FACETS
fi

if [ $BubbleTree = 'yes' ]; then
	mkdir -p $name/results/BubbleTree
fi

if [ $Delly = 'yes' ] && [ $runmode = 'MS' ] && [ $sequencing_type = 'WGS' ]; then
	mkdir -p $name/results/Delly
	mkdir -p $name/results/Chromothripsis
fi

if [ $species = 'Mouse' ]; then
	mkdir -p $name/results/Genotype
fi




if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))

echo "---- Starting ${species} Cancer Genome Analysis ----" | tee -a $name/results/QC/$name.report.txt
echo Starting pipeline using these settings: | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
echo Running sample named $name | tee -a $name/results/QC/$name.report.txt
echo Running in $runmode-mode | tee -a $name/results/QC/$name.report.txt

if [ $runmode = "MS" ] && [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Using $fastq_normal_1 and $fastq_normal_2 for normal fastqs | tee -a $name/results/QC/$name.report.txt
	echo Using $fastq_tumor_1 and $fastq_tumor_2 for tumor fastqs | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "MS" ] && [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	echo Using $bam_normal for normal bam | tee -a $name/results/QC/$name.report.txt
	echo Using $bam_tumor for tumor bam | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ] && [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ]  && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Assigning $fastq_tumor_1 and $fastq_tumor_2 as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ]  && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Assigning $fastq_normal_1 and $fastq_normal_2 as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "yes" ] && [ -z $bam_normal ] ; then
	echo Assigning $bam_tumor as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "yes" ] && [ -z $bam_tumor ] ; then
	echo Assigning $bam_normal as $types | tee -a $name/results/QC/$name.report.txt
fi
if [ $repeat_mapping = "no" ]; then
	echo Input BAMs will NOT be re-mapped | tee -a $name/results/QC/$name.report.txt
fi
echo Assuming that reads are from $species | tee -a $name/results/QC/$name.report.txt
echo Assuming that experiment is $sequencing_type | tee -a $name/results/QC/$name.report.txt
echo Reading configuration file from $config_file | tee -a $name/results/QC/$name.report.txt
echo Setting location of repository to $repository_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location of genome to $genome_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location for temporary files to $temp_dir| tee -a $name/results/QC/$name.report.txt
echo Filtering orientation bias artefacts for SNV-calling: $artefact_type | tee -a $name/results/QC/$name.report.txt
echo $filtering is setting for filtering of SNV calls | tee -a $name/results/QC/$name.report.txt

echo Quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt

echo Using GATK v$GATK | tee -a $name/results/QC/$name.report.txt
if [ $Mutect2 = "yes" ]; then
	echo Will run Mutect2 | tee -a $name/results/QC/$name.report.txt
fi
if [ $CNVKit = "yes" ]; then
	echo Will run CNVKit | tee -a $name/results/QC/$name.report.txt
fi
if [ $Delly = "yes" ]; then
	echo Will run Delly | tee -a $name/results/QC/$name.report.txt
fi
if [ $Titan = "yes" ]; then
	echo Will run Titan | tee -a $name/results/QC/$name.report.txt
fi
if [ $Absolute = "yes" ]; then
	echo Will run ABSOLUTE | tee -a $name/results/QC/$name.report.txt
fi
if [ $Facets = "yes" ]; then
	echo Will run FACETS | tee -a $name/results/QC/$name.report.txt
fi
if [ $BubbleTree = "yes" ]; then
	echo Will run BubbleTree | tee -a $name/results/QC/$name.report.txt
fi
echo Starting workflow using $threads CPU-threads and $RAM GB of RAM | tee -a $name/results/QC/$name.report.txt




#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt
# exec 2>&1 | tee $name/results/QC/$name.report.txt # print to log file and console


# RUN MUTECT2
# 1. matched tumor-normal
if [ $Mutect2 = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Running Mutect2 (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# TEMPORARY FIX
	for file in ${name}/results/Mutect2/${name}.m2.*; do
		mv "$file" "${file/.m2/.matched.m2}"
	done

	echo '---- Mutect2 Postprocessing (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	bash $repository_dir/SNV_Mutect2Postprocessing.sh \
	$name $species $config_file $filtering $artefact_type $GATK matched

	Rscript $repository_dir/SNV_SelectOutput.R $name Mutect2 $species $CGC_file $TruSight_file

	Rscript $repository_dir/SNV_Signatures.R $name $species

fi

# generate MAF file with variant effect predicted (VEP) (i.e. frameshift -> Frame_Shift_Del)
sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types

# 2. LOH analysis
if [ $Mutect2 = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Generate LOH data ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/LOH_GenerateVariantTable.R \
	$name $genome_file $repository_dir

	Rscript $repository_dir/LOH_MakePlots.R \
	$name $species $repository_dir
fi




echo '---- Generate and plot copy number data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
resolution=20000

# RUN CNVKit
# this subscript will take care of WGS/WES and SS/MS checks and runs
if [ $CNVKit = 'yes' ]; then
	echo '---- Starting CNVKit ----' | tee -a $name/results/QC/$name.report.txt

	# types is a space separated string, so it needs the ""
	sh $repository_dir/CNV_RunCNVKit.sh $name $runmode $sequencing_type $config_file $species $threads "$types"

	Rscript $repository_dir/CNV_PlotCNVKit.R $name $species $repository_dir "$types"
fi


# RUN Copywriter
if [ $sequencing_type = 'WES' ]; then

	echo '---- Run CopywriteR ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_RunCopywriter.R $name $species $threads $runmode $genome_dir $centromere_file $varregions_file $resolution $types

	echo '---- Export raw data and re-normalize using Mode ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_CopywriterGetRawData.R $name $runmode $types
	python2 $repository_dir/CNV_CopywriterGetModeCorrectionFactor.py $name
	Rscript $repository_dir/CNV_CopywriterGetModeCorrectionFactor.R $name $runmode $types

	echo '---- Plot CNV-profiles ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_PlotCopywriter.R $name $species $repository_dir
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes Copywriter $resolution $CGC_file $TruSight_file

	# clean up
	find "$name/results/Copywriter/" -type f -name "$name*SegmentsChromosome*" -exec rm -r {} +
fi

# RUN HMMCopy (bin-size 20000)
if [ $runmode = "MS" ]; then
	echo '---- Run HMMCopy (bin-size 20000) ----' | tee -a $name/results/QC/$name.report.txt

	if [ $sequencing_type = 'WES' ]; then
	echo '     ---- HMMCopy is executed on this WES data only for testing purposes, please handle the results with caution! ----' | tee -a $name/results/QC/$name.report.txt
	fi

	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file $resolution $types
fi

# RUN HMMCopy (bin-size 1000)
if [ $runmode = "MS" ] && [ $sequencing_type = 'WGS' ]; then

	echo '---- Run HMMCopy (bin-size 1000) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file 1000 $types
	Rscript $repository_dir/CNV_PlotHMMCopy.R $name $species $repository_dir $sequencing_type 1000 \
	$mapWig_file $gcWig_file $centromere_file $varregions_file $runmode $types
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes HMMCopy 1000 $CGC_file $TruSight_file
fi

# for MS + WES, PlotHMMCopy will be called but unable to find the HMMCopy/name.HMMCopy.20000.segments.txt which would be generated in the previous step (but is not because of sequencing_type=WGS)

# PLOT HMMCopy
if [ $runmode = "MS" ]; then

	echo '---- Plot HMMCopy ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_PlotHMMCopy.R $name $species $repository_dir $sequencing_type $resolution \
	$mapWig_file $gcWig_file $centromere_file $varregions_file
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes HMMCopy $resolution $CGC_file $TruSight_file
fi


echo '---- Run msisensor----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
if [ $runmode = "MS" ]; then
	msisensor msi -n $name/results/bam/$name.Normal.bam \
	-t $name/results/bam/$name.Tumor.bam \
	-o $name/results/msisensor/"$name".msisensor \
	-d $microsatellite_file -b $threads
elif [ $runmode = "SS" ]; then
	msisensor msi -t $name/results/bam/$name.$types.bam \
	-o $name/results/msisensor/"$name".$types.msisensor \
	-d $microsatellite_file -b $threads
fi

# PURITY ANALYSIS

# note: this can only be run with matched samples (which is catched at the start)
if [ $Titan = "yes" ]; then
	echo '---- Run TitanCNA ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/all_RunTitanCNA.R $name $species $repository_dir $resolution $mapWig_file $gcWig_file $exons_file $sequencing_type

	# rerun TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
	sh $repository_dir/all_RunTitanCNA.sh $name $repository_dir $threads $sequencing_type $species

	echo '		---- Mapping segments to genes ----' | tee -a $name/results/QC/$name.report.txt
	Rscript $repository_dir/LOH_MapSegmentsToGenes.R $name $species $genecode_file_genes $CGC_file $TruSight_file

	# cleanup TITAN temp dirs
	find . -maxdepth 1 -type d -name "run_ploidy*" -exec rm -r {} +
fi

# note: this can only be run with matched samples (which is catched at the start)
if [ $Facets = "yes" ]; then
	echo '---- Run FACETS ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# first generate SNP pileup for FACETS
	echo '     ---- Run snp-pileup ----' | tee -a $name/results/QC/$name.report.txt
	bash ${repository_dir}/all_RunFACETS.sh $name $species $sequencing_type $config_file

	cval=300
	ndepth=10

	echo '     ---- Run FACETS analysis ----' | tee -a $name/results/QC/$name.report.txt
	Rscript $repository_dir/all_RunFACETS.R $name $species $cval $ndepth
fi

# set CNV method (for ABSOLUTE and BubbleTree)
if [ $sequencing_type = 'WES' ]; then
CNVmethod="Copywriter"
elif [ $sequencing_type = 'WGS' ] && [ $runmode = "MS" ]; then
CNVmethod="HMMCopy"
elif [ $sequencing_type = 'WGS' ] && [ $runmode = "SS" ]; then
CNVmethod="CNVKit"
fi

if [ $Absolute = "yes" ]; then
	echo '---- Run ABSOLUTE ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# Mutect2 can be "yes" or "no"
	Rscript $repository_dir/all_RunABSOLUTE.R $name $CNVmethod $species $Mutect2 $runmode

	rm -r $name/results/ABSOLUTE/tmp/ # remove the tmp folder coming from R
fi

if [ $runmode = "MS" ] && [ $BubbleTree = 'yes' ]; then
	echo '---- Run BubbleTree ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/all_RunBubbleTree.R $name $CNVmethod
fi



rm -rf '?'

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0
