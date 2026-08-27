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
	echo "	-gatk, --GATKVersion     Set to '4.1.0.0', '4.1.3.0', '4.1.4.1' or '4.1.7.0', determining which GATK version is used. Optional. Defaults to 4.1.7.0"
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
GATK=4.1.7.0
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
	-st2|--sample_type) shift;sample_type="$1";;
	-sd|--source_dir) shift;source_dir="$1";;
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

# CUSTOM
logfile=$name.log
runmode=$sample_type
if [ $runmode = "MS" ]; then
	types="Tumor Normal"
else
	types="Tumor"
fi

# function to check if a file exists and is not empty
function CheckFile {
  file=$1

  if [ ! -s "$file" ]
  then
   echo "ERROR! File not created successfully: $file"
   exit 1
   fi
}

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

#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo '---- Checking for available reference files ----' | tee -a ${logfile}

# check if all needed files are given in the reference folder
# || exit 1 will exit if the secondary script fails and itself calls "exit 1"
$repository_dir/CheckReferenceFiles.sh $FileList  || exit 1


# CREATE SUBDIRS
if [ ! -d $temp_dir ]; then
  mkdir -p $temp_dir/
fi

if [ $Mutect2 = 'yes' ]; then
	mkdir -p $name/results/Mutect2
fi

if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))

#rerouting STDERR to report file
exec 2>> ${logfile}
# exec 2>&1 | tee ${logfile} # print to log file and console



# copy the files from SOURCE to somewhere else
source_dir=NAS/$source_dir/
cd $source_dir


#find ${name}/results/Mutect2/ -type f -name "${name}.*txt" -exec cp --parents -Rb {} ~/pipeline/ \;
find ${name}/results/Mutect2/ -type f -name "${name}.*vcf" -exec cp --parents -Rb {} ~/pipeline/ \;

# exclude some larger files we do not need
find ${name}/results/Mutect2/ -type f \( -iname "${name}.*txt" -and ! -iname "*.Positions.txt" \) -exec cp --parents -Rb {} ~/pipeline/ \;  # exclude Positions files

cd ~/pipeline/

echo '---- Running VEP/MAF ----' | tee -a ${logfile}
# generate MAF file with variant effect predicted (VEP) (i.e. frameshift -> Frame_Shift_Del)
echo sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types  | tee -a ${logfile}
sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types

rm -rf '?'

mv ${logfile} $name/${logfile}

echo '---- Finished analysis of sample '$name' ----' | tee -a ${logfile}
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a ${logfile}

exit 0
