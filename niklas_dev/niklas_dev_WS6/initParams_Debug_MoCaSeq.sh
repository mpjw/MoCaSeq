name=hPDAC02_LivMet-1
tumor=hPDAC02_LivMet-1
normal=hPDAC02_Normal
sequencing_type="WES"

types="Tumor Normal"
runmode="MS"

species=Human
repeat_mapping=yes
sequencing_type=WES
runmode=MS

fastq_normal_1=${name}/fastq/${name}.Normal.R1.fastq.gz
fastq_normal_2=${name}/fastq/${name}.Normal.R2.fastq.gz
fastq_tumor_1=${name}/fastq/${name}.Tumor.R1.fastq.gz
fastq_tumor_2=${name}/fastq/${name}.Tumor.R2.fastq.gz
bam_normal=${name}/results/bam/${name}.Normal.bam
bam_tumor=${name}/results/bam/${name}.Tumor.bam

# fastq_normal_1=/var/pipeline/raw/${normal}_R1.fastq.gz
# fastq_normal_2=/var/pipeline/raw/${normal}_R2.fastq.gz
# fastq_tumor_1=/var/pipeline/raw/${tumor}_R1.fastq.gz
# fastq_tumor_2=/var/pipeline/raw/${tumor}_R2.fastq.gz

format="tif"
resolution=20000

# default parameters
quality_control=yes
threads=60
RAM=64
temp_dir=/var/pipeline/temp
artefact_type=yes
filtering=soft
phred=
Mutect2=yes
Titan=yes
Absolute=yes
Facets=yes
BubbleTree=yes
Delly=no
GATK=4.1.7.0
test=no
memstats=0
config_file=/opt/MoCaSeq/config.sh


# set some species specific arguments
if [ $species = 'Mouse' ]; then
	echo 'Species set to Mouse'
	chromosomes=19
	echo "Species $species_lowercase"
elif [ $species = 'Human' ]; then
	echo 'Species set to Human'
	chromosomes=22
else echo "Invalid species input (${species}). Choose Mouse or Human"; sleep 2
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
else echo 'Invalid combination of input files. Either use -tf/-tr/-nf/-nr OR -tb/-nb'; #sleep 2
fi
#
# # CHECK IF ALL NEEDED FILE EXIST
# if [ "$fastq_normal_1" != "" ] && [ ! -f "$fastq_normal_1" ]; then
# echo "ERROR, File not found: $fastq_normal_1"
# sleep 2
# fi
#
# if [ "$fastq_normal_2" != "" ] && [ ! -f "$fastq_normal_2" ]; then
# echo "ERROR, File not found: $fastq_normal_2"
# sleep 2
# fi
#
# if [ "$fastq_tumor_1" != "" ] && [ ! -f "$fastq_tumor_1" ]; then
# echo "ERROR, File not found: $fastq_tumor_1"
# sleep 2
# fi
#
# if [ "$fastq_tumor_2" != "" ] && [ ! -f "$fastq_tumor_2" ]; then
# echo "ERROR, File not found: $fastq_tumor_2"
# sleep 2
# fi
#
# if [ "$bam_normal" != "" ] && [ ! -f "$bam_normal" ]; then
# echo "ERROR, File not found: $bam_normal"
# sleep 2
# fi
#
# if [ "$bam_tumor" != "" ] && [ ! -f "$bam_tumor" ]; then
# echo "ERROR, File not found: $bam_tumor"
# sleep 2
# fi


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

if [ $Mutect2 = 'yes' ] && [ $artefact_type != 'none' ]; then
	quality_control=yes
	echo 'PARAMETER CHANGE: Quality control was set to "yes", due to mutect2="yes" and artifact!="none".'
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


echo '---- Checking for available reference files ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

# Niklas: removed the "-z ...txt" since it should be "! -f" to work and therefore was not used anyways
if [ ! grep -Fxq "DONE" $genome_dir/GetReferenceData.txt ]; then
	echo '---- Reference files not found - Files will be downloaded ----' | tee -a $name/results/QC/${name}.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/${name}.report.txt
  rm -rf $genome_dir
	if [ $species = 'Mouse' ]; then
	sh $repository_dir/Preparation_GetReferenceDataMouse.sh $config_file $temp_dir
	elif [ $species = 'Human' ]; then
	sh $repository_dir/Preparation_GetReferenceDataHuman.sh $config_file $temp_dir
	fi
else
	echo '---- Reference files found! ----' | tee -a ${name}/results/QC/${name}.report.txt
fi

# check if all needed files are given in the reference folder
# || sleep 2 will exit if the secondary script fails and itself calls "sleep 2"
$repository_dir/CheckReferenceFiles.sh $FileList  || sleep 2


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
echo Assuming $artefact_type-artefacts for SNV-calling | tee -a $name/results/QC/$name.report.txt
echo $filtering is setting for filtering of SNV calls | tee -a $name/results/QC/$name.report.txt

echo Quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt

echo Using GATK v$GATK | tee -a $name/results/QC/$name.report.txt
if [ $Mutect2 = "yes" ]; then
	echo Will run Mutect2 | tee -a $name/results/QC/$name.report.txt
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
