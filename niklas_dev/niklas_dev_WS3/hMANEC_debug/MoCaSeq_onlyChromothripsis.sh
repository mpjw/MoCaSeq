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


#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

# log memory and cpu usage
logstats(){
	echo -e "date \t timestamp \t memory_usage_bytes \t cumulative_cpu_nanoseconds \t cores" > $name/results/memstats.txt
	while sleep $memstats; do ($repository_dir/Meta_logstats.sh >> $name/results/memstats.txt &) ; done
}

if [ $memstats -gt 0 ]; then
	logstats &
fi


Delly=yes
rm -r $name/results/Chromothripsis
rm -r $name/results/Chromothripsis_lowCutoffs
mkdir -p $name/results/Chromothripsis

if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))


#rerouting STDERR to report file
exec 2>> $temp_dir/$name.Chromothripsis.report.txt



# if [ $sequencing_type = 'WGS' ] && [ $Delly = 'yes' ] && [ $runmode = "MS" ]; then
# 	echo '---- Optional for WGS: Run Delly ----' | tee -a $name/results/QC/$name.report.txt
# 	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
#
# 	delly call \
# 	-o $name/results/Delly/$name.pre.bcf \
# 	-g $genome_file \
# 	$name/results/bam/$name.Tumor.bam \
# 	$name/results/bam/$name.Normal.bam
#
# 	delly filter \
# 	-f somatic -o $name/results/Delly/$name.bcf \
# 	-s $genome_dir/Samples.tsv $name/results/Delly/$name.pre.bcf
#
# 	bcftools view $name/results/Delly/$name.pre.bcf \
# 	> $name/results/Delly/$name.pre.vcf
# fi


if [ $sequencing_type = 'WGS' ] && [ $Delly = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Optional for WGS: Infer chromothripsis ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	format="tif"
	resolution=20000

	echo '---- Preparing input files and calculating coverage ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# this coverage is for control.coverage <= CoverageFilter (so lower = more stringend)
	coverage=$(sh $repository_dir/Chromothripsis_GetCoverage.sh $name)

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Delly/$name.pre.vcf \
	CHROM POS CHR2 POS2 END PE SR PRECISE IMPRECISE \
	GEN[Tumor].DR GEN[Tumor].CN GEN[Tumor].RR GEN[Tumor].RV \
	GEN[Normal].DR GEN[Normal].CN GEN[Normal].RR GEN[Normal].RV \
	MAPQ CT \
	> /var/pipeline/${name}/results/Delly/${name}.breakpoints.tab

	Rscript $repository_dir/Chromothripsis_Delly_Annotate-Filter.R \
	-n $name \
	-i $name/results/Delly/$name.breakpoints.tab \
	-o T \
	-c $coverage \
	-v 0.1 \
	-d 3000
	# outputs: $name/results/Delly/$name.breakpoints.filtered.tab

	for chr in $( seq $chromosomes ); do
	if [ $(Rscript $repository_dir/Chromothripsis_RearrangementCounter.R -i $name/results/Delly/$name.breakpoints.filtered.tab -c $chr) -ge 4 ]; then
		echo 'Analysing Chromosome '$chr
		mkdir -p $name'/results/Chromothripsis/Chr'$chr
		echo '---- Hallmark: Clustering of breakpoints for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectBreakpointClustering.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Regularity of oscillating copy number states for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_SimulateCopyNumberStates.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-o $species_lowercase -c $chr -n $name -s 1000 -a 1000 -f $format -v 0

		echo '---- Hallmark: Interspersed loss and retention of heterozygosity for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotLOHPattern.R \
		-s $name/results/HMMCopy/$name.HMMCopy.$resolution.segments.txt \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-v $name/results/LOH/$name.VariantsForLOH.txt \
		-o $species_lowercase -c $chr -n $name -f $format

		echo '---- Hallmark: Randomness of DNA fragment joins and segment order for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectRandomJoins.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Ability to walk the derivative chromosome for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_WalkDerivativeChromosome.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Visualisation: Copy number profile combined with complex rearrangements for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotRearrangementGraph.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-c $chr -n $name -f $format
	else

		# check if it did run successfully and print something to the output folder
		checkBreakpoints=$(cut -f 1 $name/results/Delly/$name.breakpoints.filtered.tab | sed -n 2p)
		if [[ $checkBreakspoints == "NA" ]]; then
			echo 'No breakpoints found for chromosome '$chr'.' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
		else
			echo 'There are too few rearrangements in chromosome '$chr'.' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
		fi

	fi
	done

fi
mv $name/results/Chromothripsis $name/results/Chromothripsis_lowCutoffs


mkdir -p $name/results/Chromothripsis
if [ $sequencing_type = 'WGS' ] && [ $Delly = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Optional for WGS: Infer chromothripsis ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	format="tif"
	resolution=20000

	echo '---- Preparing input files and calculating coverage ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# this coverage is for control.coverage <= CoverageFilter (so lower = more stringend)
	coverage=$(sh $repository_dir/Chromothripsis_GetCoverage.sh $name)

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Delly/$name.pre.vcf \
	CHROM POS CHR2 POS2 END PE SR PRECISE IMPRECISE \
	GEN[Tumor].DR GEN[Tumor].CN GEN[Tumor].RR GEN[Tumor].RV \
	GEN[Normal].DR GEN[Normal].CN GEN[Normal].RR GEN[Normal].RV \
	MAPQ CT \
	> /var/pipeline/${name}/results/Delly/${name}.breakpoints.tab

	Rscript $repository_dir/Chromothripsis_Delly_Annotate-Filter.R \
	-n $name \
	-i $name/results/Delly/$name.breakpoints.tab \
	-o T \
	-c $coverage \
	-v 0.2 \
	-d 6000
	# outputs: $name/results/Delly/$name.breakpoints.filtered.tab

	for chr in $( seq $chromosomes ); do
	if [ $(Rscript $repository_dir/Chromothripsis_RearrangementCounter.R -i $name/results/Delly/$name.breakpoints.filtered.tab -c $chr) -ge 4 ]; then
		echo 'Analysing Chromosome '$chr
		mkdir -p $name'/results/Chromothripsis/Chr'$chr
		echo '---- Hallmark: Clustering of breakpoints for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectBreakpointClustering.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Regularity of oscillating copy number states for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_SimulateCopyNumberStates.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-o $species_lowercase -c $chr -n $name -s 1000 -a 1000 -f $format -v 0

		echo '---- Hallmark: Interspersed loss and retention of heterozygosity for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotLOHPattern.R \
		-s $name/results/HMMCopy/$name.HMMCopy.$resolution.segments.txt \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-v $name/results/LOH/$name.VariantsForLOH.txt \
		-o $species_lowercase -c $chr -n $name -f $format

		echo '---- Hallmark: Randomness of DNA fragment joins and segment order for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectRandomJoins.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Ability to walk the derivative chromosome for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_WalkDerivativeChromosome.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Visualisation: Copy number profile combined with complex rearrangements for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotRearrangementGraph.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-c $chr -n $name -f $format
	else

		# check if it did run successfully and print something to the output folder
		checkBreakpoints=$(cut -f 1 $name/results/Delly/$name.breakpoints.filtered.tab | sed -n 2p)
		if [[ $checkBreakspoints == "NA" ]]; then
			echo 'No breakpoints found for chromosome '$chr'.' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
		else
			echo 'There are too few rearrangements in chromosome '$chr'.' | tee -a $name/results/Chromothripsis/$name.Chromothripsis.log
		fi

	fi
	done

fi

rm -rf '?'

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0
