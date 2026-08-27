#!/bin/bash

##########################################################################################
##
## all_RunTitanCNA.sh
##
## Loops over TitanCNA and TitanCNASolution.
##
##########################################################################################

# bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh BCNHE6-T1 Human /opt/MoCaSeq/config.sh WGS
# bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh hPDAC03_LivMet-1 Human /opt/MoCaSeq/config.sh WES "Tumor Normal"
# name=hPDAC02_PPT-1
# species=Human
# config_file=/opt/MoCaSeq/config.sh
# sequencing_type="WGS"
# temp_dir=/var/pipeline/temp/
# types="Tumor Normal"

name=$1
species=$2
config_file=$3
sequencing_type=$4
types=$5

. $config_file

sclust="/var/pipeline/Sclust/bin/Sclust"

Sclust_resfolder=$name/results/Sclust/
mkdir -p ${Sclust_resfolder}

# set chromosome variables
if [ $species = 'Human' ]; then
	chromosomes=22
	genomeBuild=hg38
elif [ $species = 'Mouse' ]; then
	chromosomes=19
	genomeBuild=mm10
fi
chromString=$(eval echo {1..$chromosomes} "X Y") # without chr
chromStringChr=$(set -f; printf 'chr%s ' $chromString; echo) # with chr

# set some sequencing type specific variables
if [ $sequencing_type = 'WES' ]; then
	genomepartition=1
elif [ $sequencing_type = 'WGS' ]; then
	genomepartition=2
fi


# THIS IS NOT ALWAYS NEEDED
echo "COPYING BAM FILES"
SOURCEDIR=raw/
cd $SOURCEDIR
find ${name}/results/bam/ -type f -name "*bam*" -exec cp --parents -Rb {} /var/pipeline/ \;
cd /var/pipeline/


# PREPROCESSING
echo "Sclust Preprocessing I (Removing contigs and \"chr\" in BAM file)"
for type in $types;
do
	bam=${name}/results/bam/${name}.${type}.bam &&
	outbam=$temp_dir/${name}.${type}.sclust.bam &&
	tmpbam=$temp_dir/${name}.${type}.tmp.bam &&

	# reduce to valid chromosomes (no contigs etc.)
	samtools view -o $tmpbam $bam $chromString &&

	# add "chr" to BAM file
	time samtools reheader -c 'perl -pe "s/^(@SQ.*)(SN:)([[:digit:]]+|X|Y|MT)(.*|\$)/\$1\$2chr\$3\$4/"' $tmpbam > $outbam &&

	# remove temp file
	rm $tmpbam &&

	# index new bam file
	samtools index $outbam  & PIDS="$PIDS $!"

done

wait $PIDS
PIDS=""

echo "Processing BAMs by chromosome"
# process each chromosome in parallel
#/var/pipeline/Sclust/bin/Sclust bamprocess -t /var/pipeline/temp/hPDAC02_PPT-1.Tumor.sclust.bam -n /var/pipeline/temp/hPDAC02_PPT-1.Normal.sclust.bam -o $name/results/Sclust/splitChroms/hPDAC02_PPT-1 -part 1 -build hg38 -r chr21

Sclust_tmpfolder=$temp_dir/${name}_Sclust/
mkdir -p ${Sclust_tmpfolder}/splitChroms

for chrom in $chromStringChr; do
	$sclust bamprocess \
	-t $temp_dir/${name}.Tumor.sclust.bam \
	-n $temp_dir/${name}.Normal.sclust.bam \
	-o ${Sclust_tmpfolder}/splitChroms/${name} \
	-part ${genomepartition} \
	-build ${genomeBuild} \
	-r ${chrom}  & PIDS="$PIDS $!"
done
wait $PIDS
PIDS=""

# merge
$sclust bamprocess \
	-i ${Sclust_tmpfolder}/splitChroms/${name} \
	-o ${Sclust_resfolder}/${name}

# get results
$sclust cn \
	-rc ${Sclust_resfolder}/${name}_rcount.txt \
	-snp ${Sclust_resfolder}/${name}_snps.txt \
	-vcf $temp_dir/HEADER.vcf \
	-o ${Sclust_resfolder}/${name}

# clean up
rm -r ${Sclust_tmpfolder}/
find ${temp_dir} -maxdepth 1 -type f -name "${name}*sclust*" -exec rm -r {} +
rm -r /var/pipeline/${name}/results/bam/

#bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh hPDAC05_LivMet-1 Human /opt/MoCaSeq/config.sh WES "Tumor Normal"
#bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh hPDAC05_PPT-1 Human /opt/MoCaSeq/config.sh WES "Tumor Normal"
