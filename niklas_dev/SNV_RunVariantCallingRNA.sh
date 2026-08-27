#!/bin/bash

##########################################################################################
##
## SNV_RunVariantCallingRNA.sh
##
## Variant-Calling pipeline from single-sample RNA using Strelka.
##
##########################################################################################

name=$1
type=$2
bam_file=$3
species=$4
threads=$5
genome_file=$6

strelka_dir=~/packages/strelka-2.9.10

mkdir $name
mkdir $name/results/
mkdir $name/results/bam
mkdir $name/results/StrelkaRNA

cp $bam_file "$name/results/bam/"$name"RNA."$type".bam"
cp ${bam_file}.bai "$name/results/bam/"$name"RNA."$type".bam.bai"

python2 $strelka_dir/bin/configureStrelkaGermlineWorkflow.py \
--bam $name"/results/bam/"$name"RNA."$type".bam" \
--ref $genome_file --runDir $name/results/StrelkaRNA/Strelka-$type \
--rna

python2 $name/results/StrelkaRNA/Strelka-$type/runWorkflow.py -m local -j $threads