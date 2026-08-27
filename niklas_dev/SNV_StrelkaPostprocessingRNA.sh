#!/bin/bash

##########################################################################################
##
## SNV_StrelkaPostprocessing.sh
##
## Postprocessing for Strelka.
##
##########################################################################################

name=$1
type=$2
species=$3
config_file=$4
filtering=$5

. $config_file

gunzip $name/results/StrelkaRNA/Strelka-$type/results/variants/variants.vcf.gz

cat $name/results/StrelkaRNA/Strelka-$type/results/variants/variants.vcf | java -jar "$snpeff_dir"/SnpSift.jar filter "( ( FILTER = 'PASS' ) )" > $name/results/StrelkaRNA/$name.$type.str.filtered.vcf

java -jar $GATK_dir/gatk.jar SelectVariants --max-indel-size 10 -V $name/results/StrelkaRNA/$name.$type.str.filtered.vcf -O $name/results/StrelkaRNA/$name.$type.str.postprocessed.vcf

bgzip $name/results/StrelkaRNA/$name.$type.str.postprocessed.vcf

bcftools norm -m -any $name/results/StrelkaRNA/$name.$type.str.postprocessed.vcf.gz -O z -o $name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf.gz

gunzip $name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf.gz

awk '{gsub(/^chr/,""); print}' $name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf > $name/results/StrelkaRNA/$name.$type.StrelkaRNA.fixed.vcf

mv $name/results/StrelkaRNA/$name.$type.StrelkaRNA.fixed.vcf $name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf

if [ $species = 'Human' ]; then

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar annotate \
	$dbsnp_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann1.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	annotate $cosmiccoding_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann1.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann2.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $cosmicnoncoding_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann2.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann3.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $clinvar_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann3.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann4.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadexome_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann4.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann5.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadgenome_file $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann5.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann6.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	DbNSFP -db $dbnsfp_file -v $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann6.vcf \
	-f MetaLR_pred,MetaSVM_pred,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PROVEAN_pred \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann7.vcf

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf.stats \
	$name/results/StrelkaRNA/$name.$type.StrelkaRNA.ann7.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf

	cat $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.one.vcf \
	 CHROM POS REF ALT "GEN[SAMPLE1].AD[0]" \
	 "GEN[SAMPLE1].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
	 ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	 dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF CNT_Coding \
	 CNT_NonCoding CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred \
	 dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_PROVEAN_pred \
	 > $name/results/StrelkaRNA/$name.$type.StrelkaRNA.txt

elif [ $species = 'Mouse' ]; then

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf.stats \
	$name/results/StrelkaRNA/$name.$type.StrelkaRNA.vcf \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf

	cat $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/StrelkaRNA/$name.$type.StrelkaRNA.annotated.one.vcf \
	CHROM POS REF ALT "GEN[SAMPLE1].AD[0]" \
	 "GEN[SAMPLE1].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
	ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	> $name/results/StrelkaRNA/$name.$type.StrelkaRNA.txt

fi

rm $name/results/StrelkaRNA/$name.Strelka.txt
rm $name/results/StrelkaRNA/$name.Tumor.str.filtered.vcf
rm $name/results/StrelkaRNA/$name.Tumor.str.postprocessed.vcf.gz
rm $name/results/StrelkaRNA/$name.Tumor.str.postprocessed.vcf.idx
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann1.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann2.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann3.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann4.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann5.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann6.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.ann7.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.annotated.one.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.annotated.vcf
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.annotated.vcf.stats
rm $name/results/StrelkaRNA/$name.Tumor.StrelkaRNA.annotated.vcf.stats.genes.txt
