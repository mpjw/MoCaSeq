#!/bin/bash

##########################################################################################
##
## SNV_Rescue.sh
##
## "Rescue SNVs/Indels from other metastasis of the same patient.
##
##########################################################################################

#parallel --eta -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human /media/rad/SSD1/DNA/configadapted.sh' ::: Mel15OP1 Mel15OP2 ::: Mel15OP1 Mel15OP2 ::: StrelkaRNA ::: Mutect2 StrelkaRNA

# cd $name"/results/bam/"
# samtools view -H $name"RNA.Tumor.bam" > header.sam
# sed 's/chr//g' header.sam > header_corrected.sam
# samtools reheader header_corrected.sam $name"RNA.Tumor.bam" > $name"RNA.Tumor.fixed.bam"
# mv $name"RNA.Tumor.fixed.bam" $name"RNA.Tumor.bam"
# rm $name"RNA.Tumor.bam.bai"
# samtools index $name"RNA.Tumor.bam"
# rm header.sam
# rm header_corrected.sam

name_rescue_to=$1
name_rescue_to_type=$2
name_rescue_to_method=$3
name_rescue_from=$4
name_rescue_from_type=$5
name_rescue_from_method=$6
species=$7
config_file=$8
cohortdatabase_file=$9

. $config_file

GetBaseCountsMultiSample=/var/pipeline/temp/packages/GetBaseCountsMultiSample-1.2.3/./GetBaseCountsMultiSample

mkdir $name_rescue_to/results/rescued/

if [ $species = 'Human' ]; then
chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,X,Y
elif [ $species = 'Mouse' ]; then
chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y
fi

echo 'Rescuing '$name_rescue_to':'$name_rescue_from_method' from '$name_rescue_from':'$name_rescue_from_method

if [ $name_rescue_to_method = 'Mutect2' ] && [ $name_rescue_from_method = 'Mutect2' ]; then
cp $name_rescue_from/results/$name_rescue_from_method/$name_rescue_from.$name_rescue_from_method.vcf $name_rescue_to/results/rescued/
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_method.txt $name_rescue_to/results/rescued/
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_method.vcf $name_rescue_to/results/rescued/
fi

if [ $name_rescue_to_method = 'Mutect2' ] && [ $name_rescue_from_method = 'StrelkaRNA' ]; then
cp $name_rescue_from/results/$name_rescue_from_method/$name_rescue_from.$name_rescue_from_type.$name_rescue_from_method.vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_method.txt $name_rescue_to/results/rescued/
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_method.vcf $name_rescue_to/results/rescued/
fi

if [ $name_rescue_to_method = 'StrelkaRNA' ] && [ $name_rescue_from_method = 'Mutect2' ]; then
cp $name_rescue_from/results/$name_rescue_from_method/$name_rescue_from.$name_rescue_from_method.vcf $name_rescue_to/results/rescued/
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_type.$name_rescue_to_method.txt $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.txt
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_type.$name_rescue_to_method.vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf
fi

if [ $name_rescue_to_method = 'StrelkaRNA' ] && [ $name_rescue_from_method = 'StrelkaRNA' ]; then
cp $name_rescue_from/results/$name_rescue_from_method/$name_rescue_from.$name_rescue_from_type.$name_rescue_from_method.vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_type.$name_rescue_to_method.txt $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.txt
cp $name_rescue_to/results/$name_rescue_to_method/$name_rescue_to.$name_rescue_to_type.$name_rescue_to_method.vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf
fi

bgzip $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz
bcftools filter $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz -r $chromosomes -O z -o $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.filtered.vcf.gz
rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz.tbi
mv $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.filtered.vcf.gz $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz
tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz
bcftools isec -C -c both -O z -w 1 \
-o $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.cleaned.vcf.gz \
$name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz \
$cohortdatabase_file
mv $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.cleaned.vcf.gz $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz
gunzip $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz
rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz.tbi

bgzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf
tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
bcftools filter $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz -r $chromosomes -O z -o $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.filtered.vcf.gz
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz.tbi
mv $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.filtered.vcf.gz $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
bcftools isec -C -c both -O z -w 1 \
-o $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.cleaned.vcf.gz \
$name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz \
$cohortdatabase_file
mv $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.cleaned.vcf.gz $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
gunzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz.tbi

if [ $name_rescue_to_method = 'Mutect2' ]; then
${GetBaseCountsMultiSample} \
--fasta $genome_file --maq 1 \
--bam Tumor:$name_rescue_to/results/bam/$name_rescue_to.Tumor.bam \
--bam Normal:$name_rescue_to/results/bam/$name_rescue_to.Normal.bam \
--vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf \
--output $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf \
--thread 8
fi

if [ $name_rescue_to_method = 'StrelkaRNA' ]; then
${GetBaseCountsMultiSample} \
--fasta $genome_file --maq 1 \
--bam "Tumor:"$name_rescue_to"/results/bam/"$name_rescue_to"RNA."$name_rescue_to_type".bam" \
--bam Normal:$name_rescue_to/results/bam/$name_rescue_to.Normal.bam \
--vcf $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf \
--output $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf \
--thread 8
fi

bgzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf
bgzip $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
bgzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf

tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
tabix -p vcf $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz

bcftools isec -C -c none -O z -w 1 \
-o $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.vcf.gz \
$name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz \
$name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz

bcftools norm -m -any $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz \
-O z -o $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.fixed.vcf.gz

bcftools norm -m -any $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.vcf.gz \
-O z -o $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.fixed.vcf.gz

mv $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.fixed.vcf.gz $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz
mv $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.fixed.vcf.gz $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.vcf.gz

gunzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz
gunzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique.vcf.gz
gunzip $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz
gunzip $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf.gz

rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.vcf.gz.tbi
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf.gz.tbi

if [ $species = 'Human' ]; then
for annotate in $name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method $name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique; do

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar annotate \
$dbsnp_file $name_rescue_to/results/rescued/$annotate.vcf \
> $name_rescue_to/results/rescued/$annotate.ann1.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
annotate $cosmiccoding_file $name_rescue_to/results/rescued/$annotate.ann1.vcf \
> $name_rescue_to/results/rescued/$annotate.ann2.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
annotate $cosmicnoncoding_file $name_rescue_to/results/rescued/$annotate.ann2.vcf \
> $name_rescue_to/results/rescued/$annotate.ann3.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
annotate $clinvar_file $name_rescue_to/results/rescued/$annotate.ann3.vcf \
> $name_rescue_to/results/rescued/$annotate.ann4.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
annotate $gnomadexome_file $name_rescue_to/results/rescued/$annotate.ann4.vcf \
> $name_rescue_to/results/rescued/$annotate.ann5.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
annotate $gnomadgenome_file $name_rescue_to/results/rescued/$annotate.ann5.vcf \
> $name_rescue_to/results/rescued/$annotate.ann6.vcf

java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
DbNSFP -db $dbnsfp_file -v $name_rescue_to/results/rescued/$annotate.ann6.vcf \
-f MetaLR_pred,MetaSVM_pred,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PROVEAN_pred \
> $name_rescue_to/results/rescued/$annotate.ann7.vcf

java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
-csvStats $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats \
$name_rescue_to/results/rescued/$annotate.ann7.vcf \
> $name_rescue_to/results/rescued/$annotate.annotated.vcf

cat $name_rescue_to/results/rescued/$annotate.annotated.vcf \
| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
> $name_rescue_to/results/rescued/$annotate.annotated.one.vcf

java -jar $snpeff_dir/SnpSift.jar extractFields \
$name_rescue_to/results/rescued/$annotate.annotated.one.vcf \
 CHROM POS REF ALT "GEN[Tumor].VF" "GEN[Tumor].RD" "GEN[Tumor].AD" \
 "GEN[Normal].RD" "GEN[Normal].AD" ANN[*].GENE  ANN[*].EFFECT \
 ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
 dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF CNT_Coding \
 CNT_NonCoding CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred \
 dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_PROVEAN_pred \
 > $name_rescue_to/results/rescued/$annotate.txt

rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.txt
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.txt
rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann1.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann2.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann3.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann4.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann5.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann6.vcf
rm -f $name_rescue_to/results/rescued/$annotate.ann7.vcf
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf
rm -f $name_rescue_to/results/rescued/$annotate.annotated.one.vcf
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats.genes.txt
done

elif [ $species = 'Mouse' ]; then
for annotate in $name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method $name_rescue_to.$name_rescue_to_method.$name_rescue_from.$name_rescue_from_method.unique; do

java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
-csvStats $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats \
$name_rescue_to/results/rescued/$annotate.vcf \
> $name_rescue_to/results/rescued/$annotate.annotated.vcf

cat $name_rescue_to/results/rescued/$annotate.annotated.vcf \
| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
> $name_rescue_to/results/rescued/$annotate.annotated.one.vcf

java -jar $snpeff_dir/SnpSift.jar extractFields \
$name_rescue_to/results/rescued/$annotate.annotated.one.vcf \
CHROM POS REF ALT "GEN[Tumor].VF" "GEN[Tumor].RD" "GEN[Tumor].AD" \
"GEN[Normal].RD" "GEN[Normal].AD" ANN[*].GENE  ANN[*].EFFECT \
ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
 > $name_rescue_to/results/rescued/$annotate.txt

rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.txt
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.txt
rm -f $name_rescue_to/results/rescued/$name_rescue_from.$name_rescue_from_method.vcf
rm -f $name_rescue_to/results/rescued/$name_rescue_to.$name_rescue_to_method.vcf
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats.genes.txt
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf.stats
rm -f $name_rescue_to/results/rescued/$annotate.annotated.one.vcf
rm -f $name_rescue_to/results/rescued/$annotate.annotated.vcf
done
fi

rm -f snpEff_summary.html
