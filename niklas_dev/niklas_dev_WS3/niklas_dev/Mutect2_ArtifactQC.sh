# to remove
$name/results/Mutect2/$name.$type.m2.filt.vcf
type=Tumor

java -jar $GATK_dir/gatk.jar LearnReadOrientationModel \
--input $name/results/Mutect2/$name.m2.f1r2.tar.gz \
--output $name/results/Mutect2/$name.$type.m2.artifact-priors.tar.gz

# no artifact filtering
java -jar $GATK_dir/gatk.jar FilterMutectCalls \
--variant $name/results/Mutect2/$name.$type.m2.vcf \
--output $name/results/Mutect2/$name.$type.m2.filt.vcf \
--reference $genome_file


# artifact filtering
java -jar $GATK_dir/gatk.jar FilterMutectCalls \
--variant $name/results/Mutect2/$name.$type.m2.vcf \
--output $name/results/Mutect2/$name.$type.m2.filt.AM.vcf \
--reference $genome_file \
-ob-priors $name/results/Mutect2/$name.$type.m2.artifact-priors.tar.gz



cat $name/results/Mutect2/$name.$type.m2.filt.vcf \
| java -jar $snpeff_dir/SnpSift.jar filter \
"(FILTER = 'PASS')" \
> $name/results/Mutect2/$name.$type.testWithArtifacts.vcf

cat $name/results/Mutect2/$name.$type.m2.filt.AM.vcf \
| java -jar $snpeff_dir/SnpSift.jar filter \
"(FILTER = 'PASS')" \
> $name/results/Mutect2/$name.$type.testNoArtifacts.vcf



# check the difference
wc -l $name/results/Mutect2/$name.$type.testWithArtifacts.vcf
wc -l $name/results/Mutect2/$name.$type.testNoArtifacts.vcf

bedtools subtract -a $name/results/Mutect2/$name.$type.testWithArtifacts.vcf -b $name/results/Mutect2/$name.$type.testNoArtifacts.vcf > diff.vcf

grep "1.*1439921" $name/results/Mutect2/$name.$type.testWithArtifacts.vcf $name/results/Mutect2/$name.$type.testNoArtifacts.vcf
grep "1.*1439921" $name/results/Mutect2/$name.$type.m2.filt.vcf $name/results/Mutect2/$name.$type.m2.filt.AM.vcf


grep "^[^#;]" $name/results/Mutect2/$name.$type.m2.filt.AM.vcf | cut -f 7 |sort |uniq -c
grep "^[^#;]" $name/results/Mutect2/$name.$type.testWithArtifacts.vcf | cut -f 7 |sort |uniq -c
