#!/bin/bash

type=Tumor
species=Human
threads=24
genome_file=/media/rad/SSD1/Genomes/GRCh38.p12/GRCh38.primary_assembly.genome.fa
filtering=none
repository_dir=/media/rad/SSD1/DNA/repository/
config_file=/media/rad/SSD1/DNA/configadapted.sh

#CALLING ALL SAMPLES AGAIN USING THE NOT-Deduplicated BAM FILES
parallel --eta -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RunVariantCallingRNA.sh \
{} Tumor /media/rad/HDD2/ImmuNeo/{}.Tumor_Aligned.bam Human 24 /media/rad/SSD1/Genomes/GRCh38.p12/GRCh38.primary_assembly.genome.fa;
sh /media/rad/SSD1/DNA/repository/SNV_StrelkaPostprocessingRNA.sh {} Tumor Human /media/rad/SSD1/DNA/configadapted_413.sh none' ::: 1MULDR_T1 1MULDR_T2 1MULDR_T4 4ATKRF_T1 4GV9FM_T1 5WCVKM_T1 64EMZ9_T2 7NBDEC_T1 8LKUL5_T1 AGKLTM_T1 ATE46U_T1 ATE46U_T2 D67ZLM_T1 E7QJNL_T1 EYVJQP_T1 GXL1B7_T1 K72C1G_T1 LFNUX6_T1 LFNUX6_T2 LXH9EV_T1 NHW4KJ_T1 NVDER5_T1 NVDER5_T2 Q1PB42_T1 Q1PB42_T2 UA4J9Y_T1 ZFQ9G4_T1 MDR1JB_T1 EYVJQP_T1 G9VHGQ_T1 SF6ZTY_T1 KK4HL9_T1 8J4F45_T1

#RESCUE UNPAIRED
for name in 9YW2AD_T1 LRE6DV_T1 M218BR_T1 NLSTH8_T1 XVM4XC_T1; do
parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: $name ::: $name ::: Mutect2 ::: Mutect2
done

for name in 4ATKRF_T1 4GV9FM_T1 5WCVKM_T1 7NBDEC_T1 8LKUL5_T1 AGKLTM_T1 D67ZLM_T1 E7QJNL_T1 EYVJQP_T1 GXL1B7_T1 K72C1G_T1 LXH9EV_T1 MDR1JB_T1 NHW4KJ_T1 UA4J9Y_T1 ZFQ9G4_T1; do
parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: $name ::: $name ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA
done

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: 1MULDR_T1 1MULDR_T2 1MULDR_T4 ::: 1MULDR_T1 1MULDR_T2 1MULDR_T4 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: ATE46U_T1 ATE46U_T2 ::: ATE46U_T1 ATE46U_T2 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: LFNUX6_T1 LFNUX6_T2 ::: LFNUX6_T1 LFNUX6_T2 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: NVDER5_T1 NVDER5_T2 ::: NVDER5_T1 NVDER5_T2 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: 42D9U7_T1 42D9U7_T2 ::: 42D9U7_T1 42D9U7_T2 ::: Mutect2 ::: Mutect2

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: 64EMZ9_T1 64EMZ9_T2 ::: 64EMZ9_T1 64EMZ9_T2 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: Q1PB42_T1 Q1PB42_T2 ::: Q1PB42_T1 Q1PB42_T2 ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA

for name in GXL1B7_T1 K72C1G_T1 LXH9EV_T1 MDR1JB_T1 NHW4KJ_T1 UA4J9Y_T1 ZFQ9G4_T1; 
do rm -rf $name/results/rescued; done
for name in GXL1B7_T1 K72C1G_T1 LXH9EV_T1 MDR1JB_T1 NHW4KJ_T1 UA4J9Y_T1 ZFQ9G4_T1; do
parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: $name ::: $name ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA
done

for name in LRE6DV_T1 MPSKVJ_T1; do parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: $name ::: $name ::: Mutect2 ::: Mutect2
done

for name in 4GV9FM 8J4F45_T1 KK4HL9_T1 G9VHGQ_T1 SF6ZTY_T1; 
do rm -rf $name/results/rescued; done
for name in 4GV9FM 8J4F45_T1 KK4HL9_T1 G9VHGQ_T1 SF6ZTY_T1; do
parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz' ::: $name ::: $name ::: Mutect2 StrelkaRNA ::: Mutect2 StrelkaRNA
don


##rerun samples 09.09.2021 because RNA or Exome files were too short 
for name in AGKLTM_T1 MDR1JB_T1; do
for tool1 in Mutect2 StrelkaRNA; do
for tool2 in Mutect2 StrelkaRNA; do
sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh $name Tumor $tool1 $name Tumor $tool2 Human /media/rad/SSD1/DNA/configadapted.sh /media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz
done
done
done


for name1 in LFNUX6_T1 LFNUX6_T2; do
for name2 in LFNUX6_T1 LFNUX6_T2; do
for tool1 in Mutect2 StrelkaRNA; do
for tool2 in Mutect2 StrelkaRNA; do
sh /media/rad/SSD1/DNA/repository/SNV_RescueVCF.sh $name1 Tumor $tool1 $name2 Tumor $tool2 Human /media/rad/SSD1/DNA/configadapted.sh /media/rad/SSD1/Genomes/GRCh38.p12/00-all.vcf.gz
done
done
done
done

####ARCHIVE####

#Calling
parallel --eta -j 2 'sh /media/rad/SSD1/DNA/repository/SNV_RunVariantCallingRNA.sh \
{} Tumor /media/rad/HDD2/immuNeo_remap/alignment/{}.Tumor/{}.Tumor_Aligned.deduped.bam Human 24 /media/rad/SSD1/Genomes/GRCh38.p12/GRCh38.primary_assembly.genome.fa;
sh /media/rad/SSD1/DNA/repository/SNV_StrelkaPostprocessingRNA.sh {} Tumor Human /media/rad/SSD1/DNA/configadapted.sh none' ::: 1MULDR_T1 1MULDR_T2 1MULDR_T4 5WCVKM_T1 64EMZ9_T2 7NBDEC_T1 8LKUL5_T1 AGKLTM_T1 ATE46U_T1 ATE46U_T2 D67ZLM_T1 E7QJNL_T1 GXL1B7_T1 K72C1G_T1 LFNUX6_T1 LFNUX6_T2 LXH9EV_T1 NHW4KJ_T1 NVDER5_T1 NVDER5_T2 Q1PB42_T2 UA4J9Y_T1

#Calling DTA2
parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RunVariantCallingRNA.sh \
{} Tumor /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=immuneo/RNA/RNA_DTA_2/alignment/{}.Tumor/{}.Tumor_Aligned.deduped.bam Human 24 /media/rad/SSD1/Genomes/GRCh38.p12/GRCh38.primary_assembly.genome.fa;
sh /media/rad/SSD1/DNA/repository/SNV_StrelkaPostprocessingRNA.sh {} Tumor Human /media/rad/SSD1/DNA/configadapted.sh none' ::: 4GV9FM_T1 EYVJQP_T1 ZFQ9G4_T1 &

parallel -j 1 'sh /media/rad/SSD1/DNA/repository/SNV_RunVariantCallingRNA.sh \
{} Tumor /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=immuneo/RNA/RNA_DTA_2/alignment/{}.Tumor/{}.Tumor_Aligned.deduped.bam Human 24 /media/rad/SSD1/Genomes/GRCh38.p12/GRCh38.primary_assembly.genome.fa;
sh /media/rad/SSD1/DNA/repository/SNV_StrelkaPostprocessingRNA.sh {} Tumor Human /media/rad/SSD1/DNA/configadapted.sh none' ::: Q1PB42_T1 4ATKRF_T1 &