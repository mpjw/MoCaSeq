# RUN WITHIN DOCKER
# working_directory=/media/rad/HDD1/hPDAC/
# ref_directory=/media/rad/SSD1/MoCaSeq_ref/
# script_directory=/media/rad/HDD1/MoCaSeq/
# sudo docker run \
# -it --entrypoint=/bin/bash \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref/ \
# -v ${script_directory}:/opt/MoCaSeq \
# -v ${working_directory}/temp/:/var/pipeline/temp/ \
# mocaseq-human2
# config_file=/opt/MoCaSeq/config.sh
# cohortdatabase_file=ref/GRCh38.p12/00-common_all.vcf.gz
# scriptDir=/opt/MoCaSeq/niklas_dev/


config_file=/media/rad/HDD1/MoCaSeq/config_notDocker.sh
cohortdatabase_file=/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz
scriptDir=/media/rad/HDD1/MoCaSeq/niklas_dev/


# SNV_RescueVCF.sh ARGUMENTS:
# name_rescue_to=hPDAC02_LivMet-1
# name_rescue_to_type=Tumor
# name_rescue_to_method=Mutect2
# name_rescue_from=hPDAC02_PPT-1
# name_rescue_from_type=Tumor
# name_rescue_from_method=Mutect2
# species=Human
# config_file=/media/rad/HDD1/MoCaSeq/config_notDocker.sh
# cohortdatabase_file=/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz

#sudo rm -r hPDAC02_LivMet-1/results/rescued/

#sudo sh /media/rad/HDD1/MoCaSeq/niklas_dev/SNV_RescueVCF.sh hPDAC02_LivMet-1 Tumor Mutect2 hPDAC02_LivMet-1 Tumor Mutect2 Human /media/rad/HDD1/MoCaSeq/config_notDocker.sh /media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz

# sample1=hPDAC02_LivMet-1
# sample2=hPDAC02_PPT-1
# type=Tumor
# method1=Mutect2
# method2=Mutect2
# sudo parallel -j 1 "echo {1} x {3} x {2} x {4}" ::: ${sample1} ${sample2} ::: ${sample1} ${sample2} ::: ${method1} ::: ${method2}
# sudo parallel -j 1 "echo {2}" ::: ${sample1} ${sample2} ::: ${sample1} ${sample2} ::: ${method1} ::: ${method2}

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC02_LivMet-1 hPDAC02_PPT-1 ::: hPDAC02_LivMet-1 hPDAC02_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC03_LivMet-1 hPDAC03_PPT-1 ::: hPDAC03_LivMet-1 hPDAC03_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC05_LivMet-1 hPDAC05_PPT-1 ::: hPDAC05_LivMet-1 hPDAC05_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC07_LivMet-1 hPDAC07_PPT-1 ::: hPDAC07_LivMet-1 hPDAC07_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC09_LivMet-1 hPDAC09_PPT-1 ::: hPDAC09_LivMet-1 hPDAC09_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC10_LivMet-1 hPDAC10_PPT-1 ::: hPDAC10_LivMet-1 hPDAC10_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC12_LivMet-1 hPDAC12_PPT-1 ::: hPDAC12_LivMet-1 hPDAC12_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC14_LivMet-1 hPDAC14_PPT-1 ::: hPDAC14_LivMet-1 hPDAC14_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: hPDAC17_LivMet-1 hPDAC17_PPT-1 ::: hPDAC17_LivMet-1 hPDAC17_PPT-1 ::: Mutect2 ::: Mutect2

echo "Done"






















# OLD

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC03_LivMet-1 hPDAC03_PPT-1 ::: hPDAC03_LivMet-1 hPDAC03_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC05_LivMet-1 hPDAC05_PPT-1 ::: hPDAC05_LivMet-1 hPDAC05_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC07_LivMet-1 hPDAC07_PPT-1 ::: hPDAC07_LivMet-1 hPDAC07_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC09_LivMet-1 hPDAC09_PPT-1 ::: hPDAC09_LivMet-1 hPDAC09_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC10_LivMet-1 hPDAC10_PPT-1 ::: hPDAC10_LivMet-1 hPDAC10_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC12_LivMet-1 hPDAC12_PPT-1 ::: hPDAC12_LivMet-1 hPDAC12_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC14_LivMet-1 hPDAC14_PPT-1 ::: hPDAC14_LivMet-1 hPDAC14_PPT-1 ::: Mutect2 ::: Mutect2

sudo parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh \
{1} Tumor {3} {2} Tumor {4} Human \
/media/rad/SSD1/DNA/configadapted.sh \
/media/rad/SSD1/Genomes/GRCh38.p12/00-common_all.vcf.gz" ::: hPDAC17_LivMet-1 hPDAC17_PPT-1 ::: hPDAC17_LivMet-1 hPDAC17_PPT-1 ::: Mutect2 ::: Mutect2
