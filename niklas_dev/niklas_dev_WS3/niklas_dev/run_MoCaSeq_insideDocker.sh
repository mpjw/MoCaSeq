#!/bin/bash
#sudo chmod 666 /var/run/docker.sock

working_directory=/media/rad/HDD1/hPDAC_ProbesV7
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
fastq_directory=/media/rad/HDD1/WES_ILSe-18855/data_combined/
script_directory=/media/rad/HDD1/MoCaSeq/

sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref \
-v ${script_directory}:/opt/MoCaSeq \
-v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
mocaseq-human2

# as root user
sudo docker run \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref \
-v ${script_directory}:/opt/MoCaSeq \
mocaseq-human2




# name=hPDAC02_LivMet-1
# species=Human
# config_file=/opt/MoCaSeq/config.sh
# filtering="none"
# artefact_type="no"
# GATK=GATK=4.1.7.0
# . $config_file






# hMANEC debug
working_directory=/media/rad/HDD1/hMANEC_2
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
fastq_directory=/media/rad/Elements/X204SC20113232-Z02-F002/raw_data/
script_directory=/media/rad/HDD1/MoCaSeq/

working_directory=/media/rad/SSD1/MoCaSeq_runs/hMANEC_2
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
script_directory=/media/rad/HDD1/MoCaSeq/
fastq_directory=/media/rad/Elements/X204SC20113232-Z02-F002/raw_data/

sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref \
-v ${script_directory}:/opt/MoCaSeq \
-v ${fastq_directory}:/var/pipeline/fastqs \
-v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
mocaseq-human2

# set this and run initParams_Debug script
name=3SH1V3
name=YHJGYY
sequencing_type="WGS"
types="Tumor Normal"
runmode="MS"
species=Human
repeat_mapping=yes

fastq_normal_1=fastqs/X_1_${name}_N1_D1/X_1_${name}_N1_D1_FKDL202622646-1a_HKKJJDSXY_L4_1.fq.gz
fastq_normal_2=fastqs/X_1_${name}_N1_D1/X_1_${name}_N1_D1_FKDL202622646-1a_HKKJJDSXY_L4_2.fq.gz
fastq_tumor_1=fastqs/X_1_${name}_T1_D1/X_1_${name}_T1_D1_FKDL202622645-1a_HKKJJDSXY_L4_1.fq.gz
fastq_tumor_2=fastqs/X_1_${name}_T1_D1/X_1_${name}_T1_D1_FKDL202622645-1a_HKKJJDSXY_L4_2.fq.gz










sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref \
-v ${script_directory}:/opt/MoCaSeq \
-v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
mocaseq-human2




while read name normal tumor; do
  #echo Running "$name with NORMAL=$normal and TUMOR=$tumor"

  #/opt/MoCaSeq/MoCaSeq.sh --name $name \
  /opt/MoCaSeq/MoCaSeq_redoMutect2.sh --name $name \
  -tf "/var/pipeline/data_combined/"$tumor"_R1.fastq.gz" \
  -tr "/var/pipeline/data_combined/"$tumor"_R2.fastq.gz" \
  -nf "/var/pipeline/data_combined/"$normal"_R1.fastq.gz" \
  -nr "/var/pipeline/data_combined/"$normal"_R2.fastq.gz" \
  --species Human \
  --repeat_mapping yes \
  --sequencing_type WES \
  --quality_control yes \
  --threads 60 \
  --RAM 32 \
  --Mutect2 yes \
  --Titan yes \
  --Absolute yes \
  --Facets yes \
  --BubbleTree yes

done < Ilse18855_WESBatch05_FolderNamesFastQ_todo.txt
