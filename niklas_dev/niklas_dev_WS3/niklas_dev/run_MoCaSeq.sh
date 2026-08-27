#chmod -R +x MoCaSeq/

working_directory=/media/rad/HDD1/WES_ILSe-18855
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
fastq_directory=/media/rad/HDD1/WES_ILSe-18855/data_combined/
script_directory=/media/rad/HDD1/MoCaSeq/

# sudo docker run \
# --user $(id -u):$(id -g) \
# -it --entrypoint=/bin/bash \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref \
# -v ${script_directory}:/opt/MoCaSeq \
# mocaseq-human

# sudo docker run \
# -it --entrypoint=/bin/bash \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref \
# -v ${script_directory}:/opt/MoCaSeq \
# mocaseq-human


while read name normal tumor; do
  echo Running "$name with NORMAL=$normal and TUMOR=$tumor"

  sudo docker run \
  --user $(id -u):$(id -g) \
  -v ${working_directory}:/var/pipeline/ \
  -v ${ref_directory}:/var/pipeline/ref \
  -v ${script_directory}:/opt/MoCaSeq \
  mocaseq-human \
  --name $name \
  -tf "/var/pipeline/data_combined/${tumor}_R1.fastq.gz" \
  -tr "/var/pipeline/data_combined/${tumor}_R2.fastq.gz" \
  -nf "/var/pipeline/data_combined/${normal}_R1.fastq.gz" \
  -nr "/var/pipeline/data_combined/${normal}_R2.fastq.gz" \
  --species Human \
  --repeat_mapping yes \
  --sequencing_type WES \
  --quality_control yes \
  --threads 40 \
  --RAM 32 \
  --Mutect2 yes \
  --Titan yes \
  --Absolute yes \
  --Facets yes \
  --BubbleTree yes

done < test.txt

#sampletab.tsv
