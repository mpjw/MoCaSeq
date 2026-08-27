sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v /media/rad/HDD2/hPDAC_WES:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v /media/nas/fastq/Studies/AGRad_hPDAC/WES/hPDAC_ProbesV7/:/var/pipeline/raw/ \
mocaseq-human

name=hPDAC02_PPT-1

cd raw/
bam_path=${name}/results/bam/
find $bam_path -type f -name "*.ba.*" -exec cp --parents -Rb {} /var/pipeline/ \;
cd /var/pipeline

bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh $name Human /opt/MoCaSeq/config.sh WES "Tumor Normal"
