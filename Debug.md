# Debugging Info for MoCaSeq pipeline

## Docker
Use this command for quick-and-dirty interactive debugging check in docker
```sh
ref_directory=/fast/SSD2/MoCaSeq_ref/GRCh38.p12/
docker run -it --rm \
--user $(id -u):$(id -g) \
-v $PWD:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq/ \
-v ${bam_directory}:/var/pipeline/raw/ \
--entrypoint=/bin/bash \
mocaseq2

# when mounting directories to docker, beware that the directory you are
# mounting into another directory must have compatible ownership/group
#
# another solution is to mount directories with incompatible permissions on 
# different mount points in the container 
docker run -it --rm \
--user $(id -u):$(id -g) \
-v $PWD:/var/pipeline/ \
-v ${ref_directory}:/fast/ref/ \
--entrypoint=/bin/bash \
mocaseq2

java -Xmx40G -Dpicard.useLegacyParser=false -jar /opt/picard-2.20.0/picard.jar SamToFastq -INPUT /var/pipeline/raw/CDS-bT1SFq.cram -FASTQ /var/pipeline/raw/CDS-bT1SFq.R1.fastq.gz -SECOND_END_FASTQ /var/pipeline/raw/CDS-bT1SFq.R2.fastq.gz -INCLUDE_NON_PF_READS true -VALIDATION_STRINGENCY LENIENT 2>&1 > CDS-bT1SFq.picard.log

# full pipeline call

docker run \
--user $(id -u):$(id -g) \
-v ${working_directory}:/var/pipeline/ \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq/ \
-v ${bam_directory}:/var/pipeline/raw/ \
--entrypoint /opt/MoCaSeq/niklas_dev/entrypoints/entrypoint_COMPASS_WGS.sh \
mocaseq2 \
-tb /var/pipeline/raw/PCSI_wgs_bam_PCSI_0465_Pa_P_526.bam \
-nb /var/pipeline/raw/PCSI_wgs_bam_PCSI_0465_Si_R.bam \
--name PCSI_0465_Pa_P \
--species Human \
--repeat_mapping yes \
--sequencing_type WGS \
--quality_control no \
--threads 40 \
--RAM 120 \
--GATKVersion 4.1.7.0 \
--filtering soft \
--artefact yes \
--Mutect2 yes \
--CNVKit yes \
--Delly no \
--BubbleTree no \
--Absolute no \
--Facets no \
--Titan no \
--para yes; mv PCSI_0465_Pa_P done_8i_1sample/
```
