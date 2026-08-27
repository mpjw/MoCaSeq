#sampletable=/home/rad/Downloads/temp1/AGRad_mPDAC_1000CLs/redoMAF_SampleSheet.tsv
#working_directory=/home/rad/Downloads/temp1/AGRad_mPDAC_WES_1000CLs_newMAFs/

sampletable=/media/rad/HDD1/data/cBioportal/AGRad_mPDAC_1000CLs_all/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/AGRad_mPDAC_WES_1000CLs_newMAFs_all/

sampletable=/media/rad/HDD1/data/cBioportal/Miguel_PhD_mCRC_organoids/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/Miguel_PhD_mCRC_organoids_newMAFs/

sampletable=/media/rad/HDD1/data/cBioportal/Miguel_PhD_mCRC_organoids_others/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/cBioportal/Miguel_PhD_mCRC_organoids_others/

sampletable=/media/rad/HDD1/data/cBioportal/hPDAC_MetastasisCohort/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/hPDAC_MetastasisCohort_newMAFs/
working_directory=/media/rad/RAD2-3000/hPDAC_MetastasisCohort_newMAFs/

sampletable=/media/rad/HDD1/data/cBioportal/mPDAC-CNV-Proteome/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/mPDAC-CNV-Proteome_newMAFs/
working_directory=/media/rad/RAD2-3000/mPDAC-CNV-Proteome_newMAFs/

sampletable=/media/rad/HDD1/data/cBioportal/AGReichert-Human-Organoids/redoMAF_SampleSheet.tsv
working_directory=/media/rad/HDD1/data/AGReichert-Human-Organoids_newMAFs/
working_directory=/media/rad/RAD2-3000/AGReichert-Human-Organoids_newMAFs/

# DO NOT CHANGE
network_directory=/run/user/1000/
ref_directory=/media/rad/HDD1/MoCaSeq_ref
script_directory=/home/rad/packages/MoCaSeq/
nParallel=3
mkdir -p ${working_directory}/temp ${working_directory}/ref ${working_directory}/raw ${working_directory}/NAS

sudo parallel --colsep '\t' -eta -j ${nParallel} "sudo time docker run \
--user $(id -u):$(id -g) \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq/ \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v ${network_directory}:/var/pipeline/NAS/ \
--entrypoint /opt/MoCaSeq/niklas_dev/redo_MAF/entrypoint_MoCaSeq_rerunVEP_forCBio.sh \
mocaseq-human \
--name {1} \
--species {2} \
--sample_type {3} \
--sequencing_type {4} \
--source_dir {5} \
--threads 20 \
--RAM 64 \
--Mutect2 yes" :::: ${sampletable}

rm -r ${working_directory}/temp ${working_directory}/ref ${working_directory}/raw ${working_directory}/NAS
rm -r ${working_directory}/*.log

exit 1




# SINGLE SAMPLE FOR DEBUG
tmpName=B34
tmpSpecies=Human
tmpSampleType=SS
tmpSourceDir=gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGReichert_hPDAC/done_SS/
tmpSeqType=WES

sudo time docker run \
-it --entrypoint=/bin/bash \
--user $(id -u):$(id -g) \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq/ \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v ${network_directory}:/var/pipeline/NAS/ \
--entrypoint /opt/MoCaSeq/niklas_dev/redo_MAF/entrypoint_MoCaSeq_rerunVEP_forCBio.sh \
mocaseq-human \
--name ${tmpName} \
--species ${tmpSpecies} \
--sample_type ${tmpSampleType} \
--source_dir ${tmpSourceDir} \
--sequencing_type ${tmpSeqType} \
--threads 20 \
--RAM 64 \
--Mutect2 yes

sudo time docker run \
-it --entrypoint=/bin/bash \
--user $(id -u):$(id -g) \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq/ \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v ${network_directory}:/var/pipeline/NAS/ \
mocaseq-human


# WITHOUT PARALLEL
# sudo time docker run \
# --user $(id -u):$(id -g) \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref/ \
# -v ${script_directory}:/opt/MoCaSeq \
# -v ${working_directory}/temp/:/var/pipeline/temp/ \
# -v ${network_directory}:/var/pipeline/NAS/ \
# --entrypoint /opt/MoCaSeq/niklas_dev/entrypoint_MoCaSeq_rerunVEP_forCBio.sh \
# mocaseq-human \
# --name mPDAC0622_OR7622f_LNMet-2 \
# --species "Mouse" \
# --sample_type "MS" \
# --source_dir "gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_mPDAC_WES_1000CLs" \
# --sequencing_type "WES" \
# --threads 20 \
# --RAM 64 \
# --Mutect2 yes




# WITHIN DOCKER
#
# sudo time docker run \
# --user $(id -u):$(id -g) \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref/ \
# -v ${script_directory}:/opt/MoCaSeq \
# -v ${working_directory}/temp/:/var/pipeline/temp/ \
# -v ${network_directory}:/var/pipeline/NAS/ \
# -it --entrypoint=/bin/bash \
# mocaseq-human
#
# name=mPDAC0622_OR7622f_LNMet-2
# config_file=/opt/MoCaSeq/config.sh
# species="Mouse"
# runmode="MS"
# types="Tumor Normal"
# repository_dir=/opt/MoCaSeq/repository/
#
# source_dir=NAS/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=fastq/Studies/AGRad_mPDAC_WES_1000CLs/
# cd $source_dir
#
# find ${name}/results/Mutect2/ -type f -name "${name}.*Mutect2.txt" -exec cp --parents -Rb {} ~/pipeline/ \;
# find ${name}/results/Mutect2/ -type f -name "${name}.*Mutect2.vcf" -exec cp --parents -Rb {} ~/pipeline/ \;
#
# cd ~/pipeline/
#
# sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types
