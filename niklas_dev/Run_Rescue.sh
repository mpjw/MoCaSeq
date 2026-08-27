# WS6

# RUN WITHIN DOCKER
working_directory=/media/rad/HDD2/MoCaSeq_runs/AGVarela_GITH/  # folder to output all the results
ref_directory=/media/rad/SSD1/MoCaSeq_ref/ # folder with reference files
script_directory=/home/rad/Packages/MoCaSeq/ # folder with the cloned github repo

# interactive
sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
mocaseq2

# install needed package
package_dir=/var/pipeline/temp/packages/
mkdir -p ${package_dir}

wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz -P ${package_dir}
tar -zxvf ${package_dir}/v2.5.2.tar.gz -C ${package_dir}
rm ${package_dir}/v2.5.2.tar.gz
cd ${package_dir}/bamtools-2.5.2
cmake -DCMAKE_INSTALL_PREFIX=/usr/local
make
make DESTDIR=${package_dir}/bamtools-2.5.2/stage/dir install
cd /var/pipeline

wget https://github.com/zengzheng123/GetBaseCountsMultiSample/archive/refs/tags/v1.2.3.tar.gz -P ${package_dir}
tar -zxvf ${package_dir}/v1.2.3.tar.gz -C ${package_dir}
rm ${package_dir}/v1.2.3.tar.gz
cd ${package_dir}/GetBaseCountsMultiSample-1.2.3
# !!!! MANUALLY CHANGE THIS IN THE MAKE FILE !!!!
#-I/var/pipeline/temp/packages/bamtools-2.5.2/stage/dir/my/install/dir/include/bamtools -L/var/pipeline/temp/packages/bamtools-2.5.2/stage/dir/my/install/dir/lib/
make
cd /var/pipeline
# exec: /var/pipeline/temp/packages/GetBaseCountsMultiSample-1.2.3/./GetBaseCountsMultiSample

species="Mouse"
config_file=/opt/MoCaSeq/config.sh
. $config_file

cohortdatabase_file=$genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz
scriptDir=/opt/MoCaSeq/niklas_dev/


# check commands
#parallel -j 1 "echo {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 1765-1 1765-2 1765-4 1765-6 ::: 1765-1 1765-2 1765-4 1765-6 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 1765-1 1765-2 1765-4 1765-6 ::: 1765-1 1765-2 1765-4 1765-6 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 757-1 757-2 757-5 757-8-1 ::: 757-1 757-2 757-5 757-8-1 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 1187-1 1187-2 1187-3 1187-4 ::: 1187-1 1187-2 1187-3 1187-4 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 1765-1 1765-2 1765-4 1765-6 ::: 1765-1 1765-2 1765-4 1765-6 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 2365-1 2365-2 2365-3 2365-5 ::: 2365-1 2365-2 2365-3 2365-5 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 2379-1 2379-2 2379-3 2379-5 2379-6 ::: 2379-1 2379-2 2379-3 2379-5 2379-6 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 2803-1 2803-2 2803-3 2803-5 2803-6 2803-7 ::: 2803-1 2803-2 2803-3 2803-5 2803-6 2803-7 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 2850-1 2850-2 2850-3 2850-5 2850-6 2850-7 ::: 2850-1 2850-2 2850-3 2850-5 2850-6 2850-7 ::: Mutect2 ::: Mutect2

echo "finished"

# MISSING FILES FOR
parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 788-6M 788-7 788-11 788-12 ::: 788-6M 788-7 788-11 788-12 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 799-1 799-4 799-7 799-8 ::: 799-1 799-4 799-7 799-8 ::: Mutect2 ::: Mutect2

parallel -j 1 "sh ${scriptDir}/SNV_RescueVCF.sh {1} Tumor {3} {2} Tumor {4} Human ${config_file} ${cohortdatabase_file}" ::: 2343-1 2343-2 2343-3 2343-5 2343-6 ::: 2343-1 2343-2 2343-3 2343-5 2343-6 ::: Mutect2 ::: Mutect2


echo "Done"


# DEBUG
name_rescue_to=1765-1
name_rescue_to_type=Tumor
name_rescue_to_method=Mutect2
name_rescue_from=1765-2
name_rescue_from_type=Tumor
name_rescue_from_method=Mutect2
species=Human
config_file=/opt/MoCaSeq/config.sh
cohortdatabase_file=/var/pipeline/ref/GRCm38.p6/MGP.v5.snp_and_indels.exclude_wild.vcf.gz
