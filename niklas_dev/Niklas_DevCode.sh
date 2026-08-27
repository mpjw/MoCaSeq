threadssudo docker run \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref \
-v ${script_directory}:/opt/MoCaSeq \
mocaseq-human


# INSTALL CNVKIT IN DOCKER

cd temp
wget https://bootstrap.pypa.io/get-pip.py
python3.6 get-pip.py
pip3.6 install Cython
pip3.6 install wheel nose cython numpy scipy networkx
apt-get update
apt-get install -y python-dev
apt-get install -y python3-dev
pip3.6 install --no-cache-dir pomegranate
pip3.6 install cnvkit


#INSTALL THETA2 IN DOCKER
cd temp
wget https://github.com/raphael-group/THetA/archive/v0.7.tar.gz
tar -xvzf v0.7.tar.gz
rm v0.7.tar.gz
cd THetA-0.7

mkdir bin
cp python/RunTHetA bin
ant
cp python/CreateExomeInput bin
cp matlab/runBAFGaussianModel.m bin

#apt-get install -y --no-install-recommends python-tk # needed but already in dockerfile
cd /opt
git clone https://github.com/bnpy/bnpy.git
cd bnpy/
pip2.7 install -e .
PYTHONPATH="${PYTHONPATH}:/opt/bnpy"

THetA2=/var/pipeline/temp/THetA-0.7/bin/RunTHetA





# this is for debugging
name=hPDAC05_HD_LivMet-1
bam_normal=${name}/results/bam/${name}.Normal.bam
bam_tumor=${name}/results/bam/${name}.Tumor.bam
CNVKit_folder=${name}/results/CNVKit
mkdir ${CNVKit_folder}
mkdir ${CNVKit_folder}/Chromosomes
temp_dir=/var/pipeline/temp2
genomes_dir=/var/pipeline/ref
genome_dir=$genomes_dir/GRCh38.p12
genome_file=$genome_dir/GRCh38.p12.fna
threads=40
#RefFlat=$genome_dir/GRCh38.p12.RefFlat

#samtools view -h ${bam_tumor} | head -n 1000 | samtools view -bS - > YHJGYY/results/bam/YHJGYY.SMALL.bam



































# PLOT
cnvkit.py heatmap ${CNVKit_folder}/${name}.Normal.cns -d -o ${CNVKit_folder}/${name}.Normal.Heatmap.d.pdf
cnvkit.py heatmap ${CNVKit_folder}/${name}.Tumor.cns -d -o ${CNVKit_folder}/${name}.Tumor.Heatmap.d.pdf

for chr in {1..19} X Y; do
        (
        cnvkit.py heatmap ${CNVKit_folder}/${name}.Tumor.cns -d -c $chr -o "${CNVKit_folder}/Chromosomes/${name}.Tumor.Heatmap.d.${chr}.pdf"
        cnvkit.py heatmap ${CNVKit_folder}/${name}.Normal.cns -d -c $chr -o "${CNVKit_folder}/Chromosomes/${name}.Normal.Heatmap.d.${chr}.pdf"
        ) &
done



# export for THeta2 analysis
mkdir ${name}/results/THetA2

cnvkit.py export theta \
  ${CNVKit_folder}/${name}.Tumor.cns \
  --reference ${CNVKit_folder}/Refernce.cnn \
  -v "${name}/results/Mutect2/${name}.Mutect2.vcf" \
  -i Tumor \
  -n Normal \
  -o ${name}/results/THetA2/${name}.interval_count

#RUN THETA2
$THetA2 \
  ${name}/results/THetA2/${name}.interval_count \
  -d ${name}/results/THetA2/ \
  --NUM_PROCESSES $threads

# parse the result output
Rscript "/opt/MoCaSeq/repository/all_ParseOutput.R" "theta2" ${name}/results/THetA2/${name}.n2.results






# SCLUST

# install
# cd temp/
# tar xvzf Sclust.tgz
# cd Sclust/src/
# cd /opt/Sclust/src/
# make -f makefile.ubuntu

# CODE IS NOW IN all_RunSclust.sh
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
script_directory=/media/rad/HDD2/MoCaSeq
fastq_directory=/media/rad/HDD1/hMANEC_combined/raw_fastqs/


working_directory=/media/rad/HDD2/hPDAC_WES
sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v /media/nas/fastq/Studies/AGRad_hPDAC/WES/hPDAC_ProbesV7/:/var/pipeline/raw/ \
mocaseq-human

name=hPDAC02_PPT-1
bash /opt/MoCaSeq/niklas_dev/all_RunSclust.sh $name Human /opt/MoCaSeq/config.sh WES "Tumor Normal"


working_directory=/media/rad/HDD2/niklas_temp
sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${working_directory}/temp/:/var/pipeline/temp/ \
-v /media/rad/HDD2/WES/:/var/pipeline/raw/ \
-v /media/rad/HDD2/hPDAC_WES:/var/pipeline/Sclust_vol \
mocaseq-human

types="Tumor Normal"
name=DNA396
temp_dir=temp/






# GISTIC
cd packages
wget ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz
mkdir GISTIC_2_0_23
mv GISTIC_2_0_23.tar.gz GISTIC_2_0_23/GISTIC_2_0_23.tar.gz
cd GISTIC_2_0_23
tar zxf GISTIC_2_0_23.tar.gz
cd MCR_Installer
unzip MCRInstaller.zip
./install -mode silent -agreeToLicense yes -destinationFolder "/var/pipeline/packages/GISTIC_2_0_23/"
XAPPLRESDIR=/var/pipeline/packages/GISTIC_2_0_23/v83/X11/app-defaults
