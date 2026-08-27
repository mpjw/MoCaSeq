ulimit -n 100000
working_directory=/media/rad/HDD1/hMANEC_combined_results
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
script_directory=/media/rad/SSD1/MoCaSeq/
network_directory=/run/user/1000/
results_directory=/media/rad/HDD1/hMANEC_combined/


sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${network_directory}:/var/pipeline/res_otherWS/ \
-v ${results_directory}:/var/pipeline/res/ \
mocaseq-human-temp


# define results folders for each workstations
WS3=res_otherWS/gvfs/sftp\:host\=imows3.med.tum.de/media/rad/HDD1/hMANEC_combined/
WS5=res_otherWS/gvfs/sftp\:host\=172.21.251.52/media/rad/HDD2/hMANEC_combined/
WS6=res/

# CNVKIT
Method=wgs
Threads=80
ref="/var/pipeline/ref/GRCh38.p12/GRCh38.p12.fna"
ann="/var/pipeline/ref/GRCh38.p12/GRCh38.p12.RefFlat"

# tumor is just a dummy, we will delete it later and only keep the reference
# cnvkit.py batch \
#   /var/pipeline/results1/1N7T85/results/bam/1N7T85.Tumor.bam \
# 	-n *Normal.bam \
# 	--output-reference PanelOfNormals/hMANEC.ref.cnn \
# 	--fasta "$ref" \
# 	--output-dir PanelOfNormals/ \
# 	-m "$Method" \
# 	-p "$Threads" \
# 	--annotate "$ann"

# run CNVKit for each sample

cnvkit.py batch \
${WS6}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir CNVKit/ \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS5}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir CNVKit/ \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS3}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir CNVKit/ \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"



repository_dir=/opt/MoCaSeq/niklas_dev/Cohort_Analysis/
species=Human
name=1N7T85


allSamples=$(find CNVKit/ -type f -name "*Tumor.cns" -exec basename {} .po \; | sed 's/.Tumor.cns//')

# this will find samples which were already plotted
finishedSamples=$(find CNVKit/ -type f -name "*.Chromosomes.CNV.CNVKit.5.pdf" -exec basename {} .po \; | sed -s 's/.Chromosomes.CNV.CNVKit.5.pdf//')
echo $allSamples | tr " " "\n" > tmp1.txt
echo $finishedSamples | tr " " "\n" > tmp2.txt
todoSamples=$(comm -3 <(sort tmp1.txt) <(sort tmp2.txt))
rm tmp1.txt tmp2.txt


for sample in $todoSamples; do
  echo $sample
  Rscript $repository_dir/CNV_PlotCNVKit_Cohort.R $sample $species $repository_dir ""
done








# PREPARE GISTIC SEG FILE
setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
Files=list.files(path = "CNVKit/", pattern = "\\.Tumor.cns$", full.names = T)

allDT <- data.table()
for(file in Files){
  dt <- fread(file, select = c("chromosome", "start", "end", "probes", "log2"))
  dt <- dt[chromosome %in% c(1:22,"X","Y")]
  dt[, name := gsub(".Tumor.cns", "", basename(file))]

  allDT <- rbind(allDT, dt)
}

setcolorder(allDT, c("name", "chromosome", "start", "end", "probes", "log2"))

allDT <- allDT[, .(Sample=name, Chromosome=chromosome, "Start Position"=start, "End Position"=end, "Num markers"=probes, Seq.CN=log2)]
fwrite(allDT, "GISTIC/all_GISTIC.seg", sep="\t", col.names = F, quote = F)
fwrite(allDT, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/all_GISTIC.seg", sep="\t", col.names = F, quote = F)






# THIS HAS TO BE RUN ON WS3
GISTIC_base=/home/rad/packages/GISTIC-2.0.23/
GISTIC=${GISTIC_base}/./gistic2
refgene=${GISTIC_base}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
working_directory=/media/rad/SSD1/niklas_tmp/hPDAC/

mkdir -p ${working_directory}/GISTIC/GISTIC_armpeel
mkdir -p ${working_directory}/GISTIC/GISTIC
mkdir -p ${working_directory}/GISTIC/GISTIC_brlen

segfile=${working_directory}/all_GISTIC.seg
segfile=${working_directory}/filtered_GISTIC.seg

cd $GISTIC_base

$GISTIC \
  -refgene $refgene \
  -seg $segfile \
  -b ${working_directory}/GISTIC/GISTIC/ \
  -run_broad_analysis 1 \
  -do_gene_gistic 1 \
  -twoside 1 \
  -ta 0.2 \
  -td 0.2

$GISTIC \
  -refgene $refgene \
  -seg $segfile \
  -b ${working_directory}/GISTIC/GISTIC_armpeel/ \
  -run_broad_analysis 1 \
  -do_gene_gistic 1 \
  -twoside 1 \
  -ta 0.2 \
  -td 0.2 \
  -armpeel

$GISTIC \
  -refgene $refgene \
  -seg $segfile \
  -b ${working_directory}/GISTIC/GISTIC_brlen/  \
  -run_broad_analysis 1  \
  -do_gene_gistic 1  \
  -twoside 1  \
  -ta 0.2  \
  -td 0.2  \
  -brlen 0.9

echo "DONE"

#
#
# mv /media/rad/SSD1/GISTIC/GISTIC/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC/
#
# mv /media/rad/SSD1/GISTIC/GISTIC_armpeel/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC_armpeel/
#
# mv /media/rad/SSD1/GISTIC/GISTIC_brlen/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC_brlen/
