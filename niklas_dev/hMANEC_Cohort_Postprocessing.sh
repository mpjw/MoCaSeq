ulimit -n 100000
working_directory=/media/rad/HDD1/hMANEC_combined_results
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
script_directory=/media/rad/SSD1/MoCaSeq/
network_directory=/run/user/1000/
results_directory=/media/rad/HDD1/hMANEC_combined/
normals_directory=/media/rad/HDD1/hMANEC_normals/

sudo docker run \
--user $(id -u):$(id -g) \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${network_directory}:/var/pipeline/res_otherWS/ \
-v ${results_directory}:/var/pipeline/res/ \
-v ${normals_directory}:/var/pipeline/normals/ \
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

output_dir=CNVKit/

# THIS WAS THE FIRST RUN; THIS IS NOW DEPRECATED IN FAVOR OF THE BELOW CODE WITH CORRECTION AND PLOTTING
# run CNVKit for each sample, this is with a hardcoded cutoffs cbs 1e-6 (this value is refined below)
cnvkit.py batch \
${WS6}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS5}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS3}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
--diagram \
--scatter \
-m "$Method" \
-p "$Threads"


# PLOTTING
repository_dir=/opt/MoCaSeq/niklas_dev/Cohort_Analysis/
species=Human

# this will find samples which were already plotted
allSamples=$(find ${output_dir} -type f -name "*Tumor.cns" -exec basename {} .po \; | sed 's/.Tumor.cns//')
finishedSamples=$(find ${output_dir} -type f -name "*.Chromosomes.CNV.CNVKit.5.pdf" -exec basename {} .po \; | sed -s 's/.Chromosomes.CNV.CNVKit.5.pdf//')
echo $allSamples | tr " " "\n" > tmp1.txt
echo $finishedSamples | tr " " "\n" > tmp2.txt
todoSamples=$(comm -3 <(sort tmp1.txt) <(sort tmp2.txt))
rm tmp1.txt tmp2.txt

for sample in $todoSamples; do
  echo $sample
  Rscript $repository_dir/CNV_PlotCNVKit_Cohort.R $sample $species $repository_dir ""
done





# THIS DID NOT WORK WELL! (the segments were too "shattered")
# try to run the above (not batch but the corresponding code) without some regions
# cnvkit.py access $ref -o access.filtered.bed -x CNV_germline.bed -x CNV_centromeres.bed -x CNV_haploAlts.bed
#
# cnvkit.py batch \
# /var/pipeline/res/1N7T85/results/bam/1N7T85.Tumor.bam \
# 	-n normals/*Normal.bam \
# 	--output-reference CNVKit_excludeTest/hMANEC.filtered.ref.cnn \
# 	--fasta "$ref" \
# 	--output-dir CNVKit_excludeTest/ \
# 	-m "$Method" \
# 	-p "$Threads" \
# 	--annotate "$ann" \
#   --access CNVKit_excludeTest/access.filtered.bed

# cnvkit.py batch \
#   /var/pipeline/res/1N7T85/results/bam/1N7T85.Tumor.bam \
#   -n normals/*Normal.bam \
#   --output-reference CNVKit_onlyGenic/hMANEC.onlygenic.ref.cnn \
#   --output-dir CNVKit_onlyGenic/ \
#   -m "$Method" \
#   -p "$Threads" \
#   --annotate "$ann" \
#   --targets CNVKit_onlyGenic/genic_regions.bed

# edit /etc/fuse.conf to allow "allow_other" by removing the comment symbol before that line
# mkdir -p /home/rad/mounts/WS2/ /home/rad/mounts/WS3/ /home/rad/mounts/WS5/
# sshfs -o allow_other rad@172.21.243.38:/ /home/rad/mounts/WS2/
# sshfs -o allow_other rad@imows3.med.tum.de:/ /home/rad/mounts/WS3/
# sshfs -o allow_other rad@172.21.251.52:/ /home/rad/mounts/WS5/


# REFINED SEGMENTS WITH DIFFERENT CBS CUTOFFS
network_directory=/home/rad/mounts/
sudo docker run \
-it --entrypoint=/bin/bash \
-v ${working_directory}:/var/pipeline/ \
-v ${script_directory}:/opt/MoCaSeq \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${network_directory}:/var/pipeline/res_otherWS/ \
-v ${results_directory}:/var/pipeline/res/ \
-v ${normals_directory}:/var/pipeline/normals/ \
mocaseq2

# WS2=res_otherWS/gvfs/sftp\:host\=172.21.243.38/home/rad/Documents/small_scripts/MoCaSeq_DevCode/
# WS3=res_otherWS/gvfs/sftp\:host\=imows3.med.tum.de/media/rad/HDD1/hMANEC_combined/
# WS5=res_otherWS/gvfs/sftp\:host\=172.21.251.52/media/rad/HDD2/hMANEC_combined/
WS2=res_otherWS/WS2/home/rad/Documents/small_scripts/MoCaSeq_DevCode/
WS3=res_otherWS/WS3/media/rad/HDD1/hMANEC_combined
WS5=res_otherWS/WS5/media/rad/HDD2/hMANEC_combined/
WS6=res/
Method=wgs
Threads=80
ref="/var/pipeline/ref/GRCh38.p12/GRCh38.p12.fna"
ann="/var/pipeline/ref/GRCh38.p12/GRCh38.p12.RefFlat"

apt-get -y install libpoppler-cpp-dev
Rscript -e 'install.packages("pdftools")'
Rscript -e 'install.packages("ggpubr")'
Rscript -e 'install.packages("plotly")'
Rscript -e 'install.packages("htmlwidgets")'

# get backup for modified script
cp /usr/local/lib/python3.7/dist-packages/cnvlib/batch.py /usr/local/lib/python3.7/dist-packages/cnvlib/batch_backup.py

# modify the hardcoded value, this will be changed from 1e-6 to 1e-NEW
cutoffValue=3 # select 3-5
sed -Ei "s/(.threshold.: 1e-).*}/\1$cutoffValue}/g" /usr/local/lib/python3.7/dist-packages/cnvlib/batch.py
grep "threshold" /usr/local/lib/python3.7/dist-packages/cnvlib/batch.py
#output_dir=CNVKit_refine-segmentsTest/cbs_1e$cutoffValue/

output_dir=CNVKit_refineSegments_1e-3/
mkdir -p $output_dir

#
# cnvkit.py batch \
# ${WS6}/P4GNJ9/results/bam/P4GNJ9.Tumor.bam \
# -r hMANEC.ref.cnn \
# --output-dir $output_dir \
# -m "$Method" \
# -p "$Threads"

cnvkit.py batch \
${WS6}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS5}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
-m "$Method" \
-p "$Threads"

cnvkit.py batch \
${WS3}/*/results/bam/*.Tumor.bam \
-r hMANEC.ref.cnn \
--output-dir ${output_dir} \
-m "$Method" \
-p "$Threads"

echo "DONE"


# reset script (or just close docker)
# cp /usr/local/lib/python3.7/dist-packages/cnvlib/batch_backup.py /usr/local/lib/python3.7/dist-packages/cnvlib/batch.py


# this will find samples which were already corrected ()
allSamples=$(find ${output_dir} -type f -name "*Tumor.cns" -exec basename {} .po \; | sed 's/.Tumor.cns//')
finishedSamples=$(find ${output_dir} -type f -name "*.Tumor.center-mode.cns" -exec basename {} .po \; | sed -s 's/.Tumor.center-mode.cns//')
echo $allSamples | tr " " "\n" > tmp1.txt
echo $finishedSamples | tr " " "\n" > tmp2.txt
todoSamples=$(comm -3 <(sort tmp1.txt) <(sort tmp2.txt))
rm tmp1.txt tmp2.txt


# CORRECTION AND PLOTTING
repository_dir=/opt/MoCaSeq/niklas_dev/Cohort_Analysis/
repository_dir=res_otherWS/WS2/home/rad/packages/MoCaSeq/niklas_dev/Cohort_Analysis/
species=Human
allSamples=$(find ${output_dir} -type f -name "*Tumor.cns" -exec basename {} .po \; | sed 's/.Tumor.cns//')

#for sample in $todoSamples; do
for sample in $allSamples; do
  echo ${sample}
	#echo "Applying center correction"

	# cnvkit.py call -m none $output_dir/${sample}.Tumor.cns --center biweight -o $output_dir/${sample}.Tumor.center-biweight.cns
	# cnvkit.py call -m none $output_dir/${sample}.Tumor.cnr --center biweight -o $output_dir/${sample}.Tumor.center-biweight.cnr
	# cnvkit.py call -m none $output_dir/${sample}.Tumor.cns --center mode -o $output_dir/${sample}.Tumor.center-mode.cns
	# cnvkit.py call -m none $output_dir/${sample}.Tumor.cnr --center mode -o $output_dir/${sample}.Tumor.center-mode.cnr

	#echo "Plotting"
	Rscript $repository_dir/CNV_PlotCNVKit_Cohort_RefinedSegments.R $sample $species $repository_dir "" ${output_dir}

	echo "Plotting (filtered)"
	Rscript $WS2/CNV_PlotCNVKit_Cohort_RefinedSegments_Filtered.R $sample $output_dir 5
	Rscript $WS2/CNV_PlotCNVKit_Cohort_RefinedSegments_Filtered.R $sample $output_dir 2

done

echo "DONE2"








# TEST ---> CBS is the winner!
# segMethod=haar
# segMethod=hmm
# segMethod=hmm-tumor
#
# output_dir=CNVKit_refine-segmentsTest/${segMethod}
# mkdir $output_dir
#
# cnvkit.py batch \
# ${WS6}/P4GNJ9/results/bam/P4GNJ9.Tumor.bam \
# -r hMANEC.ref.cnn \
# --output-dir $output_dir \
# -m "$Method" \
# -p "$Threads" \
# --segment-method ${segMethod}





# RETRY REMOVING SOME REGIONS
# tumor is just a dummy, we will delete it later and only keep the reference
# output_dir=CNVKit_refine-segmentsTest/cbs_1e${cutoffValue}_excludeCentromeres/

# cnvkit.py access $ref -o CNVKit_excludeTest/access.filtered.bed -x CNVKit_excludeTest/GRCh38_centromeres.bed
# cnvkit.py batch \
# ${WS6}/P4GNJ9/results/bam/P4GNJ9.Tumor.bam \
# 	-n normals/*Normal.bam \
# 	--output-reference CNVKit_excludeTest/hMANEC.filtered.ref.cnn \
# 	--fasta "$ref" \
# 	--output-dir ${output_dir} \
# 	-m "$Method" \
# 	-p "$Threads" \
# 	--annotate "$ann" \
#   --access CNVKit_excludeTest/access.filtered.bed
# now copy the hMANEC.filtered.ref.cnn to the main folder and use it for -r everywhere

# cnvkit.py batch \
# ${WS6}/P4GNJ9/results/bam/P4GNJ9.Tumor.bam \
# -r CNVKit_excludeTest/hMANEC.filtered.ref.cnn \
# --output-dir $output_dir \
# -m "$Method" \
# -p "$Threads"























# OPTIONAL:get the genes files for these
Rscript /home/rad/Documents/small_scripts/MoCaSeq_DevCode/CohortAnalysis/CNV_MapSegmentsToGenes_OneCohortFolder.R ${output_dir}






# OPTIONAL: PLOT HEATMAP
CNVKit_folder=removed_samples/
CNVKit_folder=CNVKit/
CNVKit_folder=CNVKit_refineSegments_1e-3/
allSamples=$(find ${CNVKit_folder}/ -type f -name "*Tumor.cns" -exec basename {} .po \; | sed 's/.Tumor.cns//')

for name in $allSamples; do
  echo $name
  cnvkit.py heatmap ${CNVKit_folder}/${name}.Tumor.cns -d -o ${CNVKit_folder}/${name}.Tumor.Heatmap.d.pdf

  for chr in {1..22} X Y; do
          (
          cnvkit.py heatmap ${CNVKit_folder}/${name}.Tumor.cns -d -c $chr -o "${CNVKit_folder}/Chromosomes_Heatmaps/${name}.Tumor.Heatmap.d.${chr}.pdf"
          ) &
  done
done

for chr in {1..22} X Y; do
        (
        cnvkit.py heatmap ${CNVKit_folder}/*.Tumor.cns -d -c ${chr} -o "${CNVKit_folder}/Combined_Heatmaps/hMANEC.Heatmap.d.${chr}.pdf"
        ) &
done

cnvkit.py heatmap ${CNVKit_folder}/*.Tumor.cns -d -o "${CNVKit_folder}/Combined_Heatmaps/hMANEC.Heatmap.d.genome.pdf"







# SPECIFIC LIST AND OUTPUT FOLDER
allSamples=$(cat samples_37.csv | tr '\n' ' ')
cd $CNVKit_folder

allCNS=$(printf '%s\n' "$allSamples" | sed 's/[^[:space:]]\{1,\}/&.Tumor.cns/g')
outputFolder=../CNVKit_Heatmap_37Samples/

allCNS=$(printf '%s\n' "$allSamples" | sed 's/[^[:space:]]\{1,\}/&.Tumor.center-mode.filtered.cns/g')
outputFolder=../CNVKit_Heatmap_37Samples_mode-filtered/

cnvkit.py heatmap ${allCNS} -d -o "${outputFolder}/hMANEC.Heatmap.d.genome.pdf"

for chr in {1..22} X Y; do
        (
        cnvkit.py heatmap ${allCNS} -d -c ${chr} -o "${outputFolder}/hMANEC.Heatmap.d.${chr}.pdf"
        ) &
done












# # PREPARE GISTIC SEG FILE
# setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
# Files=list.files(path = "CNVKit_refineSegments_1e-3/", pattern = "\\.Tumor.center-mode.cns$", full.names = T)
#
# allDT <- data.table()
# for(file in Files){
#   print(file)
#   dt <- fread(file, select = c("chromosome", "start", "end", "probes", "log2"))
#   dt <- dt[chromosome %in% c(1:22,"X","Y")]
#   dt[, name := gsub(".Tumor.*.cns", "", basename(file))]
#
#   allDT <- rbind(allDT, dt)
# }
#
# setcolorder(allDT, c("name", "chromosome", "start", "end", "probes", "log2"))
#
# allDT <- allDT[, .(Sample=name, Chromosome=chromosome, "Start Position"=start, "End Position"=end, "Num markers"=probes, Seq.CN=log2)]
#
# # keep samples
# patientfiles <- "/home/rad/Downloads/hMANEC_tmpDB/MAFtools/Oncoplot-final_patient order.txt"
# pats <- fread(patientfiles)
# pats <- pats[!V1 %in% c("", "Patient"), V1]
#
# allDT <- allDT[Sample %in% pats]
#
# fwrite(allDT, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/all_GISTIC.seg", sep="\t", col.names = F, quote = F)








# THIS HAS TO BE RUN ON WS3
GISTIC_base=/home/rad/packages/GISTIC-2.0.23/
GISTIC=${GISTIC_base}/./gistic2
refgene=${GISTIC_base}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
working_directory=/media/rad/HDD2/niklas_temp/GISTIC_TEST/

mkdir -p ${working_directory}/GISTIC/GISTIC_armpeel
mkdir -p ${working_directory}/GISTIC/GISTIC
mkdir -p ${working_directory}/GISTIC/GISTIC_brlen

#segfile=${working_directory}/GISTIC/all_GISTIC.seg
segfile=${working_directory}/GISTIC/2_filtered_GISTIC.seg
segfile=${working_directory}/GISTIC/3_filtered_ClinVar_GISTIC.seg
segfile=${working_directory}/GISTIC/filtered_37Samples_GISTIC.seg

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
  -b ${working_directory}/GISTIC/GISTIC_brlen/ \
  -run_broad_analysis 1 \
  -do_gene_gistic 1 \
  -twoside 1 -ta 0.2 -td 0.2 -brlen 0.9

echo "DONE"





#
# $GISTIC \
#   -refgene $refgene \
#   -seg ${working_directory}/GISTIC/all_GISTIC_filterTEST.seg \
#   -b ${working_directory}/GISTIC/GISTIC_filterTEST/ \
#   -run_broad_analysis 1 \
#   -do_gene_gistic 1 \
#   -twoside 1 -ta 0.2 -td 0.2 -brlen 0.9


mkdir -p ${working_directory}/GISTIC/GISTIC_relaxedfilter
$GISTIC \
  -refgene $refgene \
  -seg $segfile \
  -b ${working_directory}/GISTIC/GISTIC_relaxedfilter/ \
  -run_broad_analysis 1 \
  -do_gene_gistic 1 \
  -twoside 1 \
  -ta 0.1 \
  -td 0.1




#
#
# mv /media/rad/SSD1/GISTIC/GISTIC/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC/
#
# mv /media/rad/SSD1/GISTIC/GISTIC_armpeel/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC_armpeel/
#
# mv /media/rad/SSD1/GISTIC/GISTIC_brlen/* /run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=public/hMANEC/results_selected/GISTIC/GISTIC_brlen/










# COPY AND PROCESS FILES IN R

# #MAFfiles is taken from PrepareStudy.R (based on dbSheet)
#
# m2 <- copy(MAFfiles)
# m2[, m2 := gsub(".vep.maf", ".NoCommonSNPs.OnlyImpact.txt", maf)]
# for(file in m2$m2){
#   print(file)
#   outfile <- paste0("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/Mutect2/", basename(file))
#   file.copy(from = file, to = outfile)
# }
#
#
# stats <- copy(MAFfiles)
# stats[, basefile := gsub(paste0("Mutect2.*"), "", file)]
# stats[, qcfile := paste0(basefile, "/QC/", SampleID, "_data/multiqc_general_stats.txt")]
#
# resDT <- data.table()
# for(i in 1:nrow(stats)){
#   file <- stats[i, qcfile]
#
#   if(!file.exists(file)){
#     samp <- stats[i, SampleID]
#     add <- data.table(Sample=paste0(samp, ".Tumor"), aligned=NA, insertsizeMedian=NA, coverageMedian=NA, coverageMean=NA,bases30X=NA)
#   } else {
#     dt <- fread(file)
#
#     dt <- dt[!grep("\\.R[1\\|2]", Sample)]
#
#     add <- dt[, .(Sample,
#                   aligned=`Picard_mqc-generalstats-picard-PCT_PF_READS_ALIGNED`,
#                   insertsizeMedian=`Picard_mqc-generalstats-picard-summed_median`,
#                   coverageMedian=`Picard_mqc-generalstats-picard-MEDIAN_COVERAGE`,
#                   coverageMean=`Picard_mqc-generalstats-picard-MEAN_COVERAGE`,
#                   bases30X=`Picard_mqc-generalstats-picard-PCT_30X`)]
#   }
#
#   resDT <- rbind(resDT, add)
# }
# resDT <- resDT[grepl("Tumor", Sample) | grepl("Normal", Sample)]
# resDT[, Tumor_Sample_Barcode := gsub("(.*)\\..*", "\\1", Sample)]
# fwrite(resDT, "/home/rad/Downloads/hMANEC_tmpDB/MAFtools/samples_fastqc_statistics.tsv", sep="\t", col.names = T, quote = F)
