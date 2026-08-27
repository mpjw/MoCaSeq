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




# Run CNVKit for multisample for WES or WGS
# the commands for WGS/WES are identical, except for: -method, -access, -targets
AccessBed=${temp_dir}/access-CNVKit.bed
if [ $CNVKit = 'yes' ] && [ $runmode = "MS" ] && [ $sequencing_type = 'WES' ]; then
  echo '---- Running CNVKit (matched tumor-normal, WES) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

  cnvkit.py access $genome_file -o ${AccessBed}

  cnvkit.py batch \
    $bam_tumor \
    --normal $bam_normal \
    --fasta "$genome_file" \
    --output-reference ${CNVKit_folder}/Refernce.cnn \
    --output-dir ${CNVKit_folder} \
    --short-names \
    --diagram \
    --scatter \
    --annotate "$RefFlat" \
    --access "$AccessBed" \
    --targets "$exons_file" \
    --drop-low-coverage \
    -m hybrid \
    -p "$threads"

    # maybe improve this some day:
    #  --antitarget-avg-size 50000 --> CNVkit uses a cautious default off-target bin size that, in our experience, will typically include more reads than the average on-target bin. However, we encourage the user to examine the coverage statistics reported by CNVkit and specify a properly calculated off-target bin size for their samples in order to maximize copy number information.

elif [ $CNVKit = 'yes' ] && [ $runmode = "MS" ] && [ $sequencing_type = 'WGS' ]; then

  echo '---- Running CNVKit (matched tumor-normal, WGS) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

  #To speed up and/or improve the accuracy of WGS analyses: instead of analyzing the whole genome, use the “target” BED file to limit the analysis to just the genic regions. (like described in the CNVKit vignette)
  cnvkit.py batch \
    $bam_tumor \
    --normal $bam_normal \
    --fasta "$genome_file" \
    --output-reference ${CNVKit_folder}/Refernce.cnn \
    --output-dir ${CNVKit_folder} \
    --short-names \
    --diagram \
    --scatter \
    --annotate "$RefFlat" \
    --access "" \
    --targets "$genecode_file_genes_bed" \
    --drop-low-coverage \
    -m wgs \
    -p "$threads"
fi


# Run CNVKit for single sample (and always after MS) for WES or WGS
if [ $CNVKit = 'yes' ]; then
	echo '---- Running CNVKit (single-sample) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

  #for single sample the normal is "empty"
	for type in $types;
	do

    if [ $sequencing_type = 'WES' ]; then
      cnvkit.py access $genome_file -o ${AccessBed}

      cnvkit.py batch \
        $name/results/bam/$name.$type.bam \
        --normal "" \
        --fasta "$genome_file" \
        --output-reference ${CNVKit_folder}/Refernce.cnn \
        --output-dir ${CNVKit_folder} \
        --short-names \
        --diagram \
        --scatter \
        --annotate "$RefFlat" \
        --access "$AccessBed" \
        --targets "$exons_file" \
        --drop-low-coverage \
        -m hybrid \
        -p "$threads"

    elif [ $sequencing_type = 'WGS' ]; then

      cnvkit.py batch \
        $name/results/bam/$name.$type.bam \
        --normal "" \
        --fasta "$genome_file" \
        --output-reference ${CNVKit_folder}/Refernce.cnn \
        --output-dir ${CNVKit_folder} \
        --short-names \
        --diagram \
        --scatter \
        --annotate "$RefFlat" \
        --access "" \
        --targets "$genecode_file_genes_bed" \
        --drop-low-coverage \
        -m wgs \
        -p "$threads"
    fi
	done
fi
































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
