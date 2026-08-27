# this should be global some day
mkdir CohortAnalysis/
mkdir CohortAnalysis/MutationSignatures/

CancerSignDir=CohortAnalysis/MutationSignatures/CANCERSIGN
mkdir ${CancerSignDir}

# THIS IS CANCERSIGN

tmpdir="temp_cancersign" # this is a different folder than MoCaSeq temp, since that one could be root permission
mkdir ${tmpdir}

working_directory=/media/rad/HDD1/hPDAC_ProbesV7/
sampletable=${working_directory}/sampletable.txt

cancersignPath=/home/rad/packages/CANCERSIGN/


# iterate over all folders and generate the simple mutation txt file for cancersign
while read line ; do
    set $line
    name=$3
    echo "Converting $name" #name=hPDAC03_LivMet-1

    # # FOR RESCUE COMBINED FILE
    awk '{print $1"\tchr"$4"\t"$5"\t"$6"\t"$7}' ${name}/results/rescued/${name}_Mutect2_combined.tsv > ${tmpdir}/cancersign_temp.tsv
    sed -e '1s/FileID/sample_id/' -e '1s/chrChrom/chromosome/' -e '1s/GenomicPos/position/' -e '1s/Ref/reference/' -e '1s/Alt/mutated_to/' ${tmpdir}/cancersign_temp.tsv > ${tmpdir}/cancersign_${name}.tsv

    # FOR RAW MUTATIONS
    # awk -v sampleID="${name}" '{print sampleID"\tchr"$1"\t"$2"\t"$3"\t"$4}' ${name}/results/Mutect2/${name}.Mutect2.txt > ${tmpdir}/cancersign_temp.tsv
    # sed -e "1s/${name}/sample_id/" -e '1s/chrCHROM/chromosome/' -e '1s/POS/position/' -e '1s/REF/reference/' -e '1s/ALT/mutated_to/' ${tmpdir}/cancersign_temp.tsv > ${tmpdir}/cancersign_${name}.tsv
    # rm ${tmpdir}/cancersign_temp.tsv

done < ${sampletable}

# combine all sample files to a single file
rm -f ${tmpdir}/cancersign_input.tsv # remove existing file
cat ${tmpdir}/cancersign_* > ${tmpdir}/cancersign_input.tsv



# Run CancerSign
# A) mutation signature indentification

# generate config file
rm -f ${tmpdir}/cancersign_infer_config.txt # remove existing file
echo "input_file = ${tmpdir}/cancersign_input.tsv" >> ${tmpdir}/cancersign_infer_config.txt
echo "output_dir = CohortAnalysis/MutationSignatures/CANCERSIGN/" >> ${tmpdir}/cancersign_infer_config.txt
echo "infer_3mer_signatures = yes" >> ${tmpdir}/cancersign_infer_config.txt

# start
bash ${cancersignPath}/cancersign --config ${tmpdir}/cancersign_infer_config.txt


# B) plotting
Rscript ${cancersignPath}/src/Plot3merSignatures.R ${CancerSignDir}/infered_3mer_signatures/ 4

# C) clustering
#bash ${cancersignExec} --config config_cluster.txt






rm -r ${tmpdir}
