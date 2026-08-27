# THIS IS RESCUE ANNOTATION
working_directory=/media/rad/HDD1/hPDAC_ProbesV7/
ref_directory=/media/rad/SSD1/MoCaSeq_ref/
script_directory=/media/rad/HDD1/MoCaSeq/

sampletable=${working_directory}/sampletable.txt
sampletable=${working_directory}/sampletable_small.txt

nParallel=6

sudo parallel --colsep '\t' -eta -j ${nParallel} "sudo time docker run \
-v ${working_directory}:/var/pipeline/ \
-v ${ref_directory}:/var/pipeline/ref/ \
-v ${script_directory}:/opt/MoCaSeq \
--entrypoint /opt/MoCaSeq/entrypoint_custom.sh \
mocaseq-human2 {3}" :::: ${sampletable}


exit 1

# do not quote this variable! it only works like this
saveLocation=/run/user/1000/gvfs/smb-share\:server\=imostorage.med.tum.de\,share\=fastq/Studies/AGRad_hPDAC/WES/tmp/

# copy specific files to server
while read line ; do
    set $line
    name=$3
    echo ${name}

    #ls ${name}/results/Mutect2/${name}.Mutect2.txt
    #ls ${name}/results/rescued/${name}_Mutect2_combined.tsv
    #ls ${name}/results/VariantAnnotation/${name}_Mutect2_variant_annotation.txt

    cp ${name}/results/Mutect2/${name}.Mutect2.txt ${saveLocation}
    cp ${name}/results/rescued/${name}_Mutect2_combined.tsv ${saveLocation}
    cp ${name}/results/VariantAnnotation/${name}_Mutect2_variant_annotation.txt ${saveLocation}

done < ${sampletable}
