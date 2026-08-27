#!/usr/bin/Rscript

##########################################################################################
##
## Cohort_GenerateOverlay.R
##
## Generates Overlay for all samples in the current directory.
##
##########################################################################################



# set gene lookup
library(GenomicRanges)
library(splitstackshape)
genesLookupDT <- fread("/run/user/1000/gvfs/sftp:host=imows4.med.tum.de,user=rad/media/rad/SSD1/Genome_References/hsa/GRCh38/gencode/Lookups/36/hsa_GRCh38_gencode_v36_genes.tsv")
genesLookupGR <- makeGRangesFromDataFrame(genesLookupDT, keep.extra.columns = T)


#setwd("/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/hMANEC/results_selected/CNVKit/raw/")
setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/CNVKit/")
repository_dir="/home/rad/packages/MoCaSeq/niklas_dev/Cohort_Analysis/"

source(paste(repository_dir,"/Cohort_GenerateOverlayLibrary_CNVKit.R",sep=""))     

library(data.table)
library(openxlsx)
library(binr)
# for subsets
targets <- fread("/home/rad/Downloads/hMANEC_tmpDB/hMANEC_Combined_AllSelected/clinical_sample_data.txt", header = T)
targets <- targets[-1:-4]
colnames(targets) <- gsub(" ", "", colnames(targets))
targets[, Histological_classification_simple := gsub(" ", "_", Histological_classification_simple)]
targets[, Histological_classification := gsub(" ", "_", Histological_classification)]
targets[, Histological_classification := gsub("\\(", "", Histological_classification)]
targets[, Histological_classification := gsub("\\)", "", Histological_classification)]
targets


# exclude bad samples, by only keeping specific ones
keepSamples <- gsub(".Tumor.cnr", "", list.files(pattern = ".cnr"))
targets <- targets[SampleIdentifier %in% keepSamples]
Samples <- targets$SampleIdentifier


Method="CNVKit"
species="Human"
Paths=rep("",length(Samples))
Save_path="/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/CNVKit_Overlay/"

resolution=50000
SummaryStat <- "Mean"
chromzoom <- T # will save a wide pdf with more bars
savesegments <- T
groupname <- "All"
cominedName <- groupname

# debug
# AberrationCutoff=0.25
# ChromomsomesToRemove=c("X","Y")
# format="pdf"
# SummaryStat <- "Mean"
# Ylim <- "3"
# Suffix=paste0(".All.",SummaryStat,".",Ylim)


for (SummaryStat in c("Mean")){
  for (Ylim in c(1,3,5)){
    Suffix=paste0(".All.",SummaryStat,".",Ylim)
    RunOverlayAnalysis(Samples=Samples,Paths=Paths,Method=Method,species=species,resolution=50000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=Suffix,Save_path=Save_path)
  }
}

SummaryStat <- "Proportion"
Suffix=paste0(".All.",SummaryStat)
RunOverlayAnalysis(Samples=Samples,Paths=Paths,Method=Method,species=species,resolution=50000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=Suffix,Save_path=Save_path)




# subsets

groupname <- "Histological_classification_simple"
groups <- targets[, unique(Histological_classification_simple)]

groupname <- "Histological_classification"
groups <- targets[, unique(Histological_classification)]

groupname <- "matched_normal"
groups <- targets[, unique(matched_normal)]

for(group in groups){
  print(group)
  #group <- groups[2]
  
  #SubSamples <- targets[Histological_classification_simple == group, SampleIdentifier]
  SubSamples <- targets[matched_normal == group, SampleIdentifier]
  
  cominedName <- paste0(groupname, "-", group)
  
  if(length(SubSamples) <= 1){
    next()
  }
  
  Paths=rep("",length(SubSamples))
  
  chromzoom <- T # will save a wide pdf with more bars
  savesegments <- T
  for (Ylim in c(1,3,5)){
    print(Ylim)
    SummaryStat <- "Mean"
    baseSuffix <- paste(groupname,group,SummaryStat, sep=".")
    mySuffix <- paste0(".",baseSuffix,".",Ylim)

    RunOverlayAnalysis(Samples=SubSamples,Paths=Paths,Method=Method,species=species,resolution=50000,SummaryStat="Mean",AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=mySuffix,Save_path=Save_path)
  }
  
  chromzoom <- F
  savesegments <- F
  SummaryStat <- "Proportion"
  mySuffix <- paste0(".",groupname,".", group, ".",SummaryStat)
  RunOverlayAnalysis(Samples=SubSamples,Paths=Paths,Method=Method,species=species,resolution=50000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=mySuffix,Save_path=Save_path)
  
}



chromzoom <- F 
savesegments <- F

Ylim <- 2
SummaryStat <- "Mean"
Save_path <- "/home/rad/Downloads/temp3/"
mySuffix <- paste0(".all.",Ylim)

AberrationCutoff=0.25
ChromomsomesToRemove=c("X","Y")
format="pdf"
SummaryStat <- "Mean"
Ylim <- "3"
Suffix=paste0(".All.",SummaryStat,".",Ylim)


RunOverlayAnalysis(Samples=Samples,Paths=Paths,Method=Method,species=species,resolution=50000,SummaryStat="Mean",AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=mySuffix,Save_path=Save_path)





