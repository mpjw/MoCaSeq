#!/usr/bin/Rscript

##########################################################################################
##
## Cohort_GenerateOverlay.R
##
## Generates Overlay for all samples in the current directory.
##
##########################################################################################

library(data.table)
library(openxlsx)
library(ggplot2)
library(tidyr)

# Niklas debugging
if(F){
  # setwd("/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/hPDAC_ProbesV7/")
  # Samples = list.dirs(full.name=F,recursive=F)
  # Samples <- Samples[grepl("hPDAC", Samples)]
  # species="Human"
  # Save_path="temp3/"
  
  setwd("/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_mPDAC_WES_1000CLs/")
  Save_path="/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/hPDAC_ProbesV7/temp3/"
  
  #setwd("/media/nas/fastq/Studies/AGRad_mPDAC_WES_1000CLs/")
  #Save_path="/home/rad/TEMPS/"
  
  
  
  species="Mouse"
  targetsfile <- "210105_CohortOverview_Preliminary.xlsx"
  targets <- data.table(read.xlsx(targetsfile))
  
  # define mice with at least 1 metastasis sample
  samplesCollapsed <- targets[, paste0(SampleLocation_Short, collapse = ","), by=MouseID]
  withMet <- samplesCollapsed[grep("Met", V1), MouseID]
  targets[, MetastasisFound := "no"]
  targets[MouseID %in% withMet, MetastasisFound := "yes"]
  
  # exclude Blood, Cyst and Ascites
  targets <- targets[!SampleLocation_Short %in% c("Blood", "Cyst", "Ascites")]

  # only select a subtype for now
  Samples = list.dirs(full.name=F,recursive=F)
  Samples <- Samples[Samples %in% targets$FileID]
  

  
  
  Method="Copywriter"
  Paths=rep("",length(Samples))

  resolution=20000
  
  SummaryStat <- "Mean"
  SummaryStat <- "Proportion"
  Ylim <- 5
  
  AberrationCutoff=0.25
  ChromomsomesToRemove=c("X","Y")
  format="pdf"
  Suffix=paste0(".",SummaryStat,".",Ylim)
}


#repository_dir=args[1] #location of repository
repository_dir="/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/MoCaSeq/niklas_dev/"
repository_dir="/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/MoCaSeq/niklas_dev/"
source(paste(repository_dir,"/Cohort_GenerateOverlayLibrary_Custom.R",sep=""))     

#Samples = list.dirs(full.name=F,recursive=F)
#Method="HMMCopy"
#species="Mouse"
#Paths=rep("",length(Samples))
#Save_path=""




for (SummaryStat in c("Mean","Proportion"))
{
  for (Ylim in c(2,5))
  {
    RunOverlayAnalysis(Samples=Samples,Paths=Paths,Method=Method,species=species,resolution=20000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove=c("X","Y"),Ylim=Ylim,format="pdf",Suffix=paste0(".",SummaryStat,".",Ylim),Save_path=Save_path)
  }
}

 