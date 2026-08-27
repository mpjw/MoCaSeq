# copy results without sensitive information (i.e. SNPs) to NAS
setwd("/mnt/Human_Encrypted/TCGA-PAAD_Analysis/MoCaSeq")
resultsDir <- "/media/nas/fastq/Studies/AGRad_TCGA-GDC/PAAD/"

#setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/mnt/Human_Encrypted/TCGA-PAAD_Analysis/MoCaSeq")
#resultsDir <- "/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_TCGA-GDC/PAAD/"

sampleDirs <- basename(list.dirs(recursive = F))
sampleDirs <- sampleDirs[!sampleDirs %in% c("raw", "ref", "temp", "TCGA-F2-A7TX_WRONG_usedWGS", "TCGA-3A-A9IJ_noKRAS-BRAF-mutation", "done")]

#sample <- "TCGA-2J-AAB6"
for(sample in sampleDirs){
  print(sample)
  
  samplefolder <- paste0(sample, "/results/")
  
  # LOG
  logfile <- paste0(samplefolder,"/QC/",sample, ".report.txt")
  dir.create(unique(paste0(resultsDir, "/", dirname(logfile))), showWarnings = FALSE, recursive = T)
  file.copy(logfile, paste0(resultsDir, "/", logfile))
  
  # entire folders 1 by 1
  copyFolders <- paste0(samplefolder, c("BubbleTree/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)
  
  copyFolders <- paste0(samplefolder, c("CNVKit/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)
  
  copyFolders <- paste0(samplefolder, c("Copywriter/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)
  
  copyFolders <- paste0(samplefolder, c("HMMCopy/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)
  
  copyFolders <- paste0(samplefolder, c("ABSOLUTE/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)
  
  copyFolders <- paste0(samplefolder, c("FACETS/"))
  file.copy(copyFolders, paste0(resultsDir, "/", dirname(copyFolders), "/"), recursive=TRUE)

  # MUTECT2
  mutfiles <- paste0(sample, c(".Mutect2.txt", ".Mutect2.NoCommonSNPs.txt", 
                               ".Mutect2.NoCommonSNPs.OnlyImpact.txt", ".Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt", 
                               ".Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt"))
  mutfiles <- paste0(samplefolder, "/Mutect2/", mutfiles)
  
  dir.create(unique(paste0(resultsDir, "/", dirname(mutfiles))), showWarnings = FALSE, recursive = T)
  file.copy(mutfiles, paste0(resultsDir, "/", mutfiles))
  
  # LOH
  lohfiles <- list.files(path = paste0(samplefolder, "/LOH"), pattern = ".pdf", full.names = T)
  dir.create(unique(paste0(resultsDir, "/", dirname(lohfiles))), showWarnings = FALSE, recursive = T)
  file.copy(lohfiles, paste0(resultsDir, "/", lohfiles))
  
}


