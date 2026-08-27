suppressPackageStartupMessages(library("data.table"))
  
args <- commandArgs(TRUE)

name = args[1]
samplesheet.file = args[2]
manifest.file = args[3]
#samplesheet.file <- "/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD2/GDC_download/TCGA_PAAD_BAMS/gdc_sample_sheet.2021-11-15.tsv "
#name="TCGA-S4-A8RO"
#manifest.file <- "gdc_manifest.2021-11-11.txt"

  samplesheet <- fread(samplesheet.file)
manifest <- fread(manifest.file)

# get paired samples for a specific name
sampleinfo <- samplesheet[`Case ID` == name]

if(nrow(sampleinfo) > 2){
  
  # only keep blood
  if(all(c("Blood Derived Normal", "Solid Tissue Normal") %in% sampleinfo[, `Sample Type`])){
    sampleinfo <- sampleinfo[`Sample Type` != "Solid Tissue Normal"]
  }
  
  # only keep primary
  if(all(c("Primary Tumor", "Metastatic") %in% sampleinfo[, `Sample Type`])){
    sampleinfo <- sampleinfo[`Sample Type` != "Metastatic"]
  }
  
}

# generate manifest
submanifest <- manifest[id %in% sampleinfo$`File ID`]

outfile <- paste0(name, "/results/QC/", name, "_download_manifest.txt")
fwrite(submanifest, outfile, sep="\t", col.names = T)

