# PREPARE GISTIC SEG FILE

# 1. selected patients, modified CNVKit (1e-3), no region filtering
setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
Files=list.files(path = "CNVKit_refineSegments_1e-3/", pattern = "\\.Tumor.center-mode.cns$", full.names = T)

allDT <- data.table()
for(file in Files){
  print(file)
  dt <- fread(file, select = c("chromosome", "start", "end", "probes", "log2"))
  dt <- dt[chromosome %in% c(1:22,"X","Y")]
  dt[, name := gsub(".Tumor.*.cns", "", basename(file))]
  
  allDT <- rbind(allDT, dt)
}

setcolorder(allDT, c("name", "chromosome", "start", "end", "probes", "log2"))

allDT <- allDT[, .(Sample=name, Chromosome=chromosome, "Start Position"=start, "End Position"=end, "Num markers"=probes, Seq.CN=log2)]

# keep samples
patientfiles <- "/home/rad/Downloads/hMANEC_tmpDB/MAFtools/Oncoplot-final_patient order.txt"
pats <- fread(patientfiles)
pats <- pats[!V1 %in% c("", "Patient"), V1]

allDT <- allDT[Sample %in% pats]

fwrite(allDT, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/1_all_GISTIC.seg", sep="\t", col.names = F, quote = F)




# 2. selected patients, modified CNVKit (1e-3), region filtering
setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
Files=list.files(path = "CNVKit_refineSegments_1e-3/", pattern = "\\.Tumor.center-mode.filtered.cns$", full.names = T)

allDT <- data.table()
for(file in Files){
  print(file)
  dt <- fread(file, select = c("chromosome", "start", "end", "probes", "log2"))
  dt <- dt[chromosome %in% c(1:22,"X","Y")]
  dt[, name := gsub(".Tumor.*.cns", "", basename(file))]
  
  allDT <- rbind(allDT, dt)
}

setcolorder(allDT, c("name", "chromosome", "start", "end", "probes", "log2"))

allDT <- allDT[, .(Sample=name, Chromosome=chromosome, "Start Position"=start, "End Position"=end, "Num markers"=probes, Seq.CN=log2)]

# keep samples
patientfiles <- "/home/rad/Downloads/hMANEC_tmpDB/MAFtools/Oncoplot-final_patient order.txt"
pats <- fread(patientfiles)
pats <- pats[!V1 %in% c("", "Patient"), V1]

allDT <- allDT[Sample %in% pats]

fwrite(allDT, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/2_filtered_GISTIC.seg", sep="\t", col.names = F, quote = F)



# NOT ENOUGH DATA,. GISTIC WIL BREAK  
# # 3. selected patients, modified CNVKit (1e-3), region filtering + ClinVar filtering
# germCNV <- fread("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/ClinVar_CNVs.txt")
# germCNV <- germCNV[clinSign == "Benign"]
# germCNV <- germCNV[, .(chr=`#chrom`, start=chromStart, end=chromEnd)]
# germCNV[, class := "ClinGenCNV"]
# germCNV[, chr := gsub("chr","", chr)]
# germCNV.gr <- makeGRangesFromDataFrame(germCNV[, .(chr, start, end, excludeID=.I)])
# 
# setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
# Files=list.files(path = "CNVKit_refineSegments_1e-3/", pattern = "\\.Tumor.center-mode.filtered.cns$", full.names = T)
# 
# allDT <- data.table()
# for(file in Files){
#   print(file)
#   dt <- fread(file, select = c("chromosome", "start", "end", "probes", "log2"))
#   dt <- dt[chromosome %in% c(1:22,"X","Y")]
#   dt[, name := gsub(".Tumor.*.cns", "", basename(file))]
#   
#   gr <- makeGRangesFromDataFrame(dt)
#   hits <- findOverlaps(gr, germCNV.gr)
#   dt[unique(queryHits(hits)), hitClinVarCNV := "yes"]
#   dt <- dt[is.na(hitClinVarCNV)]
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
# fwrite(allDT, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/3_filtered_ClinVar_GISTIC.seg", sep="\t", col.names = F, quote = F)
# 
# 













# 
# 
# keepSamples1 <- cbioAnnoDT[histological_classification_simple == "NEC", Tumor_Sample_Barcode]
# tmpDT1 <- allDT[Sample %in% keepSamples1]
# fwrite(tmpDT1, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/NEC_GISTIC.seg", sep="\t", col.names = F, quote = F)
# 
# keepSamples2 <- cbioAnnoDT[histological_classification_simple == "Intestinal_type", Tumor_Sample_Barcode]
# tmpDT2 <- allDT[Sample %in% keepSamples2]
# fwrite(tmpDT2, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/Intestinal_GISTIC.seg", sep="\t", col.names = F, quote = F)
# 
# keepSamples3 <- cbioAnnoDT[histological_classification_simple == "Diffuse", Tumor_Sample_Barcode]
# tmpDT3 <- allDT[Sample %in% keepSamples3]
# fwrite(tmpDT3, "/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD2/niklas_temp/GISTIC_TEST/GISTIC/Diffuse_GISTIC.seg", sep="\t", col.names = F, quote = F)
# 
# 
# cbioAnnoDT$Tumor_Sample_Barcode[!cbioAnnoDT$Tumor_Sample_Barcode %in% allDT[, unique(Sample)]]
# 
# 
# 
# 
# 
# 
# MAFfiles[CallingMethod == "Mutect2", file2 := gsub(".Mutect2.NoCommonSNPs.OnlyImpact.txt", ".Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt", file)]
# MAFfiles[CallingMethod == "Tumor.Mutect2", file2 := file]
# 
# MAFfiles
# 
# ""
# 
# for(file in MAFfiles$file2){
#   if(!file.exists(file)){
#     print("MISSING")
#     print(file)
#   }
#   file.copy(file, paste0("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/Mutect2_AllTumorOnly/", basename(file)))
# }














