library(data.table)

summarizePurities <- function(name){
  
  ppDT <- data.table()
  
  # 1. ABSOLUTE
  abs.rawfile <- paste0(name, "/results/ABSOLUTE/", name, ".ABSOLUTE.results.tsv")
  
  if(file.exists(abs.rawfile)){
    print("ABSOLUTE: adding to summary")
    abs.raw <- fread(abs.rawfile)
    
    abs.tmp <- abs.raw[, .(sample=name, tool="ABSOLUTE", purity, ploidy, modelPriority=.I, note=paste0("model", .I))]
    ppDT <- rbind(ppDT, abs.tmp)
    
  } else {
    print("ABSOLUTE: no result file found")
  }
  
  
  # 2. FACETS
  fc.rawfile <- paste0(name, "/results/FACETS/", name, ".FACETS.txt")
  
  if(file.exists(fc.rawfile)){
    print("FACETS: adding to summary")
    fc.raw <- fread(fc.rawfile)
    
    
    fc.tmp <- fc.raw[, .(sample=name, tool="FACETS", purity=V2, ploidy=V4, modelPriority=1, note=NA)]
    ppDT <- rbind(ppDT, fc.tmp)
    
  } else {
    print("FACETS: no result file found")
  }
  
  # 3. TITAN
  dirname <- paste0(name, "/results/Titan")
  dirExists <- dir.exists(dirname)
  if(dirExists){
    print("TITAN: adding to summary")
    
    for(runploidy in 2:4){
      #print(runploidy)
      for(cluster in 1:5){
        #print(cluster)
        clustername <- paste0("cluster0", cluster)
        
        ti.rawfile <- paste0(name, "/results/Titan/run_ploidy",runploidy, "/", name, "_", clustername, ".params.txt")
        
        
        if(file.exists(fc.rawfile)){
          
          ti.raw <- fread(ti.rawfile, nrows = 2)
          ti.raw[, tempID := 1]
          ti.raw <- spread(ti.raw, V1, V2)
          names(ti.raw) <- c("tmpID", "ploidy", "normpurity")
          
          noteText <- paste0("ploidy",runploidy, "-", clustername)
          
          ti.tmp <- ti.raw[, .(sample=name, tool="TITAN", purity=1-as.numeric(normpurity), ploidy, modelPriority=NA, note=noteText)]
          ppDT <- rbind(ppDT, ti.tmp)
        }
      }
    }
  } else {
    print("TITAN: no result file found")
  }
  
  
  
  # BUBBLETREE
  bt.rawfile <- paste0(name, "/results/BubbleTree/", name, ".Bubbletree.txt")
  
  if(file.exists(bt.rawfile)){
    print("BubbleTree: adding to summary")
    bt.raw <- fread(bt.rawfile, header = F, sep = ";")
    names(bt.raw) <- c("purity", "ploidy", "deviation")
    bt.raw[, purity := gsub("Purity: ", "", purity)]
    bt.raw[, ploidy := gsub("Ploidy: ", "", ploidy)]
    
    # catch empty values
    bt.raw[purity == "Purity:", purity := NA]
    bt.raw[ploidy == " NA", ploidy := NA]
    
    bt.raw <- data.table(cSplit(bt.raw, "purity", ", ", direction = "long"))
    bt.raw[, note := paste0("subclone", .I)]
    
    # set invalid purities to NA 
    #bt.raw[purity > 1, purity := NA]
    
    bt.tmp <- bt.raw[, .(sample=name, tool="BUBBLETREE", purity, ploidy, modelPriority=1, note)]
    ppDT <- rbind(ppDT, bt.tmp)
  } else {
    print("BubbleTree: no result file found")
  }
  
  # SCLUST
  sc.rawfile <- paste0(name, "/results/Sclust/", name, "_cn_summary.txt")
  
  if(file.exists(sc.rawfile)){
    print("Sclust: adding to summary")
    sc.raw <- fread(sc.rawfile)
    
    sc.tmp <- sc.raw[, .(sample=name, tool="SCLUST", purity, ploidy, modelPriority=1, note=NA)]
    ppDT <- rbind(ppDT, sc.tmp)
  } else {
    print("Sclust: no result file found")
  }
  
  ppDT[, purity := round(as.numeric(purity), digits = 2)]
  ppDT[, ploidy := round(as.numeric(ploidy), digits = 2)]
  
  return(ppDT)
}


setwd("/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_hPDAC/WES/hPDAC_ProbesV7")

#name <- "hPDAC03_LivMet-1"
#ppDT <- summarizePurities(name)

allsamples <- list.dirs(recursive = F)
allsamples <- allsamples[grep("hPDAC", allsamples)]
allsamples <- basename(allsamples)
CohortPurity <- data.table()
for(sample in allsamples){
  print(sample)
  ppDT <- summarizePurities(sample)
  CohortPurity <- rbind(CohortPurity, ppDT)
}


fwrite(CohortPurity, "CohortAnalysis/Cohort_PurityPloidy.tsv", sep="\t", col.names = T, quote = F)

plotDT <- copy(CohortPurity)
plotDT <- plotDT[!is.na(purity)]
plotDT[, barname := paste0(sample, note)]
plotDT[, barname := str_replace_all(barname, "NA", "")]
plotDT[is.na(note), note := modelPriority]

ggplot(plotDT, aes(tool, purity, group=note, fill=tool)) + 
  geom_bar(stat="identity", position="dodge", colour="black", size=0.2) +
  facet_wrap(~sample) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CohortAnalysis/Cohort_Purity.pdf", last_plot(), device = "pdf", width = 16, height = 9)

ggplot(plotDT, aes(tool, ploidy, group=note, fill=tool)) + 
  geom_bar(stat="identity", position="dodge", colour="black", size=0.2) +
  facet_wrap(~sample) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CohortAnalysis/Cohort_Ploidy.pdf", last_plot(), device = "pdf", width = 16, height = 9)





