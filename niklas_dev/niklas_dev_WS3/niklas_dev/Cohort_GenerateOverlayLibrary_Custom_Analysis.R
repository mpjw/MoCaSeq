# Analysis: combine muliple PPT/Met samples per mouse per genotype and compare the difference

# keepVal = string or vector to define the values to keep in a specific column (hardcoded MouseGenotype_Short for now)
# keyword = keyword to grep for in colnames which will be combined by mean (e.g. PPT or Met)
RunAnalysis.MeanDif <- function(keepVal, keyword, Samples){
  #keepVal <- "PK"
  #keyword <- "Met"
  #keyword <- "PPT"
  
  # filter samples
  keepSamples <- targets[MouseGenotype_Short %in% keepVal & MetastasisFound == "yes", FileID]
  
  #keepSamples <- targets[MouseGenotype_Short %in% keepVal]
  #keepSamples <- keepSamples[grep("LungMet", SampleLocation_Short), FileID]
  
  #keepSamples <- keepSamples[grep("PPT", keepSamples)]
  #keepSamples <- keepSamples[grep("Met", keepSamples)]
  
  
  SamplesSub <- Samples[Samples %in% keepSamples]
  
  if(length(SamplesSub) == 0){
    print("No samples found for selection")
    return(NULL)
  }
  
  # 1) ORIGINAL
  suppressWarnings(rm(OverlayMat))
  BinGenome(species=species,Method=Method,resolution=resolution,ChromomsomesToRemove=ChromomsomesToRemove)
  OMat <- SetOverlayMat2(SamplesSub)
  suppressWarnings(rm(OverlayMat))
  OMat <- fillOverlayMat2(OMat,SamplesSub,Paths=Paths,AberrationCutoff=AberrationCutoff)
  tempsave.OverlayMat <- copy(OMat) # we will use this for the second part
  
  OMat <- combineOverlayMat2(OMat, SummaryStat=SummaryStat)
  OMat <- ReformatOverlayMat2(OMat)
  OM.original <- copy(OMat)
  
  newSuffix <- paste0(Suffix, ".",keepVal,".original")
  PlotOverlay2(OM.original, format=format,Ylim=Ylim,Suffix=newSuffix,Save_path=Save_path)
  #PlotOverlay2(OM.original, Ylim=Ylim)
  
  
  
  # DEBUG
  #newSuffix <- paste0(Suffix, ".",keepVal,".PPTonly")
  #newSuffix <- paste0(Suffix, ".",keepVal,".METonly")
  #PlotOverlay2(OM.original, format=format,Ylim=Ylim,Suffix=newSuffix,Save_path=Save_path)
  
  
  #PPTdata <- data.table(OM.original)
  #PPTdata[PPTdata[, which.min(Del)], ]
  
  #METdata <- data.table(OM.original)
  #METdata[PPTdata[, which.min(Del)], ]
  

  
  
  # 2) NEW (mean of paired PPT/Met samples)
  overlayDT <- data.table(tempsave.OverlayMat)
  
  subtargets <- targets[FileID %in% SamplesSub] # get subset of annotation
  
  # a) combine multiple PPT samples for a single mouse into a single (using mean)
  mouseIDs <- subtargets[, unique(MouseID)]
  for(mouseID in mouseIDs){
    
    sampleIDs <- subtargets[MouseID == mouseID, FileID]
    
    # the number of samples belonging to the keyword group (e.g. PPT, Met)
    samplesFound <- sampleIDs[grep(keyword, sampleIDs)]
    nKeywords <- length(samplesFound)
    if(nKeywords <= 1){
      next()
    }
    
    # GAIN
    sampleColnames <- paste0(samplesFound, "Gain")
    
    # get new column using mean on the found columns
    newCol <- rowMeans(overlayDT[, sampleColnames, with=F])
    
    # add new column and delete the other ones
    overlayDT[, NEW := newCol]
    setnames(overlayDT, "NEW", paste0(mouseID, "Gain"))
    overlayDT[, (sampleColnames) := NULL]
    
    
    # DEL
    sampleColnames <- paste0(samplesFound, "Del")
    
    # get new column using mean on the found columns
    newCol <- rowMeans(overlayDT[, sampleColnames, with=F])
    
    # # check for large differences
    # rs <- rowSums(overlayDT[, sampleColnames, with=F])
    # oldCol <- rs/length(sampleColnames)
    # difCol <- oldCol - newCol
    # highVal1 <- difCol[difCol > 3]
    # highVal2 <- difCol[difCol < -3]
    # if(length(highVal1) != 0 | length(highVal2) != 0){
    #   print(mouseID)
    # }

    
    # add new column and delete the other ones
    overlayDT[, NEW := newCol]
    setnames(overlayDT, "NEW", paste0(mouseID, "Del"))
    overlayDT[, (sampleColnames) := NULL]
  }

  OM.new <- data.frame(overlayDT)
  OM.new <- combineOverlayMat2(OM.new,SummaryStat=SummaryStat)
  OM.new <- ReformatOverlayMat2(OM.new)
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".collapsed", keyword)
  PlotOverlay2(OM.new, format=format,Ylim=Ylim,Suffix=newSuffix,Save_path=Save_path)
  
  # COMPARE DIFFERENCE
  dtOriginal <- data.table(OM.original[, c("Chromosome", "Start", "Stop", "Gain", "Del")])
  setnames(dtOriginal, c("Gain", "Del"), c("GainOriginal", "DelOriginal"))
  
  dtNew <-  data.table(OM.new[, c("Chromosome", "Start", "Stop", "Gain", "Del")])
  setnames(dtNew, c("Gain", "Del"), c("GainNew", "DelNew"))
  
  dt <- merge(dtOriginal, dtNew, by=c("Chromosome", "Start", "Stop"), all=T)
  # Original was samples combined with sum, so 2xPPT 
  # New is the mean of those 2xPPTs, so it is much less than 2x
  # therefore GainOriginal-GainNew show the amount of "too much gain" which is due to 2x. (negative score means the )
  # If one would take this bars and sum it with the mean bars, it would be the 2x plot
  
  # how can this be MORE after mean?   1    8580001    8600000
  # --> this is due to the formula: sum/length (taking away 1 sample but not marginally increasing the value to to equal mean, will increase/decrease the overall value)
  
  dt[, GainDif := GainOriginal-GainNew]
  dt[, DelDif := DelOriginal-DelNew]
  
  
  # # DEBUG 
  # dt[Chromosome == 5 & Start == 698100001]
  # dt[Chromosome == 4 & Start == 540760001]
  # # Chromosome 4, DelOriginal=-0.35 und DelNew=-0.18 --> damit DelDif = -0.17
  # # Chromosome 5, DelOriginal=-0.05 und DelNew=-0.09 --> damit DelDif = 0.03 (weil DelNew durch meine Korrektur größer wurde)
  # 
  # dtOriginal2 <- data.table(OM.original)
  # dtNew2 <-  data.table(OM.new)
  # tmp1 <- dtOriginal2[Chromosome == 5 & Start == 698100001, -c("Chromosome", "Start", "Stop", "Gain", "Del")]
  # tmp2 <- dtNew2[Chromosome == 5 & Start == 698100001, -c("Chromosome", "Start", "Stop", "Gain", "Del")]
  # # beides rowsums = 2 (also 2 Proben mit dem hit)
  # # aber Original hat 34 Proben und New nur 22 (die anderen wurden gemerged), damit 2/34=0.05(*-1) und 2/22=0.09(*-1) 
  
  
  DifMat <- data.frame(dt[, .(Chromosome, Start, Stop, Gain=GainDif, Del=DelDif)])
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".collapsed", keyword,".MergeDif")
  PlotOverlay2(DifMat, Ylim=Ylim, format="pdf", Suffix=newSuffix,Save_path=Save_path)
  #PlotOverlay2(DifMat, Ylim=Ylim)

  Ymax <- max(dt$GainDif, abs(dt$DelDif))
  Ymax <- round(Ymax, digits=2)+0.05
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".collapsed", keyword,".MergeDifZoom")
  PlotOverlay2(DifMat, Ylim=Ymax, format="pdf", Suffix=newSuffix,Save_path=Save_path, YTicks = 0.02)
  #PlotOverlay2(Ylim=Ymax, YTicks = 0.02)
}

genotypes <- targets[, unique(MouseGenotype_Short)]

for(genotype in genotypes){
  #RunAnalysis.MeanDif(genotypes[1], "PPT")
  #RunAnalysis.MeanDif(genotypes[1], "Met")
  print(genotype)  
  RunAnalysis.MeanDif(genotype, "PPT", Samples)
  RunAnalysis.MeanDif(genotype, "Met", Samples)
}









# Analysis: Only differences plotted
# All Metastasen Proben für z.B. die PK Kohorte. Dabei soll mit dem Mean für multiple metastasen Zelllinien einer Maus korrigiert werden 
# und nicht jede Zelllinie eigenständig gewertet werden, da sonst Mäuse mit vielen Zelllinien das Bild biasen.
# Alle PPT Samples für z.B. die PK Kohorte. Gleich wie die Metastasen, diesmal eben halt nur mit den Primärtumoren.
# Die Differenz der CNV events des PPT gegen seine gematchte Metastase und dann alle Mäuse mergen und in einem Plot darstellen.
MergeByKeyword <- function(overlayDT, keyword, mouseIDs,subtargets){
  for(mouseID in mouseIDs){
    
    #overlayDT <- copy(tmpDT)
    
    sampleIDs <- subtargets[MouseID == mouseID, FileID]

    # the number of samples belonging to the keyword group (e.g. PPT, Met)
    samplesFound <- sampleIDs[grep(keyword, sampleIDs)]
    
    # GAIN
    sampleColnames <- paste0(samplesFound, "Gain")
    
    # get new column using mean on the found columns
    newCol <- rowMeans(overlayDT[, sampleColnames, with=F])
    
    # add new column and delete the other ones
    overlayDT[, NEW := newCol]
    setnames(overlayDT, "NEW", paste(mouseID, keyword, "Gain", sep="."))
    
    overlayDT[, (sampleColnames) := NULL]
    
    # DEL
    sampleColnames <- paste0(samplesFound,"Del")
    
    # get new column using mean on the found columns
    newCol <- rowMeans(overlayDT[, sampleColnames, with=F])
    
    # add new column and delete the other ones
    overlayDT[, NEW := newCol]
    setnames(overlayDT, "NEW", paste(mouseID, keyword, "Del", sep="."))
    overlayDT[, (sampleColnames) := NULL]
  }
  return(overlayDT)
}
RunAnalysis.DifEvents <- function(keepVal, Samples){
  
  #keep <- targets[MouseGenotype_Short == "PK" & MetastasisFound == "yes", FileID]
  keep <- targets[MouseGenotype_Short == keepVal & MetastasisFound == "yes", FileID]
  SamplesSub <- Samples[Samples %in% keep]
  
  subtargets <- targets[FileID %in% SamplesSub] # get subset of annotation
  mouseIDs <- subtargets[, unique(MouseID)]

  if(nrow(subtargets) == 0){
    return(NULL)
  }
  
  # 1) ORIGINAL
  BinGenome(species=species,Method=Method,resolution=resolution,ChromomsomesToRemove=ChromomsomesToRemove)
  OMat <- SetOverlayMat2(SamplesSub)
  suppressWarnings(rm(OverlayMat))
  OMat <- fillOverlayMat2(OMat,SamplesSub,Paths=Paths,AberrationCutoff=AberrationCutoff)
  
  OMat2 <- combineOverlayMat2(OMat, SummaryStat=SummaryStat)
  OMat2 <- ReformatOverlayMat2(OMat2)

  #PlotOverlay2(OMat2, Ylim=0.4, YTicks = 0.1)

  # 2) combine all PPT and Met samples into 1 column (for gain and del independently)
  tmpDT <- data.table(OMat)

  mergedDT <- MergeByKeyword(tmpDT, "Met", mouseIDs, subtargets)
  mergedDT <- MergeByKeyword(mergedDT, "PPT", mouseIDs, subtargets)
  
  # save the raw data
  outdataDT <- copy(mergedDT)
  
  tmp <- data.table(colnames(outdataDT[, -1:-3]))
  tmp <- cSplit(tmp, "V1", ".", drop = F)
  tmp[, newcolnames := paste(V1_1, V1_3, V1_2, sep=".")]
  setorder(tmp, newcolnames)
  setcolorder(outdataDT, tmp$V1)
  setcolorder(outdataDT, c("Chromosome", "Start", "Stop"))
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"), ".tsv")
  outfile <- paste0(Save_path, "CopyNumberDataMatrix1",newSuffix)
  fwrite(outdataDT, outfile, col.names = T, sep = "\t")

  neworder <- c(colnames(outdataDT[, 1:3]), sort(colnames(outdataDT[, -1:-3])))
  setcolorder(outdataDT, neworder)  
  outfile <- paste0(Save_path, "CopyNumberDataMatrix2",newSuffix)
  fwrite(outdataDT, outfile, col.names = T, sep = "\t")

  mergedDT[, tmpID := paste0("ID-",.I)]
  
  # 3) differential CNA events between PPT and Met
  # get absolute (gain always positive, del always negative) difference between the PPT and Met (combined) sample
  
  absDT <- mergedDT[, .(Chromosome, Start, Stop)]
  relDT <- mergedDT[, .(Chromosome, Start, Stop)]
  
  for(mouseID in mouseIDs){
    subdt <- mergedDT[, grep(mouseID, colnames(mergedDT)), with=F]
    dtCols <- colnames(subdt)
    
    # GAIN
    sampleColnames <- dtCols[grep("Gain", dtCols)]
    metCol <- sampleColnames[grep("Met", sampleColnames)]
    pptCol <- sampleColnames[grep("PPT", sampleColnames)]
    
    gaindif <- subdt[, pptCol, with=F] - subdt[, metCol, with=F][[1]]
    
    relDT[, GainDif := gaindif] # add to main result dt
    setnames(relDT, "GainDif", paste("S",mouseID, "Gain", sep="."))
    
    gaindifAbs <- abs(gaindif) # absolute for event plot 
    absDT[, GainDif := gaindifAbs] # add to main result dt
    setnames(absDT, "GainDif", paste("S",mouseID, "Gain", sep="."))
    
    # Del
    sampleColnames <- dtCols[grep("Del", dtCols)]
    metCol <- sampleColnames[grep("Met", sampleColnames)]
    pptCol <- sampleColnames[grep("PPT", sampleColnames)]
    
    deldif <- abs(subdt[, pptCol, with=F]) - abs(subdt[, metCol, with=F])[[1]]
    
    relDT[, DelDif := deldif] # add to main result dt
    setnames(relDT, "DelDif", paste("S",mouseID, "Del", sep="."))
    
    deldifAbs <- abs(deldif)*(-1)
    absDT[, DelDif := deldifAbs] # add to main result dt
    setnames(absDT, "DelDif", paste("S",mouseID, "Del", sep="."))
    
    # README
    # positive number: larger gain/del in PPT (larger PPT - smaller met = positive value)
    # negative number: larger gain/del in met (smaller PPT - larger met = negative value)
    
  }
  
  baseDF <- data.frame(absDT)
  difMat <- combineOverlayMat2(baseDF, SummaryStat=SummaryStat)
  difMat <- ReformatOverlayMat2(difMat)
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".absDifEvents")
  PlotOverlay2(OMatrix=difMat, Ylim=Ylim, format="pdf", Suffix=newSuffix,Save_path=Save_path)
  PlotOverlay2(OMatrix=difMat, Ylim=1.2, format="", Suffix="",Save_path="")
  
  # with zoom
  if(SummaryStat == "Mean"){
    Ymax <- max(difMat$Gain, abs(difMat$Del))
    Ymax <- round(Ymax, digits=2)+0.05
    newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".absDifEventsZoom")
    PlotOverlay2(difMat, Ylim=Ymax, format="pdf", Suffix=newSuffix,Save_path=Save_path, YTicks = 0.2)
  }
  
  if(SummaryStat == "Mean"){
    relDF <- data.frame(relDT)
    relDF <- combineOverlayMat2(relDF, SummaryStat=SummaryStat)
    relDF <- ReformatOverlayMat2(relDF)
  } else if(SummaryStat == "Proportion"){
    relDF <- data.frame(relDT)
    relDF <- combineOverlayMatForProp(relDF, SummaryStat=SummaryStat)
    relDF <- ReformatOverlayMat2(relDF)
  }
  
  # 
  # tmpAbsDT <- data.table(difMat)
  # tmpAbsDT[Chromosome == 11 & Del < -1]
  # 
  # tmpRelDT <- data.table(relDF)
  # tmpRelDT[Chromosome == 11 & Del < -0.2]
  # 
  # 
  # 
  # rowindex2 <- tmpAbsDT[Chromosome == 11 & Start == 1596900001, which=T]
  # tmp1 <- tmpAbsDT[(rowindex2+1)]
  # tmp2 <- tmpRelDT[(rowindex2+1)]
  # 
  # rowMeans(tmp1[, grep("\\.Del", colnames(tmpAbsDT)), with=F])
  # tmp1
  # tmp2
  # rowMeans(tmp2[, grep("\\.Del", colnames(tmpAbsDT)), with=F])
  # 
  # rowMeans(tmp2[, grep("\\.Del", colnames(tmpAbsDT)), with=F][, -2])
  # median(tmp2[, grep("\\.Del", colnames(tmpAbsDT)), with=F])
  # 
  # tmp3 <- tmpRelDT[(rowindex2+1)]
  # myvec <-   tmp3[, grep("\\.Del", colnames(tmpAbsDT)), with=F]
  # rowMeans(myvec)
  # rowMeans(abs(myvec))
  # 
  # mean(c(2,1,-3,-4))
  # mean(2,2)
  # mean(-3,-3) 
  # # erwartung abs = 5
  # # aber mean ist = 2.5
  # mean(abs(c(2,1,-3,-4))) 
  # 
  # posvals <- myvec[, myvec > 0, with=F]
  # negvals <- myvec[, myvec < 0, with=F]
  # zerovals <- myvec[, myvec == 0, with=F]
  # rowMeans(posvals) - abs(rowMeans(negvals))
  # mean(c(-0.6112,0,0,0,0,0,0))
  # 
  # 
  # 
  # relDT[Chromosome == 11 & S.1712.Del == 2.4094]
  #   rowindex <- relDT[Chromosome == 11 & Start == 71200001, which=T]
  # relDT[(rowindex-1):(rowindex+1)]
  # absDT[(rowindex-1):(rowindex+1)]
  # 
  
  
  
  plotDT <- data.table(relDF)
  plotDT <- plotDT[, .(Chromosome, Start, Stop, Gain, Del)]
  
  # get chromosome axis labels 
  chromends <- plotDT[, max(Stop), by=Chromosome]$V1
  chromends <- c(0, chromends)
  chromnames <- plotDT[, paste0("chr", unique(Chromosome))]
  chromLabelPos <- diff(chromends)/2
  names(chromLabelPos) <- chromnames
  chromLabelPos <- chromends+c(chromLabelPos, 0) # with pseuso 0 after last chromosome

  # convert wide to long for ggplot facets  
  plotDT <- data.table(gather(plotDT, "CNAevent", "CNAval", Gain:Del))
  plotDT[CNAevent == "Gain", CNAcolor := "#4682B4"]
  plotDT[CNAevent == "Del", CNAcolor := "pink"]
  
  Ymax <- max(abs(plotDT$CNAval))
  Ymax <- round(Ymax, digits=2)+0.05
  newYLim <- c(-Ymax,Ymax)
  
  p1 <- ggplot(plotDT, aes(Start, CNAval, color=CNAevent)) + 
    geom_vline(xintercept = chromends, linetype=3) +
    geom_hline(yintercept = 0, linetype=1) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Differential CNA event (positive = PPT-specific, negative = Met-specific)") +
    scale_color_manual(values=c("#4682B4", "#B22222")) +
    scale_x_continuous(breaks = chromLabelPos, labels = names(chromLabelPos)) +
    facet_wrap(~CNAevent, ncol = 1) +
    geom_line() + 
    coord_cartesian(ylim=newYLim) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype=2, color="darkgrey")

  p2 <- ggplot(plotDT, aes(Start, CNAval, color=CNAevent)) + 
    geom_vline(xintercept = chromends, linetype=3) +
    geom_hline(yintercept = 0, linetype=1) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Differential CNA event (positive = PPT-specific, negative = Met-specific)") +
    scale_color_manual(values=c("#4682B4", "#B22222")) +
    scale_x_continuous(breaks = chromLabelPos, labels = names(chromLabelPos)) +
    geom_line() + 
    coord_cartesian(ylim=newYLim) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype=2, color="darkgrey")
  
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".relDifEventsZoom.pdf")
  outfile <- paste0(Save_path, "CopyNumberOverlay",newSuffix)
  ggsave(outfile, p1, device = "pdf", width = 16, height = 9)
  
  newSuffix <- paste0(Suffix, ".",paste0(keepVal, collapse = "_"),".relDifEventsZoom2.pdf")
  outfile <- paste0(Save_path, "CopyNumberOverlay",newSuffix)
  ggsave(outfile, p2, device = "pdf", width = 16, height = 9)
}



SummaryStat <- "Mean"
Suffix=paste0(".",SummaryStat,".",Ylim)
#RunAnalysis.DifEvents("PK", Samples)
for(genotype in genotypes){
  print(genotype)  
  RunAnalysis.DifEvents(genotype, Samples)
}
  

SummaryStat <- "Proportion"
Suffix=paste0(".",SummaryStat,".",Ylim)
#RunAnalysis.DifEvents("PK", Samples)
for(genotype in genotypes){
  print(genotype)  
  RunAnalysis.DifEvents(genotype, Samples)
}
















