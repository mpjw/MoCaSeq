library(ggpubr)
library(scales)
library(dplyr)
library(tidyr)
library(data.table)
library(zoo)
library(GenomicRanges)
library(ggplot2)

method <- cnaMethod <- "Copywriter"
segmentsFrom <- "CNA"
cutoff.CNA <- 0.1
smoothCNA <- T
plotQC <- T

LoadCNA <- function(name, method="Copywriter", smoothCNA=T, MinDifToMerge=0.05, MinSegSize1=5000000, MinSegSize2=10000000, MaxDistToNext=10000000, plotQC=T){
  
  # MinSegSize1 is applied before merging. Here we are a bit less stringed than after merging (MinSegSize2)
  # MaxDistToNext defines the size of the "gap to be jumped" while merging neighboring segments
  # MinDifToMerge defines the "up/low shift" between neighboring segments to merge (values below the treshold are regarded equal and merged to 1 single segment using mean)
  
  # load data
  if(method == "Copywriter"){
    seg.file <- paste0(name,"/results/Copywriter/",name,".Copywriter.segments.Mode.txt")
  } else if(method == "HMMCopy"){
    seg.file <- paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt")
  } else {
    stop("other methods not implemented yet")
  }
  
  cnaSegs.raw <- fread(seg.file)
  cnaSegs.raw <- cnaSegs.raw[Chrom %in% 1:22]
  cnaSegs.raw[, cnasegID := .I]
  cnaSegs.raw[, width := round(End-Start)]
  
  if(smoothCNA){
    cnaSegs <- copy(cnaSegs.raw)
    setorder(cnaSegs, Chrom, Start)
    cnaSegs <- cnaSegs[width > MinSegSize1]
    
    # merge similar neighboring segments
    cnaSegs[, DifNextSeg := abs(Mean - data.table::shift(Mean)), by = Chrom]
    cnaSegs[, DistNextSeg := abs(Start - data.table::shift(End)), by = Chrom]
    cnaSegs[, mergeID := as.numeric(.I)]
    cnaSegs[DifNextSeg < MinDifToMerge & DistNextSeg < MaxDistToNext, mergeID := NA]
    cnaSegs[, mergeID := na.locf(mergeID)]
    
    cnaSegs.merged <- unique(cnaSegs[, .(Chrom,Start=min(Start), End=max(End), Mean=mean(Mean), width=(max(End)-min(Start))), by=mergeID])
    cnaSegs.merged[, mergeID := NULL]
    cnaSegs.merged[, width := round(End-Start)]
    cnaSegs.merged <- cnaSegs.merged[width >= MinSegSize2]
    
    if(plotQC){
      pCNAsmoothing <- ggplot() + 
        facet_wrap(~Chrom, scales = "free_x") +
        geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=Mean, yend=Mean), color="darkgreen", size=3, alpha=0.5) +
        geom_segment(data=cnaSegs.merged, aes(x=Start, xend=End, y=Mean, yend=Mean), color="red", size=0.8, linetype=1) +
        geom_segment(data=cnaSegs.raw[!cnasegID %in% cnaSegs$cnasegID], aes(x=Start, xend=End, y=Mean, yend=Mean), color="darkgrey", size=0.8)
    }
    
    cnaSegs <- copy(cnaSegs.merged)
    
  } else {
    cnaSegs <- copy(cnaSegs.raw)
  }
  return(list(cnaSegs=cnaSegs, pCNAsmoothing=pCNAsmoothing))
}
EstimatePurity <- function(name, cnaSegs, lohSegs, outdir, segmentsFrom="CNA", cutoff.CNA=-0.1, plotQC=F, plotSteps=T, cutoff.width=5000000, cutoff.LOH=0.55){
  
  lohSegs[, Chrom := as.numeric(Chrom)]
  cnaSegs[, Chrom := as.numeric(Chrom)]
  
  # filtering cutoffs
  #cutoff.LOH <- 0.55
  #cutoff.width <- 1000000
  #cutoff.width <- 5000000
  
  # apply size filter 
  lohSegs <- lohSegs[width > cutoff.width]
  cnaSegs <- cnaSegs[width > cutoff.width]
  lohSegs.gr <- makeGRangesFromDataFrame(lohSegs)
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs[, .(chr=Chrom, start=Start, end=End)])
  
  # segments either from CNA or from LOH segmentation 
  if(segmentsFrom == "CNA"){
    delSegs <- copy(cnaSegs)
    setnames(delSegs, "Mean", "CNAlog2")
    delSegs[, CNA := 2^CNAlog2]
    delSegs[, segID := .I]
    setcolorder(delSegs, "segID")
    delSegs.gr <- makeGRangesFromDataFrame(delSegs)
    
    hits <- findOverlaps(delSegs.gr, lohSegs.gr)
    # CNA segments without LOH segment are discarded
    delTmp <- cbind(delSegs[queryHits(hits)], lohSegs[subjectHits(hits), .(lohUpSeg=upSeg, start.other=start, end.other=end)]) 
    
    delSegs <- unique(delTmp[delTmp[, .I[lohUpSeg == max(lohUpSeg)], by=segID]$V1])
    #delSegs <- unique(delTmp[, lohUpSeg := max(lohUpSeg), by=segID]) # deprecated
    
  } else if(segmentsFrom == "LOH"){
    
    delSegs <- lohSegs[, .(Chrom, Start=start, End=end, lohUpSeg=upSeg)]
    delSegs[, segID := .I]
    setcolorder(delSegs, "segID")
    delSegs.gr <- makeGRangesFromDataFrame(delSegs) 
    
    hits <- findOverlaps(delSegs.gr, cnaSegs.gr)
    delTmp <- cbind(delSegs[queryHits(hits)], cnaSegs[subjectHits(hits), .(CNAlog2=Mean, start.other=Start, end.other=End)]) # LOH segments without CNA segment are discarded
    delTmp[, CNA := 2^CNAlog2]
    
    delSegs <- unique(delTmp[delTmp[, .I[abs(CNAlog2) == max(abs(CNAlog2))], by=segID]$V1])
  } else {
    stop("wrong segmentsFrom")
  }
  
  # reduce the segment size to the overlapping intersect range (using pseudo chromosome to ensure 1:1 mapping)
  regionIntersect <- as.data.table(GenomicRanges::intersect(makeGRangesFromDataFrame(delSegs[,.(Chrom=.I, Start, End)]), 
                                                            makeGRangesFromDataFrame(delSegs[,.(Chrom=.I, Start=start.other, End=end.other)])))
  
  delSegs <- cbind(delSegs[, -c("Start", "End", "start.other", "end.other")], regionIntersect[, .(Start=start, End=end)])
  delSegs[, width := round(End-Start)]
  
  if(plotQC){
    # these are QC plots, they are independent of the final segments used
    chromDummy <- data.table(Chrom=names(chrom.sizes), End=chrom.sizes)
    chromDummyStarts <- merge(lohSegs[, min(start), by=Chrom], cnaSegs[, min(Start), by=Chrom], by="Chrom", all=T)
    chromDummyStarts[is.na(chromDummyStarts)] <- 10000000000000 # pseudo large value for NA
    chromDummyStarts[V1.x <= V1.y, Start  := V1.x]
    chromDummyStarts[V1.x > V1.y, Start  := V1.y]
    chromDummy <- merge(chromDummy[, .(Chrom=as.numeric(Chrom), End)], chromDummyStarts[, .(Chrom=as.numeric(Chrom), Start)], by="Chrom")
    
    borders <- unique(rbind(lohSegs[, .(Chrom,pos=start)], lohSegs[, .(Chrom,pos=end)], cnaSegs[, .(Chrom,pos=Start)], cnaSegs[, .(Chrom,pos=End)]))
    p4raw <- ggplot() +
      geom_vline(data=borders, aes(xintercept=pos), alpha=0.3) +
      geom_segment(data=lohSegs, aes(x=start, xend=end, y=upSeg, yend=upSeg), color="red", size=0.8) +
      geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=Mean, yend=Mean), color="blue", size=0.8) +
      geom_hline(yintercept = c(0, 0.5), alpha=0.3) +
      facet_wrap(~factor(Chrom, levels=c(1:22)), scales = "free_x") +
      geom_blank(data=chromDummy, aes(x=Start, xend=End))
    
    lohtmp1 <- lohSegs[upSeg >= cutoff.LOH]
    cnatmp1 <- cnaSegs[Mean <= cutoff.CNA]
    borders <- unique(rbind(lohtmp1[, .(Chrom,pos=start)], lohtmp1[, .(Chrom,pos=end)], cnatmp1[, .(Chrom,pos=Start)], cnatmp1[, .(Chrom,pos=End)]))
    p4filtered <- ggplot() +
      geom_vline(data=borders, aes(xintercept=pos), alpha=0.3) +
      geom_segment(data=lohtmp1, aes(x=start, xend=end, y=upSeg, yend=upSeg), color="red", size=0.8) +
      geom_segment(data=cnatmp1, aes(x=Start, xend=End, y=Mean, yend=Mean), color="blue", size=0.8) +
      geom_hline(yintercept = c(0, 0.5), alpha=0.3) +
      facet_wrap(~factor(Chrom, levels=c(1:22)), scales = "free_x") +
      geom_blank(data=chromDummy, aes(x=Start, xend=End))
    
    p4 <- ggarrange(p4raw,p4filtered)
    ggsave(paste0(outdir, "/",name,"_QC1_CNA-LOH-segments.png"), p4, width = 16, height = 9)
  }
  
  
  # FILTERING
  delSegs[, filterstatus := "used"]
  delSegs[lohUpSeg <= cutoff.LOH, filterstatus := "removed"]
  delSegs[CNAlog2 >= cutoff.CNA, filterstatus := "removed"]
  delSegs[width <= cutoff.width, filterstatus := "removed"]
  
  delSegsRaw <- copy(delSegs) # save all for qc plotting
  delSegs <- delSegs[filterstatus == "used"]
  delSegs <- delSegs[, -c("filterstatus")]
  
  if(nrow(delSegs) <= 1){
    print("Not enough LOH segments found. You could try and reduce the cutoffs for LOH and CNA, but probably your purity is <10% anyways.")
    
    parametersStopped <- data.table(name = c("sample", "genomic contamination (stroma content)", "genomic purity (tumor content)", "segments based on",
                                      "cutoff value (minimum) for LOH segments", "cutoff value (minimum) for CNA segments",
                                      "cutoff length (minimum) for LOH/CNA segment widths", "number of segments used for analysis",
                                      "number of chromosomes with segments", "range of CNA values",
                                      "range of LOH values", "range of segment widths"), 
                             value = c(name, 1, 0, segmentsFrom, cutoff.LOH, cutoff.CNA, format(cutoff.width, scientific=F),
                                       0, 0, 0, 0, 0),
                             comment = c("",">80%","<10%","","","", "", "", "", "", "", ""))
    
    fwrite(parametersStopped, paste0(outdir, "/",name,"_results.txt"), sep="\t", col.names = T, quote = F)
    
    return(NULL)
  }
  
  if(plotQC){
    pSegments <- ggplot() + 
      facet_wrap(~Chrom, scales = "free_x") +
      geom_segment(data=delSegs, aes(x=Start, xend=End, y=lohUpSeg, yend=lohUpSeg), color="red", size=1) +
      geom_segment(data=delSegs, aes(x=Start, xend=End, y=CNAlog2, yend=CNAlog2), color="blue", size=1) +
      geom_blank(data=chromDummy[Chrom %in% delSegs$Chrom], aes(x=Start, xend=End)) +
      geom_hline(yintercept = c(0, 0.5), alpha=0.3)
    
    ggsave(paste0(outdir, "/",name,"_QC2_purity-segments.png"), pSegments, width = 16, height = 9)
  }
  
  # init the summary statistics and results output (cont and pur will be replaced later)
  nChroms <- delSegs[, length(unique(Chrom))]
  cnaRange <- paste0(delSegs[, round(range(CNAlog2), digits=2)], collapse = " - ")
  lohRange <- paste0(delSegs[, round(range(lohUpSeg), digits=2)], collapse = " - ")
  widthRange <- paste0(delSegs[, range(width)], collapse = " - ")
  
  parameters <- data.table(name = c("sample", "genomic contamination (stroma content)", "genomic purity (tumor content)", "segments based on",
                                    "cutoff value (minimum) for LOH segments", "cutoff value (minimum) for CNA segments",
                                    "cutoff length (minimum) for LOH/CNA segment widths", "number of segments used for analysis",
                                    "number of chromosomes with segments", "range of CNA values",
                                    "range of LOH values", "range of segment widths"), 
                           value = c(name, as.numeric(NA), as.numeric(NA), segmentsFrom, cutoff.LOH, cutoff.CNA, format(cutoff.width, scientific=F),
                                     nrow(delSegs), nChroms, cnaRange, lohRange, widthRange),
                           comment = c("",">80%","<10%","","","", "", "", "", "", "", ""))
  
  my_binwidth <- 0.05
  
  delSegs[, gCont := (2*CNA) - (2*CNA*lohUpSeg)]
  delSegs[, gContRounded := round(gCont, digits = 2)]
  delSegs[, gPur := 1-gCont]
  delSegs[, gPurRounded := round(gPur, digits = 2)]
  
  values <- delSegs$gPur
  modes <- GetModes2(values,adjust=1,2,min(values)-1,max(values)+1, bandwidth = my_binwidth)
  modeIndex <- modes[, which.max(modes)] # select the smaller contamination value
  gPur <- round(modes[modeIndex, modes], digits = 4)
  gCont <- 1-gPur
  
  # values <- delSegs$gCont
  # modes <- GetModes2(values,adjust=1,2,min(values)-1,max(values)+1, bandwidth = my_binwidth)
  # modeIndex <- modes[, which.min(modes)] # select the smaller contamination value
  # gCont <- round(modes[modeIndex, modes], digits = 4)
  # gPur <- 1-gCont

  pMain <- ggplot(delSegs, aes(gPur)) + 
    geom_histogram(binwidth = my_binwidth, aes(fill=factor(gContRounded))) + 
    geom_density(bw=my_binwidth, color="red", aes(y = 0.1 * ..count..), show.legend = FALSE) +
    scale_fill_viridis_d(limits = factor(seq(0,1,0.01)), breaks=seq(0,1,0.1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), legend.position = "none") +
    guides(fill=guide_legend(title="Genomic purity")) +
    geom_vline(xintercept = gPur, linetype=2) +
    scale_x_continuous(breaks = breaks_pretty(n = 20), limits = c(0.01,0.99)) +
    scale_y_continuous(breaks = breaks_pretty(n = 10)) +
    ylab("N segments") + 
    annotate("label", x=gPur, y=-0.3, label=gPur)
  pMain
  
  ggsave(paste0(outdir, "/",name,"_histogram.png"), pMain, width = 16, height = 9)
  
  # this is optional: redo the analysis with 10% steps using the formula to show the segments behavior 
  if(plotSteps){
    
    steps10 <- data.table()
    valueDensData10 <- data.table()
    valueData10 <- data.table()
    formulaData10 <- data.table()
    outlierData10 <- data.table()
    
    seq10 <- round(seq(0,1,0.1), digits = 2) # we need to round, because float point comparison is a problem for %in% (only for seq, not for c())
    seq10 <- sort(unique(c(seq10, gCont)))
    
    for(ContTheor in seq10){
      values <- delSegs[, (lohUpSeg - ((0.5*(1/CNA))*ContTheor)) / (1 - ((1/CNA)*ContTheor))] # formula applied to each segment
      
      # DEBUG
      formulaData10tmp <- copy(delSegs)
      formulaData10tmp[, gCont := ContTheor]
      formulaData10tmp[, gPur := 1-ContTheor]
      formulaData10tmp[, formulaResult := values]
      formulaData10tmp[, upperTerm := lohUpSeg - ((0.5*(1/CNA))*ContTheor)]
      formulaData10tmp[, lowerTerm := (1 - ((1/CNA)*ContTheor))]
      formulaData10tmp[, lowerTermRight := (1/CNA)*ContTheor]
      setorder(formulaData10tmp, CNAlog2)
      formulaData10 <- rbind(formulaData10, formulaData10tmp)
      
      ContTheor <- as.character(ContTheor)
      outliers <- values[!values %between% c(-2,2)]
      values <- values[values %between% c(-2,2)]
      
      if(length(values) <= 1){
        values <- c(0,0)
      }
      
      # save the outliers into DT
      if(length(outliers) > 0){
        outlierData10Tmp <- data.table(gCont=ContTheor, gPur=as.character(round(1-as.numeric(ContTheor), digits=4)), formulavalue=outliers)
        outlierData10 <- rbind(outlierData10, outlierData10Tmp)
      }
      
      # save the raw formula values into DT
      valueDataTmp <- data.table(gCont=ContTheor, gPur=as.character(round(1-as.numeric(ContTheor), digits=4)), formulavalue=values)
      valueData10 <- rbind(valueData10, valueDataTmp)
      
      # get the density values for plotting and value extraction (mode)
      valueDensDataTmp <- data.table(gCont=ContTheor, gPur=as.character(round(1-as.numeric(ContTheor), digits=4)), X=density(values, bw=my_binwidth)$x, Y=density(values, bw=my_binwidth)$y)
      valueDensData10 <- rbind(valueDensData10, valueDensDataTmp)
      
      modes <- GetModes2(values,adjust=1,2,min(values)-1,max(values)+1, bandwidth = my_binwidth)
      modeIndex <- modes[1:2, which(abs(modes - 1) == min(abs(modes - 1), na.rm=T))]     # from the top 2 modes, select the one closer to 1
      valueMode <- modes[modeIndex, modes]
      valueHeight <- modes[modeIndex, density]
      
      distTo1 <- abs(1-valueMode)
      steps10 <- rbind(steps10, data.table(gCont=ContTheor, gPur=as.character(round(1-as.numeric(ContTheor), digits=4)), Mode=valueMode, modeHeight=valueHeight, ModeDistTo1=distTo1))
    }

    # QC plot
    if(plotQC){
      qcData <- copy(formulaData10)
      qcData <- qcData[formulaResult %between% c(-2,2)]
      setorder(qcData, formulaResult)
      qcData <- data.table(gather(qcData, name, value, c("CNAlog2","lohUpSeg", "width")))
      
      p3a <- ggplot(qcData[name=="CNAlog2"], aes(formulaResult, value)) + 
        geom_point() +
        facet_wrap(~gPur, ncol = 1) +
        coord_cartesian(ylim = c(NA,0.1), xlim = c(-2,2)) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 0) +
        theme(axis.title.x=element_blank()) +
        ylab("CNA [log2]")
      
      p3b <- ggplot(qcData[name=="lohUpSeg"], aes(formulaResult, value)) + 
        geom_point() +
        facet_wrap(~gPur, ncol = 1) +
        coord_cartesian(ylim = c(0,1), xlim = c(-2,2)) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 0) +
        theme(axis.title.x=element_blank()) +
        ylab("LOH")
      
      p3c <- ggplot(qcData[name=="width"], aes(formulaResult, log10(value))) + 
        geom_point() +
        facet_wrap(~gPur, ncol = 1) +
        coord_cartesian(ylim = c(5,NA), xlim = c(-2,2)) +
        geom_vline(xintercept = 1) +
        theme(axis.title.x=element_blank()) +
        ylab("segment size [log10]")
      
      p3 <- ggarrange(p3a, p3b, p3c, ncol = 3)
      
      ggsave(paste0(outdir, "/",name,"_QC3_segment-stats.png"), p3, width = 16, height = 9)
      
      p5a <- ggplot(delSegsRaw, aes(log10(width),CNAlog2, color=lohUpSeg)) +
        geom_point(aes(size=filterstatus), alpha=0.8) +
        geom_vline(xintercept = log10(cutoff.width)) +
        geom_hline(yintercept = cutoff.CNA)
      
      ggsave(paste0(outdir, "/",name,"_QC4a_segment-filtering-v1.png"), p5a, width = 16, height = 9)
      
      p5b <- ggplot(delSegsRaw, aes(log10(width),CNAlog2, color=lohUpSeg)) +
        geom_point(aes(shape=filterstatus),size=10) +
        scale_shape_manual(values=c("\u2620", "\u2665")) +
        geom_vline(xintercept = log10(cutoff.width)) +
        geom_hline(yintercept = cutoff.CNA)
      
      ggsave(paste0(outdir, "/",name,"_QC4b_segment-filtering-v2.png"), p5b, width = 16, height = 9)
    }

    # calculate the scale factor (to reduce the Y altitude for the red density line, so it will fit the histogram)
    nValues <- valueData10[, .N, by=gCont]
    valueDensData10 <- merge(valueDensData10, nValues)
    valueDensData10[, scaleFactorValue := N/10] # specific density line scale values
    
    # show the outliers at the plot borders
    if(nrow(outlierData10) > 0){
      outlierData10[formulavalue > 2, formulavalue := 2]
      outlierData10[formulavalue < -2, formulavalue := -2]
    }

    gContTmp <- gCont
    gPur.order <- sort(valueData10[,unique(gPur)], decreasing = T)
    #gPur.order <- as.character(sort(as.numeric(valueData10[,unique(gPur)]), decreasing = T))
    pSub10 <- ggplot() +
      geom_vline(xintercept = 1) +
      geom_histogram(data=valueData10, binwidth = my_binwidth, aes(formulavalue, fill=gPur), color="white", size=0.3) +
      geom_line(data=valueDensData10, aes(X,Y*scaleFactorValue), col = 2, size=0.5) +
      scale_fill_viridis_d(direction = -1) +
      facet_wrap(~factor(gPur, levels=gPur.order)) +
      theme_bw() +
      theme(axis.title.x=element_blank(), legend.position = "none", panel.grid.major = element_blank()) +
      guides(fill=guide_legend(title="contamination \n (theoretical)")) +
      scale_x_continuous(breaks = breaks_pretty(n = 15), limits = c(-2.1,2.1)) +
      scale_y_continuous(breaks = breaks_pretty(n = 10)) +
      geom_rect(data = valueData10[gCont == gContTmp], fill = NA, colour = "darkgreen", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,size=1) +
      ylab("N segments")
    
    if(nrow(outlierData10) > 0){
      pSub10 <- pSub10 + geom_histogram(data=outlierData10, binwidth = my_binwidth, aes(formulavalue),alpha=0.3)
    }
    
    ggsave(paste0(outdir, "/",name,"_histogram2.png"), pSub10, width = 16, height = 9)
    fwrite(formulaData10, paste0(outdir, "/",name,"_segment-values.txt"), sep="\t", col.names = T)
    
    
    # additional (optional but pretty) plots
    pData <- ggplot_build(pSub10)
    pDataNames <- data.table(pData$layout$layout)
    pDataNames[, gPur := as.character(`factor(gPur, levels = gPur.order)`)]
    pData <- data.table(pData$data[[2]])
    pData <- pData[count > 0]
    pData <- merge(pData, pDataNames[, .(PANEL, gPur)], by="PANEL", sort=F)
    pData[, gPurName := paste0(as.numeric(gPur)*100, "%")]
    pData[, gCont := as.character(1-as.numeric(gPur))]
    
    bestHitIndex <- pData[, which(unique(gCont) == as.character(gContTmp))]
    pData[, gPur := factor(gPur, levels = rev(gPur.order))]
    
    # VERTICAL BARS
    pData[, PANEL := 12-as.numeric(PANEL)+1] # invert PANEL order
    pData[, countP := count / max(count)]
    pData[, lowerEnd := as.numeric(PANEL)-(countP/2)+0.02]
    pData[, upperEnd := as.numeric(PANEL)+(countP/2)-0.02]
    
    ggplot(pData, aes(x=gPur,y=x)) + 
      coord_flip() + 
      theme_bw() +
      geom_hline(yintercept = 1) + 
      geom_hline(yintercept = c(1.01, 0.99), linetype=2, alpha=0.2) +
      geom_vline(xintercept = as.character(gPur), alpha=0.5) +
      geom_point(size=0) + 
      scale_color_viridis_d(direction = -1) +
      scale_y_continuous(breaks = breaks_pretty(n = 15), limits = c(-2.1,2.1)) +
      geom_segment(aes(x=lowerEnd, xend=upperEnd, y=x, yend=x, color=gPur), size=4)
    ggsave(paste0(outdir, "/",name,"_bars.png"), last_plot(), width = 16, height = 9)
    
    # VERTICAL BARS POLAR
    ggplot(pData, aes(x=gPur,y=x)) +
      geom_hline(yintercept = 1) +
      geom_point(size=0) +
      scale_color_viridis_d(direction = -1) +
      scale_size(range=c(3, 15)) +
      theme_bw() +
      scale_y_continuous(breaks = breaks_pretty(n = 15), limits = c(-2.1,2.1)) +
      geom_segment(aes(x=lowerEnd, xend=upperEnd, y=x, yend=x, color=gPur), size=3) +
      coord_polar() +
      geom_vline(xintercept = bestHitIndex+1, alpha=0.5)
    ggsave(paste0(outdir, "/",name,"_bars_polar.png"), last_plot(), width = 16, height = 9)
  }
  
  parameters[2,2] <- gCont
  parameters[3,2] <- gPur
  parameters[2,3] <- paste0(round(gCont, digits = 2)*100, "%")
  parameters[3,3] <- paste0(round(gPur, digits = 2)*100, "%")
  
  fwrite(parameters, paste0(outdir, "/",name,"_results.txt"), sep="\t", col.names = T, quote = F)
  
  return(gCont)
}
CopyNeutralPurity <- function(name, lohSegs){
  
  lohSegs.gr <- makeGRangesFromDataFrame(lohSegs)
  
  # load data
  if(method == "Copywriter"){
    seg.file <- paste0(name,"/results/Copywriter/",name,".Copywriter.segments.Mode.txt")
    cnaSegs <- fread(seg.file)
    cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs[, .(chr=Chrom, start=Start, end=End)])
  } else if(method == "HMMCopy"){
    seg.file <- paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt")
    cnaSegs <- fread(seg.file)
    cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs[, .(chr=Chrom, start=Start, end=End)])
  } else {
    stop("other methods not implemented yet")
  }
  cnaSegs <- cnaSegs[Chrom %in% c(1:22)]
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  
  hits <- findOverlaps(lohSegs.gr, cnaSegs.gr)
  tmpSegs <- cbind(lohSegs[queryHits(hits)], cnaSegs[subjectHits(hits), Mean])
  setnames(tmpSegs, "V2", "CNAlog2")
  tmpSegs <- unique(tmpSegs[, CNAlog2 := mean(CNAlog2), by=segID])
  
  combinedSegs <- merge(lohSegs, tmpSegs[, .(segID2, CNAlog2)], by="segID2", all.x=T)
  
  combinedSegs[is.na(CNAlog2), CNAlog2 := 0]
  combinedSegs[, status := "neutral"]
  combinedSegs[CNAlog2 > 0.1, status := "amp"]
  combinedSegs[CNAlog2 < -0.1, status := "del"]
  neutralSegs <- combinedSegs[status == "neutral"]
  
  maxValue <- max(neutralSegs$segDif)
  return(maxValue)
}
GetModes <- function(x,adjust,signifi,from,to) {  
  # https://stackoverflow.com/questions/27418461/calculate-the-modes-in-a-multimodal-distribution-in-r
  # Jeffrey Evans in Peak of the kernel density estimation 
  den <- density(x, kernel=c("gaussian"),adjust=adjust,from=from,to=to)
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  den.sign <- sign(s.1$y)
  a<-c(1,1+which(diff(den.sign)!=0))
  b<-rle(den.sign)$values
  df<-data.frame(a,b)
  df = df[which(df$b %in% -1),]
  modes<-s.1$x[df$a]
  density<-s.0$y[df$a]
  df2<-data.frame(modes,density)
  
  df2 <- df2[df2$density > 0.05,]
  df2 <- data.table(df2)
  setorder(df2, -density)
  return(df2)
}
GetModes2 <- function(x,adjust,signifi,from,to,bandwidth) {  
  den <- density(x, kernel=c("gaussian"),adjust=adjust,from=from,to=to, bw=bandwidth)
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  den.sign <- sign(s.1$y)
  a<-c(1,1+which(diff(den.sign)!=0))
  b<-rle(den.sign)$values
  df<-data.frame(a,b)
  df = df[which(df$b %in% -1),]
  modes<-s.1$x[df$a]
  density<-s.0$y[df$a]
  df2<-data.frame(modes,density)
  
  df2 <- df2[df2$density > 0.05,]
  df2 <- data.table(df2)
  setorder(df2, -density)
  return(df2)
}
GetModesSimple <- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  return(modes)
}

if(species=="Human"){
  chrom.sizes = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
                  138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
                  83257441,80373285,58617616,64444167,46709983,50818468)
  names(chrom.sizes) = c(1:22)
} else if(species=="Mouse"){
  chrom.sizes = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,
                  130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
  names(chrom.sizes) = c(1:19)
}

# approach1: ~70 seconds
# approach2: ~1.5 seconds





# now we can correct CNA, LOH, SNV and SNP data

# QC: given a bam file and a data.table with mutations (chrom, pos, pos), will return MAPQ and BAQ values for the mutation position
Calculate_MAPQ_BAQ <- function(inputbam, inputSNV, tempfolder){
  # tempfolder <- "/home/rad/Downloads/temp4/" # for samtools and jvarkit calculations 
  # inputbam <- copy(tbam)
  # inputSNV <- copy(SNVs)
  
  # first check if index was found, else try to copy or stop (we could execute samtools index here as well)
  indexname1 <- paste0(inputbam, ".bai") # this was what we need
  indexname2 <- gsub("\\.bam$", ".bai", inputbam) # sometimes this is the filename
  if(!file.exists(indexname1)){
    if(file.exists(indexname2)){
      file.copy(indexname2, indexname1)
    } else {
      stop("No index file for BAM found, please rerun samtools index first")
    }
  }
  
  inputbed <- unique(inputSNV[, .(CHROM, POS, POS)])
  
  outbam.file <- paste0(tempfolder,"/",name,"_snv-regions.bam")
  regionString <- paste0(inputbed[, paste0(CHROM, ":", POS, "-", POS)], collapse = " ") # like this, since -L is slow and gives weird results
  command1 <- paste0("samtools view -bh -F 4 -F 1024 ",inputbam," ",regionString," > ",outbam.file)
  system(command1)
  
  outbam.file2 <- paste0(tempfolder,"/",name,"_snv-regions-sorted.bam")
  command2 <- paste0("samtools sort -f ",outbam.file ," ", outbam.file2)
  system(command2)
  
  command3 <- paste0("samtools index ",outbam.file2)
  system(command3)
  
  outstats.file <- paste0(tempfolder,"/",name,"_snv-regions.stats")
  command4 <- paste0("java -jar /home/rad/packages/jvarkit/dist/sam2tsv.jar -R /media/rad/HDD1/MoCaSeq_ref/GRCh38.p12/GRCh38.p12.fna ",outbam.file," > ",outstats.file)
  system(command4)
  
  # read and process the QC data
  readstats <- fread(outstats.file)
  
  file.remove(outbam.file, outbam.file2, paste0(outbam.file2, ".bai"), outstats.file)
  
  colnames(readstats) <- gsub("-","",colnames(readstats))
  setnames(readstats, "#ReadName", "ReadName")
  readstats <- unique(readstats) 
  readstats[, baseID := .I]
  readstats[READQUAL != ".", READQUALNUM := utf8ToInt(READQUAL), by=baseID]
  readstats[, READQUALNUM := READQUALNUM - 33]
  readstats[, mutID2 := paste(CHROM, REFPOS1, REFBASE, READBASE, sep="-")]
  readstats[, CHROM := as.character(CHROM)]
  
  # assign mutation positions
  subs.positions <- inputSNV[type %in% c("SNV", "MNV"), POS]
  indel.positions <- inputSNV[type %in% c("del", "ins", "dup"), POS]
  
  # only get MAPQ for these mutations
  if(length(indel.positions) >= 1){
    readstats.indels <- readstats[REFPOS1 %in% indel.positions]
    readstats.indels <- unique(readstats.indels[, .(POS=as.numeric(REFPOS1), MAPQ=mean(MAPQ, na.rm=T), BAQ=NA), by=.(CHROM, REFPOS1)])
    readstats.indels[, REFPOS1 := NULL]
  } else {
    readstats.indels <- data.table()
  }
  
  # only single/multi nucleotide variants
  readstats <- readstats[REFPOS1 %in% subs.positions]
  
  # remove all reads with REF=ALT
  readstats <- readstats[READBASE != REFBASE]
  
  inputSNV[, mutID2 := paste(CHROM, POS, substr(REF, 1, 1), substr(ALT, 1, 1), sep="-")] # this way we include multi substitutions and insertions 
  readstats <- merge(readstats, inputSNV[, .(mutID, mutID2, REF)], by="mutID2", all.x=T)
  
  # here we remove mutations which are at the given locations but show a different base substitution (e.g. C>T but there is also a C>A which is removed)
  readstats <- readstats[!is.na(REF)]
  
  # reduce to single information per mutation
  readstats <- unique(readstats[, .(MAPQ=mean(MAPQ, na.rm=T), BAQ=mean(READQUALNUM, na.rm=T)), by=.(mutID)])
  
  # first bind the simple mutations onto the output using mutation ID, then indels only based on CHROM/POS
  readstats <- merge(inputSNV[POS %in% subs.positions], readstats, by=c("mutID"))
  readstats.indels <- merge(inputSNV[POS %in% indel.positions], readstats.indels, by=c("CHROM", "POS"), all.x=T)
  
  # merging by CHROM/POS is dangerous, since there can be a SNV at the same position as an indel, we catch it like this
  readstats.indels <- readstats.indels[!grepl(">", HGVS_C)]
  
  outputSNV <- rbind(readstats, readstats.indels)
  
  # some mutations are only found in the m2 bam from Mutect2, here we flag these
  
  errorSNV <- inputSNV[!inputSNV$mutID %in% outputSNV$mutID]
  errorSNV[, MAPQ := 0]
  errorSNV[, BAQ := 0]
  outputSNV <- rbind(outputSNV, errorSNV)
  
  if(nrow(inputSNV) != nrow(outputSNV)){
    stop("error A")
    # test1 <- merge(inputSNV[POS %in% subs.positions], readstats, by=c("mutID"))
    # test2 <- merge(inputSNV[POS %in% indel.positions], readstats.indels, by=c("CHROM", "POS"), all.x=T)
    # outputSNV <- rbind(test1, test2)
    # nrow(outputSNV)
    # nrow(inputSNV)
  }
  
  setcolorder(outputSNV, colnames(inputSNV))
  outputSNV[, mutID2 := NULL]
  
  outputSNV[, MAPQ := round(MAPQ, digits = 1)]
  outputSNV[, BAQ := round(BAQ, digits = 2)]
  
  return(outputSNV)
}

# DEBUG 
if(F){
  cutoff.LOH <- 0.55
  cutoff.width <- 5000000
  cutoff.CNA <- 0.5 # this is different then used for purity analysis
}
ProcessSample <- function(name, purityDT, LOHfile, CNAfile, SNVfile=NULL, cutoff.LOH=0.55, cutoff.width=5000000, cutoff.CNA=0.5){
  
  sample.gCont <- purityDT[sample == name, gCont]
  
  # CNA AND LOH PREPROCESSING
  # LOHfile <- paste0(LOHpath, "/", name, "_LOH_Segments.tsv")
  lohSegs <- fread(LOHfile)
  # 
  # CNAfile <- paste0(outdir, "/",name,"_CNA_",cnaMethod,"_smooth.tsv")
  cnaSegs <- fread(CNAfile)
  
  # combined LOH and CNA
  lohSegs[, Chrom := as.numeric(Chrom)]
  lohSegs[, width := end-start]
  cnaSegs[, Chrom := as.numeric(Chrom)]
  
  # apply size filter 
  lohSegs <- lohSegs[width > cutoff.width]
  cnaSegs <- cnaSegs[width > cutoff.width]
  lohSegs <- lohSegs[upSeg > cutoff.LOH]
  cnaSegs <- cnaSegs[Mean < cutoff.CNA]
  lohSegsRaw <- copy(lohSegs)
  cnaSegsRaw <- copy(cnaSegs)
  lohSegs.gr <- makeGRangesFromDataFrame(lohSegs)
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  
  # CNA
  setnames(cnaSegsRaw, "Mean", "CNAlog2")
  cnaSegsRaw[, CNA := 2^CNAlog2]
  cnaSegsRaw[, segID := .I]
  hits <- findOverlaps(cnaSegs.gr, lohSegs.gr)
  cnaTmp <- cbind(cnaSegsRaw[queryHits(hits)], lohSegsRaw[subjectHits(hits), .(lohUpSeg=upSeg, lohLowSeg=lowSeg, start.other=start, end.other=end)]) 
  cnaSegs <- unique(cnaTmp[cnaTmp[, .I[lohUpSeg == max(lohUpSeg)], by=segID]$V1])

  # LOH
  lohSegsRaw <- lohSegsRaw[, .(Chrom, Start=start, End=end, lohUpSeg=upSeg, lohLowSeg=lowSeg)]
  lohSegsRaw[, segID := .I]
  hits <- findOverlaps(lohSegs.gr, cnaSegs.gr)
  lohTmp <- cbind(lohSegsRaw[queryHits(hits)], cnaSegsRaw[subjectHits(hits), .(CNAlog2, CNA, start.other=Start, end.other=End)])
  lohSegs <- unique(lohTmp[lohTmp[, .I[abs(CNAlog2) == max(abs(CNAlog2))], by=segID]$V1])

  lohSegs.gr <- makeGRangesFromDataFrame(lohSegs)
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  
  
  # SNVs
  #SNVfile <- paste0(name,"/results/Mutect2/",name,".Mutect2.NoCommonSNPs.OnlyImpact.txt")
  SNVs <- fread(SNVfile)
  SNVs[, AF := NULL]
  
  SNVs[, mutID := paste(CHROM, POS, REF, ALT, sep="-")]
  colnames(SNVs) <- gsub("(GEN|ANN)\\[(.*)\\]\\.", "\\2", colnames(SNVs))
  colnames(SNVs) <- gsub("\\[|\\]", "", colnames(SNVs))
  colnames(SNVs) <- gsub("\\*", "", colnames(SNVs))
  setnames(SNVs, c("TumorAD0","TumorAD1","NormalAD0","NormalAD1"), c("TumorRefC","TumorAltC","NormalRefC","NormalAltC"))
  SNVs[, TumorTotalC := TumorRefC+TumorAltC]
  SNVs[, gCont := sample.gCont]
  SNVs <- SNVs[, -c("dbNSFP_MetaLR_pred", "dbNSFP_MetaSVM_pred", "dbNSFP_PROVEAN_pred",
                    "dbNSFP_SIFT_pred", "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_Polyphen2_HVAR_pred",
                    "CLNSIG", "CLNREVSTAT", "CLNDN")]
  
  # apply hard filtering flag (filter later)
  SNVs[, filtered := "no"]
  SNVs[TumorAF < 0.1, filtered := "yes"]
  SNVs[TumorRefC + TumorAltC < 20, filtered := "yes"]
  SNVs[NormalRefC + NormalAltC < 20, filtered := "yes"]
  SNVs[TumorAltC < 3, filtered := "yes"]
  SNVs[NormalAltC != 0, filtered := "yes"]
  
  SNVs.gr <- makeGRangesFromDataFrame(SNVs[, .(chr=CHROM, start=POS, end=POS)])
  
  # categorize mutations
  SNVs[grepl(">", HGVS_C), type  := "SNV"]
  SNVs[nchar(REF) > 1 & nchar(ALT) > 1 & nchar(REF) == nchar(ALT), type := "MNV"]
  SNVs[grepl("ins", HGVS_C) & !grepl("del", HGVS_C), type := "ins"]
  SNVs[grepl("del", HGVS_C) & !grepl("ins", HGVS_C), type := "del"]
  SNVs[grepl("dup", HGVS_C), type := "dup"]
  
  SNVs[, LOH := 0.5] # defaults
  SNVs[, CNA := 0] # defaults
  SNVs[, CNAlog2 := 0]
  hits1 <- findOverlaps(SNVs.gr, lohSegs.gr)
  SNVs[queryHits(hits1), LOH := lohSegs[subjectHits(hits1), round(lohUpSeg, digits = 2)]]
  hits2 <- findOverlaps(SNVs.gr, cnaSegs.gr)
  SNVs[queryHits(hits2), CNAlog2 := cnaSegs[subjectHits(hits2), CNAlog2]]
  SNVs[, CNAlog2 := round(CNAlog2, digits = 1)]
  SNVs[, CNA := 2^CNAlog2]
  
  # CORRECT SNV
  # only do this if the purity is above 20%
  if(sample.gCont <= 0.8){
    SNVs[, TumorAF.corrected := TumorAltC / (TumorTotalC - (TumorTotalC / CNA * sample.gCont))]
    SNVs[, TumorAF.corrected := round(TumorAF.corrected, digits = 3)]
  } else {
    SNVs[, TumorAF.corrected := as.numeric(NA)]
  }
  
  SNVs[, sample := name]
  
  # CORRECT CNA
  cnaSegs[, CNA.corrected := (CNA - sample.gCont) / (1 - sample.gCont)]
  setcolorder(cnaSegs, c("Chrom", "Start", "End", "width"))
  
  # CORRECT LOH
  lohSegs[, lohUpSeg.corrected := (lohUpSeg - ((0.5*(1/CNA))*sample.gCont)) / (1 - ((1/CNA)*sample.gCont))]
  lohSegs[, lohLowSeg.corrected := (lohLowSeg - ((0.5*(1/CNA))*sample.gCont)) / (1 - ((1/CNA)*sample.gCont))]
  
  # POST FILTERING OF READS (or not)
  # tempfolder <- "/home/rad/Downloads/temp4/" # for samtools and jvarkit calculations 
  # tbam <- paste0(getwd(), "/",name,"/results/bam/",name,".Tumor.bam")
  # SNVs <- Calculate_MAPQ_BAQ(tbam, SNVs, tempfolder)
  # 
  # SNVs[, filtered2 := "no"]
  # SNVs[MAPQ < 60, filtered2 := "yes"]
  # SNVs[BAQ < 30, filtered2 := "yes"]
  
  return(list(SNV=SNVs, CNA=cnaSegs, LOH=lohSegs))
}


if(F){
  
  sampleSNVs <- copy(SNVs)
  
  # PLOT SNVs multie sample
  sampleSNVs.filtered <- sampleSNVs[filtered == "no"]
  sampleSNVs.filtered <- sampleSNVs.filtered[!is.na(TumorAF.corrected)]
  SNV.plotdata <- data.table(gather(sampleSNVs.filtered[, .(mutID, TumorAF, TumorAF.corrected, sample)], name, value, TumorAF:TumorAF.corrected))
  # ggplot(SNV.plotdata, aes(value, fill=name)) + 
  #   geom_density(alpha=0.5, bw=0.01) + 
  #   facet_wrap(~sample) +
  #   scale_x_continuous(breaks = pretty_breaks(n=10)) +
  #   geom_vline(xintercept = c(0, 0.5, 0.75,1), linetype=2, alpha=0.5) +
  #   coord_cartesian(xlim = c(-0.1, 1.5))
  
  # PLOT SNVs single sample
  SNV.plotdata <- gather(sampleSNVs.filtered[, .(mutID, TumorAF, TumorAF.corrected)], name, value, -mutID)
  
  #SNVs <- sampleSNVs[sample == "hPDAC02_LivMet-1" & GENE %in% keepGenes]
  
  SNV.plotdata2a <- sampleSNVs.filtered[, .(CHROM, POS, TumorAF, TumorAF.corrected)]
  SNV.plotdata2b <- gather(SNV.plotdata2a, name, value, TumorAF:TumorAF.corrected)
  
  # ggplot() + 
  #   geom_point(data = SNV.plotdata2b, aes(x=POS, y=value, color=name), size=5) +
  #   geom_segment(data = SNV.plotdata2a, aes(x=POS, xend=POS, y=TumorAF, yend=TumorAF.corrected), size = 0.5, arrow = arrow(length = unit(0.3, "cm")))
  
  SNV.plotdata3a <- sampleSNVs.filtered[, .(GENE, POS, ID=paste(GENE, CHROM, POS, REF, ALT, sep="-"), TumorAF, TumorAF.corrected)]
  SNV.plotdata3a <- SNV.plotdata3a[!is.na(TumorAF.corrected)]
  setorder(SNV.plotdata3a, TumorAF.corrected)
  SNV.plotdata3a[, ID := factor(ID, levels = SNV.plotdata3a[, unique(ID)])]
  SNV.plotdata3b <- gather(SNV.plotdata3a, name, value, TumorAF:TumorAF.corrected)
  
  #  geom_segment(data = SNV.plotdata3a, aes(y=ID, yend=ID, x=TumorAF, xend=TumorAF.corrected), size = 0.3, arrow = arrow(length = unit(0.2, "cm"))) +
  ggplot() + 
    geom_point(data=SNV.plotdata3b, aes(x=value, y=ID, color=name), size=7) +
    scale_y_discrete(labels= SNV.plotdata3b$GENE) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    geom_segment(data = SNV.plotdata3a, aes(y=ID, yend=ID, x=TumorAF, xend=TumorAF.corrected), size = 0.2, arrow = arrow(length = unit(0.2, "cm"))) +
    geom_vline(xintercept = c(0,1))
  
  #ggsave("/home/rad/Downloads/temp4/SNV_corrected.pdf", last_plot(), width = 16, height = 9)
  
  ggplot() + 
    geom_point(data=SNV.plotdata3b, aes(x=ID, y=value, color=name), size=3) +
    geom_segment(data = SNV.plotdata3a, aes(y=TumorAF, yend=TumorAF.corrected, x=ID, xend=ID), size = 0.2, arrow = arrow(length = unit(0.2, "cm"))) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    scale_x_discrete(labels= SNV.plotdata3b$GENE)
  
  
  CNA.plotdata <- data.table(gather(cnaSegs[, .(Chrom, Start, End, CNA, CNA.corrected)], name, value, CNA:CNA.corrected))
  ggplot(CNA.plotdata) + 
    geom_histogram(aes(value, fill=name), bins=100) + 
    scale_x_continuous(breaks = pretty_breaks(n=30)) +
    facet_wrap(~name, ncol = 1) +
    geom_vline(xintercept = c(0.5, 0.75, 1), linetype=2, alpha=0.5)
  
  ggsave("/home/rad/Downloads/temp4/CNA_corrected_histogram.pdf", last_plot(), width = 16, height = 9)
  
  cnaBins <- fread("hPDAC02_LivMet-1/results/Copywriter/hPDAC02_LivMet-1.Copywriter.log2RR.Mode.txt")
  cnaBins <- fread("hPDAC05_LivMet-1/results/Copywriter/hPDAC05_LivMet-1.Copywriter.log2RR.Mode.txt")
  cnaBins <- fread("B25_NFM_IP_P52/results/Copywriter/B25_NFM_IP_P52.Copywriter.log2RR.Mode.txt")
  cnaBins <- cnaBins[log2Ratio %between% c(-5,5)]
  cnaSegs[, Mid := Start+((End-Start)/2)]
  cnaBins <- cnaBins[Chrom %in% c(1:22)]
  
  # ggplot() +
  #   geom_hline(yintercept = 1) +
  #   geom_point(data=cnaBins, aes(Start, 2^log2Ratio), shape=".") +
  #   geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=CNA, yend=CNA), color="red") +
  #   geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=CNA.corrected, yend=CNA.corrected), color="blue") +
  #   geom_segment(data=cnaSegs, aes(y=CNA, yend=CNA.corrected, x=Mid, xend=Mid), size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
  #   facet_wrap(~Chrom, scales = "free_x", nrow = 1)
  
  
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-0.25, -0.5, -1, 0.25, 0.5, 1), linetype=2) +
    geom_point(data=cnaBins, aes(Start, log2Ratio), shape=".", alpha=0.1) +
    geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=log2(CNA), yend=log2(CNA)), color="red", size=1) +
    geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=log2(CNA.corrected), yend=log2(CNA.corrected)), color="blue", size=1) +
    facet_wrap(~factor(Chrom, levels=c(1:22)), scales = "free_x", nrow=1) +
    theme_minimal() +
    geom_segment(data=cnaSegs, aes(y=log2(CNA), yend=log2(CNA.corrected), x=Mid, xend=Mid), size = 0.5, arrow = arrow(length = unit(0.2, "cm")), color="black")
  
  #ggsave("/home/rad/Downloads/temp4/CNA_corrected_segments.pdf", last_plot(), width = 16, height = 9)
  
  lohTmp <- copy(lohSegs)
  lohTmp[, Mid := Start+((End-Start)/2)]
  LOH.plotdata <- data.table(gather(lohTmp[, .(Chrom, Start, End, lohUpSeg, lohUpSeg.corrected)], name, value, lohUpSeg:lohUpSeg.corrected))
  ggplot(LOH.plotdata) + 
    geom_histogram(aes(value, fill=name), bins=100) + 
    scale_x_continuous(breaks = pretty_breaks(n=30)) +
    facet_wrap(~name, ncol = 1) +
    geom_vline(xintercept = c(0.5, 0.75, 1), linetype=2, alpha=0.5)
  
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_segment(data=lohTmp, aes(x=Start, xend=End, y=lohUpSeg, yend=lohUpSeg), color="red", size=1) +
    geom_segment(data=lohTmp, aes(x=Start, xend=End, y=lohUpSeg.corrected, yend=lohUpSeg.corrected), color="blue", size=1) +
    facet_wrap(~factor(Chrom, levels=c(1:22)), scales = "free_x", nrow=1) +
    theme_minimal() +
    geom_segment(data=lohTmp[abs(lohUpSeg-lohUpSeg.corrected) > 0.01], aes(y=lohUpSeg, yend=lohUpSeg.corrected, x=Mid, xend=Mid), size = 0.5, arrow = arrow(length = unit(0.2, "cm")), color="black")
  
  
  # allele specific copy number?!
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  hits <- findOverlaps(cnaSegs.gr, lohSegs.gr)
  cnaSegs[, segID := paste(Chrom, Start, End, sep="-")]

  tmp <- cbind(cnaSegs[queryHits(hits)], lohSegs[subjectHits(hits), .(lohUpSeg.corrected, lohLowSeg.corrected)]) 
  tmp <- unique(tmp[tmp[, .I[lohUpSeg.corrected == max(lohUpSeg.corrected)], by=segID]$V1])
  tmp <- merge(cnaSegs, tmp[, .(segID, lohUpSeg.corrected, lohLowSeg.corrected)], all.x=T)
  
  tmp[, CNA.corrected.A := CNA.corrected * lohUpSeg.corrected]
  tmp[, CNA.corrected.B := CNA.corrected * lohLowSeg.corrected]
  
  tmp[lohUpSeg.corrected != 0.5 & lohUpSeg.corrected != 0.5]
  CNA.plotdata <- data.table(gather(tmp[, .(Chrom, Start, End, CNA, CNA.corrected, CNA.corrected.A, CNA.corrected.B)], name, value, CNA.corrected:CNA.corrected.B))
  
  tmp2 <- CNA.plotdata[name %in% c("CNA.corrected.A", "CNA.corrected.B")]
  tmp2[, name := "CNA.corrected.combined"]
  CNA.plotdata <- rbind(CNA.plotdata, tmp2)
  
  ggplot(CNA.plotdata) + 
    geom_histogram(aes(value, fill=name), bins=50) + 
    scale_x_continuous(breaks = pretty_breaks(n=30)) +
    facet_wrap(~name, ncol = 1) +
    geom_vline(xintercept = c(0.25, 0.33, 0.5, 0.66, 0.75, 1), linetype=2, alpha=0.5)
  
  ggsave("/home/rad/Downloads/temp4/CNA_corrected_alleles.pdf", last_plot(), width = 16, height = 9)
  #ggsave("/home/rad/Downloads/temp4/CNA_corrected_alleles_raw.pdf", last_plot(), width = 16, height = 9)
  
  tmp[CNA.corrected.A %between% c(0.45, 0.55)]
  
}


# correct CNA for ALL segments without filtering and merging with LOH?!
ProcessCNA.nofilters <- function(name, purityDT){
  
  sample.gCont <- purityDT[sample == name, gCont]

  CNAfile <- paste0(outdir, "/",name,"_CNA_",cnaMethod,"_smooth.tsv")
  cnaSegs <- fread(CNAfile)
  cnaSegs[, Chrom := as.numeric(Chrom)]

  cnaSegsRaw <- copy(cnaSegs)
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  
  # CNA
  setnames(cnaSegsRaw, "Mean", "CNAlog2")
  cnaSegsRaw[, CNA := 2^CNAlog2]
  cnaSegsRaw[, segID := .I]
  cnaSegs <- copy(cnaSegsRaw)
  cnaSegs.gr <- makeGRangesFromDataFrame(cnaSegs)
  
  # CORRECT CNA
  cnaSegs[, CNA.corrected := (CNA - sample.gCont) / (1 - sample.gCont)]
  setcolorder(cnaSegs, c("Chrom", "Start", "End", "width"))

  return(cnaSegs)
}

if(F){
  cnaSegs <- ProcessCNA.nofilters("hPDAC05_LivMet-1", purityDT)
  CNA.plotdata <- data.table(gather(cnaSegs[, .(Chrom, Start, End, CNA, CNA.corrected)], name, value, CNA:CNA.corrected))
  ggplot(CNA.plotdata) + 
    geom_histogram(aes(value, fill=name), bins=100) + 
    scale_x_continuous(breaks = pretty_breaks(n=30)) +
    facet_wrap(~name, ncol = 1) +
    geom_vline(xintercept = c(0.5, 0.75, 1), linetype=2, alpha=0.5)
  ggsave("/home/rad/Downloads/temp4/CNA_corrected_histogram_raw.pdf", last_plot(), width = 16, height = 9)
  
  cnaBins <- fread("hPDAC05_LivMet-1/results/Copywriter/hPDAC05_LivMet-1.Copywriter.log2RR.Mode.txt")
  cnaBins <- cnaBins[log2Ratio %between% c(-5,5)]
  cnaSegs[, Mid := Start+((End-Start)/2)]
  cnaBins <- cnaBins[Chrom %in% c(1:22)]

  ggplot() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-0.25, -0.5, -1, 0.25, 0.5, 1), linetype=2) +
    geom_point(data=cnaBins, aes(Start, log2Ratio), shape=".", alpha=0.1) +
    geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=log2(CNA), yend=log2(CNA)), color="red", size=1) +
    geom_segment(data=cnaSegs, aes(x=Start, xend=End, y=log2(CNA.corrected), yend=log2(CNA.corrected)), color="blue", size=1) +
    facet_wrap(~factor(Chrom, levels=c(1:22)), scales = "free_x", nrow=1) +
    theme_minimal() +
    geom_segment(data=cnaSegs, aes(y=log2(CNA), yend=log2(CNA.corrected), x=Mid, xend=Mid), size = 0.5, arrow = arrow(length = unit(0.2, "cm")), color="black")
  ggsave("/home/rad/Downloads/temp4/CNA_corrected_segments_raw.pdf", last_plot(), width = 16, height = 9)
  
}






# DEBUG: which of the 4 approaches is the best?
if(F){
  
  dir.create(paste0(outdir, "/segmentsFrom-CNA_cutoff01"))
  dir.create(paste0(outdir, "/segmentsFrom-CNA_cutoff-01"))
  dir.create(paste0(outdir, "/segmentsFrom-LOH_cutoff01"))
  dir.create(paste0(outdir, "/segmentsFrom-LOH_cutoff-01"))
  
  for(name in allsamples){
    print(name)
    
    #nres <- CopyNeutralPurity(name, lohSegs) # totally deprecated, unsmoothed segments etc.
    
    # Load LOH
    LOHfile <- paste0(LOHpath, "/", name, "_LOH_Segments.tsv")
    lohSegs <- fread(LOHfile)
    lohSegs[, segID2 := paste0(Chrom,"-",segID)]
    
    # Load CNA
    cnaRes <- LoadCNA(name, "Copywriter", smoothCNA=T)
    cnaSegs <- cnaRes$cnaSegs
    ggsave(paste0(outdir, "/",name,"_QC_CNAsmoothing.png"), cnaRes$pCNAsmoothing, width = 16, height = 9)
    
    # from CNV segments
    res1 <- EstimatePurity(name, lohSegs, paste0(outdir, "segmentsFrom-CNA_cutoff-01"), segmentsFrom = "CNA", cutoff.CNA = -0.1)
    res2 <- EstimatePurity(name, lohSegs, paste0(outdir, "segmentsFrom-CNA_cutoff01"), segmentsFrom = "CNA", cutoff.CNA = 0.1)
    
    # from LOH segments
    res3 <- EstimatePurity(name, lohSegs, paste0(outdir, "segmentsFrom-LOH_cutoff-01"), segmentsFrom = "LOH", cutoff.CNA = -0.1)
    res4 <- EstimatePurity(name, lohSegs, paste0(outdir, "segmentsFrom-LOH_cutoff01"), segmentsFrom = "LOH", cutoff.CNA = 0.1)
    
    # DEBUG comparison, for this the function needs to return the plot instead of the data
    if(is.null(res1)){
      res1 <- ggplot()
    }
    if(is.null(res2)){
      res2 <- ggplot()
    }
    if(is.null(res3)){
      res3 <- ggplot()
    }
    if(is.null(res4)){
      res4 <- ggplot()
    }
    
    resAll <- ggarrange(res1+ggtitle("segmentsFrom=CNA, cutoff=-0.1"),
                        res2+ggtitle("segmentsFrom=CNA, cutoff=+0.1"),
                        res3+ggtitle("segmentsFrom=LOH, cutoff=-0.1"),
                        res4+ggtitle("segmentsFrom=LOH, cutoff=+0.1"))
    ggsave(paste0(outdir, "/", name, "_histograms_comparison.png"),resAll, width = 16, height = 9)
    
  }
}