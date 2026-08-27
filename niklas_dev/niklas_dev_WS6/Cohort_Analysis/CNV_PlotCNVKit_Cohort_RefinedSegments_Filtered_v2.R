suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(pdftools))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))
args = commandArgs(trailingOnly=TRUE)

# DEBUG
# setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")
# name <- "P4GNJ9"
# ylim <- 5
# workdir <- "CNVKit_refine-segmentsTest/cbs_1e3"

# INPUT
name <- args[1] #used for naming in- and output files
workdir <- args[2]
ylim <- as.numeric(args[3])

setwd(workdir)
SegmentsFile <- paste0(name,".Tumor.center-mode.cns")
CountsFile <- paste0(name,".Tumor.center-mode.cnr")

Output <- paste0(name, ".Chromosomes.CNV.CNVKit.",ylim,".filtered.pdf")

# LOAD AND PREPARE LOOKUPS
GetMappability <- function(binDT){
  binDT[, ID := .I]
  segDT <- unique(binDT[, .(chr=chromosome, start=start,end=end, ID)])
  a <- makeGRangesFromDataFrame(segDT)
  hits <- findOverlaps(a, mappab.gr)
  
  mappab.subset <- mappab[subjectHits(hits)] # get all the overlapping region (which have mappability = 1)
  mappab.subset[, ID := segDT[queryHits(hits), ID]]
  mappab.subset[, width := end-start] # calculate the width of each mappa region
  
  # some regions overlap at the borders, here we want to calcualte the % of overlap to get the true width 
  overlaps <- pintersect(a[queryHits(hits)], mappab.gr[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(mappab.gr[subjectHits(hits)])
  mappab.subset[, pOverlap := percentOverlap]
  mappab.subset[, width := round(width * pOverlap)]
  
  segment.maprange <- mappab.subset[, sum(width), by=ID] # for each segment, sum the amount of mappa regions (1-regions)
  
  out <- copy(binDT)
  out <- merge(out, segment.maprange, by="ID", all.x=T, sort=F)
  setnames(out, "V1", "map1width")
  out[, width := start-end] # calculate the percentage of mappa=1 position (assuming everything else is 0)
  out[, mappa := map1width / width]
  out[mappa > 1, mappa := 1] # rounding error from width and percantage
  
  out[, mappAbility := round(mappa, digits=2)]
  
  out <- unique(out[, .(ID, mappAbility)])
  return(out)
}
#load("/media/rad/HDD1/Lookups/hg38_cents-haplo-mappability.RData")
load("../hg38_cents-haplo-mappability.RData")

germCNV <- fread("../ClinVar_CNVs.txt")
germCNV <- germCNV[clinSign == "Benign"]
germCNV <- germCNV[, .(chr=`#chrom`, start=chromStart, end=chromEnd)]
germCNV[, class := "ClinGenCNV"]
germCNV[, chr := gsub("chr","", chr)]
germCNV.gr <- makeGRangesFromDataFrame(germCNV[, .(chr, start, end, excludeID=.I)])

custom <- data.table(chr=1, start=1, end=1)
custom[, class := "Custom"]
custom.gr <- makeGRangesFromDataFrame(custom[, .(chr, start, end, excludeID=.I)])

#removeGR <- list(cents.gr, haplo.gr, custom.gr)
removeGR <- list(cents.gr, haplo.gr)
removeGR <- do.call(c, as(removeGR, "GRangesList"))


# LOAD DATA
segDT <- fread(SegmentsFile)
countsDT <- fread(CountsFile)

# LOAD DATA
#countsDT <- data.table(Counts[[4]])
countsDT[, gene := NULL]
countsDT[, ID := .I]
countsGR <- makeGRangesFromDataFrame(countsDT[, .(chr=chromosome, start, end)])

#segDT <- data.table(Segments[[2]])
segDT[, gene := NULL]
segDT[, ID := .I]

# SHRINK OUTLIERS
countsDT[, outlier := F]
countsDT[, plotValue := log2]
countsDT[log2 >= ylim, plotValue := ylim] 
countsDT[log2 <= -ylim, plotValue := -ylim]

segDT[, outlier := F]
segDT[, plotValue := log2]
segDT[log2 >= ylim, plotValue := ylim]
segDT[log2 <= -ylim, plotValue := -ylim]

# FILTER REGIONS
# remove centromeres and alternative haplotype regions from the data, data was generated in Cohort_GenerateOverlay_FIlterNoise
hits <- findOverlaps(countsGR,cents.gr)
countsDT[unique(queryHits(hits)), hitCent := "yes"]

hits <- findOverlaps(countsGR,haplo.gr)
countsDT[unique(queryHits(hits)), hitHaploAlt := "yes"]

hits <- findOverlaps(countsGR, germCNV.gr)
countsDT[unique(queryHits(hits)), hitClinVarCNV := "yes"]

hits <- findOverlaps(countsGR,custom.gr)
countsDT[unique(queryHits(hits)), hitCustom := "yes"]

# input: x[, .(chr=chromosome, start,end, ID)])
mapvals <- GetMappability(countsDT)
mapvals[is.na(mappAbility), mappAbility := 0]
countsDT <- merge(countsDT, mapvals, all.x=T, sort=F, by="ID")
countsDT[, mappAbility := abs(mappAbility)]

countsDT[, noise := "valid SNPs"]
mappabilityCutoff <- 0.5

countsDT[mappAbility < mappabilityCutoff , noise := "low mappability"]
countsDT[!is.na(hitCustom), noise := "custom removed"]
countsDT[!is.na(hitClinVarCNV), noise := "ClinVar germline CNVs"]
countsDT[!is.na(hitHaploAlt), noise := "alternative haplotype"]
countsDT[!is.na(hitCent), noise := "centromere"]

# set color vector
noiseColors <- c("#ed0cc1", "#0ced7d", "#0cc7ed", "#f3d61f", "#4c1ff3", "#000000")
names(noiseColors) <- c("low mappability", "custom removed", "ClinVar germline CNVs", "alternative haplotype", "centromere", "valid SNPs")


# REMOVE BY FILTER

# segments with low mappability can be removed like this
mapvals <- GetMappability(segDT)
mapvals[is.na(mappAbility), mappAbility := 0]
segDT <- merge(segDT, mapvals, all.x=T, sort=F, by="ID")
segDT[, mappAbility := abs(mappAbility)]

# set new data
plotsegDT <- segDT[mappAbility > mappabilityCutoff]
segGR <- makeGRangesFromDataFrame(plotsegDT[, .(chr=chromosome, start, end)])

# this will "cut" out overlapping noise regions out of the existing segments (and then merge the information back in)
tmpDT <- plotsegDT[, .(chr=chromosome, start, end)]
tmpDT[, start := start+1]
tmpDT[, end := end-1]
tmpGR <- makeGRangesFromDataFrame(tmpDT)
tmpGR <- setdiff(tmpGR, removeGR, ignore.strand=T)
tmpDT <- as.data.table(tmpGR)
hits <- findOverlaps(tmpGR, segGR)
plotsegDT <- cbind(tmpDT[queryHits(hits), -("seqnames")], plotsegDT[subjectHits(hits), -c("start", "end")])

removedSegs <- segDT[mappAbility < mappabilityCutoff] # we will use this to further filter counts
removedSegsGR <- makeGRangesListFromDataFrame(removedSegs[, .(chr=chromosome, start, end)])


plotcountsDT <- copy(countsDT)

# additionally, filter all SNPS if the entire region is of low mappability
plotcountsGR <- makeGRangesListFromDataFrame(plotcountsDT[, .(chr=chromosome, start, end)])
hits <- findOverlaps(plotcountsGR, removedSegsGR)
plotcountsDT[queryHits(hits), noise := "low mappability"]

plotcountsDT <- plotcountsDT[noise %in% c("valid SNPs", "ClinVar germline CNVs")]

#save(plotcountsDT, plotsegDT, file = "../test.rdata")

PlotChrom <- function(myChrom){
  pBase <- ggplot() +
    coord_cartesian(ylim=c(-ylim,ylim)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype=1, color="black") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = pretty_breaks(n = 10), labels=function(x)x/1000000)
  
  pFiltered <- pBase + 
    geom_point(data=plotcountsDT[chromosome == myChrom], aes(start, plotValue), size=0.001, alpha=0.6) +
    geom_segment(data=plotsegDT[chromosome == myChrom], aes(x=start, xend=end, y=plotValue, yend=plotValue), color="red", size=1.5)
  
  # pRaw <- pBase  +
  #   geom_point(data=countsDT[chromosome == myChrom], aes(start, plotValue), size=0.001, alpha=0.6) +
  #   geom_segment(data=segDT[chromosome == myChrom], aes(x=start, xend=end, y=plotValue, yend=plotValue), color="red", size=1.5) +
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  pRawNoiseMarked <- pBase +
    geom_point(data=countsDT[chromosome == myChrom], aes(start, plotValue, color=noise), size=0.001, alpha=0.6) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(values = noiseColors) +
    theme(legend.position="bottom", legend.title = element_blank()) +
    geom_segment(data=segDT[chromosome == myChrom], aes(x=start, xend=end, y=plotValue, yend=plotValue), color="red", size=1.5) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  #pOut <- ggarrange(pRaw, pRawNoiseMarked, pFiltered, ncol=1)
  pOut <- ggarrange(pRawNoiseMarked, pFiltered, ncol=1)
  
  return(pOut)
}
#PlotChrom(1)

tmpdirname <- paste0(name,"_Chromosomes")
dir.create(tmpdirname)
filesVec <- c()
for(myChrom in segDT[, unique(chromosome)]){
  pOut <- PlotChrom(myChrom)
  pOut <- annotate_figure(pOut, top = text_grob(paste0("chr", myChrom), color = "black", face = "bold", size = 12))

  tmpOut <- gsub(".pdf", paste0(".", myChrom, ".pdf"), Output)
  tmpOut <- paste0(tmpdirname, "/", tmpOut)
  
  filesVec <- c(filesVec, tmpOut)
  ggsave(tmpOut, pOut, width=16, height=9, device="pdf")
}

filesString <- paste0(filesVec, collapse = " ")
pdf_combine(filesVec, output = Output)

#unlink(tmpdirname, recursive = T)

# as PNG
for(myChrom in segDT[, unique(chromosome)]){
  
  pBase <- ggplot() +
    coord_cartesian(ylim=c(-ylim,ylim)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype=1, color="black") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = pretty_breaks(n = 10), labels=function(x)x/1000000)
  
  pFiltered <- pBase + 
    geom_point(data=plotcountsDT[chromosome == myChrom], aes(start, plotValue), size=0.001, alpha=0.6) +
    geom_segment(data=plotsegDT[chromosome == myChrom], aes(x=start, xend=end, y=plotValue, yend=plotValue), color="red", size=1.5)
  

  tmpOut <- gsub(".pdf", paste0(".", myChrom, ".png"), Output)
  tmpOut <- paste0(tmpdirname, "/", tmpOut)
  
  ggsave(tmpOut, pFiltered, width=16, height=9, device="png")
}




# full genome version
# convert genomic coordinates to a single x axis for plotting
Genomic2AxisCoordinates <- function(GenomeBins){
  GenomeBins[, start := as.numeric(start)]
  GenomeBins[, end := as.numeric(end)]
  Len <- GenomeBins[chromosome == 1, max(end)]
  LabelPos = Len/2
  Lentmp = Len
  chromosomes = as.character(GenomeBins[chromosome != 1, unique(chromosome)])
  GenomeBins[chromosome==1, plotstart := start]
  GenomeBins[chromosome==1, plotend := end]
  for(ichromosome in chromosomes){
    Lentmp <- Lentmp + GenomeBins[chromosome==ichromosome,max(end)]
    GenomeBins[chromosome==ichromosome, plotstart  := start+max(Len)]
    GenomeBins[chromosome==ichromosome, plotend  := end+max(Len)]
    labelpos = (max(GenomeBins[chromosome==ichromosome,"end"])/2) + (max(Len)/2)
    LabelPos = c(LabelPos,labelpos)
    Len = c(Len,Lentmp)
  }
  return(GenomeBins)
}

a <- plotcountsDT[, .(chromosome, start, end, plotValue,data="counts")]
b <- plotsegDT[, .(chromosome, start, end, plotValue,data="segs")]
genomeDT <- rbind(a,b)
genomeDT <- Genomic2AxisCoordinates(genomeDT)

plotcountsDT.genome <- genomeDT[data=="counts"]
plotsegDT.genome <- genomeDT[data=="segs"]

# get chromosome axis labels  
chromends <- genomeDT[, max(plotend), by=chromosome]$V1
chromends <- c(0, chromends)
chromnames <- genomeDT[, paste0("chr", unique(chromosome))]
chromLabelPos <- diff(chromends)/2
names(chromLabelPos) <- chromnames
chromLabelPos <- chromends+c(chromLabelPos, 0) # with pseuso 0 after last chromosome


pGenome <- ggplot() +
coord_cartesian(ylim=c(-ylim,ylim)) +
theme_classic() +
geom_hline(yintercept = 0, linetype=1, color="black") +
theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
axis.text.x = element_text(angle=45, hjust=1)) +
geom_point(data=plotcountsDT.genome, aes(plotstart, plotValue), shape=".", alpha=0.3) +
geom_segment(data=plotsegDT.genome, aes(x=plotstart, xend=plotend, y=plotValue, yend=plotValue), color="red", size=1.5) +
  scale_x_continuous(breaks = chromLabelPos, labels = names(chromLabelPos)) +
  geom_vline(xintercept = chromends, linetype=3, color="grey")

Output2 <- paste0(name, ".CNV.CNVKit.",ylim,".filtered.pdf")
ggsave(Output2,pGenome, width=16, height=9)

# as interactive segments with Plotly
pBase <- ggplot() +
  coord_cartesian(ylim=c(-ylim,ylim)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype=1, color="black") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks = pretty_breaks(n = 10), labels=function(x)x/1000000)

pFiltered <- pBase + 
  geom_segment(data=plotsegDT, aes(x=start, xend=end, y=plotValue, yend=plotValue), color="red", size=1.5) +
  facet_wrap(~chromosome, scales = "free")

fig <- ggplotly(pFiltered)

tmpOut <- gsub(".pdf", ".segments.html", Output)
tmpOut <- paste0(tmpdirname, "/", tmpOut)
htmlwidgets::saveWidget(as_widget(fig), tmpOut)








