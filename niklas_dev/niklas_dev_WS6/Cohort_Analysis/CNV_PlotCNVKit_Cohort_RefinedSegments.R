#!/usr/bin/Rscript

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
species=args[2]
repository_dir=args[3]  #location of repository
types=args[4]
workdir=args[5]

# debug
# name <- "P4GNJ9.Tumor.antitargetcoverage.cnn"
# species <- "Human"
# repository_dir <- "/home/rad/packages/MoCaSeq/niklas_dev/Cohort_Analysis/"
# types <- ""
# workdir <- "CNVKit_refine-segmentsTest/cbs_1e3/"

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

baseWD <- getwd()
#setwd("/run/user/1000/gvfs/sftp:host=172.21.251.53,user=rad/media/rad/HDD1/hMANEC_combined_results/")

# set the first directory
setwd(workdir)

chrom.sizes = DefineChromSizes(species)
if (species=="Human"){
  chromosomes=22
} else if (species=="Mouse"){
  chromosomes=19
}

normalization=""

PlotFunction <- function(SegmentsFile, CountsFile, SampleName){
  ylim <- 0
  for (y_axis in c("CNV_5","CNV_2")){
    Segments = ProcessSegmentData(segmentdata=SegmentsFile,chrom.sizes,method="CNVKit")
    Counts = ProcessCountData(countdata=CountsFile,chrom.sizes,method="CNVKit")
    
    plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=SampleName,method="CNV",toolname="CNVKit",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")
    
    for (i in 1:chromosomes){
      plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=SampleName,chromosome=i,method="CNV",toolname="CNVKit",normalization=normalization,y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
    }
    
    command <- paste("pdfunite ",SampleName,"_Chromosomes/",SampleName,".Chr?.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf ",SampleName,"_Chromosomes/",SampleName,".Chr??.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf ",SampleName,".Chromosomes.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf",sep="")
    system(command)
  }
}

system(paste("mkdir -p ",name,"_Chromosomes",sep=""))

# first run the matched sample
SegmentsFile = paste(name,".Tumor.center-mode.cns",sep="")
CountsFile = paste(name,".Tumor.center-mode.cnr",sep="")

if(!file.exists(SegmentsFile)){
  stop(paste0("Segments file not found: ", SegmentsFile))
}
if(!file.exists(CountsFile)){
  stop(paste0("Segments file not found: ", CountsFile))
}

PlotFunction(SegmentsFile, CountsFile, name)

system(paste("rm -r ",name,"_Chromosomes",sep=""))












