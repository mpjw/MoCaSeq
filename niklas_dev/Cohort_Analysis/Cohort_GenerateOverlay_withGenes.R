ProcessSegmentData2 = function(segmentdata="",chrom.sizes=chrom.sizes,method="")
{
  Outlist = list()
  FirstPosition = c()
  if(method!="aCGH")
  {
    segmentdata = read.table(segmentdata,header=T,sep="\t")
  }
  SetVariableNamesSegments(method)
  cnSeg <- data.frame()
  borders <- c()
  last = 0
  for(i in names(chrom.sizes))
  {
    cur <- segmentdata[segmentdata[,SegmentChromosome]==i,c(SegmentChromosome,SegmentStart,SegmentStop,SegmentMean)]
    #cur[cur[,SegmentMean]<=(-5),SegmentMean]=-4.9
    #cur[cur[,SegmentMean]>=(5),SegmentMean]=4.9
    cnSeg <- rbind(cnSeg,data.frame(Chromosome=cur[,SegmentChromosome],start=cur[,SegmentStart]+last,
                                    stop=cur[,SegmentStop]+last,copy=cur[,SegmentMean]))
    borders <- c(borders,last)
    last = last + chrom.sizes[i]
  }
  Outlist[["CN"]] = cnSeg
  Outlist[["NonProcessed"]] = segmentdata
  return(Outlist)
}


setwd("/run/user/1000/gvfs/sftp:host=imows3.med.tum.de,user=rad/media/rad/HDD1/hMANEC/results_selected/CNVKit/raw/")
name="1N7T85"
species="Human"
repository_dir="/home/rad/packages/MoCaSeq/repository/"
CNSfile <- paste0(name,".Tumor.cns")
CNRfile <- paste0(name,".Tumor.cnr")
source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

chrom.sizes = DefineChromSizes(species)

if (species=="Human"){
  chromosomes=22
} else if (species=="Mouse"){
  chromosomes=19
}

normalization=""
y_axis <- "CNV_5"
ylim <- 5

Segments = CNSfile
Counts = CNRfile
Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="CNVKit")
Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="CNVKit")

plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="CNVKit",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3)

  

Segments2 = ProcessSegmentData2(segmentdata=CNSfile,chrom.sizes,method="CNVKit")[[1]]
myylim <- 5
hits <- Segments2[!Segments2$copy %between% c(-myylim, myylim),]
hits <- data.table(hits)

hits[copy < (-myylim), ypos := -5]
hits[copy > (myylim), ypos := 5]

segDT <- fread(CNSfile)
tmp <- segDT[log2 %in% hits$copy]
tmp[gene == "-", gene := "none"]

hits <- merge(hits, tmp[, .(log2, gene)], by.x="copy", by.y="log2")
hits[, copy := round(copy, digits = 1)]
hits[, label := paste0(gene, " (", copy, ")")]



#tmp <- segDT[log2 < -5]



plotGlobalRatioProfile2 = function(cn=cn,ChromBorders=ChromBorders,cnSeg="",samplename="",method="",toolname="",normalization="",y_axis="",Transparency=20, Cex=0.2,outformat="", ylim <- c(-5,5), outfolder="")
{
  SetYAxis(y_axis)
  #ylim <- c(-10,10)
  
  # ylow <- ylim[1]
  # yhigh <- ylim[2]
  # ylim <- c(ylow,yhigh)
  # ypos <- seq(ylow,yhigh,1)
  # ylabels <- seq(ylow,yhigh,1)
  # y_output <- paste0(".",yhigh,".")
  # yborder <- yhigh
  
  
  SetVariableNamesSegments(toolname)
  SetVariableNames(toolname)
  ChromNamePos = FindCordsChromNames(ChromBorders)
  ChromBorders=c(0,ChromBorders)
  YaxisPosition = min(ChromBorders)-130000000
  LastEntry = length(ChromBorders)
  if(normalization != "")
  {
    normalization=paste(".",normalization,sep="")
  }
  if(outformat=="emf")
  {
    emf(paste(samplename,".",method,".",toolname,normalization,y_output,"emf",sep=""), bg="white", width=25,height=10)
  }
  if(outformat=="pdf")
  {
    pdf(paste(samplename,".",method,".",toolname,normalization,y_output,"pdf",sep=""),width=25,height=10)
  }
  pdf("~/Downloads/hMANEC_CNVKit_CNV_withGenes.pdf",width=25,height=10)
  par(mar=c(4, 4, 0, 0))
  #if(toolname=="CNVKit"){cn <- cn[c(TRUE,rep(FALSE,9)), ]}
  plot(cn$position,cn$copy,pch=".",
       ylim=ylim,xaxt="n",yaxt="n",bty="n",col=paste("#000000",Transparency,sep=""),yaxs="i",ylab="",xlab="",yaxt="n")
  segments(min(ChromBorders),0,ChromBorders[LastEntry],col="#A9A9A9",lwd=2)
  axis(2,las=1,pos=YaxisPosition, outer=T, at=ypos,labels=ylabels)
  axis(1,las=1, labels=rep("",length(ChromBorders)),at=ChromBorders)
  mtext(side=1,line=1,at=ChromNamePos,names(ChromNamePos))
  if(method=="CNV")
  {
    segments(cnSeg$start,cnSeg$copy,cnSeg$stop,cnSeg$copy,col="#FF4D4D",lwd=2)
  }
  ChromBordersReduced = ChromBorders[-c(1,LastEntry)]
  if(method == "CNV")
  {
    segments(ChromBordersReduced,-yborder,ChromBordersReduced,yborder,lty=3,col="grey40",lwd=0.9)
  }
  if(method =="LOH" || method=="LOH_Germline" || method=="LOH_raw")
  {
    segments(ChromBordersReduced,0,ChromBordersReduced,1,lty=3,col="grey40",lwd=0.9)
  }
  
  text(hits$start, hits$ypos, hits$label, cex=0.8)
  
  garbage <- dev.off()
}