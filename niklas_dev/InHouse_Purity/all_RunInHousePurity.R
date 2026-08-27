suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(BubbleTree))

args <- commandArgs(TRUE)

name <- args[1]
species <- args[2]
repository_dir <- args[3]7
seqType <- args[4]

cnaMethod <- "Copywriter"
smoothCNA <- T

#species="Human"
#seqType="WES"
#repository_dir <- "/home/rad/packages/MoCaSeq/repository/"

source(paste0(repository_dir, "../niklas_dev/InHouse_Purity/Calculate_Purity_Custom.R"))

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

# ORIGINAL DATA
setwd("/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGReichert_hPDAC/done_MS/")
name <- "B25_NFM_IP_P52"

allnames <- list.dirs(recursive = F)
allnames <- basename(allnames)
allnames <- allnames[!basename(allnames) %in% c("mixed_bams", "raw", "ref", "temp")]

#for(name in allnames){

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# WORKFLOW FOR PURITY ANALYSIS #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# within mocaseq folder structure
outdir.loh <- paste0(name, "/results/LOH/")
outdir.cna <- paste0(name, "/results/", cnaMethod, "/")
outdir.purity <- paste0(name, "/results/InhousePurity/")


segfile.loh <- paste0(outdir.loh, "/",name, ".LOH.Segments.tsv")
segfile.cna <- paste0(outdir.cna, "/",name,"_CNA_",cnaMethod,"_smooth.tsv")
file.snv <- paste0(name,"/results/Mutect2/",name,".Mutect2.NoCommonSNPs.OnlyImpact.txt") # to correct mutation AF after predicting the purity

# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 0. INIT #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
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


# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 1. LOH segmentation #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
outdir.lohChroms <- paste0(outdir.loh, "/",name,"_Chromosomes/")

# this is now deprecated and will be removed soon, the more stringend cutoffs can be regarded default and there will always only be 1 version in the result folders
# res.loh <- Get_LOH_segments(name, seqType, species)
# fwrite(res.loh$lohSegs, paste0(outdir.loh, "/",name, ".LOH.Segments.tsv"), sep="\t", col.names = T)
# fwrite(res.loh$lohDT, paste0(outdir.loh, "/",name, ".LOH.Segments.Variants.tsv"), sep="\t", col.names = T)
# ggsave(paste0(outdir.loh, "/",name, ".LOH.LOH.Segments.pdf"), res.loh$plotCleanXContinous, width = 16, height = 9)
# ggsave(paste0(outdir.loh, "/",name, ".LOH.LOH.Segments.QC.pdf"), res.loh$plot, width = 16, height = 9)
# 
# if(!dir.exists(outdir.lohChroms)){
#   dir.create(outdir.lohChroms)
# }
# lapply(seq_along(res.loh$plotChromosomes),function(x) ggsave(filename=paste0(outdir.lohChroms, "/",name, ".Chr",x,".LOH.LOH.Segments.png"), plot=res.loh$plotChromosomes[[x]], width=16, height=9, bg="white"))

# also with more stringend parameters

res.loh <- Get_LOH_segments(name, seqType, species, minSegSize=5000000, minSNPs=100)
fwrite(res.loh$lohSegs, segfile.loh, sep="\t", col.names = T)
fwrite(res.loh$lohDT, paste0(outdir.loh, "/",name, ".LOH.Segments.Variants.tsv"), sep="\t", col.names = T)
ggsave(paste0(outdir.loh, "/",name, ".LOH.LOH.Segments.pdf"), res.loh$plotCleanXContinous, width = 16, height = 9)
ggsave(paste0(outdir.loh, "/",name, ".LOH.LOH.Segments.QC.pdf"), res.loh$plot, width = 16, height = 9)

lapply(seq_along(res.loh$plotChromosomes),function(x) ggsave(filename=paste0(outdir.lohChroms, "/",name, ".Chr",x,".LOH.LOH.Segments.png"), plot=res.loh$plotChromosomes[[x]], width=16, height=9, bg="white"))

# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 2. CNA smoothing #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
cnaRes <- LoadCNA(name, cnaMethod, smoothCNA=smoothCNA)
cnaSegs <- cnaRes$cnaSegs
fwrite(cnaSegs, segfile.cna, sep="\t")
ggsave(paste0(outdir.cna, "/",name,"_QC0_CNAsmoothing.png"), cnaRes$pCNAsmoothing, width = 16, height = 9)

# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 3. Purity estimation #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
lohSegs <- res.loh$lohSegs
lohSegs[, segID2 := paste0(Chrom,"-",segID)]
if(!dir.exists(outdir.purity)){dir.create(outdir.purity)}
res <- EstimatePurity(name, cnaSegs, lohSegs, outdir.purity, segmentsFrom = "CNA", cutoff.CNA = 0.1, plotQC=T, cutoff.LOH = 0.5)

purity.resfile <- paste0(outdir.purity, "/", name,"_results.txt")

#}


# run purity on pre-calculated data
if(F){
  cnaSegs <- fread(segfile.cna)
  lohSegs <- fread(segfile.loh)
  EstimatePurity(name, cnaSegs, lohSegs, outdir.purity, segmentsFrom = "CNA", cutoff.CNA = -0.1, plotQC=T, cutoff.LOH = 0.5, plotSteps=F)
}

# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 4. Correct SNVs, LOH and CNAs by purity #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# read data
purityDT <- fread(purity.resfile)
purityDT <- data.table(t(purityDT[1:3,2]))
names(purityDT) <- c("sample", "gCont", "gPur")
purityDT[, gCont := as.numeric(gCont)]
purityDT[, gPur := as.numeric(gPur)]

# purityDT[, sampleGroup := gsub("(hPDAC\\d*)_.*", "\\1", sample)]
# ggplot(purityDT, aes(gPur, factor(sample, levels = purityDT$sample))) +
#   geom_bar(stat="identity", aes(fill=factor(gPur))) +
#   scale_fill_viridis_d(direction = -1)

# correct data
resList <- ProcessSample(name, purityDT, LOHfile = segfile.loh, CNAfile = segfile.cna, SNVfile = file.snv)
SNVs <- resList$SNV
lohSegs <- resList$LOH
cnaSegs <- resList$CNA

# filter to relevant genes
keepGenes <- fread(paste0(name,"/results/Mutect2/",name,".Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt"))
keepGenes <- keepGenes[, unique(`ANN[*].GENE`)]



