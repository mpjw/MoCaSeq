#!/usr/bin/Rscript

##########################################################################################
##
## LOH_GenerateVariantTable.R
##
## Calculate datapoints needed for plotting LOH while transforming them to B-allele frequencies.
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
genome_fasta=args[2] #location of the reference genome fasta
repository_dir=args[3] #location of repository

source(paste(repository_dir,"/LOH_Library.R",sep=""))

#read input files
normal = read.table(paste(name,"/results/Mutect2/",name,".Normal.Mutect2.Positions.txt", sep=""), header=T, sep="\t")

#adjust column names
colnames(normal) = c("Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount", "MapQ", "BaseQ")

#filter reads
normal=LOH_FilterReads(normal)


normal[,"UniquePos"]=paste(normal[,"Chrom"], normal[,"Pos"], normal[,"Ref"], normal[,"Alt"], sep="_")
normal = normal[,c("UniquePos", "Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount")]
colnames(normal) = c("UniquePos", "Chrom", "Pos", "Ref", "Alt", "Normal_Freq", "Normal_RefCount", "Normal_AltCount")
normal[,"Chrom"]=as.character(normal[,"Chrom"])

#write out table used for plotting of germline variants
write.table(normal,file=paste(name,"/results/LOH/",name,".VariantsForLOHGermline_onlyNormal.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
