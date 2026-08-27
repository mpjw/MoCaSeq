library(data.table)
library(splitstackshape)
args = commandArgs(trailingOnly=TRUE)

InputFile <- args[1]
#InputFile <- "AGKLTM_T1/results/rescued/AGKLTM_T1.StrelkaRNA.AGKLTM_T1.StrelkaRNA_BODY.vcf"
OutFile <- gsub(".vcf$", ".scientificNotationFix.vcf", InputFile)

dt <- fread(InputFile)

dtSplit <- cSplit(dt, "Tumor", ":")
cols <- colnames(dtSplit)[grep("Tumor_", colnames(dtSplit))] # get colnames to iterate (=subcols)
dtSplit[, (cols) := lapply(.SD, format, scientific=F), .SDcols = cols] # remove scientific notation for all subcols
dtSplit[, Tumor := gsub(" ", "",do.call(paste, c(dtSplit[,cols, with=F], sep=":")))] # merge the subcols back to original 
dtSplit[, (cols) := NULL] # remove sub columns
dt <- copy(dtSplit)

dtSplit <- cSplit(dt, "Normal", ":")
cols <- colnames(dtSplit)[grep("Normal_", colnames(dtSplit))] # get colnames to iterate (=subcols)
dtSplit[, (cols) := lapply(.SD, format, scientific=F), .SDcols = cols] # remove scientific notation for all subcols
dtSplit[, Normal := gsub(" ", "",do.call(paste, c(dtSplit[,cols, with=F], sep=":")))] # merge the subcols back to original 
dtSplit[, (cols) := NULL] # remove sub columns
dt <- copy(dtSplit)

fwrite(dt, OutFile, sep="\t", col.names = T)
   
