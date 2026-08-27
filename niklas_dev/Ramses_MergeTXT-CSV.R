suppressPackageStartupMessages(library(data.table))

args = commandArgs(TRUE)
txtfile=args[1]
csvfile=args[2]
#txtfile <- "2365-3/results/Ramses/2365-3.Ramses.Mutations.Annotated.txt"
#csvfile <- "2365-3/results/Ramses/2365-3.Ramses.Mutations.csv"

txt <- fread(txtfile)
csv <- fread(csvfile)

txt[, id := paste(CHROM, POS, REF, ALT, sep="-")]
csv[, id := paste(Chromosome, Position, Ref_base, Mut_base, sep="-")]

out <- merge(txt, csv, all=T)
out[, `GEN[Tumor].AF` := round(Observed_Frequency, digits = 2)]

out[, id := NULL]

fwrite(out, txtfile, sep="\t", col.names = T)

out2 <- out[`ANN[*].IMPACT` %in% c("MODERATE", "HIGH")]
filteredtextfile <- gsub(".txt", ".OnlyImpact.txt", txtfile)
fwrite(out2, filteredtextfile, sep="\t", col.names = T)