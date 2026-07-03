#' Genomic data and functions
#'
#' This script contains global genomic variables and genomic functions relevent
#' to the MoCaSeq pipeline.
#'
#' @name genome.R
NULL

# GRCh38p12 (hg38) chromosome lengths
HG38_AUTOMOSOME_SIZES <- c(
  "1" = 248956422,
  "2" = 242193529,
  "3" = 198295559,
  "4" = 190214555,
  "5" = 181538259,
  "6" = 170805979,
  "7" = 159345973,
  "8" = 145138636,
  "9" = 138394717,
  "10" = 133797422,
  "11" = 135086622,
  "12" = 133275309,
  "13" = 114364328,
  "14" = 107043718,
  "15" = 101991189,
  "16" = 90338345,
  "17" = 83257441,
  "18" = 80373285,
  "19" = 58617616,
  "20" = 64444167,
  "21" = 46709983,
  "22" = 50818468
)

HG38_SEX_CHROMOSOME_SIZES <- c(
  "X" = 156040895,
  "Y" = 57227415
)

# GRCm38p6 (mm10) chromosome lengths
MM10_AUTOMOSOME_SIZES <- c(
  "1" = 195471971,
  "2" = 182113224,
  "3" = 160039680,
  "4" = 156508116,
  "5" = 151834684,
  "6" = 149736546,
  "7" = 145441459,
  "8" = 129401213,
  "9" = 124595110,
  "10" = 130694993,
  "11" = 122082543,
  "12" = 120129022,
  "13" = 120421639,
  "14" = 124902244,
  "15" = 104043685,
  "16" = 98207768,
  "17" = 94987271,
  "18" = 90702639,
  "19" = 61431566
)

MM10_SEX_CHROMOSOME_SIZES <- c(
  "X" = 171031299,
  "Y" = 91744698
)

#' Define chromosome sizes based on species.
#'
#' @param species Character name of spieces one of "Human" or "Mouse"
#' @param only_autosomes Logical if to include only autosomes [default: TRUE]
define_chromosome_sizes <- function(species = "Human", only_autosomes = TRUE) {
  switch(
    species,
    "Human" = if (only_autosomes) {
      HG38_AUTOMOSOME_SIZES
    } else {
      c(HG38_AUTOMOSOME_SIZES, HG38_SEX_CHROMOSOME_SIZES)
    },
    "Mouse" = if (only_autosomes) {
      MM10_AUTOMOSOME_SIZES
    } else {
      c(MM10_AUTOMOSOME_SIZES, MM10_SEX_CHROMOSOME_SIZES)
    }
  )
}

#' Process count data from LOH
#'
#' @param chromosome_sizes
#' @param count_data Character path to count data
#' @param method Character
#' @note version 2 uses data.table::fread instead of read.table
process_count_data <- function(
  chromosome_sizes,
  count_data = "",
  method = ""
) {
  out_list <- list()
  count_data <- data.table::fread(count_data, header = TRUE, sep = "\t")
  SetVariableNames(method)
  out_list <- ConvertGenomicCords(
    data.frame(count_data),
    chromosome_sizes,
    Start,
    CopyNumber,
    Chromosome
  )
  return(out_list)
}

#' Convert genomic coordinates to a single continuous x axis for plotting
#'
#' @param dt_bins data.table holding genomic bins with columns
#' * Chrom
#' * start
#' * end
convert_genomic_to_continuous_axis <- function(dt_bins) {
  dt_bins[, Chrom := as.numeric(Chrom)]
  setorder(dt_bins, Chrom)
  dt_bins[, start := as.numeric(start)]
  dt_bins[, end := as.numeric(end)]
  # TODO: try leverage base::cumsum here
  Len <- dt_bins[Chrom == 1, max(end)]
  LabelPos = Len / 2
  Lentmp = Len
  chromosomes = as.character(dt_bins[Chrom != 1, unique(Chrom)])
  dt_bins[Chrom == 1, plotstart := start]
  dt_bins[Chrom == 1, plotend := end]
  for (chromosome in chromosomes) {
    Lentmp <- Lentmp + dt_bins[Chrom == chromosome, max(end)]
    dt_bins[Chrom == chromosome, plotstart := start + max(Len)]
    dt_bins[Chrom == chromosome, plotend := end + max(Len)]
    labelpos = (max(dt_bins[Chrom == chromosome, "end"]) / 2) +
      (max(Len) / 2)
    LabelPos = c(LabelPos, labelpos)
    Len = c(Len, Lentmp)
  }
  return(dt_bins)
}


# QC: given a bam file and a data.table with mutations (chrom, pos, pos), will return MAPQ and BAQ values for the mutation position
calculate_mapq_baq <- function(inputbam, inputSNV, tempfolder) {
  # tempfolder <- "/home/rad/Downloads/temp4/" # for samtools and jvarkit calculations
  # inputbam <- copy(tbam)
  # inputSNV <- copy(SNVs)

  # first check if index was found, else try to copy or stop (we could execute samtools index here as well)
  indexname1 <- paste0(inputbam, ".bai") # this was what we need
  indexname2 <- gsub("\\.bam$", ".bai", inputbam) # sometimes this is the filename
  if (!file.exists(indexname1)) {
    if (file.exists(indexname2)) {
      file.copy(indexname2, indexname1)
    } else {
      stop("No index file for BAM found, please rerun samtools index first")
    }
  }

  inputbed <- unique(inputSNV[, .(CHROM, POS, POS)])

  outbam.file <- paste0(tempfolder, "/", name, "_snv-regions.bam")
  regionString <- paste0(
    inputbed[, paste0(CHROM, ":", POS, "-", POS)],
    collapse = " "
  ) # like this, since -L is slow and gives weird results
  command1 <- paste0(
    "samtools view -bh -F 4 -F 1024 ",
    inputbam,
    " ",
    regionString,
    " > ",
    outbam.file
  )
  system(command1)

  outbam.file2 <- paste0(tempfolder, "/", name, "_snv-regions-sorted.bam")
  command2 <- paste0("samtools sort -f ", outbam.file, " ", outbam.file2)
  system(command2)

  command3 <- paste0("samtools index ", outbam.file2)
  system(command3)

  outstats.file <- paste0(tempfolder, "/", name, "_snv-regions.stats")
  command4 <- paste0(
    "java -jar /home/rad/packages/jvarkit/dist/sam2tsv.jar -R /fast/GRCh38.p12/GRCh38.p12.fna ",
    outbam.file,
    " > ",
    outstats.file
  )
  system(command4)

  # read and process the QC data
  readstats <- fread(outstats.file)

  file.remove(
    outbam.file,
    outbam.file2,
    paste0(outbam.file2, ".bai"),
    outstats.file
  )

  colnames(readstats) <- gsub("-", "", colnames(readstats))
  setnames(readstats, "#ReadName", "ReadName")
  readstats <- unique(readstats)
  readstats[, baseID := .I]
  readstats[READQUAL != ".", READQUALNUM := utf8ToInt(READQUAL), by = baseID]
  readstats[, READQUALNUM := READQUALNUM - 33]
  readstats[, mutID2 := paste(CHROM, REFPOS1, REFBASE, READBASE, sep = "-")]
  readstats[, CHROM := as.character(CHROM)]

  # assign mutation positions
  subs.positions <- inputSNV[type %in% c("SNV", "MNV"), POS]
  indel.positions <- inputSNV[type %in% c("del", "ins", "dup"), POS]

  # only get MAPQ for these mutations
  if (length(indel.positions) >= 1) {
    readstats.indels <- readstats[REFPOS1 %in% indel.positions]
    readstats.indels <- unique(readstats.indels[,
      .(POS = as.numeric(REFPOS1), MAPQ = mean(MAPQ, na.rm = T), BAQ = NA),
      by = .(CHROM, REFPOS1)
    ])
    readstats.indels[, REFPOS1 := NULL]
  } else {
    readstats.indels <- data.table()
  }

  # only single/multi nucleotide variants
  readstats <- readstats[REFPOS1 %in% subs.positions]

  # remove all reads with REF=ALT
  readstats <- readstats[READBASE != REFBASE]

  inputSNV[,
    mutID2 := paste(CHROM, POS, substr(REF, 1, 1), substr(ALT, 1, 1), sep = "-")
  ] # this way we include multi substitutions and insertions
  readstats <- merge(
    readstats,
    inputSNV[, .(mutID, mutID2, REF)],
    by = "mutID2",
    all.x = T
  )

  # here we remove mutations which are at the given locations but show a different base substitution (e.g. C>T but there is also a C>A which is removed)
  readstats <- readstats[!is.na(REF)]

  # reduce to single information per mutation
  readstats <- unique(readstats[,
    .(MAPQ = mean(MAPQ, na.rm = T), BAQ = mean(READQUALNUM, na.rm = T)),
    by = .(mutID)
  ])

  # first bind the simple mutations onto the output using mutation ID, then indels only based on CHROM/POS
  readstats <- merge(
    inputSNV[POS %in% subs.positions],
    readstats,
    by = c("mutID")
  )
  readstats.indels <- merge(
    inputSNV[POS %in% indel.positions],
    readstats.indels,
    by = c("CHROM", "POS"),
    all.x = T
  )

  # merging by CHROM/POS is dangerous, since there can be a SNV at the same position as an indel, we catch it like this
  readstats.indels <- readstats.indels[!grepl(">", HGVS_C)]

  outputSNV <- rbind(readstats, readstats.indels)

  # some mutations are only found in the m2 bam from Mutect2, here we flag these

  errorSNV <- inputSNV[!inputSNV$mutID %in% outputSNV$mutID]
  errorSNV[, MAPQ := 0]
  errorSNV[, BAQ := 0]
  outputSNV <- rbind(outputSNV, errorSNV)

  if (nrow(inputSNV) != nrow(outputSNV)) {
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
