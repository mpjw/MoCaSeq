#' IO code for MoCaSeq
#'
#' Functions for building paths, and reading or writing data related to the
#' MoCaSeq pipeline

SAMPLE_TYPES <- c("Tumor", "matched", "Normal")

NUCLEOTIDE_VARIANT_TOOLS <- c("Mutect2", "Strelka")

COPY_NUMBER_CALLERS <- c("CNVKit", "Copywriter", "HMMCopy")

MOCASEQ_TOOLS <- c(NUCLEOTIDE_VARIANT_TOOLS, COPY_NUMBER_CALLERS, "LOH")

#' Build path to SNP file in MoCaSeq results
#'
#' MoCaSeq uses MuTect2 for mutation calling, this function creates the standard
#' path to the (postprocessed) SNP file given the sample name and sample type.
#' Please note, cancer genomics requires discerning germline and somatic
#' nucleotide variants. We refer to germline variants as single nucleotide
#' polymorphisms (SNPs) and somatic variants as single nucleotide variants
#' (SNVs).
#'
#' @param sample_name Character, name of sample
#' @param sample_type Character, one of "Tumor", "Normal", "matched" indicating
#'  type of sample
#' @param postprocessing Character, describing postprocessing applied to
#'  mutations E.g. "Positions", "NoCommonSNPs"
get_mocaseq_snv_file <- function(
  sample_name,
  sample_type,
  mutation_caller = "Mutect2",
  postprocessing = "NoCommonSNPs",
  pipeline_version = "bash"
) {
  if (pipeline_version == "bash") {
    switch(
      mutation_caller,
      Mutect2 = {
        mutect2_file_prefix <- switch(
          sample_type,
          Tumor = ,
          Normal = paste0(
            c(sample_name, sample_type, "Mutect2"),
            collapse = "."
          ),
          matched = paste0(c(sample_name, "Mutect2"), collapse = "."),
        )
        mutect2_file_suffix <- switch(
          postprocessing,
          Positions = ,
          NoCommonSNPs = paste0(c(postprocessing, "txt"), collapse = "."),
          OnlyImpact = "NoCommonSNPs.OnlyImpact.txt",
          TruSight = "NoCommonSNPs.OnlyImpact.TruSight.txt",
          CGC = "NoCommonSNPs.OnlyImpact.CGC.txt",
          "txt"
        )
        paste0(mutect2_file_prefix, ".", mutect2_file_suffix)
      },
      Strelka = {
        file.path("Strelka-Tumor", "results", "variants", "variants.vcf.gz")
      }
    )
  } else if (pipeline_version == "nextflow") {
    switch(
      mutation_caller,
      Mutect2 = {
        mutect2_file_prefix <- paste0(
          c(sample_name, sample_type, "Mutect2"),
          collapse = "."
        )
        mutect2_file_suffix <- switch(
          postprocessing,
          Positions = ,
          NoCommonSNPs = paste0(c(postprocessing, "txt"), collapse = "."),
          OnlyImpact = "NoCommonSNPs.OnlyImpact.txt",
          TruSight = "NoCommonSNPs.OnlyImpact.TruSight.txt",
          CGC = "NoCommonSNPs.OnlyImpact.CGC.txt",
          "txt"
        )
        paste0(mutect2_file_prefix, ".", mutect2_file_suffix)
      },
      Strelka = {
        file.path("Strelka-Tumor", "results", "variants", "variants.vcf.gz")
      }
    )
  } else {
    stop(paste0("Pipeline version '", pipeline_version, "' not recognized!"))
  }
}

#' Detect MoCaSeq pipeline version
#'
#' @param results_path Character path to results for a given sample. I.e. this
#' path should end with sample_id/results
#' @return Character 'bash', or 'nextflow' indicating MoCaSeq pipeline version
detect_mocaseq_version <- function(results_path) {
  bash_qc_exists <- dir.exists(file.path(results_path, "QC"))
  matched_files <- Sys.glob(file.path(results_path, "*", "*.matched.*"))
  if (bash_qc_exists || length(matched_files) == 0) {
    "bash"
  } else {
    "nextflow"
  }
}

#' Detect CNV caller in MoCaSeq results
#'
#' MoCaSeq employs CNVKit, HMMCopy and Copywriter for calling copy number
#' variations (CNVs). It stores any results under \code{sample_name/results/} in
#' a subfolder with the respective \code{tool_name}. This function exploits this
#' structure by searching all available tools for a CNV caller.
#'
#' @param results_path Character, path to results for a certain sample.
#' @return Character name of CNV caller one of "CNVKit", "HMMCopy" or
#'  "CopyWriter".
detect_mocaseq_cnv_caller <- function(results_path) {
  stopifnot(dir.exists(results_path))
  result_tools <- basename(list.dirs(results_path, recursive = FALSE))
  if ("CNVKit" %in% result_tools) {
    "CNVKit"
  } else if ("HMMCopy" %in% result_tools) {
    "HMMCopy"
  } else if ("Copywriter" %in% result_tools) {
    "Copywriter"
  } else {
    stop(paste0(
      "Cannot detect CNV caller from MoCaSeq output! ",
      "No CNV caller found in: ",
      results_path
    ))
  }
}

#' Build path to copy number segments file.
#'
#' Copy number segment (cns) files produced by MoCaSeq can be either from
#' CNVKit, HMMcopy, or CopyWriter. This function resolves possible cns files.
#'
#' @param sample_name Character, name of sample
#' @param sample_type Character, type of sample ("Tumor", "Normal", or
#'  "matched")
#' @param data_path Character, path to results (default: current working dir).
#' @param cnv_caller Character, name of CNV caller ("CNVKit", "HMMCopy" or
#'  "Copywriter")
#' @param result_type Character, type of CNV result ("segments", "log2RR")
#' @param seg_size Integer, segment size of HMMCopy output, only used in case
#'  of \code{cnv_caller = "HMMCopy"}.
get_mocaseq_cnv_file <- function(
  sample_name,
  sample_type,
  cnv_caller,
  result_type = "segments",
  seg_size = 20000
) {
  stopifnot(sample_type %in% c("matched", "Tumor", "Normal"))
  if (result_type == "segments") {
    switch(
      cnv_caller,
      CNVKit = file.path(sample_type, paste0(sample_name, ".cns")),
      HMMCopy = paste0(
        c(sample_name, "HMMCopy", seg_size, "segments", "txt"),
        collapse = "."
      ),
      Copywriter = paste0(
        c(sample_name, "Copywriter", "segments", "Mode", "txt"),
        collapse = "."
      ),
      stop("Unknown CNV caller: '", cnv_caller, "'")
    )
  } else if (result_type == "ratios") {
    switch(
      cnv_caller,
      CNVKit = file.path(sample_type, paste0(sample_name, ".cnr")),
      HMMCopy = paste0(
        c(sample_name, "HMMCopy", seg_size, "log2RR", "txt"),
        collapse = "."
      ),
      Copywriter = paste0(
        c(sample_name, "Copywriter", "log2RR", "Mode", "txt"),
        collapse = "."
      ),
      stop("Unknown CNV caller: '", cnv_caller, "'")
    )
  } else {
    stop("Unknown CNV result: '", result_type, "'")
  }
}

#' Build path to LOH variants from MoCaSeq
#'
#' MoCaSeq produces loss of heterozygosity data. This function returns the name
#' of the file containing these variants, given a sample name and which type of
#' variant desired (i.e. germline or somatic).
#'
#' @param sample_name Character sample name from MoCaSeq run
#' @param variant_type Character variant type one of "germline" or "somatic"
get_mocaseq_loh_file <- function(
  sample_name,
  result_type = "variants_for_LOH",
  variant_type = "somatic"
) {
  switch(
    result_type,
    segments = paste0(c(sample_name, "LOH", "Segments", "tsv"), collapse = "."),
    segments_flex = paste0(
      c(sample_name, "LOH", "Segments", "flex", "tsv"),
      collapse = "."
    ),
    variants = paste0(c(sample_name, "LOH", "Variants", "tsv"), collapse = "."),
    variants = paste0(
      c(sample_name, "LOH", "Variants", "flex", "tsv"),
      collapse = "."
    ),
    variants_for_LOH = {
      switch(
        variant_type,
        mixed = ,
        somatic = paste0(sample_name, ".VariantsForLOH.txt"),
        germline = paste0(sample_name, ".VariantsForLOHGermline.txt")
      )
    },
    stop("Unknown type of LOH result:, '", result_type, "'!")
  )
}


#' Construct path for MoCaSeq result file
#'
#' @param sample_name Character name of sample from MoCaSeq run
#' @param sample_type Character sample type, one of Tumor, matched or Normal
#' @param tool_name Character name of tool from MoCaSeq pipeline
#' @param base_path Character of path to directory with MoCaSeq output
#'
#' TODO: expand/recycle parameters other than sample_name if multiple sample
#' names are provided
#' @returns Character path to described file, or NULL
#' @export
get_mocaseq_path <- function(
  sample_name,
  sample_type,
  tool_name,
  base_path = ".",
  variant_type = "mixed",
  verbose = FALSE,
  ...
) {
  stopifnot(sample_type %in% SAMPLE_TYPES)
  stopifnot(variant_type %in% c("mixed", "germline", "somatic"))
  # paste0("Unknown sample type: '", sample_type, "'! Expected:", SAMPLE_TYPES)
  stopifnot(tool_name %in% MOCASEQ_TOOLS)
  results_path <- file.path(base_path, sample_name, "results")

  if (!dir.exists(results_path)) {
    warning("MoCaSeq results not found at: ", results_path)
    return(NULL)
  }

  if (!exists("pipeline_version") || pipeline_version == NULL) {
    pipeline_version <- detect_mocaseq_version(results_path)
    if (verbose) print(paste("detected pipeline version", pipeline_version))
  }

  file_name <- switch(
    tool_name,
    Mutect2 = get_mocaseq_snv_file(
      sample_name,
      sample_type,
      mutation_caller = tool_name,
      pipeline_version = pipeline_version,
      ...
    ),
    LOH = get_mocaseq_loh_file(sample_name, ...),
    CNVKit = ,
    HMMCopy = ,
    Copywriter = get_mocaseq_cnv_file(
      sample_name,
      sample_type,
      cnv_caller = tool_name,
      ...
    )
  )

  if (verbose) {
    print(paste("file name built:", file_name))
  }
  file_path <- file.path(results_path, tool_name, file_name)
  if (file.exists(file_path)) {
    file_path
  } else {
    warning(paste(file_path, "does not exist"))
  }
}

#' Parse available CNV caller results form MoCaSeq folder
get_mocaseq_cna_callers <- function(
  sample_id,
  result_dir_path = "."
) {
  sample_path <- file.path(result_dir_path, sample_id)
  if (!dir.exists(sample_path)) {
    stop(paste0("Cannot find MoCaSeq results at", sample_path))
  }

  check_cna_caller_exists <- function(cnv_caller) {
    #  & file.exists(get_mocaseq_cnv_file(cnv_caller))
    if (dir.exists(file.path(sample_path, "results", cnv_caller))) {
      cnv_caller
    }
  }

  unlist(lapply(X = COPY_NUMBER_CALLERS, FUN = check_cna_caller_exists))
}
