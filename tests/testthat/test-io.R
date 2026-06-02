library(testthat)

source("../../R/io.R")

test_that("get_mocaseq_cnv_file resolves paths correctly", {
  sample_id <- "PCSI_0509_Pa_P"

  # copy number segment files
  result_type <- "segments"
  expect_equal(
    get_mocaseq_cnv_file(sample_id, "matched", "HMMCopy", result_type),
    paste0(sample_id, ".HMMCopy.20000.segments.txt")
  )

  expect_equal(
    get_mocaseq_cnv_file(sample_id, "matched", "Copywriter", result_type),
    paste0(sample_id, ".Copywriter.segments.Mode.txt")
  )

  expect_equal(
    get_mocaseq_cnv_file(sample_id, "Tumor", "CNVKit", result_type),
    file.path("Tumor", paste0(sample_id, ".cns"))
  )

  expect_error(
    get_mocaseq_cnv_file(sample_id, "foo", "Copywriter", result_type)
  )

  # copy number ratio files
  result_type <- "ratios"
  expect_equal(
    get_mocaseq_cnv_file(sample_id, "matched", "HMMCopy", result_type),
    paste0(sample_id, ".HMMCopy.20000.log2RR.txt")
  )

  expect_equal(
    get_mocaseq_cnv_file(sample_id, "Normal", "Copywriter", result_type),
    paste0(sample_id, ".Copywriter.log2RR.Mode.txt")
  )

  expect_equal(
    get_mocaseq_cnv_file(sample_id, "Tumor", "CNVKit", result_type),
    file.path("Tumor", paste0(sample_id, ".cnr"))
  )

  expect_error(
    get_mocaseq_cnv_file(sample_id, "foo", "CNVKit", result_type)
  )
})
