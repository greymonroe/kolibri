test_that("read_bed12 reads a valid BED12 file", {
  path <- system.file("extdata", "alignVA2.bed12", package = "kolibri")
  skip_if(!file.exists(path), "Example BED12 file not found")

  bed <- read_bed12(path)

  expected_cols <- c("chrom","start","end","name","score","strand",
                     "thickStart","thickEnd","itemRgb","blockCount",
                     "blockSizes","blockStarts",
                     "is_split","aligned_length","nameN")
  expect_true(all(expected_cols %in% names(bed)))
  expect_s3_class(bed, "data.table")
  expect_true(nrow(bed) > 0)
})

test_that("read_bed12 adds is_split column correctly", {
  path <- system.file("extdata", "alignVA2.bed12", package = "kolibri")
  skip_if(!file.exists(path), "Example BED12 file not found")

  bed <- read_bed12(path)
  expect_type(bed$is_split, "logical")
  expect_true(all(bed$is_split == (bed$blockCount > 1L)))
})

test_that("read_bed12 adds aligned_length as end - start", {
  path <- system.file("extdata", "alignVA2.bed12", package = "kolibri")
  skip_if(!file.exists(path), "Example BED12 file not found")

  bed <- read_bed12(path)
  expect_true(all(bed$aligned_length == pmax(0L, bed$end - bed$start)))
})

test_that("read_bed12 errors on missing file", {
  expect_error(read_bed12("/nonexistent/path/file.bed12"))
})
