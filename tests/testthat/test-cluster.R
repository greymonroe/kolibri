library(data.table)

make_test_dt <- function() {
  data.table(
    chrom = c("Chr1","Chr1","Chr1","Chr2","Chr2"),
    start = c(100L, 200L, 5000L, 100L, 200L),
    end   = c(300L, 400L, 5200L, 300L, 400L),
    name  = c("r1","r2","r3","r4","r5"),
    score = rep(60L, 5),
    aligned_length = rep(200L, 5)
  )
}

test_that("cluster_reads_sweep adds required columns", {
  dt <- make_test_dt()
  res <- cluster_reads_sweep(dt, verbose = FALSE)

  expect_true(all(c("cluster_id","cluster_n","chrom_cluster_id") %in% names(res)))
})

test_that("cluster_reads_sweep clusters overlapping intervals", {
  dt <- make_test_dt()
  res <- cluster_reads_sweep(dt, delta = 0L, verbose = FALSE)

  # Chr1 rows 1 and 2 (100-300 and 200-400) should be in the same cluster
  chr1 <- res[chrom == "Chr1"]
  ids  <- chr1[start %in% c(100L, 200L), cluster_id]
  expect_equal(length(unique(ids)), 1L)
})

test_that("cluster_reads_sweep separates distant intervals", {
  dt <- make_test_dt()
  res <- cluster_reads_sweep(dt, delta = 0L, verbose = FALSE)

  # Chr1 rows 1/2 (100-400) and row 3 (5000-5200) should be in different clusters
  chr1   <- res[chrom == "Chr1"]
  near   <- unique(chr1[start < 1000, cluster_id])
  far    <- unique(chr1[start > 1000, cluster_id])
  expect_false(any(near %in% far))
})

test_that("cluster_reads_sweep cluster_n counts reads per cluster", {
  dt <- make_test_dt()
  res <- cluster_reads_sweep(dt, delta = 100L, verbose = FALSE)

  expect_true(all(res$cluster_n >= 1L))
  # The two nearby Chr1 reads should form a cluster of size 2
  chr1 <- res[chrom == "Chr1" & cluster_n > 0]
  expect_true(any(chr1$cluster_n == 2L))
})

test_that("cluster_reads_sweep chrom_cluster_id has chrom prefix", {
  dt <- make_test_dt()
  res <- cluster_reads_sweep(dt, verbose = FALSE)
  labeled <- res[cluster_id > 0]
  if (nrow(labeled)) {
    expect_true(all(grepl("^Chr", labeled$chrom_cluster_id)))
  }
})
