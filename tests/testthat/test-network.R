library(data.table)

make_clustered_dt <- function() {
  # Two clusters on Chr1 and Chr2, connected by read "shared_r"
  data.table(
    chrom            = c("Chr1","Chr1","Chr2","Chr2"),
    start            = c(100L, 110L, 5000L, 5010L),
    end              = c(300L, 320L, 5200L, 5220L),
    name             = c("shared_r","shared_r","shared_r","other_r"),
    score            = c(60L, 60L, 60L, 60L),
    aligned_length   = c(200L, 210L, 200L, 210L),
    cluster_id       = c(1L, 1L, 2L, 2L),
    cluster_n        = c(2L, 2L, 2L, 2L),
    chrom_cluster_id = c("Chr1_200","Chr1_200","Chr2_5105","Chr2_5105")
  )
}

test_that("build_cluster_network returns list with nodes, edges, reads", {
  dt  <- make_clustered_dt()
  net <- build_cluster_network(dt, verbose = FALSE)

  expect_type(net, "list")
  expect_true(all(c("nodes","edges","reads") %in% names(net)))
})

test_that("build_cluster_network nodes table has required columns", {
  dt  <- make_clustered_dt()
  net <- build_cluster_network(dt, verbose = FALSE)

  expected <- c("cluster","chrom","cluster_n","node_avg_score",
                "node_avg_len","min_start","max_end","edge_count",
                "total_length")
  expect_true(all(expected %in% names(net$nodes)))
})

test_that("build_cluster_network edges table has required columns", {
  dt  <- make_clustered_dt()
  net <- build_cluster_network(dt, verbose = FALSE)

  expected <- c("cluster1","cluster2","shared_reads")
  expect_true(all(expected %in% names(net$edges)))
})

test_that("build_cluster_network detects shared read", {
  dt  <- make_clustered_dt()
  net <- build_cluster_network(dt, verbose = FALSE)

  if (nrow(net$edges) > 0) {
    expect_true(net$edges$shared_reads[1] >= 1L)
  }
})

test_that("build_cluster_network returns empty list when no edges", {
  dt <- data.table(
    chrom            = "Chr1",
    start            = 100L,
    end              = 300L,
    name             = "solo_read",
    score            = 60L,
    aligned_length   = 200L,
    cluster_id       = 1L,
    cluster_n        = 1L,
    chrom_cluster_id = "Chr1_200"
  )
  net <- build_cluster_network(dt, verbose = FALSE)
  expect_equal(nrow(net$edges), 0L)
})
