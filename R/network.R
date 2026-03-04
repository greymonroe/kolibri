#' Find reads shared between two clusters
#'
#' @title Shared reads between two clusters
#' @description Returns the read names (or rows) that appear in both cluster
#'   \code{a} and cluster \code{b}.
#'
#' @param dt A \code{data.table} produced by one of the
#'   \code{cluster_reads_*()} functions.
#' @param a Identifier for the first cluster (character if
#'   \code{id_type = "chrom_cluster_id"}, integer if
#'   \code{id_type = "cluster_id"}).
#' @param b Identifier for the second cluster.
#' @param id_type Which cluster ID column to use: \code{"chrom_cluster_id"}
#'   (default) or \code{"cluster_id"}.
#' @param chrom_a Optional chromosome for cluster \code{a} (needed when
#'   \code{id_type = "cluster_id"} and IDs are not globally unique).
#' @param chrom_b Optional chromosome for cluster \code{b}.
#' @param read_col Name of the read name column. Default \code{"name"}.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param return Either \code{"reads"} (character vector of shared read names)
#'   or \code{"rows"} (rows from \code{dt} for those reads, tagged with
#'   \code{cluster_side}).
#'
#' @return A character vector (if \code{return = "reads"}) or a
#'   \code{data.table} (if \code{return = "rows"}).
#'
#' @examples
#' \dontrun{
#' shared <- shared_reads_between(clusters, a = "Chr1_12000", b = "Chr3_45000")
#' }
#'
#' @family network
#' @importFrom data.table is.data.table rbindlist
#' @export
shared_reads_between <- function(dt, a, b,
                                 id_type = c("chrom_cluster_id","cluster_id"),
                                 chrom_a = NULL, chrom_b = NULL,
                                 read_col = "name",
                                 chrom_col = "chrom",
                                 return = c("reads","rows")) {
  stopifnot(data.table::is.data.table(dt))
  id_type <- match.arg(id_type)
  return  <- match.arg(return)

  if (!all(read_col %in% names(dt))) stop("`read_col` not found in dt.")
  if (id_type == "chrom_cluster_id" && !"chrom_cluster_id" %in% names(dt))
    stop("`chrom_cluster_id` column not found.")
  if (id_type == "cluster_id" && !"cluster_id" %in% names(dt))
    stop("`cluster_id` column not found in dt.")

  if (id_type == "chrom_cluster_id") {
    A <- dt[get("chrom_cluster_id") == a]
    B <- dt[get("chrom_cluster_id") == b]
  } else {
    A <- dt[cluster_id == a]
    B <- dt[cluster_id == b]
    if (!is.null(chrom_a)) A <- A[get(chrom_col) == chrom_a]
    if (!is.null(chrom_b)) B <- B[get(chrom_col) == chrom_b]
    if (!nrow(A))
      stop(sprintf("No rows for cluster_id=%s%s.", as.character(a),
                   if (!is.null(chrom_a)) paste0(" (chrom=", chrom_a, ")") else ""))
    if (!nrow(B))
      stop(sprintf("No rows for cluster_id=%s%s.", as.character(b),
                   if (!is.null(chrom_b)) paste0(" (chrom=", chrom_b, ")") else ""))
    if (is.null(chrom_a) && length(unique(A[[chrom_col]])) > 1L)
      stop("cluster_id 'a' appears on multiple chromosomes; supply chrom_a.")
    if (is.null(chrom_b) && length(unique(B[[chrom_col]])) > 1L)
      stop("cluster_id 'b' appears on multiple chromosomes; supply chrom_b.")
  }

  reads_A <- unique(A[[read_col]])
  reads_B <- unique(B[[read_col]])
  shared  <- intersect(reads_A, reads_B)

  if (return == "reads") {
    return(shared)
  } else {
    A2 <- A[get(read_col) %in% shared][, cluster_side := "A"]
    B2 <- B[get(read_col) %in% shared][, cluster_side := "B"]
    return(data.table::rbindlist(list(A2, B2), use.names = TRUE, fill = TRUE))
  }
}


#' Build a shared-read cluster network
#'
#' @title Build cluster network
#' @description Takes a clustered read table and builds a network where nodes
#'   are clusters and edges connect clusters that share at least
#'   \code{min_shared} read names. An iterative filter ensures that the
#'   returned reads appear in exactly two of the kept clusters.
#'
#' @param dt A \code{data.table} from \code{cluster_reads_sweep()},
#'   \code{cluster_reads_dbscan()}, or \code{cluster_reads_kde()}.
#' @param cluster_col Name of the cluster label column. Default
#'   \code{"chrom_cluster_id"}.
#' @param read_col Name of the read name column. Default \code{"name"}.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param start_col Name of the start coordinate column. Default \code{"start"}.
#' @param end_col Name of the end coordinate column. Default \code{"end"}.
#' @param id_col Name of the numeric cluster ID column (used for noise
#'   filtering). Default \code{"cluster_id"}.
#' @param include_noise Logical. If \code{FALSE} (default), rows with
#'   \code{cluster_id == 0} are excluded.
#' @param restrict_within_chrom Logical. If \code{TRUE}, only edges between
#'   clusters on the same chromosome are kept.
#' @param min_shared Integer. Minimum number of shared reads for an edge.
#'   Default \code{1L}.
#' @param max_iter Integer. Maximum number of filtering iterations. Default
#'   \code{5L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{nodes}{A \code{data.table} with columns \code{cluster},
#'       \code{chrom}, \code{cluster_n}, \code{node_avg_score},
#'       \code{node_avg_len}, \code{min_start}, \code{max_end},
#'       \code{edge_count}, \code{total_length}.}
#'     \item{edges}{A \code{data.table} with columns \code{cluster1},
#'       \code{cluster2}, \code{shared_reads}, \code{shared_avg_score},
#'       \code{shared_avg_len}.}
#'     \item{reads}{Rows from \code{dt} for reads that appear in exactly two
#'       of the kept clusters.}
#'   }
#'
#' @examples
#' \dontrun{
#' clusters <- cluster_reads_sweep(split_reads, delta = 1000)
#' net <- build_cluster_network(clusters, min_shared = 2)
#' }
#'
#' @family network
#' @importFrom data.table is.data.table copy rbindlist setnames setkey first
#'   uniqueN
#' @export
build_cluster_network <- function(dt,
                                  cluster_col = "chrom_cluster_id",
                                  read_col    = "name",
                                  chrom_col   = "chrom",
                                  start_col   = "start",
                                  end_col     = "end",
                                  id_col      = "cluster_id",
                                  include_noise = FALSE,
                                  restrict_within_chrom = FALSE,
                                  min_shared = 1L,
                                  max_iter = 5L,
                                  verbose = TRUE) {
  stopifnot(data.table::is.data.table(dt))
  req_cols <- c(cluster_col, read_col, chrom_col, start_col, end_col,
                "score", "aligned_length")
  miss <- setdiff(req_cols, names(dt))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

  x <- data.table::copy(dt)

  if (!include_noise && id_col %in% names(x)) {
    x <- x[x[[id_col]] > 0L]
  }

  empty_result <- function() list(
    nodes = data.table::data.table(
      cluster = character(), chrom = character(),
      cluster_n = integer(), node_avg_score = numeric(),
      node_avg_len = numeric(), min_start = integer(), max_end = integer(),
      edge_count = integer(), total_length = integer()),
    edges = data.table::data.table(
      cluster1 = character(), cluster2 = character(),
      shared_reads = integer(), shared_avg_score = numeric(),
      shared_avg_len = numeric()),
    reads = x[0]
  )

  if (!nrow(x)) return(empty_result())

  build_cr <- function(xx) {
    unique(xx[, .(cluster = .SD[[1]], read = .SD[[2]], chrom = .SD[[3]]),
              .SDcols = c(cluster_col, read_col, chrom_col)])
  }

  compute_edges <- function(xx, crx) {
    read_global <- xx[, .(
      read_avg_score = mean(score, na.rm = TRUE),
      read_avg_len   = mean(aligned_length, na.rm = TRUE)
    ), by = read_col]
    data.table::setnames(read_global, read_col, "read")

    data.table::setkey(crx, read)
    pairs <- crx[crx, allow.cartesian = TRUE, nomatch = 0,
                 .(read, c1 = i.cluster, ch1 = i.chrom,
                   c2 = x.cluster, ch2 = x.chrom)]
    pairs <- pairs[c1 != c2]
    if (!nrow(pairs)) return(data.table::data.table(
      cluster1=character(), cluster2=character(),
      shared_reads=integer(), shared_avg_score=numeric(),
      shared_avg_len=numeric()))
    ord <- pairs$c1 > pairs$c2
    pairs[, `:=`(
      cluster1 = ifelse(ord, c2, c1),
      cluster2 = ifelse(ord, c1, c2),
      chrom1   = ifelse(ord, ch2, ch1),
      chrom2   = ifelse(ord, ch1, ch2)
    )][, c("c1","c2","ch1","ch2") := NULL]

    if (restrict_within_chrom) pairs <- pairs[chrom1 == chrom2]
    if (!nrow(pairs)) return(data.table::data.table(
      cluster1=character(), cluster2=character(),
      shared_reads=integer(), shared_avg_score=numeric(),
      shared_avg_len=numeric()))

    pairs <- read_global[pairs, on = .(read)]
    edges <- pairs[, .(
      shared_reads     = .N,
      shared_avg_score = mean(read_avg_score, na.rm = TRUE),
      shared_avg_len   = mean(read_avg_len,   na.rm = TRUE)
    ), by = .(cluster1, cluster2)]

    if (min_shared > 1L) edges <- edges[shared_reads >= min_shared]
    edges[order(-shared_reads, cluster1, cluster2)]
  }

  compute_nodes <- function(xx, keep_clusters) {
    if (!length(keep_clusters)) return(data.table::data.table(
      cluster=character(), chrom=character(), cluster_n=integer(),
      node_avg_score=numeric(), node_avg_len=numeric(),
      min_start=integer(), max_end=integer(), edge_count=integer(),
      total_length=integer()))
    yy <- xx[get(cluster_col) %chin% keep_clusters]

    cr_detail <- yy[, .(
      read_cluster_score = mean(score, na.rm = TRUE),
      read_cluster_len   = mean(aligned_length, na.rm = TRUE)
    ), by = .(cluster = yy[[cluster_col]], read = yy[[read_col]])]

    cl_sizes <- unique(yy[, .(cluster = .SD[[1]], read = .SD[[2]], chrom = .SD[[3]]),
                          .SDcols = c(cluster_col, read_col, chrom_col)])[
                            , .(cluster_n = .N, chrom = data.table::first(chrom)),
                            by = cluster]

    spans <- yy[, .(
      min_start = suppressWarnings(min(as.numeric(get(start_col)), na.rm = TRUE)),
      max_end   = suppressWarnings(max(as.numeric(get(end_col)),   na.rm = TRUE))
    ), by = .(cluster = get(cluster_col))]

    nodes <- cl_sizes[cr_detail, on = .(cluster), nomatch = 0][
      , .(
        chrom          = data.table::first(chrom),
        cluster_n      = data.table::first(cluster_n),
        node_avg_score = mean(read_cluster_score, na.rm = TRUE),
        node_avg_len   = mean(read_cluster_len,   na.rm = TRUE)
      ), by = cluster][spans, on = .(cluster)]

    nodes
  }

  cr0    <- build_cr(x)
  edges0 <- compute_edges(x, cr0)
  keep   <- unique(c(edges0$cluster1, edges0$cluster2))

  if (!length(keep)) {
    if (verbose) message("No edges after initial pass; returning empty network.")
    return(empty_result())
  }

  iter <- 1L
  read_set <- NULL
  repeat {
    x_keep <- x[get(cluster_col) %chin% keep]
    if (!nrow(x_keep)) break

    crk <- build_cr(x_keep)
    r_counts <- crk[, .(n_clusters = .N), by = read]
    reads2   <- r_counts[n_clusters == 2L, read]

    if (!length(reads2)) {
      if (verbose) message("No reads with exactly 2 clusters; returning empty.")
      return(empty_result())
    }

    xf <- x_keep[get(read_col) %chin% reads2]
    crf    <- build_cr(xf)
    edgesf <- compute_edges(xf, crf)
    keep2  <- unique(c(edgesf$cluster1, edgesf$cluster2))

    crf_keep <- crf[cluster %chin% keep2]
    r2_again <- crf_keep[, .(n_clusters = .N), by = read]
    reads2b  <- r2_again[n_clusters == 2L, read]
    xf2      <- xf[get(read_col) %chin% reads2b]

    if (setequal(keep2, keep) && setequal(reads2b, read_set)) {
      x_final     <- xf2
      edges_final <- edgesf
      keep_final  <- keep2
      break
    }

    keep    <- keep2
    read_set <- reads2b
    iter <- iter + 1L
    if (iter > max_iter) {
      if (verbose) message("Reached max_iter; using current filtered network.")
      x_final     <- xf2
      edges_final <- edgesf
      keep_final  <- keep2
      break
    }
  }

  if (!nrow(edges_final)) {
    if (verbose) message("No edges after filtering to reads in exactly two clusters.")
    return(empty_result())
  }

  nodes_final <- compute_nodes(x_final, keep_final)

  deg_tbl <- data.table::rbindlist(list(
    data.table::data.table(cluster = edges_final$cluster1),
    data.table::data.table(cluster = edges_final$cluster2)
  ))[, .(edge_count = .N), by = cluster]
  nodes_final <- deg_tbl[nodes_final, on = .(cluster)]
  nodes_final[is.na(edge_count), edge_count := 0L]
  nodes_final[, total_length := as.integer(max_end - min_start)]

  nodes_final <- nodes_final[order(-cluster_n, cluster)]
  edges_final <- edges_final[order(-shared_reads, cluster1, cluster2)]
  x_final     <- x_final[order(get(chrom_col), get(start_col), get(end_col))]

  if (verbose) {
    message(sprintf(
      "Network built in %d iteration(s). Kept %d clusters, %d edges, %d reads.",
      iter, length(keep_final), nrow(edges_final),
      data.table::uniqueN(x_final[[read_col]])))
  }

  list(nodes = nodes_final[], edges = edges_final[], reads = x_final[])
}
