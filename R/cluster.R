# ---- Internal helpers --------------------------------------------------------

#' Split clusters by maximum consecutive gap
#'
#' @description After an initial cluster assignment, re-splits any cluster whose
#'   consecutive sorted intervals are separated by more than \code{max_gap} bp.
#'
#' @param dt A \code{data.table} with at least \code{id_col}, \code{start}, and
#'   \code{end} columns.
#' @param id_col Name of the cluster ID column. Default \code{"cluster_id"}.
#' @param start Name of the start coordinate column. Default \code{"start"}.
#' @param end Name of the end coordinate column. Default \code{"end"}.
#' @param max_gap Maximum allowed consecutive gap in bp. \code{Inf} (default)
#'   disables splitting.
#' @param chrom_col Optional chromosome column name for informative messages.
#' @param verbose Logical. Emit messages for each split. Default \code{TRUE}.
#' @param return_notes Logical. If \code{TRUE}, return a list with elements
#'   \code{dt} and \code{split_notes}.
#'
#' @return The modified \code{data.table} with updated \code{id_col} values, or
#'   a named list if \code{return_notes = TRUE}.
#'
#' @keywords internal
#' @noRd
.split_by_max_gap <- function(dt,
                              id_col   = "cluster_id",
                              start    = "start",
                              end      = "end",
                              max_gap  = Inf,
                              chrom_col = NULL,
                              verbose  = TRUE,
                              return_notes = FALSE) {
  if (!is.finite(max_gap)) {
    out <- data.table::copy(dt)
    attr(out, "split_notes") <- data.table::data.table()
    return(out)
  }
  x <- data.table::copy(dt)
  n <- nrow(x)
  new_ids <- integer(n)
  next_id <- 1L

  s   <- suppressWarnings(as.numeric(x[[start]]))
  e   <- suppressWarnings(as.numeric(x[[end]]))
  cid <- x[[id_col]]

  notes <- list()

  for (g in sort(unique(cid[cid > 0L]))) {
    idx  <- which(cid == g)
    ord  <- order(s[idx], e[idx], na.last = TRUE)
    idxo <- idx[ord]

    prev_end   <- -Inf
    cur_global <- NA_integer_
    prev_subid <- NA_integer_

    for (ii in seq_along(idxo)) {
      si <- s[idxo[ii]]; ei <- e[idxo[ii]]
      if (is.na(si) || is.na(ei)) next

      if (is.na(cur_global) || si > prev_end + max_gap) {
        if (!is.na(cur_global)) {
          gap_bp <- as.numeric(si - prev_end)
          note <- data.table::data.table(
            original_cluster    = g,
            new_subcluster      = next_id,
            previous_subcluster = prev_subid,
            split_at_start      = si,
            previous_end        = prev_end,
            gap_bp              = gap_bp,
            threshold_bp        = max_gap,
            chrom               = if (!is.null(chrom_col) && chrom_col %in% names(x))
              x[[chrom_col]][idxo[ii]] else NA_character_
          )
          notes[[length(notes) + 1L]] <- note
          if (isTRUE(verbose)) {
            msg_chr <- if (!is.null(chrom_col) && chrom_col %in% names(x))
              paste0(" on ", note$chrom) else ""
            message(sprintf("Split cluster %s%s: gap %d bp > max_gap %d; new subcluster %d",
                            g, msg_chr, gap_bp, max_gap, next_id))
          }
        }
        cur_global <- next_id
        prev_subid <- cur_global
        next_id    <- next_id + 1L
      }

      new_ids[idxo[ii]] <- cur_global
      prev_end <- max(prev_end, ei, na.rm = TRUE)
    }
  }

  out <- data.table::copy(x)
  out[[id_col]] <- ifelse(new_ids > 0L, new_ids, 0L)

  notes_dt <- if (length(notes)) data.table::rbindlist(notes) else data.table::data.table()
  attr(out, "split_notes") <- notes_dt

  if (isTRUE(return_notes)) {
    return(list(dt = out, split_notes = notes_dt))
  }
  out
}


#' Finalize cluster assignments
#'
#' @description Adds \code{cluster_n} (number of reads per cluster per
#'   chromosome) and optionally warns about clusters that are closer than a
#'   given distance.
#'
#' @param dt A \code{data.table} with \code{cluster_id} and \code{chrom_col}
#'   columns.
#' @param start Name of start column. Default \code{"start"}.
#' @param end Name of end column. Default \code{"end"}.
#' @param chrom_col Name of chromosome column. Default \code{"chrom"}.
#' @param warn_if_clusters_closer_than Numeric. Emit warnings when cluster
#'   centres are closer than this many bp. \code{NA} disables.
#' @param verbose Logical. Print warning table. Default \code{TRUE}.
#'
#' @return The input \code{data.table} with \code{cluster_n} added.
#'
#' @keywords internal
#' @noRd
.finalize_clusters <- function(dt,
                               start = "start",
                               end = "end",
                               chrom_col = "chrom",
                               warn_if_clusters_closer_than = NA_real_,
                               verbose = TRUE) {
  stopifnot(data.table::is.data.table(dt))
  x <- data.table::copy(dt)

  sizes <- x[cluster_id > 0L, .N, by = c(chrom_col, "cluster_id")]
  x[, cluster_n := 0L]
  if (nrow(sizes)) {
    x[sizes, on = c(chrom_col, "cluster_id"), cluster_n := i.N]
  }

  warn_list <- list()
  if (is.finite(warn_if_clusters_closer_than)) {
    chroms <- unique(x[[chrom_col]])
    for (cc in chroms) {
      sub <- x[(x[[chrom_col]] == cc) & cluster_id > 0L]
      if (!nrow(sub)) next

      centers <- sub[, {
        s <- suppressWarnings(as.numeric(.SD[[1]]))
        e <- suppressWarnings(as.numeric(.SD[[2]]))
        mid <- (s + e) / 2
        .(center = stats::median(mid, na.rm = TRUE))
      }, by = cluster_id, .SDcols = c(start, end)][is.finite(center)][order(center)]

      if (nrow(centers) >= 2L) {
        d <- diff(centers$center)
        hits <- which(is.finite(d) & d < warn_if_clusters_closer_than)
        if (length(hits)) {
          warn_list[[as.character(cc)]] <- data.table::data.table(
            chrom     = cc,
            cluster_a = centers$cluster_id[hits],
            cluster_b = centers$cluster_id[hits + 1L],
            dist_bp   = as.numeric(d[hits])
          )
        }
      }
    }
  }

  warn_tbl <- if (length(warn_list)) data.table::rbindlist(warn_list) else data.table::data.table()
  attr(x, "cluster_warnings") <- warn_tbl
  if (verbose && nrow(warn_tbl)) {
    message(sprintf(
      "[cluster warning] %d pair(s) closer than %s bp:",
      nrow(warn_tbl), format(warn_if_clusters_closer_than, big.mark = ",")
    ))
    print(warn_tbl)
  }
  x
}


#' Apply a clustering function per chromosome
#'
#' @description Applies a clustering function to each chromosome subset and
#'   optionally re-scopes cluster IDs globally.
#'
#' @param dt A \code{data.table}.
#' @param chrom_col Name of chromosome column.
#' @param f A function that takes a per-chromosome \code{data.table} and
#'   returns one with a \code{cluster_id} column.
#' @param id_scope Either \code{"per_chrom"} (IDs restart at 1 per chrom) or
#'   \code{"global"} (IDs are unique across all chroms).
#'
#' @return A \code{data.table} with results from all chromosomes combined.
#'
#' @keywords internal
#' @noRd
.apply_by_chrom <- function(dt, chrom_col, f, id_scope = c("per_chrom","global")) {
  id_scope <- match.arg(id_scope)
  x <- data.table::copy(dt)
  chroms <- unique(x[[chrom_col]])
  out_list <- vector("list", length(chroms))
  names(out_list) <- as.character(chroms)
  next_id <- 1L
  for (i in seq_along(chroms)) {
    cc <- chroms[i]
    sub <- x[x[[chrom_col]] == cc]
    res <- f(sub)
    if (id_scope == "global" && any(res$cluster_id > 0L, na.rm = TRUE)) {
      max_here <- max(res$cluster_id, na.rm = TRUE)
      res[, cluster_id := data.table::fifelse(cluster_id > 0L,
                                              cluster_id + (next_id - 1L), 0L)]
      next_id <- next_id + max_here
    }
    out_list[[i]] <- res
  }
  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
}


#' Add chromosome-scoped cluster ID labels
#'
#' @description Creates a human-readable \code{chrom_cluster_id} column of the
#'   form \code{"Chr1_12345678"} (chromosome + cluster centre position).
#'
#' @param dt A \code{data.table} with \code{cluster_id} and \code{chrom_col}.
#' @param chrom_col Name of chromosome column. Default \code{"chrom"}.
#' @param start Name of start column. Default \code{"start"}.
#' @param end Name of end column. Default \code{"end"}.
#' @param id_col Name of cluster ID column. Default \code{"cluster_id"}.
#' @param label_col Name of output label column. Default \code{"chrom_cluster_id"}.
#' @param zero_label How to label noise (\code{cluster_id == 0}):
#'   \code{"chrom_0"} or \code{"NA"}.
#' @param pos_stat Summary statistic for cluster centre: \code{"median"} or
#'   \code{"mean"}.
#' @param pos_round Logical. Round position to nearest integer. Default
#'   \code{TRUE}.
#'
#' @return The input \code{data.table} with \code{label_col} added.
#'
#' @keywords internal
#' @noRd
.add_chrom_cluster_id <- function(dt,
                                  chrom_col  = "chrom",
                                  start      = "start",
                                  end        = "end",
                                  id_col     = "cluster_id",
                                  label_col  = "chrom_cluster_id",
                                  zero_label = c("chrom_0","NA"),
                                  pos_stat   = c("median","mean"),
                                  pos_round  = TRUE) {
  stopifnot(data.table::is.data.table(dt))
  zero_label <- match.arg(zero_label)
  pos_stat   <- match.arg(pos_stat)

  x <- data.table::copy(dt)

  centers <- x[get(id_col) > 0L,
               {
                 s <- suppressWarnings(as.numeric(.SD[[1]]))
                 e <- suppressWarnings(as.numeric(.SD[[2]]))
                 mid <- (s + e)/2
                 cval <- if (pos_stat == "median") stats::median(mid, na.rm = TRUE) else mean(mid, na.rm = TRUE)
                 if (pos_round && is.finite(cval)) cval <- as.integer(round(cval))
                 .(center = cval)
               },
               by = c(chrom_col, id_col),
               .SDcols = c(start, end)
  ]

  x[, (label_col) := NA_character_]
  if (nrow(centers)) {
    x[centers, on = c(chrom_col, id_col),
      (label_col) := paste0(get(chrom_col), "_", i.center)]
  }

  if (zero_label == "chrom_0") {
    x[get(id_col) == 0L, (label_col) := paste0(get(chrom_col), "_0")]
  } else {
    x[get(id_col) == 0L, (label_col) := NA_character_]
  }
  x
}


# ---- Public clustering functions --------------------------------------------

#' Cluster reads using a sweep algorithm
#'
#' @title Sweep-based read clustering
#' @description Clusters genomic intervals per chromosome by a greedy
#'   left-to-right sweep: two consecutive intervals are merged when the gap
#'   between them is at most \code{delta} bp. After sweeping, any cluster
#'   with an internal gap larger than \code{max_gap} is split.
#'
#' @param dt A \code{data.table} with at least \code{chrom_col}, \code{start},
#'   and \code{end} columns.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param start Name of the start coordinate column. Default \code{"start"}.
#' @param end Name of the end coordinate column. Default \code{"end"}.
#' @param delta Integer. Maximum gap (bp) for merging consecutive intervals
#'   during the sweep. Default \code{0L} (overlapping only).
#' @param max_gap Numeric. After sweeping, split any cluster whose consecutive
#'   sorted reads are separated by more than this many bp. Default \code{Inf}
#'   (no splitting).
#' @param max_span Deprecated alias for \code{max_gap}.
#' @param warn_if_clusters_closer_than Numeric. Warn when two cluster centres
#'   are closer than this many bp. Default \code{NA} (no warnings).
#' @param id_scope Character. Either \code{"global"} (cluster IDs unique across
#'   all chromosomes) or \code{"per_chrom"} (IDs restart per chromosome).
#' @param verbose Logical. Print split/warning messages. Default \code{TRUE}.
#'
#' @return The input \code{data.table} with three new columns: \code{cluster_id}
#'   (integer), \code{cluster_n} (number of reads in the cluster), and
#'   \code{chrom_cluster_id} (character label, e.g. \code{"Chr1_12345678"}).
#'
#' @examples
#' \dontrun{
#' bed <- read_bed12(system.file("extdata", "alignVA2.bed12", package = "kolibri"))
#' split_reads <- bed[nameN > 1 & score > 30]
#' clusters <- cluster_reads_sweep(split_reads, delta = 1000, max_gap = 2000)
#' }
#'
#' @family clustering
#' @importFrom data.table copy is.data.table setkeyv rbindlist fifelse
#' @importFrom stats median
#' @export
cluster_reads_sweep <- function(dt, chrom_col = "chrom",
                                start = "start", end = "end",
                                delta = 0L,
                                max_gap = Inf,
                                max_span = NULL,
                                warn_if_clusters_closer_than = NA_real_,
                                id_scope = c("global"),
                                verbose = TRUE) {
  if (!is.null(max_span)) {
    max_gap <- max_span
    if (verbose) message("Using max_span as max_gap (consecutive-gap rule).")
  }
  id_scope <- match.arg(id_scope)
  stopifnot(data.table::is.data.table(dt))

  per_chrom <- function(sub) {
    y <- data.table::copy(sub)[, .idx := .I]

    s <- suppressWarnings(as.numeric(y[[start]]))
    e <- suppressWarnings(as.numeric(y[[end]]))
    o <- order(s, e, na.last = TRUE)
    s <- s[o]; e <- e[o]
    n <- length(s)

    cl <- integer(n); cur <- 0L
    cur_end <- -Inf
    for (i in seq_len(n)) {
      si <- s[i]; ei <- e[i]
      if (is.na(si) || is.na(ei)) { cl[i] <- 0L; next }
      if (si <= cur_end + delta) {
        cur_end <- max(cur_end, ei, na.rm = TRUE)
      } else {
        cur <- cur + 1L
        cur_end <- ei
      }
      cl[i] <- cur
    }

    y[o, cluster_id := cl]
    data.table::setorder(y, .idx); y[, .idx := NULL]

    y <- .split_by_max_gap(y, "cluster_id", start, end, max_gap)
    .finalize_clusters(y, start, end, chrom_col,
                       warn_if_clusters_closer_than, verbose)
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  .add_chrom_cluster_id(res,
                        chrom_col = chrom_col,
                        start = start, end = end,
                        id_col = "cluster_id",
                        label_col = "chrom_cluster_id",
                        zero_label = "chrom_0",
                        pos_stat = "median",
                        pos_round = TRUE)
}


#' Cluster reads using DBSCAN
#'
#' @title DBSCAN-based read clustering
#' @description Clusters genomic intervals per chromosome using the DBSCAN
#'   density-based algorithm. Three distance metrics are supported: gap-based
#'   (distance = max(0, gap between intervals)), Euclidean (on midpoint and
#'   length), and L-infinity.
#'
#' @param dt A \code{data.table} with at least \code{chrom_col}, \code{start},
#'   and \code{end} columns.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param start Name of the start coordinate column. Default \code{"start"}.
#' @param end Name of the end coordinate column. Default \code{"end"}.
#' @param distance Distance metric: \code{"gap"}, \code{"euclidean"}, or
#'   \code{"linf"}.
#' @param eps DBSCAN epsilon radius. If \code{NULL}, defaults to the maximum
#'   interval length on the chromosome.
#' @param minPts DBSCAN minimum points. Default \code{3L}.
#' @param max_gap Post-clustering gap split threshold. Default \code{Inf}.
#' @param max_span Deprecated alias for \code{max_gap}.
#' @param warn_if_clusters_closer_than Numeric. Warn when cluster centres are
#'   closer than this many bp. Default \code{NA}.
#' @param id_scope \code{"per_chrom"} or \code{"global"}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return The input \code{data.table} with \code{cluster_id}, \code{cluster_n},
#'   and \code{chrom_cluster_id} columns added.
#'
#' @examples
#' \dontrun{
#' clusters <- cluster_reads_dbscan(split_reads, eps = 5000, minPts = 3)
#' }
#'
#' @family clustering
#' @importFrom data.table copy is.data.table rbindlist fifelse
#' @importFrom dbscan dbscan
#' @importFrom stats median
#' @export
cluster_reads_dbscan <- function(dt, chrom_col = "chrom",
                                 start = "start", end = "end",
                                 distance = c("gap","euclidean","linf"),
                                 eps = NULL, minPts = 3L,
                                 max_gap = Inf, max_span = NULL,
                                 warn_if_clusters_closer_than = NA_real_,
                                 id_scope = c("per_chrom","global"),
                                 verbose = TRUE) {
  if (!is.null(max_span)) {
    max_gap <- max_span
    if (verbose) message("Using max_span as max_gap (consecutive-gap rule).")
  }
  distance <- match.arg(distance); id_scope <- match.arg(id_scope)
  stopifnot(data.table::is.data.table(dt))

  per_chrom <- function(sub) {
    s <- as.numeric(sub[[start]])
    e <- as.numeric(sub[[end]]); n <- length(s)
    local_eps <- if (is.null(eps)) max(e - s, na.rm = TRUE) else eps

    if (distance == "gap") {
      D <- matrix(0, n, n)
      for (i in seq_len(n)) {
        gi <- pmax(0, pmax(s[i], s) - pmin(e[i], e))
        D[i, ] <- gi; D[, i] <- gi
      }
      fit <- dbscan::dbscan(stats::as.dist(D), eps = local_eps,
                            minPts = minPts, borderPoints = TRUE)
    } else if (distance == "euclidean") {
      X <- cbind((s + e) / 2, e - s)
      fit <- dbscan::dbscan(X, eps = local_eps, minPts = minPts)
    } else {
      D <- matrix(0, n, n)
      for (i in seq_len(n)) {
        di <- pmax(abs(s[i] - s), abs(e[i] - e))
        D[i, ] <- di; D[, i] <- di
      }
      fit <- dbscan::dbscan(stats::as.dist(D), eps = local_eps,
                            minPts = minPts, borderPoints = TRUE)
    }

    out <- data.table::copy(sub)
    out[, cluster_id := as.integer(fit$cluster)]
    out <- .split_by_max_gap(out, "cluster_id", start, end, max_gap)
    .finalize_clusters(out, start, end, chrom_col,
                       warn_if_clusters_closer_than, verbose)
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  .add_chrom_cluster_id(res,
                        chrom_col = chrom_col,
                        start = start, end = end,
                        id_col = "cluster_id",
                        label_col = "chrom_cluster_id",
                        zero_label = "chrom_0",
                        pos_stat = "median",
                        pos_round = TRUE)
}


#' Cluster reads using kernel density estimation
#'
#' @title KDE-based read clustering
#' @description Clusters genomic intervals per chromosome by fitting a weighted
#'   kernel density estimate to the coverage landscape and assigning each
#'   interval to the nearest density mode. Intervals further than
#'   \code{max_gap} from any mode are labelled as unassigned.
#'
#' @param dt A \code{data.table} with at least \code{chrom_col}, \code{start},
#'   and \code{end} columns.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param start Name of the start coordinate column. Default \code{"start"}.
#' @param end Name of the end coordinate column. Default \code{"end"}.
#' @param bw Numeric bandwidth. If \code{NULL}, computed automatically via
#'   \code{bw_method}.
#' @param bw_method Bandwidth selection method: \code{"SJ"} (default),
#'   \code{"nrd0"}, \code{"ucv"}, or \code{"bcv"}.
#' @param bw_adjust Multiplicative adjustment to the computed bandwidth.
#'   Default \code{1}.
#' @param kernel Kernel name passed to \code{\link[stats]{density}}. Default
#'   \code{"gaussian"}.
#' @param max_gap Post-clustering gap split threshold. Default \code{Inf}.
#' @param max_span Deprecated alias for \code{max_gap}.
#' @param warn_if_clusters_closer_than Numeric warning threshold. Default
#'   \code{NA}.
#' @param id_scope \code{"per_chrom"} or \code{"global"}.
#' @param verbose Logical. Default \code{TRUE}.
#' @param unassigned_label Integer label for unassigned reads. Default
#'   \code{0L}.
#'
#' @return The input \code{data.table} with \code{cluster_id}, \code{cluster_n},
#'   and \code{chrom_cluster_id} columns added.
#'
#' @examples
#' \dontrun{
#' clusters <- cluster_reads_kde(split_reads, bw_method = "SJ", bw_adjust = 1.5)
#' }
#'
#' @family clustering
#' @importFrom data.table copy is.data.table rbindlist fifelse
#' @importFrom stats density bw.SJ bw.nrd0 bw.ucv bw.bcv median
#' @export
cluster_reads_kde <- function(dt, chrom_col = "chrom",
                              start = "start", end = "end",
                              bw = NULL, bw_method = c("SJ","nrd0","ucv","bcv"),
                              bw_adjust = 1, kernel = "gaussian",
                              max_gap = Inf, max_span = NULL,
                              warn_if_clusters_closer_than = NA_real_,
                              id_scope = c("per_chrom","global"),
                              verbose = TRUE,
                              unassigned_label = 0L) {
  if (!is.null(max_span)) {
    max_gap <- max_span
    if (verbose) message("Using max_span as max_gap (consecutive-gap rule).")
  }
  bw_method <- match.arg(bw_method); id_scope <- match.arg(id_scope)
  stopifnot(data.table::is.data.table(dt))

  per_chrom <- function(sub) {
    s <- as.numeric(sub[[start]]); e <- as.numeric(sub[[end]])
    if (!length(s)) {
      z <- data.table::copy(sub)
      z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]
      return(z)
    }

    ev_pos <- c(s, e); ev_val <- c(rep(1L, length(s)), rep(-1L, length(e)))
    o <- order(ev_pos, -ev_val); ev_pos <- ev_pos[o]; ev_val <- ev_val[o]
    pos_u <- unique(ev_pos)
    if (length(pos_u) < 2L) {
      z <- data.table::copy(sub)
      z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]
      return(z)
    }
    cov <- integer(length(pos_u) - 1L); cur <- 0L; k <- 1L
    for (i in seq_along(ev_val)) {
      cur <- cur + ev_val[i]
      if (i < length(ev_val) && ev_pos[i] < ev_pos[i+1]) { cov[k] <- cur; k <- k + 1L }
    }
    seg_start <- pos_u[-length(pos_u)]
    seg_end   <- pos_u[-1L]
    seg_mid   <- (seg_start + seg_end) / 2
    seg_len   <- pmax(0, seg_end - seg_start)
    w <- cov * seg_len
    keep <- w > 0

    if (!any(keep)) {
      z <- data.table::copy(sub)
      z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]
      return(z)
    }

    x <- seg_mid[keep]; w <- w[keep] / sum(w[keep])
    local_bw <- if (is.null(bw)) switch(bw_method,
                                        SJ   = stats::bw.SJ(x, method="ste"),
                                        nrd0 = stats::bw.nrd0(x),
                                        ucv  = stats::bw.ucv(x),
                                        bcv  = stats::bw.bcv(x)) else bw

    kd <- stats::density(x, weights = w, bw = local_bw * bw_adjust,
                         kernel = kernel)
    midx <- which(diff(sign(diff(kd$y))) == -2) + 1
    modes <- kd$x[midx]

    lab <- integer(length(s))
    for (i in seq_along(s)) {
      if (!length(modes)) { lab[i] <- unassigned_label; next }
      d <- ifelse(modes >= s[i] & modes <= e[i], 0,
                  pmin(abs(modes - s[i]), abs(modes - e[i])))
      j <- which.min(d)
      lab[i] <- if (length(j)==0 || d[j] > max_gap) unassigned_label else j
    }
    out <- data.table::copy(sub)
    out[, cluster_id := as.integer(lab)]

    out <- .split_by_max_gap(out, "cluster_id", start, end, max_gap)
    out <- .finalize_clusters(out, start, end, chrom_col,
                              warn_if_clusters_closer_than, verbose)
    attr(out, "kde_modes") <- modes
    attr(out, "kde_bw")    <- kd$bw
    out
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  .add_chrom_cluster_id(res,
                        chrom_col = chrom_col,
                        start = start, end = end,
                        id_col = "cluster_id",
                        label_col = "chrom_cluster_id",
                        zero_label = "chrom_0",
                        pos_stat = "median",
                        pos_round = TRUE)
}
