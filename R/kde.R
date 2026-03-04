#' Boundary-reflected kernel density estimate
#'
#' @title Reflection-corrected KDE
#' @description Computes a kernel density estimate with boundary reflection to
#'   reduce edge bias. Points near the left and right boundaries are mirrored
#'   before estimation, and the density is evaluated only over the original
#'   data range.
#'
#' @param x Numeric vector of observations.
#' @param bw Numeric bandwidth. If \code{NULL}, computed from \code{bw_method}.
#' @param bw_method Bandwidth selection method: \code{"SJ"} (default),
#'   \code{"nrd0"}, \code{"ucv"}, or \code{"bcv"}.
#' @param kernel Kernel name. Default \code{"gaussian"}.
#' @param grid_n Number of evaluation points. Default \code{2048}.
#' @param reflect_multiplier Points within this many bandwidths of a boundary
#'   are reflected. Default \code{3}.
#'
#' @return A named list with elements \code{bw} (bandwidth used), \code{density}
#'   (the full \code{\link[stats]{density}} object), \code{modes} (x-positions
#'   of local maxima), and \code{mode_idx} (indices into the density grid).
#'
#' @examples
#' \dontrun{
#' res <- kde_reflect(c(1000, 2000, 2100, 3000, 5000))
#' plot(res$density)
#' }
#'
#' @family KDE
#' @importFrom stats density bw.SJ bw.nrd0 bw.ucv bw.bcv
#' @export
kde_reflect <- function(x,
                        bw = NULL,
                        bw_method = c("SJ","nrd0","ucv","bcv"),
                        kernel = "gaussian",
                        grid_n = 2048,
                        reflect_multiplier = 3) {
  bw_method <- match.arg(bw_method)

  if (is.null(bw)) {
    bw <- switch(
      bw_method,
      SJ   = stats::bw.SJ(x, method = "ste"),
      nrd0 = stats::bw.nrd0(x),
      ucv  = stats::bw.ucv(x),
      bcv  = stats::bw.bcv(x)
    )
  }

  L <- min(x); R <- max(x)

  left_keep  <- x < (L + reflect_multiplier * bw)
  right_keep <- x > (R - reflect_multiplier * bw)

  x_left_ref  <- 2 * L - x[left_keep]
  x_right_ref <- 2 * R - x[right_keep]

  x_ref <- c(x, x_left_ref, x_right_ref)

  kd <- stats::density(
    x_ref,
    bw     = bw,
    kernel = kernel,
    from   = L,
    to     = R,
    n      = grid_n
  )

  y  <- kd$y
  xx <- kd$x
  interior_idx <- which(diff(sign(diff(y))) == -2) + 1
  mode_idx <- interior_idx
  if (y[1] > y[2]) mode_idx <- c(1L, mode_idx)
  n <- length(y)
  if (y[n] > y[n-1]) mode_idx <- c(mode_idx, n)

  list(
    bw       = kd$bw,
    density  = kd,
    modes    = xx[mode_idx],
    mode_idx = mode_idx
  )
}


#' Estimate 1-D density modes
#'
#' @title Estimate density modes
#' @description Fits a kernel density estimate and returns the x-positions of
#'   all interior local maxima (modes).
#'
#' @param x Numeric vector of observations.
#' @param bw Numeric bandwidth. If \code{NULL}, computed from \code{bw_method}.
#' @param bw_method Bandwidth selection method. Default \code{"SJ"}.
#' @param kernel Kernel name. Default \code{"gaussian"}.
#' @param grid_n Number of evaluation grid points. Default \code{2048}.
#'
#' @return A named list with elements \code{bw}, \code{modes}, and
#'   \code{density}.
#'
#' @examples
#' \dontrun{
#' res <- estimate_modes_1d(c(1000, 2000, 2100, 5000))
#' res$modes
#' }
#'
#' @family KDE
#' @importFrom stats density bw.SJ bw.nrd0 bw.ucv bw.bcv
#' @export
estimate_modes_1d <- function(x,
                              bw = NULL,
                              bw_method = c("SJ","nrd0","ucv","bcv"),
                              kernel = "gaussian",
                              grid_n = 2048) {
  bw_method <- match.arg(bw_method)
  if (is.null(bw)) {
    bw <- switch(bw_method,
                 SJ   = stats::bw.SJ(x, method = "ste"),
                 nrd0 = stats::bw.nrd0(x),
                 ucv  = stats::bw.ucv(x),
                 bcv  = stats::bw.bcv(x))
  }
  kd <- stats::density(x, bw = bw, kernel = kernel, n = grid_n)
  y <- kd$y; idx <- which(diff(sign(diff(y))) == -2) + 1
  list(bw = kd$bw, modes = kd$x[idx], density = kd)
}


#' Plot read breakpoint KDE for a network node
#'
#' @title KDE plot of read breakpoints for a node
#' @description For a chosen node, pools all read start and end positions and
#'   overlays a reflection-corrected KDE on top of a histogram to visualise
#'   breakpoint density.
#'
#' @param net A network list from \code{\link{build_cluster_network}}.
#' @param node Character. The \code{chrom_cluster_id} of the node.
#' @param cluster_col Name of the cluster column in \code{net$reads}. Default
#'   \code{"chrom_cluster_id"}.
#' @param bw_method Bandwidth selection method. Default \code{"SJ"}.
#' @param bins Integer. Number of histogram bins. Default \code{50L}.
#' @param kernel Kernel name. Default \code{"gaussian"}.
#' @param lw KDE line width. Default \code{1}.
#'
#' @return A named list with elements \code{plot} (a \code{ggplot2} object),
#'   \code{breaks} (numeric vector of start/end positions), \code{dens_df}
#'   (density data.table), \code{bw_kde}, and \code{bw_hist}.
#'
#' @examples
#' \dontrun{
#' res <- plot_node_breaks_kde(net, node = "Chr1_12953910")
#' res$plot
#' }
#'
#' @family KDE
#' @importFrom data.table as.data.table data.table
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line labs theme_bw
#' @export
plot_node_breaks_kde <- function(net,
                                 node,
                                 cluster_col = "chrom_cluster_id",
                                 bw_method   = "SJ",
                                 bins        = 50L,
                                 kernel      = "gaussian",
                                 lw          = 1) {
  dt_reads <- data.table::as.data.table(net$reads)
  if (!all(c(cluster_col, "start", "end") %in% names(dt_reads)))
    stop("net$reads must have columns: ", cluster_col, ", start, end")

  node_rows <- dt_reads[get(cluster_col) == node]
  if (!nrow(node_rows))
    stop(sprintf("No reads found for node '%s'", node))

  start_ends <- c(node_rows$start, node_rows$end)
  res_kde    <- kde_reflect(start_ends, bw_method = bw_method, kernel = kernel)

  dt      <- data.table::data.table(x = start_ends)
  bw_hist <- diff(range(dt$x)) / bins

  dens_df <- data.table::data.table(
    x = res_kde$density$x,
    y = res_kde$density$y
  )
  dens_df$y_scaled <- dens_df$y * nrow(dt) * bw_hist

  x_peak <- dens_df[which.max(y_scaled), x]
  dens_df[, rel_x := x - x_peak]

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(binwidth = bw_hist,
                            fill = "grey90", color = "grey80") +
    ggplot2::geom_line(data = dens_df,
                       ggplot2::aes(x = x, y = y_scaled),
                       color = "steelblue", linewidth = lw) +
    ggplot2::labs(
      title = paste0("Kernel density for ", node),
      x     = "Read breaks",
      y     = "Count"
    ) +
    ggplot2::theme_bw(base_size = 6)

  list(
    plot    = p,
    breaks  = start_ends,
    dens_df = dens_df,
    bw_kde  = res_kde$density$bw,
    bw_hist = bw_hist
  )
}
