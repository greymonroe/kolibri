#' Plot cluster network with ggraph
#'
#' @title Plot shared-read cluster network
#' @description Renders a force-directed (or other layout) graph of the cluster
#'   network. Nodes are scaled by cluster size, coloured by chromosome, and
#'   optionally labelled. Low-count nodes can be filtered while keeping their
#'   direct neighbours.
#'
#' @param net A network list from \code{\link{build_cluster_network}}.
#' @param layout ggraph layout algorithm. Default \code{"fr"}
#'   (Fruchterman-Reingold). Other options: \code{"kk"}, \code{"drl"}, etc.
#' @param cap Maximum node size (cluster_n capped at this value). Default
#'   \code{30}.
#' @param alpha_var Column in \code{net$nodes} to map to node transparency.
#'   Default \code{"node_avg_score"}. Set to \code{NULL} for uniform alpha.
#' @param alpha_range Numeric vector of length 2 giving the alpha range.
#'   Default \code{c(0.3, 1.0)}.
#' @param label Logical. Whether to label nodes. Default \code{FALSE}.
#' @param label_col Column to use for node labels. Default \code{"cluster"}.
#' @param label_size Text size for labels. Default \code{2}.
#' @param label_repel Logical. Use repelling labels. Default \code{TRUE}.
#' @param seed Random seed for layout reproducibility. Default \code{1}.
#' @param min_cluster_n Minimum \code{cluster_n} (or \code{edge_count}) for a
#'   node to be retained. Default \code{1L} (keep all).
#' @param keep_neighbors_of_kept Logical. Retain direct neighbours of kept
#'   nodes even if they fall below \code{min_cluster_n}. Default \code{TRUE}.
#'
#' @return A \code{ggraph}/\code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' net <- build_cluster_network(clusters)
#' plot_cluster_network_ggraph(net, label = TRUE, min_cluster_n = 3)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes scale_size_binned scale_alpha theme_void
#'   labs annotate guide_bins
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @export
plot_cluster_network_ggraph <- function(net,
                                        layout = "fr",
                                        cap = 30,
                                        alpha_var = "node_avg_score",
                                        alpha_range = c(0.3, 1.0),
                                        label = FALSE,
                                        label_col = "cluster",
                                        label_size = 2,
                                        label_repel = TRUE,
                                        seed = 1,
                                        min_cluster_n = 1L,
                                        keep_neighbors_of_kept = TRUE) {
  stopifnot(is.list(net), all(c("nodes","edges") %in% names(net)))
  nodes <- data.table::as.data.table(net$nodes)
  edges <- data.table::as.data.table(net$edges)

  if (!nrow(nodes) || !nrow(edges)) {
    message("No nodes/edges to plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = "No edges to plot"))
  }

  if (!is.null(min_cluster_n) && min_cluster_n > 1L) {
    kept <- nodes[edge_count > min_cluster_n | cluster_n >= min_cluster_n, cluster]
    if (isTRUE(keep_neighbors_of_kept) && length(kept)) {
      nbr1 <- edges[cluster1 %in% kept, cluster2]
      nbr2 <- edges[cluster2 %in% kept, cluster1]
      kept <- unique(c(kept, nbr1, nbr2))
    }
    nodes <- nodes[cluster %in% kept]
    if (!nrow(nodes)) {
      message("No nodes remain after threshold+neighbour filter.")
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::annotate("text", x = 0, y = 0,
                                 label = "No nodes after filter"))
    }
    edges <- edges[cluster1 %in% nodes$cluster & cluster2 %in% nodes$cluster]
    if (!nrow(edges)) message("No edges remain after filter.")
  }

  if (!is.null(alpha_var) && alpha_var %in% names(nodes)) {
    rng <- range(nodes[[alpha_var]], na.rm = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
      nodes[, node_alpha := (get(alpha_var) - rng[1]) / (rng[2] - rng[1])]
    } else nodes[, node_alpha := 1]
  } else nodes[, node_alpha := 1]

  g <- igraph::graph_from_data_frame(
    d = if (nrow(edges)) edges[, .(from = cluster1, to = cluster2)] else NULL,
    directed = FALSE,
    vertices = nodes[, .(name = if ("name" %in% names(nodes)) name else cluster,
                         cluster    = cluster,
                         chrom      = chrom,
                         cluster_n  = cluster_n,
                         node_alpha = node_alpha)]
  )

  set.seed(seed)
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(color = "black", width = 0.25, alpha = 0.6,
                           show.legend = FALSE) +
    ggraph::geom_node_point(
      ggplot2::aes(size = pmin(cluster_n, cap), color = chrom,
                   alpha = node_alpha),
      stroke = 0.2
    ) +
    ggplot2::scale_size_binned(
      range       = c(0.8, 5),
      breaks      = c(1, 3, 5, 10, 20, cap),
      labels      = c("1","3","5","10","20", paste0(cap, "+")),
      limits      = c(1, cap),
      nice.breaks = FALSE,
      guide       = ggplot2::guide_bins(title = "Cluster size")
    ) +
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::theme_void(base_size = 6) +
    ggplot2::labs(color = "Chromosome")

  if (isTRUE(label)) {
    lab <- if (label_col %in% names(nodes)) label_col else "cluster"
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = .data[[lab]]),
      repel = isTRUE(label_repel),
      size  = label_size
    )
  }
  p
}


#' Plot reads for a node and all its network partners
#'
#' @title Plot reads by cluster for a node
#' @description For a chosen primary node, plots all reads belonging to it
#'   across the primary cluster and any partner clusters it is connected to in
#'   the network. Reads are arranged with one consistent y-position per read
#'   name. Facets represent individual clusters.
#'
#' @param net A network list from \code{\link{build_cluster_network}}.
#' @param node Character. The \code{chrom_cluster_id} of the primary node.
#' @param cluster_col Name of the cluster column in \code{net$reads}. Default
#'   \code{"chrom_cluster_id"}.
#' @param read_id_col Name of the read name column. Default \code{"name"}.
#' @param start_col Name of the start coordinate column. Default \code{"start"}.
#' @param end_col Name of the end coordinate column. Default \code{"end"}.
#' @param alpha_col Column to map to transparency (e.g. \code{"score"}). Default
#'   \code{"score"}.
#' @param alpha_range Numeric vector of length 2. Default \code{c(0.3, 1.0)}.
#' @param primary_color Colour for segments in the primary node. Default
#'   \code{"red2"}.
#' @param partner_color Colour for segments in partner nodes. Default
#'   \code{"grey60"}.
#' @param base_size Base font size. Default \code{6}.
#' @param ensure_strip_space Logical. Apply theme tweaks for readable strip
#'   labels. Default \code{TRUE}.
#' @param min_facet_width_mm Minimum facet width in mm (used when ggh4x is
#'   available). Default \code{18}.
#' @param strip_wrap_width Character wrap width for strip labels. Default
#'   \code{20}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' plot_node_reads_by_cluster(net, node = "Chr1_12953910")
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table setkey setorder uniqueN fifelse
#' @importFrom ggplot2 ggplot aes geom_segment facet_grid scale_color_manual
#'   scale_alpha scale_x_continuous guides theme_bw theme element_blank
#'   element_text margin coord_cartesian label_wrap_gen
#' @importFrom grid unit
#' @export
plot_node_reads_by_cluster <- function(net, node,
                                       cluster_col   = "chrom_cluster_id",
                                       read_id_col   = "name",
                                       start_col     = "start",
                                       end_col       = "end",
                                       alpha_col     = "score",
                                       alpha_range   = c(0.3, 1.0),
                                       primary_color = "red2",
                                       partner_color = "grey60",
                                       base_size     = 6,
                                       ensure_strip_space = TRUE,
                                       min_facet_width_mm = 18,
                                       strip_wrap_width   = 20) {
  stopifnot(is.list(net), all(c("edges","reads") %in% names(net)))
  edges <- data.table::as.data.table(net$edges)
  reads <- data.table::as.data.table(net$reads)
  nodes <- if ("nodes" %in% names(net)) data.table::as.data.table(net$nodes) else NULL

  partners <- unique(c(
    edges[cluster1 == node, cluster2],
    edges[cluster2 == node, cluster1]
  ))
  facets <- unique(c(node, partners))
  if (!length(facets)) stop(sprintf("Node '%s' has no partners in net$edges.", node))

  if (!all(c(cluster_col, read_id_col, start_col, end_col) %in% names(reads))) {
    stop("Expected columns missing in net$reads.")
  }
  primary_reads <- unique(reads[get(cluster_col) == node, get(read_id_col)])
  if (!length(primary_reads))
    stop(sprintf("No reads found in primary node '%s'.", node))

  sub <- reads[
    get(read_id_col) %chin% primary_reads & get(cluster_col) %chin% facets
  ][, `:=`(
    cl  = as.character(get(cluster_col)),
    rid = as.character(get(read_id_col))
  )]

  sub[, cluster_type := data.table::fifelse(cl == node, "primary", "partner")]

  partner_levels <- sort(setdiff(facets, node))
  facet_levels   <- c(node, partner_levels)
  sub[, (cluster_col) := factor(cl, levels = facet_levels)]

  if (!is.null(nodes) && all(c("cluster","cluster_n") %in% names(nodes))) {
    cl_size <- nodes[cluster %chin% facet_levels, .(cluster, cluster_n)]
  } else {
    cl_size <- sub[, .(cluster_n = data.table::uniqueN(rid)), by = .(cluster = cl)]
  }
  data.table::setkey(cl_size, cluster)

  per_read_cl <- sub[, .(
    read_min_start = suppressWarnings(min(get(start_col), na.rm = TRUE))
  ), by = .(cl, rid)]
  data.table::setkey(per_read_cl, rid, cl)

  read_clusters <- per_read_cl[, .(clusters = list(cl)), by = rid]
  data.table::setkey(read_clusters, rid)

  choose_best_secondary <- function(rid_val, current_cl) {
    cls <- read_clusters[.(rid_val), clusters][[1]]
    cls <- setdiff(cls, current_cl)
    cls <- setdiff(cls, node)
    if (length(cls) == 0L) return(list(NA_character_, 0L, Inf))
    cand <- data.table::data.table(cluster = cls)
    cand <- cl_size[cand, on = .(cluster)]
    cand[is.na(cluster_n), cluster_n := 0L]
    st <- per_read_cl[.(rid_val, cand$cluster), read_min_start]
    if (length(st) != nrow(cand)) {
      st <- per_read_cl[.(rid_val, cand$cluster), on = .(rid, cl),
                        x.read_min_start]
    }
    st[is.na(st)] <- Inf
    cand[, read_min_start := st]
    data.table::setorder(cand, -cluster_n, read_min_start, cluster)
    list(cand$cluster[1], as.integer(cand$cluster_n[1]),
         cand$read_min_start[1])
  }

  primary_start <- per_read_cl[.(rid = primary_reads, cl = node),
                               on = .(rid, cl)]
  primary_start <- primary_start[, .(rid, primary_start = read_min_start)]
  data.table::setkey(primary_start, rid)

  order_reads <- data.table::data.table(rid = primary_reads)
  tmp <- order_reads[, c("sec_cl","sec_n","sec_start") :=
                       choose_best_secondary(rid, node), by = rid]
  order_reads <- tmp
  order_reads <- primary_start[order_reads, on = .(rid)]
  order_reads[is.na(primary_start), primary_start := Inf]

  data.table::setorder(order_reads, -sec_n, sec_start, primary_start, rid)
  order_reads[, y := (.N + 1L) - seq_len(.N)]

  data.table::setkey(order_reads, rid)
  data.table::setkey(sub, rid)
  sub <- order_reads[sub]

  if (!alpha_col %in% names(sub))
    stop(sprintf("alpha_col '%s' not found in net$reads.", alpha_col))
  av  <- as.numeric(sub[[alpha_col]])
  rng <- range(av, na.rm = TRUE)
  a0  <- if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
    (av - rng[1]) / (rng[2] - rng[1])
  } else rep(1, length(av))
  a0[!is.finite(a0)] <- 0.5
  sub[, alpha_val := alpha_range[1] + a0 * (alpha_range[2] - alpha_range[1])]

  p <- ggplot2::ggplot(sub) +
    ggplot2::geom_segment(
      ggplot2::aes(
        y = y, yend = y,
        x = .data[[start_col]], xend = .data[[end_col]],
        color = cluster_type, alpha = alpha_val
      ),
      linewidth = 0.25, lineend = "round"
    ) +
    ggplot2::facet_grid(
      ~ .data[[cluster_col]],
      scales   = "free_x",
      space    = "free_x",
      labeller = ggplot2::label_wrap_gen(width = strip_wrap_width)
    ) +
    ggplot2::scale_color_manual(
      values = c(primary = primary_color, partner = partner_color)) +
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::scale_x_continuous(
      position     = "top",
      breaks       = function(lims) {
        from <- floor(lims[1] / 1000) * 1000
        to   <- ceiling(lims[2] / 1000) * 1000
        seq(from, to, by = 1000L)
      },
      labels       = NULL,
      minor_breaks = NULL,
      expand       = c(0.02, 0.02)
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.placement     = "outside",
      strip.background    = ggplot2::element_blank(),
      strip.text          = ggplot2::element_text(
        angle = 90, vjust = 0.5,
        margin = ggplot2::margin(t = 4, b = 4)),
      strip.clip          = "off",
      panel.grid          = ggplot2::element_blank(),
      panel.spacing.x     = grid::unit(6, "pt"),
      axis.title          = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.title.x.bottom = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


#' Plot all reads overlapping a node's genomic window
#'
#' @title Plot reads around a node
#' @description Fetches all reads from the original BED file that overlap a
#'   window around the node's genomic span and plots them as horizontal
#'   segments. Reads that belong to the node are coloured differently from
#'   background reads.
#'
#' @param net A network list from \code{\link{build_cluster_network}}.
#' @param node Character. The \code{chrom_cluster_id} of the node.
#' @param bed A \code{data.table} with the full BED alignment (from
#'   \code{\link{read_bed12}}).
#' @param window_kb Numeric. Window around the node span in kilobases. Default
#'   \code{20}.
#' @param primary_color Colour for node reads. Default \code{"red2"}.
#' @param other_color Colour for background reads. Default \code{"grey60"}.
#' @param alpha_col Optional column name for transparency mapping. Default
#'   \code{NULL}.
#' @param alpha_range Numeric vector of length 2. Default \code{c(0.3, 1.0)}.
#' @param base_size Base font size. Default \code{6}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' plot_reads_around_node(net, node = "Chr1_12953910", bed = bed12)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table setkey setorder
#' @importFrom ggplot2 ggplot aes geom_segment scale_color_manual scale_alpha
#'   scale_x_continuous guides theme_bw theme element_blank
#' @importFrom grid unit
#' @export
plot_reads_around_node <- function(net, node, bed,
                                   window_kb   = 20,
                                   primary_color = "red2",
                                   other_color   = "grey60",
                                   alpha_col     = NULL,
                                   alpha_range   = c(0.3, 1.0),
                                   base_size     = 6) {
  stopifnot(is.list(net), "reads" %in% names(net))
  reads_dt <- data.table::as.data.table(net$reads)
  bed_dt   <- data.table::as.data.table(bed)

  need_reads_cols <- c("chrom_cluster_id","chrom","start","end","name")
  miss_reads <- setdiff(need_reads_cols, names(reads_dt))
  if (length(miss_reads))
    stop("net$reads missing columns: ", paste(miss_reads, collapse=", "))

  need_bed_cols <- c("chrom","start","end","name")
  miss_bed <- setdiff(need_bed_cols, names(bed_dt))
  if (length(miss_bed))
    stop("bed is missing columns: ", paste(miss_bed, collapse=", "))

  chrom_node <- NULL; node_start <- NULL; node_end <- NULL
  if ("nodes" %in% names(net) &&
      all(c("cluster","min_start","max_end","chrom") %in% names(net$nodes))) {
    node_row <- data.table::as.data.table(net$nodes)[cluster == node]
    if (nrow(node_row)) {
      chrom_node <- node_row$chrom[1]
      node_start <- as.numeric(node_row$min_start[1])
      node_end   <- as.numeric(node_row$max_end[1])
    }
  }
  if (is.null(chrom_node) || any(!is.finite(c(node_start, node_end)))) {
    node_reads <- reads_dt[chrom_cluster_id == node]
    if (!nrow(node_reads))
      stop(sprintf("No reads found for node '%s' in net$reads.", node))
    chrom_node <- as.character(node_reads$chrom[1])
    node_start <- min(as.numeric(node_reads$start), na.rm = TRUE)
    node_end   <- max(as.numeric(node_reads$end),   na.rm = TRUE)
  }

  pad       <- as.integer(round(window_kb * 1000))
  win_start <- max(0L, as.integer(floor(node_start - pad)))
  win_end   <- as.integer(ceiling(node_end + pad))

  window_reads <- bed_dt[chrom == chrom_node &
                           start <= win_end & end >= win_start]
  if (!nrow(window_reads)) stop("No reads found in the specified window.")

  window_reads[, cluster_type := ifelse(nameN > 1, "node", "other")]

  window_reads[, read_min_start :=
                 suppressWarnings(min(start, na.rm = TRUE)), by = name]
  ord <- unique(window_reads[, .(name, read_min_start)])
  data.table::setorder(ord, read_min_start, name)
  ord[, y := (.N + 1L) - seq_len(.N)]
  data.table::setkey(ord, name); data.table::setkey(window_reads, name)
  window_reads <- ord[window_reads]

  if (!is.null(alpha_col)) {
    if (!alpha_col %in% names(window_reads))
      stop(sprintf("alpha_col '%s' not found.", alpha_col))
    av  <- as.numeric(window_reads[[alpha_col]])
    r   <- range(av, na.rm = TRUE)
    a01 <- if (is.finite(r[1]) && is.finite(r[2]) && diff(r) > 0) {
      (av - r[1]) / (r[2] - r[1])
    } else rep(1, length(av))
    a01[!is.finite(a01)] <- 0.5
    window_reads[, alpha_val :=
                   alpha_range[1] + a01 * (alpha_range[2] - alpha_range[1])]
  } else {
    window_reads[, alpha_val := 1]
  }

  ggplot2::ggplot(window_reads) +
    ggplot2::geom_segment(
      ggplot2::aes(y = y, yend = y, x = start, xend = end,
                   color = cluster_type, alpha = alpha_val),
      linewidth = 0.5, lineend = "round") +
    ggplot2::scale_color_manual(
      values = c(node = primary_color, other = other_color)) +
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::scale_x_continuous(
      position     = "top",
      breaks       = function(lims) {
        from <- floor(lims[1] / 1000) * 1000
        to   <- ceiling(lims[2] / 1000) * 1000
        seq(from, to, by = 1000L)
      },
      labels       = NULL,
      minor_breaks = NULL,
      expand       = c(0, 0)
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid          = ggplot2::element_blank(),
      axis.title          = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.title.x.bottom = ggplot2::element_blank(),
      strip.placement     = "outside",
      strip.background    = ggplot2::element_blank(),
      axis.ticks.length   = grid::unit(-1.5, "mm")
    )
}
