#' Derive genome extents from a BED table
#'
#' @title Genome ranges from BED data
#' @description Computes per-chromosome min/max coordinates from any BED-like
#'   \code{data.table}, returning a simple genome extent table.
#'
#' @param bed A \code{data.frame} or \code{data.table} with at least
#'   \code{chrom_col}, \code{start_col}, and \code{end_col}.
#' @param chrom_col Name of the chromosome column. Default \code{"chrom"}.
#' @param start_col Name of the start column. Default \code{"start"}.
#' @param end_col Name of the end column. Default \code{"end"}.
#'
#' @return A \code{data.table} with columns \code{chrom}, \code{start}, and
#'   \code{end}, one row per chromosome, sorted by \code{chrom}.
#'
#' @examples
#' \dontrun{
#' genome_dt <- genome_from_bed(bed12)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table is.data.table
#' @export
genome_from_bed <- function(bed,
                            chrom_col = "chrom",
                            start_col = "start",
                            end_col   = "end") {
  dt <- if (data.table::is.data.table(bed)) bed else data.table::as.data.table(bed)

  req <- c(chrom_col, start_col, end_col)
  if (!all(req %in% names(dt)))
    stop("Missing required columns: ", paste(setdiff(req, names(dt)),
                                             collapse = ", "))

  dt[, (start_col) := as.integer(get(start_col))]
  dt[, (end_col)   := as.integer(get(end_col))]

  dt[, .(
    start = min(get(start_col), na.rm = TRUE),
    end   = max(get(end_col),   na.rm = TRUE)
  ), by = .(chrom = get(chrom_col))][order(chrom)]
}


#' Build a cumulative genome coordinate map
#'
#' @title Cumulative genome map
#' @description Converts a per-chromosome genome extent table into a cumulative
#'   map suitable for plotting all chromosomes on a single linear axis.
#'
#' @param genome_dt A \code{data.table} with columns \code{chrom},
#'   \code{start}, and \code{end} (typically from \code{\link{genome_from_bed}}).
#' @param chrom_order Optional character vector giving the desired chromosome
#'   order. Defaults to the order in \code{genome_dt}.
#'
#' @return A \code{data.table} with the original columns plus \code{length},
#'   \code{cum_start}, \code{cum_end}, \code{xmid}, and \code{y_row}.
#'
#' @examples
#' \dontrun{
#' gmap <- make_genome_map(genome_from_bed(bed12))
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table
#' @export
make_genome_map <- function(genome_dt, chrom_order = NULL) {
  g <- data.table::as.data.table(genome_dt)[, .(chrom, start, end)]
  g[, length := as.numeric(end) - as.numeric(start)]
  if (is.null(chrom_order)) chrom_order <- g$chrom
  chrom_order <- unique(chrom_order[chrom_order %chin% g$chrom])
  g <- g[match(chrom_order, chrom)]

  g[, cum_start := c(0, head(cumsum(length), -1))]
  g[, cum_end   := cum_start + length]
  g[, xmid      := (cum_start + cum_end)/2]
  g[, y_row     := seq_len(.N)]
  g[]
}


#' Map genomic positions to cumulative x coordinates
#'
#' @title Map positions to global x axis
#' @description Given a \code{data.table} with \code{chrom} and a position
#'   column, adds a column \code{x_global} equal to the chromosome's cumulative
#'   start offset plus the local position.
#'
#' @param dt A \code{data.table} with at least columns \code{"chrom"} and the
#'   column named by \code{pos_col}.
#' @param genome_map A genome map from \code{\link{make_genome_map}}.
#' @param pos_col Name of the genomic position column. Default \code{"pos"}.
#'
#' @return The input \code{data.table} with \code{x_global} added (and the
#'   temporary \code{cum_start} column removed).
#'
#' @examples
#' \dontrun{
#' dt[, pos := (start + end) %/% 2L]
#' dt <- map_to_genome_x(dt, gmap)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table
#' @export
map_to_genome_x <- function(dt, genome_map, pos_col = "pos") {
  stopifnot(all(c("chrom", pos_col) %in% names(dt)))
  M <- genome_map[, .(chrom, cum_start)]
  dt[M, on = "chrom", x_global := cum_start + get(pos_col)][, cum_start := NULL]
  dt[]
}


#' Plot genome overview with inter-chromosomal arrows (linear layout)
#'
#' @title Genome arrow plot (linear/stacked layout)
#' @description Draws a stacked genome view where each chromosome is a
#'   horizontal backbone on its own y row. Primary sites are marked with
#'   downward triangles and secondary sites with upward red triangles. Curved
#'   arrows connect primary to secondary locations.
#'
#' @param genome_dt A genome extent table (e.g. from \code{\link{genome_from_bed}}).
#' @param prepped_dt A prepared data table from
#'   \code{prep_primary_secondary_bins()} with columns \code{name},
#'   \code{chrom}, \code{start}, \code{end}, \code{in_primary},
#'   \code{read_has_secondary}, \code{region_midpoint}, \code{sec_bin}.
#' @param chrom_order Optional chromosome order.
#' @param base_size Base font size. Default \code{9}.
#' @param backbone_outline Linewidth for backbone outline. Default \code{5}.
#' @param backbone_fill Linewidth for backbone fill. Default \code{3}.
#' @param arrow_len_pt Arrow head length in pt. Default \code{3}.
#' @param arrow_curvature Arrow curvature. Default \code{0.25}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' plot_genome_arrows_linear(genome_dt, prepped_dt)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table copy
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_curve
#'   scale_x_continuous scale_y_continuous coord_cartesian theme_void theme
#'   element_text margin
#' @importFrom grid arrow unit
#' @export
plot_genome_arrows_linear <- function(genome_dt,
                                      prepped_dt,
                                      chrom_order      = NULL,
                                      base_size        = 9,
                                      backbone_outline = 5,
                                      backbone_fill    = 3,
                                      arrow_len_pt     = 3,
                                      arrow_curvature  = 0.25) {
  g <- data.table::as.data.table(genome_dt)
  p <- data.table::as.data.table(prepped_dt)
  need_g <- c("chrom","start","end")
  if (!all(need_g %in% names(g))) stop("genome_dt must have: chrom, start, end")
  need_p <- c("name","chrom","start","end","in_primary","read_has_secondary",
               "region_midpoint","sec_bin")
  if (!all(need_p %in% names(p)))
    stop("prepped_dt must include: ", paste(need_p, collapse=", "))

  gmap <- make_genome_map(g, chrom_order)
  backbone <- data.table::copy(gmap)[, .(chrom, x0 = cum_start, x1 = cum_end,
                                         y = y_row)]

  prim_mid <- p[in_primary == TRUE,
                .(prim_chrom = chrom[1L],
                  prim_pos   = as.integer(round(mean((start + end) %/% 2L)))),
                by = name]
  prim_mid <- map_to_genome_x(
    prim_mid[, .(name, chrom = prim_chrom, pos = prim_pos)],
    gmap, pos_col = "pos"
  )
  prim_mid[gmap, on = "chrom", y := i.y_row]

  sec_pts <- p[!is.na(region_midpoint),
               .(name, chrom, pos = as.integer(region_midpoint), sec_bin)]
  sec_pts <- unique(sec_pts, by = c("name","chrom","pos","sec_bin"))
  sec_pts <- map_to_genome_x(sec_pts, gmap, pos_col = "pos")
  sec_pts[gmap, on = "chrom", y := i.y_row]

  prim_only <- p[in_primary == TRUE & read_has_secondary == FALSE,
                 .(chrom, pos = as.integer((start + end) %/% 2L))]
  prim_only <- unique(prim_only)
  prim_only <- map_to_genome_x(prim_only, gmap, pos_col = "pos")
  prim_only[gmap, on = "chrom", y := i.y_row]

  sec_mark <- data.table::copy(sec_pts)[, .(chrom, x_global, y)]
  sec_mark <- unique(sec_mark)

  arrows <- if (nrow(sec_pts) && nrow(prim_mid)) {
    merge(
      prim_mid[, .(name, x0 = x_global, y0 = y)],
      sec_pts[,  .(name, x1 = x_global, y1 = y, sec_bin)],
      by = "name", allow.cartesian = TRUE
    )
  } else data.table::data.table()

  y_offset_up   <-  0.14
  y_offset_down <- -0.14
  prim_only[, y_pt := y + y_offset_up]
  sec_mark[,  y_pt := y + y_offset_down]
  if (nrow(arrows)) {
    arrows[, y0 := y0 + 0.10]
    arrows[, y1 := y1 - 0.10]
  }

  a_up <- grid::arrow(type = "closed",
                      length = grid::unit(arrow_len_pt, "pt"))

  breaks <- gmap$xmid
  labs   <- gmap$chrom

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = backbone,
      ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = backbone_outline, lineend = "round", color = "black"
    ) +
    ggplot2::geom_segment(
      data = backbone,
      ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = backbone_fill, lineend = "round", color = "grey70"
    ) +
    ggplot2::geom_point(
      data = prim_only,
      ggplot2::aes(x = x_global, y = y_pt),
      shape = 25, size = 2.4, stroke = 0, fill = "black"
    ) +
    ggplot2::geom_point(
      data = sec_mark,
      ggplot2::aes(x = x_global, y = y_pt),
      shape = 24, size = 2.4, stroke = 0, fill = "#D62728"
    ) +
    { if (nrow(arrows)) ggplot2::geom_curve(
      data = arrows,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1,
                   group = interaction(name, sec_bin)),
      curvature = arrow_curvature, color = "#D62728",
      linewidth = 0.4, lineend = "round", arrow = a_up
    ) } +
    ggplot2::scale_x_continuous(
      breaks = breaks, labels = labs, expand = c(0.001, 0.001)
    ) +
    ggplot2::scale_y_continuous(NULL, breaks = backbone$y,
                                labels = rep("", length(backbone$y))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(
        size = base_size, margin = ggplot2::margin(t = 4)),
      plot.margin  = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
    )
}


#' Plot genome overview with inter-chromosomal arrows (flat layout)
#'
#' @title Genome arrow plot (flat single-row layout)
#' @description Draws a flat, single-row genome view where chromosomes are
#'   placed side-by-side with optional gaps. Curved arcs connect primary to
#'   secondary insertion sites above the backbone.
#'
#' @param genome_dt A genome extent table.
#' @param prepped_dt A prepared data table from
#'   \code{prep_primary_secondary_bins()}.
#' @param gap_bp Numeric. Visual gap between chromosomes in bp. Default
#'   \code{1e6}.
#' @param base_size Base font size. Default \code{8}.
#' @param chrom_outline_size Backbone outline linewidth. Default \code{5}.
#' @param chrom_fill_size Backbone fill linewidth. Default \code{3}.
#' @param arrow_color Colour for arc arrows. Default \code{"#C43B3B"}.
#' @param curvature_abs Arc curvature magnitude (0–1). Default \code{0.25}.
#' @param lift Vertical lift for arc endpoints above the backbone. Default
#'   \code{0.06}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' plot_genome_arrows_flat(genome_dt, prepped_dt)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table copy
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_curve
#'   scale_x_continuous coord_cartesian theme_void theme element_text
#'   element_blank expansion
#' @importFrom grid arrow unit
#' @export
plot_genome_arrows_flat <- function(
    genome_dt,
    prepped_dt,
    gap_bp = 1e6,
    base_size = 8,
    chrom_outline_size = 5,
    chrom_fill_size    = 3,
    arrow_color = "#C43B3B",
    curvature_abs = 0.25,
    lift = 0.06
) {
  g0 <- data.table::as.data.table(genome_dt)
  p0 <- data.table::as.data.table(prepped_dt)

  stopifnot(all(c("chrom","start","end") %in% names(g0)))
  need_cols <- c("chrom","start","end","in_primary","read_has_secondary",
                 "bin_start","bin_end","name","region_midpoint","sec_bin")
  if (!all(need_cols %in% names(p0)))
    stop("prepped_dt must include: ", paste(need_cols, collapse = ", "),
         ". (Use prep_primary_secondary_bins().)")

  g <- data.table::copy(g0)[order(chrom)]
  g[, chr_len := as.numeric(end - start)]
  g[, offset  := c(0, cumsum(head(chr_len + gap_bp, -1)))]
  g[, `:=`(x0 = offset + start, x1 = offset + end, y = 0)]

  offmap    <- setNames(g$offset, g$chrom)
  to_global <- function(chrom, pos) as.numeric(offmap[chrom]) + as.numeric(pos)

  prim_pts <- unique(
    p0[in_primary == TRUE & read_has_secondary == FALSE,
       .(chrom, x = as.integer((bin_start + bin_end) %/% 2L))]
  )
  if (nrow(prim_pts)) prim_pts[, gx := to_global(chrom, x)]

  sec_pts <- unique(
    p0[!in_primary & !is.na(region_midpoint),
       .(chrom, x = as.integer(region_midpoint))]
  )
  if (nrow(sec_pts)) sec_pts[, gx := to_global(chrom, x)]

  prim_mid <- p0[in_primary == TRUE,
                 .(chrom_primary = chrom[1L],
                   x_primary     = as.integer(
                     stats::median((bin_start + bin_end) %/% 2L))),
                 by = name]

  sec_rows <- p0[!in_primary & !is.na(sec_bin) & !is.na(region_midpoint)]
  sec_mid  <- sec_rows[, .(
    chrom_secondary = chrom[1L],
    x_secondary     = as.integer(region_midpoint[1L])
  ), by = name]

  arrows <- merge(prim_mid, sec_mid, by = "name", all = FALSE)
  if (nrow(arrows)) {
    arrows[, `:=`(
      x0 = to_global(chrom_primary,   x_primary),
      x1 = to_global(chrom_secondary, x_secondary),
      y0 = 0,
      y1 = lift
    )]
  }

  ar_lr <- arrows[x1 >= x0]
  ar_rl <- arrows[x1 <  x0]

  chr_labs <- g[, .(chrom, gx = (x0 + x1) / 2)]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = g,
      ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = chrom_outline_size, lineend = "round", color = "black") +
    ggplot2::geom_segment(
      data = g,
      ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = chrom_fill_size, lineend = "round", color = "grey70") +
    { if (nrow(prim_pts)) ggplot2::geom_point(
      data = prim_pts, ggplot2::aes(x = gx, y = 0),
      shape = 25, size = 2.4, stroke = 0, fill = "black"
    ) } +
    { if (nrow(sec_pts)) ggplot2::geom_point(
      data = sec_pts, ggplot2::aes(x = gx, y = 0),
      shape = 24, size = 2.4, stroke = 0, fill = "#C43B3B"
    ) } +
    { if (nrow(ar_lr)) ggplot2::geom_curve(
      data = ar_lr,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature = -curvature_abs,
      color = arrow_color, linewidth = 0.25, lineend = "round",
      arrow = grid::arrow(type = "closed",
                          length = grid::unit(0, "pt"))
    ) } +
    { if (nrow(ar_rl)) ggplot2::geom_curve(
      data = ar_rl,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature = +curvature_abs,
      color = arrow_color, linewidth = 0.25, lineend = "round",
      arrow = grid::arrow(type = "closed",
                          length = grid::unit(0, "pt"))
    ) } +
    ggplot2::scale_x_continuous(
      breaks = chr_labs$gx, labels = chr_labs$chrom,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.01, 0.3), clip = "off") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(
        margin = grid::unit(c(0, 0, 0, 0), "pt")),
      axis.ticks.x = ggplot2::element_blank()
    )

  return(p)
}
