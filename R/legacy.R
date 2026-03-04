# Legacy functions retained for backwards compatibility.
# These may be removed or substantially changed in a future version.

#' Reshape paired reads to wide format (legacy)
#'
#' @title Paired reads in wide format
#' @description Filters a BED table to reads with exactly two alignment
#'   segments on specified chromosomes and reshapes to one row per read pair
#'   with columns for each segment's chromosome, start, end, and strand.
#'
#' @param bed A \code{data.table} from \code{\link{read_bed12}}.
#' @param score_min Minimum mapping score. Default \code{60}.
#' @param length_min Minimum aligned length. Default \code{1000}.
#' @param chroms Integer vector of chromosome numbers. Default \code{1:5}.
#' @param chrom_prefix Character prefix for chromosome names. Default
#'   \code{"Chr"}.
#' @param require_exact_pairs Logical. If \code{TRUE} (default), keep only
#'   reads with exactly 2 segments.
#'
#' @return A wide-format \code{data.table} with columns \code{name},
#'   \code{chrom_1}, \code{chrom_2}, \code{start_1}, \code{start_2},
#'   \code{end_1}, \code{end_2}, \code{strand_1}, \code{strand_2}.
#'
#' @family legacy
#' @importFrom data.table as.data.table dcast rowid
#' @export
make_pairs_wide_strand <- function(
    bed, score_min = 60, length_min = 1000,
    chroms = 1:5, chrom_prefix = "Chr",
    require_exact_pairs = TRUE
) {
  stopifnot(all(c("chrom","start","end","name","score",
                  "aligned_length","strand") %in% names(bed)))
  dt <- data.table::as.data.table(bed)

  x <- dt[
    score >= score_min & aligned_length > length_min &
      chrom %chin% paste0(chrom_prefix, chroms)
  ][, chrom := as.integer(sub(paste0("^", chrom_prefix), "", chrom))]

  if (require_exact_pairs) {
    x <- x[, if (.N == 2) .SD, by = name]
  } else {
    x <- x[order(name, chrom, start)][, .SD[1:2], by = name]
  }

  x <- x[order(name, chrom, start)][,
         site := data.table::rowid(name)]
  data.table::dcast(x, name ~ site,
                    value.var = c("chrom","start","end","strand"))
}


#' Cluster read pairs into big-window bins (legacy)
#'
#' @title Big-window cluster pairs
#' @description Bins paired-read entries into coarse genomic windows and
#'   labels each pair with a cluster string of the form
#'   \code{"ChrA_X-ChrB_Y"}.
#'
#' @param dt A wide-format \code{data.table} from
#'   \code{\link{make_pairs_wide_strand}}.
#' @param tol1 Window size in bp for the first segment. Default \code{100000L}.
#' @param tol2 Window size in bp for the second segment. Default \code{tol1}.
#' @param label_prefix Chromosome label prefix. Default \code{"Chr"}.
#' @param min_support Minimum number of reads per cluster. Default \code{1L}.
#'
#' @return The input \code{data.table} with \code{cluster} and \code{clusterN}
#'   columns added.
#'
#' @family legacy
#' @importFrom data.table copy is.data.table
#' @export
cluster_pairs_bigwindows <- function(dt,
                                     tol1 = 100000L, tol2 = tol1,
                                     label_prefix = "Chr",
                                     min_support = 1L) {
  stopifnot(data.table::is.data.table(dt))
  req <- c("chrom_1","chrom_2","start_1","end_1","start_2","end_2")
  if (!all(req %in% names(dt)))
    stop("Missing required columns: ",
         paste(setdiff(req, names(dt)), collapse = ", "))

  z <- data.table::copy(dt)[, `:=`(
    mid1 = as.integer((start_1 + end_1) %/% 2L),
    mid2 = as.integer((start_2 + end_2) %/% 2L)
  )]

  half1 <- tol1 %/% 2L; half2 <- tol2 %/% 2L
  z[, bin1_idx    := (mid1 + half1) %/% tol1]
  z[, bin2_idx    := (mid2 + half2) %/% tol2]
  z[, cluster_key := paste(chrom_1, bin1_idx, chrom_2, bin2_idx, sep = "||")]

  lab <- z[, .(
    support = .N,
    chrom1  = as.integer(names(sort(table(chrom_1), decreasing = TRUE))[1]),
    chrom2  = as.integer(names(sort(table(chrom_2), decreasing = TRUE))[1]),
    X1      = as.integer(round(mean(mid1))),
    X2      = as.integer(round(mean(mid2)))
  ), by = cluster_key
  ][, cluster := sprintf("%s%d_%d-%s%d_%d",
                         label_prefix, chrom1, X1,
                         label_prefix, chrom2, X2)
  ][, .(cluster_key, cluster, support)]

  z <- lab[z, on = "cluster_key"]
  if (min_support > 1L) z[support < min_support, cluster := NA_character_]
  z[, c("mid1","mid2","bin1_idx","bin2_idx","cluster_key","support") := NULL]
  z[, clusterN := .N, by = cluster]
  z[]
}


#' Fetch reads near a genomic region (legacy)
#'
#' @title Reads near a genomic region
#' @description Identifies reads overlapping a primary window around a target
#'   position and their secondary alignments in extended windows.
#'
#' @param bed A \code{data.table} from \code{\link{read_bed12}}.
#' @param region Character region string of the form \code{"ChrX_POS"}.
#' @param window_kb Numeric. Half-width of the primary window in kb. Default
#'   \code{10}.
#' @param sec_window_kb Numeric. Half-width of secondary windows in kb. If
#'   \code{NULL}, uses \code{window_kb}.
#' @param merge_secondary Logical. Merge overlapping secondary windows. Default
#'   \code{TRUE}.
#'
#' @return A named list with elements \code{names_primary},
#'   \code{names_secondary}, \code{hits_primary}, \code{hits_secondary},
#'   \code{secondary_windows}, and \code{all_by_name}.
#'
#' @family legacy
#' @importFrom data.table as.data.table setkey setorder data.table
#' @export
reads_near_region <- function(
    bed, region, window_kb = 10, sec_window_kb = NULL, merge_secondary = TRUE
) {
  if (!inherits(bed, "data.table")) bed <- data.table::as.data.table(bed)
  toks <- strsplit(region, "_", fixed = TRUE)[[1]]
  if (length(toks) < 2L) stop("region must look like 'ChrX_POS'")
  chrom_str <- toks[1]
  pos_bp    <- suppressWarnings(as.integer(toks[length(toks)]))
  if (is.na(pos_bp)) stop("Failed to parse position from: ", region)

  w  <- as.integer(round(window_kb * 1000))
  sw <- as.integer(round((if (is.null(sec_window_kb)) window_kb
                          else sec_window_kb) * 1000))

  win_start <- max(0L, pos_bp - w)
  win_end   <- pos_bp + w

  hits_primary <- bed[chrom == chrom_str &
                        start <= win_end & end >= win_start]
  if (nrow(hits_primary) == 0L) {
    out_empty <- bed[0]
    attr(out_empty, "region")    <- region
    attr(out_empty, "window_bp") <- c(chrom = chrom_str,
                                      start = win_start, end = win_end)
    return(list(
      names_primary     = character(0),
      names_secondary   = character(0),
      hits_primary      = bed[0],
      hits_secondary    = bed[0],
      secondary_windows = data.table::data.table(
        chrom = character(), start = integer(), end = integer()),
      all_by_name       = out_empty
    ))
  }

  names_primary       <- unique(hits_primary$name)
  all_primary_segs    <- bed[name %in% names_primary]
  all_primary_segs[, in_target := (chrom == chrom_str &
                                     start <= win_end & end >= win_start)]
  secondary_segs <- all_primary_segs[in_target == FALSE]

  secondary_windows <- data.table::data.table(
    chrom = character(), start = integer(), end = integer())
  if (nrow(secondary_segs)) {
    secondary_windows <- secondary_segs[
      , .(chrom, start = pmax(0L, as.integer(start - sw)),
          end = as.integer(end + sw))]
    if (merge_secondary) {
      data.table::setorder(secondary_windows, chrom, start, end)
      secondary_windows[, grp := {
        grp <- integer(.N); cur <- 0L; last_end <- -1L
        for (i in seq_len(.N)) {
          if (i == 1L || start[i] > last_end) cur <- cur + 1L
          grp[i] <- cur
          last_end <- max(last_end, end[i])
        }
        grp
      }, by = chrom]
      secondary_windows <- secondary_windows[,
        .(start = min(start), end = max(end)),
        by = .(chrom, grp)][, grp := NULL][]
    }
  }

  hits_secondary <- bed[0]; names_secondary <- character(0)
  if (nrow(secondary_windows)) {
    data.table::setkey(secondary_windows, chrom)
    data.table::setkey(bed, chrom)
    cand <- bed[secondary_windows, nomatch = 0L, allow.cartesian = TRUE]
    hits_secondary <- cand[start <= end & end >= start,
                           .(chrom, start, end, name, score, strand,
                             thickStart, thickEnd, itemRgb, blockCount,
                             blockSizes, blockStarts, is_split,
                             aligned_length, nameN)]
    if (nrow(hits_secondary)) names_secondary <- unique(hits_secondary$name)
  }

  all_names   <- unique(c(names_primary, names_secondary))
  all_by_name <- bed[name %in% all_names]
  attr(all_by_name, "region")    <- region
  attr(all_by_name, "window_bp") <- c(chrom = chrom_str,
                                      start = win_start, end = win_end)
  list(
    names_primary     = names_primary,
    names_secondary   = setdiff(all_names, names_primary),
    hits_primary      = hits_primary[],
    hits_secondary    = hits_secondary[],
    secondary_windows = secondary_windows[],
    all_by_name       = all_by_name[]
  )
}


#' Prepare primary/secondary alignment bins for plotting (legacy)
#'
#' @title Prepare primary/secondary bins
#' @description Takes the output of \code{\link{reads_near_region}} and
#'   organises reads into primary (target window) and secondary (best off-target
#'   bin) segments for downstream plotting.
#'
#' @param cr Output from \code{\link{reads_near_region}}.
#' @param bin_kb Integer. Bin width in kb. Default \code{100L}.
#'
#' @return A \code{data.table} suitable for \code{\link{plot_primary_secondary_bins}}
#'   and the genome arrow plot functions.
#'
#' @family legacy
#' @importFrom data.table as.data.table setorder setnames
#' @importFrom scales comma
#' @export
prep_primary_secondary_bins <- function(cr, bin_kb = 100L) {
  dt <- data.table::as.data.table(cr$all_by_name)
  if (!nrow(dt)) stop("No rows in cr$all_by_name.")

  primary_names <- unique(cr$hits_primary$name)
  if (!length(primary_names)) stop("No primary hits found in cr$hits_primary.")
  dt <- dt[name %chin% primary_names]
  if (!nrow(dt)) stop("After restricting to primary names, nothing to prepare.")

  wb <- attr(dt, "window_bp")
  if (is.null(wb)) stop("Missing 'window_bp' attribute on cr$all_by_name.")
  tgt_chr   <- as.character(wb[["chrom"]])
  tgt_start <- as.integer(wb[["start"]])
  tgt_end   <- as.integer(wb[["end"]])

  dt[, in_primary := (chrom == tgt_chr & start <= tgt_end & end >= tgt_start)]

  bin_bp <- as.integer(bin_kb * 1000L)
  dt[, mid       := as.integer((start + end) %/% 2L)]
  dt[, bin_start := (mid %/% bin_bp) * bin_bp]
  dt[, bin_end   := bin_start + bin_bp]
  dt[, bin_label := sprintf("%s %s\u2013%s kb",
                            chrom,
                            scales::comma(bin_start/1000),
                            scales::comma(bin_end/1000))]

  tgt_mid    <- as.integer((tgt_start + tgt_end) %/% 2L)
  tgt_bstart <- (tgt_mid %/% bin_bp) * bin_bp
  tgt_bend   <- tgt_bstart + bin_bp
  target_bin_label <- sprintf("%s %s\u2013%s kb", tgt_chr,
                              scales::comma(tgt_bstart/1000),
                              scales::comma(tgt_bend/1000))

  sec_counts <- dt[in_primary == FALSE, .N, by = .(name, bin_label)]
  if (nrow(sec_counts)) {
    data.table::setorder(sec_counts, name, -N, bin_label)
    best_sec <- sec_counts[, .SD[1], by = name]
    data.table::setnames(best_sec, c("bin_label","N"), c("sec_bin","secN"))
  } else {
    best_sec <- data.table::data.table(
      name    = character(0),
      sec_bin = character(0),
      secN    = integer(0))
  }

  bin_mid  <- unique(dt[, .(bin_label,
                             region_midpoint = as.integer((bin_start + bin_end) %/% 2L))])
  best_sec <- merge(best_sec, bin_mid,
                    by.x = "sec_bin", by.y = "bin_label", all.x = TRUE)

  primary_meta <- dt[in_primary == TRUE,
                     .(primary_start = min(start, na.rm = TRUE)), by = name]

  tmp       <- best_sec[dt, on = "name"]
  sec_start <- tmp[!is.na(sec_bin) & bin_label == sec_bin,
                   .(sec_start = min(start, na.rm = TRUE)), by = name]

  all_reads          <- unique(dt$name)
  reads_with_sec     <- unique(best_sec$name)
  reads_primary_only <- setdiff(all_reads, reads_with_sec)

  sec_pop <- if (nrow(best_sec))
    best_sec[, .(reads = .N), by = sec_bin][order(-reads, sec_bin)]
  else data.table::data.table(sec_bin = character(0), reads = integer(0))

  prim_only_order <- primary_meta[
    name %chin% reads_primary_only][order(primary_start, name), name]

  best_sec_rank <- merge(best_sec, sec_start, by = "name", all.x = TRUE)
  best_sec_rank <- merge(best_sec_rank, primary_meta, by = "name", all.x = TRUE)
  best_sec_rank[is.na(sec_start), sec_start := primary_start]

  y_order <- c(
    prim_only_order,
    unlist(lapply(sec_pop$sec_bin, function(k) {
      best_sec_rank[sec_bin == k][order(sec_start, name), name]
    }))
  )
  y_order <- unique(y_order[y_order %chin% all_reads])

  dt        <- best_sec[dt, on = "name"]
  plot_dt   <- dt[in_primary == TRUE | (!is.na(sec_bin) & bin_label == sec_bin)]
  if (!nrow(plot_dt))
    stop("No primary or chosen-secondary segments after filtering.")

  plot_dt[, is_secondary_segment := (!in_primary) & (!is.na(sec_bin)) &
            (bin_label == sec_bin)]
  bins_with_secondary <- plot_dt[
    , .(has_secondary = any(is_secondary_segment)), by = bin_label][
      has_secondary == TRUE, bin_label]
  keep_bins <- unique(c(as.character(bins_with_secondary), target_bin_label))
  plot_dt   <- plot_dt[bin_label %chin% keep_bins]
  if (!nrow(plot_dt))
    stop("After keeping secondary/target bins, nothing remains.")

  plot_dt[, y := (length(y_order) + 1L) - match(name, y_order)]

  has_secondary_map <- setNames(rep(FALSE, length(all_reads)), all_reads)
  has_secondary_map[reads_with_sec] <- TRUE
  plot_dt[, read_has_secondary := has_secondary_map[name]]

  data.table::setorder(plot_dt, chrom, bin_start)
  facet_levels <- unique(plot_dt$bin_label)
  plot_dt[, bin_label := factor(bin_label, levels = facet_levels)]

  list_cols <- names(plot_dt)[vapply(plot_dt, is.list, TRUE)]
  if (length(list_cols)) plot_dt[, (list_cols) := NULL]

  drop_extra <- intersect(
    c("itemRgb","thickStart","thickEnd","blockCount","aligned_length",
      "nameN","score","strand"),
    names(plot_dt))
  if (length(drop_extra)) plot_dt[, (drop_extra) := NULL]

  keep_cols <- c("name","chrom","start","end",
                 "in_primary","read_has_secondary",
                 "sec_bin","secN","region_midpoint",
                 "mid","bin_start","bin_end","bin_label",
                 "is_secondary_segment","y")
  keep_cols <- intersect(keep_cols, names(plot_dt))
  plot_dt   <- plot_dt[, ..keep_cols]

  attr(plot_dt, "bin_kb")           <- bin_kb
  attr(plot_dt, "target_bin_label") <- target_bin_label
  attr(plot_dt, "facet_levels")     <- facet_levels
  plot_dt[]
}


#' Plot primary/secondary bins (legacy)
#'
#' @title Plot primary/secondary alignment bins
#' @description Plots reads prepared by \code{\link{prep_primary_secondary_bins}}
#'   as horizontal segments faceted by genomic bin, with reads coloured by
#'   whether they have a secondary alignment.
#'
#' @param prepped_dt Output from \code{\link{prep_primary_secondary_bins}}.
#' @param base_size Base font size. Default \code{6}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @family legacy
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes geom_segment facet_grid scale_x_continuous
#'   scale_color_manual guides theme_bw theme element_blank element_text
#'   element_rect margin
#' @export
plot_primary_secondary_bins <- function(prepped_dt, base_size = 6) {
  plot_dt <- data.table::as.data.table(prepped_dt)
  if (!nrow(plot_dt)) stop("Empty prepared data provided.")

  breaks_1kb <- function(lims) {
    from <- floor(lims[1] / 1000) * 1000
    to   <- ceiling(lims[2] / 1000) * 1000
    seq(from, to, by = 1000L)
  }

  ggplot2::ggplot(plot_dt) +
    ggplot2::geom_segment(
      ggplot2::aes(y = y, yend = y, x = start, xend = end,
                   color = read_has_secondary),
      linewidth = 0.5, lineend = "round") +
    ggplot2::facet_grid(~ bin_label, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_continuous(
      position     = "top",
      breaks       = breaks_1kb,
      labels       = NULL,
      minor_breaks = NULL,
      expand       = c(.01, .01)
    ) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#D62728", `FALSE` = "grey60")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.background    = ggplot2::element_rect(fill = NA, colour = NA),
      panel.grid          = ggplot2::element_blank(),
      strip.background    = ggplot2::element_blank(),
      strip.text          = ggplot2::element_text(face = "plain", angle = 90),
      strip.placement     = "outside",
      axis.title          = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.title.x.bottom = ggplot2::element_blank()
    )
}
