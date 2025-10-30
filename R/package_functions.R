# functions.R
library(data.table)
library(ggplot2)
library(polymorphology2)

# ---- I/O ----

# Minimal BED12 reader (single-block friendly)
# - Reads as plain columns (no list-cols)
# - If any blockCount > 1, just message() but DO NOT stop
# - Adds: is_split, aligned_length, nameN
read_bed12 <- function(path, key_on = c("chrom","start","end")) {
  cols <- c("chrom","start","end","name","score","strand",
            "thickStart","thickEnd","itemRgb","blockCount",
            "blockSizes","blockStarts")
  cc <- c("character","integer","integer","character","integer","character",
          "integer","integer","character","integer",
          # keep these as character to avoid list-cols
          "character","character")

  dt <- fread(path, sep = "\t", header = FALSE, fill = TRUE, quote = "",
              col.names = cols, colClasses = cc, showProgress = FALSE)

  # Soft notice only if multi-block lines exist
  n_multi <- sum(dt$blockCount > 1L, na.rm = TRUE)
  if (n_multi > 0L) {
    message(sprintf("Note: detected %d row(s) with blockCount > 1 (multi-block BED12).",
                    n_multi))
  }

  # Plain derived columns (single-block safe)
  dt[, is_split       := blockCount > 1L]
  dt[, aligned_length := pmax(0L, end - start)]
  dt[, nameN          := .N, by = name]

  if (!is.null(key_on)) data.table::setkeyv(dt, key_on)
  dt[]
}

# ---- Pairing & clustering used by your pipeline ----
make_pairs_wide_strand <- function(
    bed, score_min = 60, length_min = 1000,
    chroms = 1:5, chrom_prefix = "Chr",
    require_exact_pairs = TRUE
) {
  stopifnot(all(c("chrom","start","end","name","score","aligned_length","strand") %in% names(bed)))
  dt <- as.data.table(bed)

  x <- dt[
    score >= score_min & aligned_length > length_min &
      chrom %chin% paste0(chrom_prefix, chroms)
  ][ , chrom := as.integer(sub(paste0("^", chrom_prefix), "", chrom)) ]

  if (require_exact_pairs) {
    x <- x[, if (.N == 2) .SD, by = name]
  } else {
    x <- x[order(name, chrom, start)][, .SD[1:2], by = name]
  }

  x <- x[order(name, chrom, start)][, site := data.table::rowid(name)]

  dcast(x, name ~ site, value.var = c("chrom","start","end","strand"))
}

cluster_pairs_bigwindows <- function(dt,
                                     tol1 = 100000L, tol2 = tol1,
                                     label_prefix = "Chr",
                                     min_support = 1L) {
  stopifnot(inherits(dt, "data.table"))
  req <- c("chrom_1","chrom_2","start_1","end_1","start_2","end_2")
  if (!all(req %in% names(dt))) stop("Missing required columns: ", paste(setdiff(req, names(dt)), collapse = ", "))

  z <- copy(dt)[ , `:=`(
    mid1 = as.integer((start_1 + end_1) %/% 2L),
    mid2 = as.integer((start_2 + end_2) %/% 2L)
  )]

  half1 <- tol1 %/% 2L; half2 <- tol2 %/% 2L
  z[, bin1_idx := (mid1 + half1) %/% tol1]
  z[, bin2_idx := (mid2 + half2) %/% tol2]
  z[, cluster_key := paste(chrom_1, bin1_idx, chrom_2, bin2_idx, sep = "||")]

  lab <- z[, .(
    support = .N,
    chrom1  = as.integer(names(sort(table(chrom_1), decreasing = TRUE))[1]),
    chrom2  = as.integer(names(sort(table(chrom_2), decreasing = TRUE))[1]),
    X1      = as.integer(round(mean(mid1))),
    X2      = as.integer(round(mean(mid2)))
  ), by = cluster_key
  ][ , cluster := sprintf("%s%d_%d-%s%d_%d", label_prefix, chrom1, X1, label_prefix, chrom2, X2)
  ][ , .(cluster_key, cluster, support)]

  z <- lab[z, on = "cluster_key"]
  if (min_support > 1L) z[support < min_support, cluster := NA_character_]
  z[, c("mid1","mid2","bin1_idx","bin2_idx","cluster_key","support") := NULL]
  z[, clusterN := .N, by = cluster]
  z[]
}

# ---- Region fetch / prep / plot ----
reads_near_region <- function(
    bed, region, window_kb = 10, sec_window_kb = NULL, merge_secondary = TRUE
) {
  if (!inherits(bed, "data.table")) bed <- as.data.table(bed)
  toks <- strsplit(region, "_", fixed = TRUE)[[1]]
  if (length(toks) < 2L) stop("region must look like 'ChrX_POS'")
  chrom_str <- toks[1]
  pos_bp    <- suppressWarnings(as.integer(toks[length(toks)]))
  if (is.na(pos_bp)) stop("Failed to parse position from: ", region)

  w  <- as.integer(round(window_kb * 1000))
  sw <- as.integer(round((if (is.null(sec_window_kb)) window_kb else sec_window_kb) * 1000))

  win_start <- max(0L, pos_bp - w)
  win_end   <- pos_bp + w

  hits_primary <- bed[chrom == chrom_str & start <= win_end & end >= win_start]
  if (nrow(hits_primary) == 0L) {
    out_empty <- bed[0]
    attr(out_empty, "region")    <- region
    attr(out_empty, "window_bp") <- c(chrom = chrom_str, start = win_start, end = win_end)
    return(list(
      names_primary     = character(0),
      names_secondary   = character(0),
      hits_primary      = bed[0],
      hits_secondary    = bed[0],
      secondary_windows = data.table(chrom = character(), start = integer(), end = integer()),
      all_by_name       = out_empty
    ))
  }

  names_primary <- unique(hits_primary$name)
  all_primary_segments <- bed[name %in% names_primary]
  all_primary_segments[, in_target := (chrom == chrom_str & start <= win_end & end >= win_start)]
  secondary_segs <- all_primary_segments[in_target==F]

  secondary_windows <- data.table(chrom = character(), start = integer(), end = integer())
  if (nrow(secondary_segs)) {
    secondary_windows <- secondary_segs[
      , .(chrom, start = pmax(0L, as.integer(start - sw)), end = as.integer(end + sw))]
    if (merge_secondary) {
      setorder(secondary_windows, chrom, start, end)
      secondary_windows[, grp := {
        grp <- integer(.N); cur <- 0L; last_end <- -1L
        for (i in seq_len(.N)) {
          if (i == 1L || start[i] > last_end) cur <- cur + 1L
          grp[i] <- cur
          last_end <- max(last_end, end[i])
        }
        grp
      }, by = chrom]
      secondary_windows <- secondary_windows[, .(start = min(start), end = max(end)), by = .(chrom, grp)][, grp := NULL][]
    }
  }

  hits_secondary <- bed[0]; names_secondary <- character(0)
  if (nrow(secondary_windows)) {
    setkey(secondary_windows, chrom); setkey(bed, chrom)
    cand <- bed[secondary_windows, nomatch = 0L, allow.cartesian = TRUE]
    hits_secondary <- cand[start <= end & end >= start,
                           .(chrom,start,end,name,score,strand,thickStart,thickEnd,
                             itemRgb,blockCount,blockSizes,blockStarts,is_split,
                             aligned_length,blockEnds,nameN)]
    if (nrow(hits_secondary)) names_secondary <- unique(hits_secondary$name)
  }

  all_names   <- unique(c(names_primary, names_secondary))
  all_by_name <- bed[name %in% all_names]
  attr(all_by_name, "region")    <- region
  attr(all_by_name, "window_bp") <- c(chrom = chrom_str, start = win_start, end = win_end)
  list(
    names_primary     = names_primary,
    names_secondary   = setdiff(all_names, names_primary),
    hits_primary      = hits_primary[],
    hits_secondary    = hits_secondary[],
    secondary_windows = secondary_windows[],
    all_by_name       = all_by_name[]
  )
}

prep_primary_secondary_bins <- function(cr, bin_kb = 100L) {
  dt <- data.table::as.data.table(cr$all_by_name)
  if (!nrow(dt)) stop("No rows in cr$all_by_name.")

  # Restrict to reads whose names are in PRIMARY hits
  primary_names <- unique(cr$hits_primary$name)
  if (!length(primary_names)) stop("No primary hits found in cr$hits_primary.")
  dt <- dt[name %chin% primary_names]
  if (!nrow(dt)) stop("After restricting to primary names, nothing to prepare.")

  # Target (primary) window from reads_near_region()
  wb <- attr(dt, "window_bp")
  if (is.null(wb)) stop("Missing 'window_bp' attribute on cr$all_by_name.")
  tgt_chr   <- as.character(wb[["chrom"]])
  tgt_start <- as.integer(wb[["start"]])
  tgt_end   <- as.integer(wb[["end"]])

  # Primary flag
  dt[, in_primary := (chrom == tgt_chr & start <= tgt_end & end >= tgt_start)]

  # Midpoint bins (facets)
  bin_bp <- as.integer(bin_kb * 1000L)
  dt[, mid := as.integer((start + end) %/% 2L)]
  dt[, bin_start := (mid %/% bin_bp) * bin_bp]
  dt[, bin_end   := bin_start + bin_bp]
  dt[, bin_label := sprintf("%s %s–%s kb",
                            chrom,
                            scales::comma(bin_start/1000),
                            scales::comma(bin_end/1000))]

  # Target bin label — always keep this bin
  tgt_mid    <- as.integer((tgt_start + tgt_end) %/% 2L)
  tgt_bstart <- (tgt_mid %/% bin_bp) * bin_bp
  tgt_bend   <- tgt_bstart + bin_bp
  target_bin_label <- sprintf("%s %s–%s kb", tgt_chr,
                              scales::comma(tgt_bstart/1000),
                              scales::comma(tgt_bend/1000))

  # Best secondary bin per read (exclude primary segments)
  sec_counts <- dt[in_primary == FALSE, .N, by = .(name, bin_label)]
  if (nrow(sec_counts)) {
    data.table::setorder(sec_counts, name, -N, bin_label)
    best_sec <- sec_counts[, .SD[1], by = name]       # name, bin_label, N
    data.table::setnames(best_sec, c("bin_label","N"), c("sec_bin","secN"))
  } else {
    best_sec <- data.table::data.table(name = character(0), sec_bin = character(0), secN = integer(0))
  }

  # Numeric midpoint for each bin_label; map to each read's sec_bin -> region_midpoint
  bin_mid <- unique(dt[, .(bin_label, region_midpoint = as.integer((bin_start + bin_end) %/% 2L))])
  best_sec <- merge(best_sec, bin_mid, by.x = "sec_bin", by.y = "bin_label", all.x = TRUE)

  # Per-read primary start (for primary-only and fallback)
  primary_meta <- dt[in_primary == TRUE, .(primary_start = min(start, na.rm = TRUE)), by = name]

  # Per-read secondary start inside chosen secondary bin
  tmp <- best_sec[dt, on = "name"]
  sec_start <- tmp[!is.na(sec_bin) & bin_label == sec_bin,
                   .(sec_start = min(start, na.rm = TRUE)), by = name]

  # Reads with/without secondary among the primary-name set
  all_reads          <- unique(dt$name)
  reads_with_sec     <- unique(best_sec$name)
  reads_primary_only <- setdiff(all_reads, reads_with_sec)

  # Popularity of secondary bins (group order)
  sec_pop <- if (nrow(best_sec)) best_sec[, .(reads = .N), by = sec_bin][order(-reads, sec_bin)]
  else data.table::data.table(sec_bin = character(0), reads = integer(0))

  # ---- Build Y order ----
  prim_only_order <- primary_meta[name %chin% reads_primary_only][order(primary_start, name), name]

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

  # Keep ONLY: primary segments OR segments in the read’s chosen secondary bin
  dt <- best_sec[dt, on = "name"]  # adds sec_bin/secN/region_midpoint (NA for primary-only)
  plot_dt <- dt[in_primary == TRUE | (!is.na(sec_bin) & bin_label == sec_bin)]
  if (!nrow(plot_dt)) stop("No primary or chosen-secondary segments to prepare after filtering to primary names.")

  # Identify bins with any secondary segment; keep those + the target bin
  plot_dt[, is_secondary_segment := (!in_primary) & (!is.na(sec_bin)) & (bin_label == sec_bin)]
  bins_with_secondary <- plot_dt[, .(has_secondary = any(is_secondary_segment)), by = bin_label][
    has_secondary == TRUE, bin_label]
  keep_bins <- unique(c(as.character(bins_with_secondary), target_bin_label))
  plot_dt <- plot_dt[bin_label %chin% keep_bins]
  if (!nrow(plot_dt)) stop("After keeping secondary bins and target bin, nothing remains to prepare.")

  # Final Y mapping (primary-only first at TOP)
  plot_dt[, y := (length(y_order) + 1L) - match(name, y_order)]

  # Color by read: red if the read has any secondary, gray if primary-only
  has_secondary_map <- setNames(rep(FALSE, length(all_reads)), all_reads)
  has_secondary_map[reads_with_sec] <- TRUE
  plot_dt[, read_has_secondary := has_secondary_map[name]]

  # Stable facet order
  data.table::setorder(plot_dt, chrom, bin_start)
  facet_levels <- unique(plot_dt$bin_label)
  plot_dt[, bin_label := factor(bin_label, levels = facet_levels)]

  # ---- SLIM DOWN: drop list columns so this is a clean data.table ----
  list_cols <- names(plot_dt)[vapply(plot_dt, is.list, TRUE)]
  if (length(list_cols)) plot_dt[, (list_cols) := NULL]

  # Also drop noisy BED12 columns you don’t use for plotting/ranking (optional)
  drop_extra <- intersect(
    c("itemRgb","thickStart","thickEnd","blockCount","aligned_length","nameN","score","strand"),
    names(plot_dt)
  )
  if (length(drop_extra)) plot_dt[, (drop_extra) := NULL]

  # Keep useful columns only (order them nicely)
  keep_cols <- c("name","chrom","start","end",
                 "in_primary","read_has_secondary",
                 "sec_bin","secN","region_midpoint",
                 "mid","bin_start","bin_end","bin_label",
                 "is_secondary_segment","y")
  keep_cols <- intersect(keep_cols, names(plot_dt))
  plot_dt <- plot_dt[, ..keep_cols]

  # Return with attrs
  attr(plot_dt, "bin_kb")           <- bin_kb
  attr(plot_dt, "target_bin_label") <- target_bin_label
  attr(plot_dt, "facet_levels")     <- facet_levels
  plot_dt[]
}

plot_primary_secondary_bins <- function(prepped_dt, base_size = 6) {
  plot_dt <- as.data.table(prepped_dt)
  if (!nrow(plot_dt)) stop("Empty prepared data provided.")

  breaks_1kb <- function(lims) {
    from <- floor(lims[1] / 1000) * 1000
    to   <- ceiling(lims[2] / 1000) * 1000
    seq(from, to, by = 1000L)
  }

  ggplot(plot_dt) +
    geom_segment(aes(y = y, yend = y, x = start, xend = end,
                     color = read_has_secondary),
                 linewidth = 0.5, lineend = "round") +
    facet_grid(~ bin_label, scales = "free_x", space = "free_x") +
    scale_x_continuous(
      position = "top",
      breaks = breaks_1kb,
      labels = NULL,
      minor_breaks = NULL,
      expand = c(.01, .01)
    ) +
    scale_color_manual(values = c(`TRUE` = "#D62728", `FALSE` = "grey60")) +
    guides(color = "none") +
    theme_bw(base_size = base_size) +
    theme(
      panel.background    = element_rect(fill = NA, colour = NA),
      panel.grid          = element_blank(),
      strip.background    = element_blank(),
      strip.text          = element_text(face = "plain", angle = 90),
      strip.placement     = "outside",
      axis.title          = element_blank(),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank()
    )
}



# 1) Genome ranges per chromosome (min start / max end)
genome_from_bed <- function(bed,
                            chrom_col = "chrom",
                            start_col = "start",
                            end_col   = "end") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")
  dt <- if (data.table::is.data.table(bed)) bed else data.table::as.data.table(bed)

  req <- c(chrom_col, start_col, end_col)
  if (!all(req %in% names(dt))) {
    stop("Missing required columns: ", paste(setdiff(req, names(dt)), collapse = ", "))
  }

  dt[, (start_col) := as.integer(get(start_col))]
  dt[, (end_col)   := as.integer(get(end_col))]

  start_col="start"
  end_col="end"
  dt[, .(
    start = min(get(start_col), na.rm = TRUE),
    end   = max(get(end_col),   na.rm = TRUE)
  ), by = .(chrom = get(chrom_col))][order(chrom)]
}


# Build a cumulative map so each chromosome occupies a contiguous x-range
# genome_dt must have: chrom, start, end (start<=end; typically start=0/min)
make_genome_map <- function(genome_dt, chrom_order = NULL) {
  library(data.table)
  g <- as.data.table(genome_dt)[, .(chrom, start, end)]
  g[, length := as.numeric(end) - as.numeric(start)]
  if (is.null(chrom_order)) chrom_order <- g$chrom
  chrom_order <- unique(chrom_order[chrom_order %chin% g$chrom])
  g <- g[match(chrom_order, chrom)]

  g[, cum_start := c(0, head(cumsum(length), -1))]
  g[, cum_end   := cum_start + length]
  g[, xmid      := (cum_start + cum_end)/2]
  g[, y_row     := seq_len(.N)]  # one row per chromosome (top to bottom)
  g[]
}

# Map a data.table with (chrom, pos) to global x = cum_start + pos
# `pos_col` is a single column name containing base-pair coordinates.
map_to_genome_x <- function(dt, genome_map, pos_col = "pos") {
  library(data.table)
  stopifnot(all(c("chrom", pos_col) %in% names(dt)))
  M <- genome_map[, .(chrom, cum_start)]
  dt[M, on="chrom", x_global := cum_start + get(pos_col)][, cum_start := NULL]
  dt[]
}


# Plot inter-chrom arrows:
# - Chromosome backbones on separate y rows (no facet)
# - Primary-only sites as black ▼
# - Secondary sites as red ▲
# - Curved arrows primary -> secondary across rows (inter- & intra-chrom)
plot_genome_arrows_linear <- function(genome_dt,
                                      prepped_dt,
                                      chrom_order      = NULL,
                                      base_size        = 9,
                                      backbone_outline = 5,
                                      backbone_fill    = 3,
                                      arrow_len_pt     = 3,
                                      arrow_curvature  = 0.25) {
  library(data.table)
  library(ggplot2)

  # --- checks
  g <- as.data.table(genome_dt)
  p <- as.data.table(prepped_dt)
  need_g <- c("chrom","start","end")
  if (!all(need_g %in% names(g))) stop("genome_dt must have: chrom, start, end")
  need_p <- c("name","chrom","start","end","in_primary","read_has_secondary","region_midpoint","sec_bin")
  if (!all(need_p %in% names(p))) stop(
    "prepped_dt must include: ",
    paste(need_p, collapse=", "),
    " (from prep_primary_secondary_bins())"
  )

  # --- genome mapping (cumulative x + chrom rows)
  gmap <- make_genome_map(g, chrom_order)
  # backbone data with y rows
  backbone <- copy(gmap)[, .(chrom, x0 = cum_start, x1 = cum_end, y = y_row)]

  # --- primary midpoint per read (on target chrom)
  prim_mid <- p[in_primary == TRUE,
                .(prim_chrom = first(chrom),
                  prim_pos   = as.integer(round(mean((start + end) %/% 2L)))),
                by = name]
  # map to global x and y row
  prim_mid <- map_to_genome_x(
    prim_mid[, .(name, chrom = prim_chrom, pos = prim_pos)],
    gmap, pos_col = "pos"
  )
  prim_mid[gmap, on="chrom", y := i.y_row]

  # --- secondary midpoints per read (chosen bin); may be multiple per read
  sec_pts <- p[!is.na(region_midpoint),
               .(name, chrom, pos = as.integer(region_midpoint), sec_bin)]
  sec_pts <- unique(sec_pts, by = c("name","chrom","pos","sec_bin"))
  sec_pts <- map_to_genome_x(sec_pts, gmap, pos_col = "pos")
  sec_pts[gmap, on="chrom", y := i.y_row]

  # --- primary-only points (reads with no secondary), using bin midpoint for stability
  prim_only <- p[in_primary == TRUE & read_has_secondary == FALSE,
                 .(chrom, pos = as.integer((start + end) %/% 2L))]
  prim_only <- unique(prim_only)
  prim_only <- map_to_genome_x(prim_only, gmap, pos_col = "pos")
  prim_only[gmap, on="chrom", y := i.y_row]

  # --- markers for secondary (one per sec point)
  sec_mark <- copy(sec_pts)[, .(chrom, x_global, y)]
  sec_mark <- unique(sec_mark)

  # --- arrows: join per read primary -> all secondary bins for that read
  arrows <- if (nrow(sec_pts) && nrow(prim_mid)) {
    merge(
      prim_mid[, .(name, x0 = x_global, y0 = y)],
      sec_pts[,  .(name, x1 = x_global, y1 = y, sec_bin)],
      by = "name", allow.cartesian = TRUE
    )
  } else data.table()

  # adjust small vertical offsets so glyphs/arrows are visible against the backbone
  y_offset_up   <-  0.14
  y_offset_down <- -0.14
  # primary-only ▼ on the backbone row (slightly above)
  prim_only[, y_pt := y + y_offset_up]
  # secondary ▲ slightly below
  sec_mark[,  y_pt := y + y_offset_down]
  # curve ends slightly above/below their rows
  if (nrow(arrows)) {
    arrows[, y0 := y0 + 0.10]
    arrows[, y1 := y1 - 0.10]
  }

  a_up <- grid::arrow(type = "closed", length = grid::unit(arrow_len_pt, "pt"))

  # x-axis breaks at chromosome midpoints with labels = chrom
  breaks <- gmap$xmid
  labs   <- gmap$chrom

  ggplot() +
    # chromosome backbones (black outline then grey fill for rounded look)
    geom_segment(
      data = backbone,
      aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = backbone_outline, lineend = "round", color = "black"
    ) +
    geom_segment(
      data = backbone,
      aes(x = x0, xend = x1, y = y, yend = y),
      linewidth = backbone_fill, lineend = "round", color = "grey70"
    ) +
    # primary-only markers (black ▼)
    geom_point(
      data = prim_only,
      aes(x = x_global, y = y_pt),
      shape = 25, size = 2.4, stroke = 0, fill = "black"
    ) +
    # secondary markers (red ▲)
    geom_point(
      data = sec_mark,
      aes(x = x_global, y = y_pt),
      shape = 24, size = 2.4, stroke = 0, fill = "#D62728"
    ) +
    # curved arrows (primary -> secondary across chromosomes)
    { if (nrow(arrows)) geom_curve(
      data = arrows,
      aes(x = x0, y = y0, xend = x1, yend = y1, group = interaction(name, sec_bin)),
      curvature = arrow_curvature, color = "#D62728",
      linewidth = 0.4, lineend = "round", arrow = a_up
    ) } +
    scale_x_continuous(
      breaks = breaks, labels = labs, expand = c(0.001, 0.001)
    ) +
    scale_y_continuous(NULL, breaks = backbone$y, labels = rep("", length(backbone$y))) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = base_size, margin = margin(t = 4)),
      plot.margin = margin(t = 6, r = 8, b = 6, l = 8)
    )
}



plot_genome_arrows_flat <- function(genome_dt,
                                    prepped_dt,
                                    chrom_order      = NULL,
                                    gap_bp           = 5e6,   # << gap between chromosomes
                                    sep_lines        = TRUE,  # vertical separators in the gaps
                                    base_size        = 9,
                                    backbone_outline = 5,
                                    backbone_fill    = 3,
                                    marker_offset    = 0.08,
                                    arc_height       = 0.8,
                                    arrow_len_pt     = 3,
                                    curvature        = 0.25) {

  stopifnot(all(c("chrom","start","end") %in% names(genome_dt)))
  need <- c("name","chrom","start","end","in_primary",
            "read_has_secondary","region_midpoint")
  stopifnot(all(need %in% names(prepped_dt)))

  library(data.table); library(ggplot2)

  # ----- Linear genome with gaps -----
  g <- as.data.table(genome_dt)[, .(chrom, start = as.numeric(start), end = as.numeric(end))]
  g[, length := end - start]

  if (is.null(chrom_order)) chrom_order <- g$chrom
  chrom_order <- unique(chrom_order[chrom_order %chin% g$chrom])
  g <- g[match(chrom_order, chrom)]

  # cumulative positions with gaps
  gap_bp <- as.numeric(gap_bp)
  n <- nrow(g)
  g[, cum_start := 0]
  if (n > 1) {
    for (i in 2:n) {
      g$cum_start[i] <- g$cum_start[i-1] + g$length[i-1] + gap_bp
    }
  }
  g[, cum_end := cum_start + length]
  g[, xmid   := (cum_start + cum_end)/2]

  # backbones
  backbone <- g[, .(chrom, x0 = cum_start, x1 = cum_end, y = 0)]

  # optional separators (center of each gap)
  seps <- NULL
  if (sep_lines && n > 1) {
    seps <- data.table(x = g$cum_end[-n] + gap_bp/2)
  }

  # helper: map chromosome bp -> global x with gaps
  map_to_global <- function(dt, pos_col) {
    x <- as.data.table(dt)
    x[g[, .(chrom, cum_start)], on = "chrom",
      x_global := cum_start + get(pos_col)]
    x[]
  }

  # ----- Primary / Secondary coordinates -----
  p <- as.data.table(prepped_dt)

  prim_mid <- p[in_primary == TRUE,
                .(chrom = chrom[1L],
                  pos   = as.integer(round(mean((start + end) %/% 2L)))),
                by = name]
  prim_mid <- map_to_global(prim_mid, "pos")

  sec_pts <- p[!is.na(region_midpoint),
               .(name, chrom, pos = as.integer(region_midpoint))]
  sec_pts <- unique(sec_pts, by = c("name","chrom","pos"))
  sec_pts <- map_to_global(sec_pts, "pos")

  prim_only <- p[in_primary == TRUE & read_has_secondary == FALSE,
                 .(chrom, pos = as.integer((start + end) %/% 2L))]
  prim_only <- unique(prim_only)
  prim_only <- map_to_global(prim_only, "pos")
  prim_only[, y := 0 + marker_offset]

  sec_unique <- unique(sec_pts[, .(x_global)])
  sec_unique[, y := 0 + marker_offset]

  arrows <- if (nrow(sec_pts) && nrow(prim_mid)) {
    merge(prim_mid[, .(name, x0 = x_global)],
          sec_pts[,  .(name, x1 = x_global)],
          by = "name", allow.cartesian = TRUE)
  } else data.table()
  if (nrow(arrows)) {
    arrows[, `:=`(y0 = 0 + 0.02, y1 = 0 + 0.01)]
  }

  backbone   <- backbone[is.finite(x0) & is.finite(x1)]
  prim_only  <- prim_only[is.finite(x_global)]
  sec_unique <- sec_unique[is.finite(x_global)]
  arrows     <- arrows[is.finite(x0) & is.finite(x1)]

  x_breaks <- g$xmid
  x_labs   <- g$chrom
  y_lim    <- c(-marker_offset*2, arc_height)

  ggplot() +
    # optional gap separators
    { if (!is.null(seps)) geom_segment(
      data = seps, aes(x = x, xend = x, y = -marker_offset*1.5, yend = arc_height*0.25),
      linewidth = 0.2, color = "grey80", inherit.aes = FALSE
    ) } +
    # chromosome backbones
    geom_segment(data = backbone,
                 aes(x = x0, xend = x1, y = 0, yend = 0),
                 linewidth = backbone_outline, lineend = "round", color = "black") +
    geom_segment(data = backbone,
                 aes(x = x0, xend = x1, y = 0, yend = 0),
                 linewidth = backbone_fill, lineend = "round", color = "grey70") +
    # markers
    geom_point(data = prim_only, aes(x = x_global, y = y),
               shape = 25, size = 2.4, stroke = 0, fill = "black") +
    geom_point(data = sec_unique, aes(x = x_global, y = y),
               shape = 24, size = 2.4, stroke = 0, fill = "#D62728") +
    # arrows
    { if (nrow(arrows)) geom_curve(
      data = arrows,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature = curvature, color = "#C43B3B", linewidth = 0.5, lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(arrow_len_pt, "pt"))
    ) } +
    scale_x_continuous(breaks = x_breaks, labels = x_labs, expand = c(0.002, 0.002)) +
    scale_y_continuous(NULL, limits = y_lim, breaks = NULL, labels = NULL, expand = c(0.1, 0.1)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = base_size) +
    theme(
      axis.text.x = element_text(margin = ggplot2::margin(t = 4)),
      plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
    )
}


plot_genome_arrows_flat <- function(
    genome_dt,
    prepped_dt,
    gap_bp = 1e6,              # visual gap inserted between chromosomes
    base_size = 8,
    chrom_outline_size = 5,
    chrom_fill_size    = 3,
    arrow_color = "#C43B3B",
    curvature_abs = 0.25,      # 0..1; larger = more arch
    lift = 0.06                 # vertical lift so arcs sit above backbone
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Please install 'ggplot2'.")
  library(data.table); library(ggplot2)

  g0 <- as.data.table(genome_dt)
  p0 <- as.data.table(prepped_dt)

  # --- checks ---
  stopifnot(all(c("chrom","start","end") %in% names(g0)))
  need_cols <- c("chrom","start","end","in_primary","read_has_secondary",
                 "bin_start","bin_end","name","region_midpoint","sec_bin")
  if (!all(need_cols %in% names(p0))) {
    stop("prepped_dt must include: ", paste(need_cols, collapse = ", "),
         ". (Use prep_primary_secondary_bins().)")
  }

  # ---- place chromosomes on one long x-axis with gaps ----
  g <- copy(g0)[order(chrom)]
  g[, chr_len := as.numeric(end - start)]
  g[, offset  := c(0, cumsum(head(chr_len + gap_bp, -1)))]
  g[, `:=`(x0 = offset + start, x1 = offset + end, y = 0)]

  # map function: (chrom, pos) -> global x
  offmap <- setNames(g$offset, g$chrom)
  to_global <- function(chrom, pos) as.numeric(offmap[chrom]) + as.numeric(pos)

  # --- markers ---
  # primary-only markers (down triangles)
  prim_pts <- unique(
    p0[in_primary == TRUE & read_has_secondary == FALSE,
       .(chrom, x = as.integer((bin_start + bin_end) %/% 2L))]
  )
  if (nrow(prim_pts)) prim_pts[, gx := to_global(chrom, x)]

  # secondary markers (up triangles) at region_midpoint
  sec_pts <- unique(
    p0[!in_primary & !is.na(region_midpoint),
       .(chrom, x = as.integer(region_midpoint))]
  )
  if (nrow(sec_pts)) sec_pts[, gx := to_global(chrom, x)]

  # --- arrows (primary -> chosen secondary for each read with a sec_bin) ---
  # primary midpoint per read (on the target chromosome)
  prim_mid <- p0[in_primary == TRUE,
                 .(chrom_primary = chrom[1L],
                   x_primary     = as.integer(median((bin_start + bin_end) %/% 2L))),
                 by = name]

  # secondary midpoint per read (use region_midpoint for the chosen sec_bin)
  sec_rows <- p0[!in_primary & !is.na(sec_bin) & !is.na(region_midpoint)]
  sec_mid  <- sec_rows[, .(
    chrom_secondary = chrom[1L],
    x_secondary     = as.integer(region_midpoint[1L])
  ), by = name]

  arrows <- merge(prim_mid, sec_mid, by = "name", all = FALSE)
  if (nrow(arrows)) {
    arrows[, `:=`(
      x0 = to_global(chrom_primary,  x_primary),
      x1 = to_global(chrom_secondary, x_secondary),
      y0 = lift,
      y1 = lift
    )]
  }

  # split by direction so both sets bend upward
  ar_lr <- arrows[x1 >= x0]
  ar_rl <- arrows[x1 <  x0]

  # chromosome label positions
  chr_labs <- g[, .(chrom, gx = (x0 + x1) / 2)]

  # ---- plot ----
  p <- ggplot() +
    # chromosome backbone: outline then fill to get rounded track look
    geom_segment(data = g,
                 aes(x = x0, xend = x1, y = y, yend = y),
                 linewidth = chrom_outline_size, lineend = "round", color = "black") +
    geom_segment(data = g,
                 aes(x = x0, xend = x1, y = y, yend = y),
                 linewidth = chrom_fill_size, lineend = "round", color = "grey70") +

    # primary-only markers: black down triangles
    { if (nrow(prim_pts)) geom_point(
      data = prim_pts, aes(x = gx, y = 0),
      shape = 25, size = 2.4, stroke = 0, fill = "black"
    ) } +

    # secondary markers: red up triangles
    { if (nrow(sec_pts)) geom_point(
      data = sec_pts, aes(x = gx, y = 0),
      shape = 24, size = 2.4, stroke = 0, fill = "#C43B3B"
    ) } +

    # upward-curving arrows (two layers with opposite curvature signs)
    { if (nrow(ar_lr)) geom_curve(
      data = ar_lr,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature =  -curvature_abs,
      color = "#C43B3B", linewidth = 0.5, lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(3, "pt"))
    ) } +
    { if (nrow(ar_rl)) geom_curve(
      data = ar_rl,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature = curvature_abs,
      color = "#C43B3B", linewidth = 0.5, lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(3, "pt"))
    ) } +

    # chromosome labels (x axis look)
    scale_x_continuous(
      breaks = chr_labs$gx, labels = chr_labs$chrom,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(ylim = c(-0.12, 0.5), clip = "off") +
    theme_void(base_size = base_size) +
    theme(
      axis.text.x  = element_text(),
      axis.ticks.x = element_blank()
    )

  return(p)
}


plot_genome_arrows_flat <- function(
    genome_dt,
    prepped_dt,
    gap_bp = 1e6,              # visual gap inserted between chromosomes
    base_size = 8,
    chrom_outline_size = 5,
    chrom_fill_size    = 3,
    arrow_color = "#C43B3B",
    curvature_abs = 0.25,      # 0..1; larger = more arch
    lift = 0.06                 # vertical lift so arcs sit above backbone
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Please install 'ggplot2'.")
  library(data.table); library(ggplot2)

  g0 <- as.data.table(genome_dt)
  p0 <- as.data.table(prepped_dt)

  # --- checks ---
  stopifnot(all(c("chrom","start","end") %in% names(g0)))
  need_cols <- c("chrom","start","end","in_primary","read_has_secondary",
                 "bin_start","bin_end","name","region_midpoint","sec_bin")
  if (!all(need_cols %in% names(p0))) {
    stop("prepped_dt must include: ", paste(need_cols, collapse = ", "),
         ". (Use prep_primary_secondary_bins().)")
  }

  # ---- place chromosomes on one long x-axis with gaps ----
  g <- copy(g0)[order(chrom)]
  g[, chr_len := as.numeric(end - start)]
  g[, offset  := c(0, cumsum(head(chr_len + gap_bp, -1)))]
  g[, `:=`(x0 = offset + start, x1 = offset + end, y = 0)]

  # map function: (chrom, pos) -> global x
  offmap <- setNames(g$offset, g$chrom)
  to_global <- function(chrom, pos) as.numeric(offmap[chrom]) + as.numeric(pos)

  # --- markers ---
  # primary-only markers (down triangles)
  prim_pts <- unique(
    p0[in_primary == TRUE & read_has_secondary == FALSE,
       .(chrom, x = as.integer((bin_start + bin_end) %/% 2L))]
  )
  if (nrow(prim_pts)) prim_pts[, gx := to_global(chrom, x)]

  # secondary markers (up triangles) at region_midpoint
  sec_pts <- unique(
    p0[!in_primary & !is.na(region_midpoint),
       .(chrom, x = as.integer(region_midpoint))]
  )
  if (nrow(sec_pts)) sec_pts[, gx := to_global(chrom, x)]

  # --- arrows (primary -> chosen secondary for each read with a sec_bin) ---
  # primary midpoint per read (on the target chromosome)
  prim_mid <- p0[in_primary == TRUE,
                 .(chrom_primary = chrom[1L],
                   x_primary     = as.integer(median((bin_start + bin_end) %/% 2L))),
                 by = name]

  # secondary midpoint per read (use region_midpoint for the chosen sec_bin)
  sec_rows <- p0[!in_primary & !is.na(sec_bin) & !is.na(region_midpoint)]
  sec_mid  <- sec_rows[, .(
    chrom_secondary = chrom[1L],
    x_secondary     = as.integer(region_midpoint[1L])
  ), by = name]

  arrows <- merge(prim_mid, sec_mid, by = "name", all = FALSE)
  if (nrow(arrows)) {
    arrows[, `:=`(
      x0 = to_global(chrom_primary,  x_primary),
      x1 = to_global(chrom_secondary, x_secondary),
      y0 = 0,
      y1 = lift
    )]
  }

  # split by direction so both sets bend upward
  ar_lr <- arrows[x1 >= x0]
  ar_rl <- arrows[x1 <  x0]

  # chromosome label positions
  chr_labs <- g[, .(chrom, gx = (x0 + x1) / 2)]

  # ---- plot ----
  p <- ggplot() +
    # chromosome backbone: outline then fill to get rounded track look
    geom_segment(data = g,
                 aes(x = x0, xend = x1, y = y, yend = y),
                 linewidth = chrom_outline_size, lineend = "round", color = "black") +
    geom_segment(data = g,
                 aes(x = x0, xend = x1, y = y, yend = y),
                 linewidth = chrom_fill_size, lineend = "round", color = "grey70") +

    # primary-only markers: black down triangles
    { if (nrow(prim_pts)) geom_point(
      data = prim_pts, aes(x = gx, y = 0),
      shape = 25, size = 2.4, stroke = 0, fill = "black"
    ) } +

    # secondary markers: red up triangles
    { if (nrow(sec_pts)) geom_point(
      data = sec_pts, aes(x = gx, y = 0),
      shape = 24, size = 2.4, stroke = 0, fill = "#C43B3B"
    ) } +

    # upward-curving arrows (two layers with opposite curvature signs)
    { if (nrow(ar_lr)) geom_curve(
      data = ar_lr,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature =  -curvature_abs,
      color = "#C43B3B", linewidth = 0.25, lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(0, "pt"))
    ) } +
    { if (nrow(ar_rl)) geom_curve(
      data = ar_rl,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      curvature = +curvature_abs,
      color = "#C43B3B", linewidth = 0.25, lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(0, "pt"))
    ) } +

    # chromosome labels (x axis look)
    scale_x_continuous(
      breaks = chr_labs$gx, labels = chr_labs$chrom,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(ylim = c(-0.01, 0.3), clip = "off") +
    theme_void(base_size = base_size) +
    theme(
      axis.text.x  = element_text(margin = unit(c(0, 0, 0, 0), "pt")),
      axis.ticks.x = element_blank()
    )

  return(p)
}

library(data.table)
library(ggplot2)
library(grid)  # for unit()

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

  # ---- 1) Define facets: primary node + its partners
  partners <- unique(c(
    edges[cluster1 == node, cluster2],
    edges[cluster2 == node, cluster1]
  ))
  facets <- unique(c(node, partners))
  if (!length(facets)) stop(sprintf("Node '%s' has no partners in net$edges.", node))

  # ---- 2) Keep only reads that are in the PRIMARY node
  if (!all(c(cluster_col, read_id_col, start_col, end_col) %in% names(reads))) {
    stop("Expected columns missing in net$reads.")
  }
  primary_reads <- unique(reads[get(cluster_col) == node, get(read_id_col)])
  if (!length(primary_reads)) {
    stop(sprintf("No reads found in primary node '%s'.", node))
  }

  # All segments for those primary reads, but only in (primary + partner) facets
  sub <- reads[
    get(read_id_col) %chin% primary_reads & get(cluster_col) %chin% facets
  ][
    , `:=`(
      cl  = as.character(get(cluster_col)),
      rid = as.character(get(read_id_col))
    )
  ]

  # Label primary vs partner for color
  sub[, cluster_type := data.table::fifelse(cl == node, "primary", "partner")]

  # Facet order: primary first, then partners alphabetically
  partner_levels <- sort(setdiff(facets, node))
  facet_levels <- c(node, partner_levels)
  sub[, (cluster_col) := factor(cl, levels = facet_levels)]

  # ---- 3) Cluster sizes (for ranking the "best" secondary cluster)
  if (!is.null(nodes) && all(c("cluster","cluster_n") %in% names(nodes))) {
    cl_size <- nodes[cluster %chin% facet_levels, .(cluster, cluster_n)]
  } else {
    # Fallback: unique reads per cluster within 'sub'
    cl_size <- sub[, .(cluster_n = data.table::uniqueN(rid)), by = .(cluster = cl)]
  }
  data.table::setkey(cl_size, cluster)

  # ---- 4) Per-(cluster, read) min start (so we can compare within clusters)
  per_read_cl <- sub[, .(
    read_min_start = suppressWarnings(min(get(start_col), na.rm = TRUE))
  ), by = .(cl, rid)]
  data.table::setkey(per_read_cl, rid, cl)

  # ---- 5) For each read, choose its "best" secondary cluster
  # Best secondary = among partner facets that the read hits: largest cluster_n (ties: earlier start, then name)
  read_clusters <- per_read_cl[, .(clusters = list(cl)), by = rid]
  data.table::setkey(read_clusters, rid)

  choose_best_secondary <- function(rid_val, current_cl) {
    cls <- read_clusters[.(rid_val), clusters][[1]]
    cls <- setdiff(cls, current_cl)                 # exclude the current facet
    cls <- setdiff(cls, node)                       # exclude primary; "secondary" means partner
    if (length(cls) == 0L) return(list(NA_character_, 0L, Inf))

    cand <- data.table::data.table(cluster = cls)
    cand <- cl_size[cand, on = .(cluster)]
    cand[is.na(cluster_n), cluster_n := 0L]

    # Tie-break by earlier start for THIS read in that cluster
    # (if missing, treat as +Inf)
    st <- per_read_cl[.(rid_val, cand$cluster), read_min_start]
    if (length(st) != nrow(cand)) {
      # align safely
      st <- per_read_cl[.(rid_val, cand$cluster), on = .(rid, cl), x.read_min_start]
    }
    st[is.na(st)] <- Inf

    cand[, read_min_start := st]
    data.table::setorder(cand, -cluster_n, read_min_start, cluster)
    list(cand$cluster[1], as.integer(cand$cluster_n[1]), cand$read_min_start[1])
  }

  # Build one row per READ (global ordering keys)
  # Also keep primary-start for tie-breaking
  primary_start <- per_read_cl[.(rid = primary_reads, cl = node), on = .(rid, cl)]
  primary_start <- primary_start[, .(rid, primary_start = read_min_start)]
  data.table::setkey(primary_start, rid)

  order_reads <- data.table::data.table(rid = primary_reads)
  tmp <- order_reads[, c("sec_cl","sec_n","sec_start") := choose_best_secondary(rid, node), by = rid]
  order_reads <- tmp
  order_reads <- primary_start[order_reads, on = .(rid)]
  order_reads[is.na(primary_start), primary_start := Inf]

  # Global ordering of READS: by sec_n desc → sec_start asc → primary_start asc → rid
  data.table::setorder(order_reads, -sec_n, sec_start, primary_start, rid)

  # Assign a single y per read (top row = smallest index)
  order_reads[, y := (.N + 1L) - seq_len(.N)]

  # Join that y back to EVERY segment for that read across all facets
  data.table::setkey(order_reads, rid)
  data.table::setkey(sub, rid)
  sub <- order_reads[sub]

  # ---- 6) Alpha mapping per segment
  if (!alpha_col %in% names(sub)) {
    stop(sprintf("alpha_col '%s' not found in net$reads.", alpha_col))
  }
  av <- as.numeric(sub[[alpha_col]])
  rng <- range(av, na.rm = TRUE)
  a0 <- if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
    (av - rng[1]) / (rng[2] - rng[1])
  } else {
    rep(1, length(av))
  }
  a0[!is.finite(a0)] <- 0.5
  sub[, alpha_val := alpha_range[1] + a0 * (alpha_range[2] - alpha_range[1])]

  # ---- 7) Plot ---------------------------------------------------------------
  breaks_1kb <- function(lims) {
    from <- floor(lims[1] / 1000) * 1000
    to   <- ceiling(lims[2] / 1000) * 1000
    seq(from, to, by = 1000L)
  }

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
      scales = "free_x",
      space  = "free_x",
      labeller = ggplot2::label_wrap_gen(width = strip_wrap_width)
    ) +
    ggplot2::scale_color_manual(values = c(primary = primary_color, partner = partner_color)) +
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::scale_x_continuous(
      position = "top",
      breaks = function(lims) {
        from <- floor(lims[1] / 1000) * 1000
        to   <- ceiling(lims[2] / 1000) * 1000
        seq(from, to, by = 1000L)
      },
      labels = NULL,
      minor_breaks = NULL,
      expand = c(0.02, 0.02)  # small padding to avoid panel/strip crowding
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.placement     = "outside",
      strip.background    = ggplot2::element_blank(),
      strip.text          = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                  margin = ggplot2::margin(t = 4, b = 4)),
      strip.clip          = "off",                # <- prevent clipping
      panel.grid          = ggplot2::element_blank(),
      panel.spacing.x     = grid::unit(6, "pt"),  # <- a bit more room between facets
      axis.title          = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.title.x.bottom = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")        # <- allow strip to render fully

  if (isTRUE(ensure_strip_space)) {
    # If ggh4x is installed, enforce a minimum width for each facet
    # if (requireNamespace("ggh4x", quietly = TRUE)) {
    #   ncols <- nlevels(sub[[cluster_col]])
    #   p <- force_panelsizes(
    #     p,
    #     cols_widths = grid::unit(rep(min_facet_width_mm, ncols), "mm"),
    #     respect = TRUE
    #   )
    # }
    # If ggh4x isn’t available, the theme tweaks above still help a lot.
  }

  p
}



plot_reads_around_node <- function(net, node, bed,
                                   window_kb   = 20,
                                   primary_color = "red2",
                                   other_color   = "grey60",
                                   alpha_col     = NULL,        # e.g., "score"
                                   alpha_range   = c(0.3, 1.0),
                                   base_size     = 6) {
  stopifnot(is.list(net), "reads" %in% names(net))
  reads_dt <- as.data.table(net$reads)
  bed_dt   <- as.data.table(bed)

  # ---- Required columns ----
  need_reads_cols <- c("chrom_cluster_id","chrom","start","end","name")
  miss_reads <- setdiff(need_reads_cols, names(reads_dt))
  if (length(miss_reads)) stop("net$reads missing columns: ", paste(miss_reads, collapse=", "))

  need_bed_cols <- c("chrom","start","end","name")
  miss_bed <- setdiff(need_bed_cols, names(bed_dt))
  if (length(miss_bed)) stop("bed is missing columns: ", paste(miss_bed, collapse=", "))

  # ---- Node span (from net$nodes if present, else from net$reads) ----
  chrom_node <- NULL; node_start <- NULL; node_end <- NULL
  if ("nodes" %in% names(net) && all(c("cluster","min_start","max_end","chrom") %in% names(net$nodes))) {
    node_row <- as.data.table(net$nodes)[cluster == node]
    if (nrow(node_row)) {
      chrom_node <- node_row$chrom[1]
      node_start <- as.numeric(node_row$min_start[1])
      node_end   <- as.numeric(node_row$max_end[1])
    }
  }
  if (is.null(chrom_node) || any(!is.finite(c(node_start, node_end)))) {
    node_reads <- reads_dt[chrom_cluster_id == node]
    if (!nrow(node_reads)) stop(sprintf("No reads found for node '%s' in net$reads.", node))
    chrom_node <- as.character(node_reads$chrom[1])
    node_start <- min(as.numeric(node_reads$start), na.rm = TRUE)
    node_end   <- max(as.numeric(node_reads$end),   na.rm = TRUE)
  }

  pad <- as.integer(round(window_kb * 1000))
  win_start <- max(0L, as.integer(floor(node_start - pad)))
  win_end   <- as.integer(ceiling(node_end + pad))

  # Node read names (from net$reads)
  node_read_names <- unique(reads_dt[chrom_cluster_id == node]$name)

  # Reads in window (from BED)
  window_reads <- bed_dt[chrom == chrom_node & start <= win_end & end >= win_start]
  if (!nrow(window_reads)) stop("No reads found in the specified window.")

  # Label node vs other
  window_reads[, cluster_type := ifelse(nameN>1, "node", "other")]

  # ---- Order reads: smallest start at TOP ----
  window_reads[, read_min_start := suppressWarnings(min(start, na.rm = TRUE)), by = name]
  ord <- unique(window_reads[, .(name, read_min_start)])
  setorder(ord, read_min_start, name)
  ord[, y := (.N + 1L) - seq_len(.N)]
  setkey(ord, name); setkey(window_reads, name)
  window_reads <- ord[window_reads]  # adds y by name

  # ---- Optional alpha mapping ----
  if (!is.null(alpha_col)) {
    if (!alpha_col %in% names(window_reads)) {
      stop(sprintf("alpha_col '%s' not found in BED/window data.", alpha_col))
    }
    av  <- as.numeric(window_reads[[alpha_col]])
    r   <- range(av, na.rm = TRUE)
    a01 <- if (is.finite(r[1]) && is.finite(r[2]) && diff(r) > 0) {
      (av - r[1]) / (r[2] - r[1])
    } else rep(1, length(av))
    a01[!is.finite(a01)] <- 0.5
    window_reads[, alpha_val := alpha_range[1] + a01 * (alpha_range[2] - alpha_range[1])]
  } else {
    window_reads[, alpha_val := 1]
  }

  # ---- Plot ----
  breaks_1kb <- function(lims) {
    from <- floor(lims[1] / 1000) * 1000
    to   <- ceiling(lims[2] / 1000) * 1000
    seq(from, to, by = 1000L)
  }

  ggplot(window_reads) +
    geom_segment(aes(y = y, yend = y,
                     x = start, xend = end,
                     color = cluster_type, alpha = alpha_val),
                 linewidth = 0.5, lineend = "round") +
    scale_color_manual(values = c(node = primary_color, other = other_color)) +
    scale_alpha(range = alpha_range, guide = "none") +
    scale_x_continuous(
      position = "top",
      breaks = breaks_1kb,
      labels = NULL,
      minor_breaks = NULL,
      expand = c(0, 0)
    ) +
    guides(color = "none") +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid          = element_blank(),
      axis.title          = element_blank(),
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      legend.position     = "none",
      axis.text.x.bottom  = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      strip.placement     = "outside",
      strip.background    = element_blank(),
      axis.ticks.length   = grid::unit(-1.5, "mm")  # ticks into panel
    )
}


library(data.table)
library(dbscan)

# --- split clusters by a maximum allowed consecutive gap (in bp) --------------
# Adds:
#  - verbose messages for each split (gap > max_gap)
#  - an attribute "split_notes" on the returned data.table (data.table of splits)
#  - optional return of notes via return_notes = TRUE
.split_by_max_gap <- function(dt,
                              id_col   = "cluster_id",
                              start    = "start",
                              end      = "end",
                              max_gap  = Inf,
                              chrom_col = NULL,       # optional: include chrom in notes if available
                              verbose  = TRUE,
                              return_notes = FALSE) {
  if (!is.finite(max_gap)) {
    out <- data.table::copy(dt)
    attr(out, "split_notes") <- data.table::data.table()  # empty notes
    return(out)
  }
  x <- data.table::copy(dt)
  n <- nrow(x)
  new_ids <- integer(n)          # keep 0 for noise
  next_id <- 1L

  s   <- suppressWarnings(as.numeric(x[[start]]))
  e   <- suppressWarnings(as.numeric(x[[end]]))
  cid <- x[[id_col]]

  # collect split notes
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

      # need new subcluster?
      if (is.na(cur_global) || si > prev_end + max_gap) {
        # if not the very first interval in this original cluster, record a split
        if (!is.na(cur_global)) {
          gap_bp <- as.numeric(si - prev_end)
          note <- data.table::data.table(
            original_cluster = g,
            new_subcluster   = next_id,                # about to assign
            previous_subcluster = prev_subid,
            split_at_start   = si,
            previous_end     = prev_end,
            gap_bp           = gap_bp,
            threshold_bp     = max_gap,
            chrom            = if (!is.null(chrom_col) && chrom_col %in% names(x))
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
        # start a new subcluster
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

  # attach notes
  notes_dt <- if (length(notes)) data.table::rbindlist(notes) else data.table::data.table()
  attr(out, "split_notes") <- notes_dt

  if (isTRUE(return_notes)) {
    return(list(dt = out, split_notes = notes_dt))
  }
  out
}


# --- add cluster_n and optional "too-close" warnings (per chrom) ---------------
.finalize_clusters <- function(dt,
                               start = "start",
                               end = "end",
                               chrom_col = "chrom",
                               warn_if_clusters_closer_than = NA_real_,
                               verbose = TRUE) {
  stopifnot(data.table::is.data.table(dt))
  x <- copy(dt)

  # cluster_n (size per cluster, per chromosome); noise 0 gets 0
  sizes <- x[cluster_id > 0L, .N, by = c(chrom_col, "cluster_id")]
  x[, cluster_n := 0L]
  if (nrow(sizes)) {
    x[sizes, on = c(chrom_col, "cluster_id"), cluster_n := i.N]
  }

  # Optional warnings: cluster centers per chromosome
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
          warn_list[[as.character(cc)]] <- data.table(
            chrom     = cc,
            cluster_a = centers$cluster_id[hits],
            cluster_b = centers$cluster_id[hits + 1L],
            dist_bp   = as.numeric(d[hits])
          )
        }
      }
    }
  }

  warn_tbl <- if (length(warn_list)) rbindlist(warn_list) else data.table()
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

# --- apply a function per chromosome, combine, and scope IDs -------------------
# robust per-chrom apply (no `get()`)
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
      res[, cluster_id := fifelse(cluster_id > 0L, cluster_id + (next_id - 1L), 0L)]
      next_id <- next_id + max_here
    }
    out_list[[i]] <- res
  }
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

# label clusters as "chrom_pos" where pos = center of cluster (midpoint median/mean)
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


cluster_reads_sweep <- function(dt, chrom_col = "chrom",
                                start = "start", end = "end",
                                delta = 0L,
                                max_gap = Inf,  # ← what you meant by "span"
                                max_span = NULL,  # deprecated alias → treated as max_gap
                                warn_if_clusters_closer_than = NA_real_,
                                id_scope = c("global"),
                                verbose = TRUE) {
  if (!is.null(max_span)) { max_gap <- max_span; if (verbose) message("Using max_span as max_gap (consecutive-gap rule).") }
  id_scope <- match.arg(id_scope)
  stopifnot(data.table::is.data.table(dt))

  per_chrom <- function(sub) {
    y <- copy(sub)[, .idx := .I]

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
      # ONLY the overlap / ≤ delta rule; no total-width cap here
      if (si <= cur_end + delta) {
        cur_end <- max(cur_end, ei, na.rm = TRUE)
      } else {
        cur <- cur + 1L
        cur_end <- ei
      }
      cl[i] <- cur
    }

    y[o, cluster_id := cl]
    setorder(y, .idx); y[, .idx := NULL]

    # split by consecutive gap > max_gap
    y <- .split_by_max_gap(y, "cluster_id", start, end, max_gap)

    # finalize per chrom
    .finalize_clusters(y, start, end, chrom_col,
                       warn_if_clusters_closer_than, verbose)
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  res <- .add_chrom_cluster_id(res,
                               chrom_col = chrom_col,
                               start = start, end = end,
                               id_col = "cluster_id",
                               label_col = "chrom_cluster_id",
                               zero_label = "chrom_0",     # or "NA"
                               pos_stat = "median",
                               pos_round = TRUE)
  res
}

cluster_reads_dbscan <- function(dt, chrom_col = "chrom",
                                 start = "start", end = "end",
                                 distance = c("gap","euclidean","linf"),
                                 eps = NULL, minPts = 3L,
                                 max_gap = Inf, max_span = NULL,
                                 warn_if_clusters_closer_than = NA_real_,
                                 id_scope = c("per_chrom","global"),
                                 verbose = TRUE) {
  if (!is.null(max_span)) { max_gap <- max_span; if (verbose) message("Using max_span as max_gap (consecutive-gap rule).") }
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
      fit <- dbscan(as.dist(D), eps = local_eps, minPts = minPts, borderPoints = TRUE)
    } else if (distance == "euclidean") {
      X <- cbind((s + e) / 2, e - s)
      fit <- dbscan(X, eps = local_eps, minPts = minPts)
    } else { # linf
      D <- matrix(0, n, n)
      for (i in seq_len(n)) {
        di <- pmax(abs(s[i] - s), abs(e[i] - e))
        D[i, ] <- di; D[, i] <- di
      }
      fit <- dbscan(as.dist(D), eps = local_eps, minPts = minPts, borderPoints = TRUE)
    }

    out <- copy(sub)
    out[, cluster_id := as.integer(fit$cluster)]
    out <- .split_by_max_gap(out, "cluster_id", start, end, max_gap)
    .finalize_clusters(out, start, end, chrom_col,
                       warn_if_clusters_closer_than, verbose)
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  res <- .add_chrom_cluster_id(res,
                               chrom_col = chrom_col,
                               start = start, end = end,
                               id_col = "cluster_id",
                               label_col = "chrom_cluster_id",
                               zero_label = "chrom_0",     # or "NA"
                               pos_stat = "median",
                               pos_round = TRUE)
  res
}

cluster_reads_kde <- function(dt, chrom_col = "chrom",
                              start = "start", end = "end",
                              bw = NULL, bw_method = c("SJ","nrd0","ucv","bcv"),
                              bw_adjust = 1, kernel = "gaussian",
                              max_gap = Inf, max_span = NULL,
                              warn_if_clusters_closer_than = NA_real_,
                              id_scope = c("per_chrom","global"),
                              verbose = TRUE,
                              unassigned_label = 0L) {
  if (!is.null(max_span)) { max_gap <- max_span; if (verbose) message("Using max_span as max_gap (consecutive-gap rule).") }
  bw_method <- match.arg(bw_method); id_scope <- match.arg(id_scope)
  stopifnot(data.table::is.data.table(dt))

  per_chrom <- function(sub) {
    s <- as.numeric(sub[[start]]); e <- as.numeric(sub[[end]])
    if (!length(s)) { z <- copy(sub); z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]; return(z) }

    # compress coverage into segments
    ev_pos <- c(s, e); ev_val <- c(rep(1L, length(s)), rep(-1L, length(e)))
    o <- order(ev_pos, -ev_val); ev_pos <- ev_pos[o]; ev_val <- ev_val[o]
    pos_u <- unique(ev_pos)
    if (length(pos_u) < 2L) {
      z <- copy(sub); z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]; return(z)
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
      z <- copy(sub); z[, `:=`(cluster_id = unassigned_label, cluster_n = 0L)]
      return(z)
    }

    x <- seg_mid[keep]; w <- w[keep] / sum(w[keep])
    local_bw <- if (is.null(bw)) switch(bw_method,
                                        SJ   = bw.SJ(x, method="ste"),
                                        nrd0 = bw.nrd0(x),
                                        ucv  = bw.ucv(x),
                                        bcv  = bw.bcv(x)) else bw

    kd <- density(x, weights = w, bw = local_bw * bw_adjust, kernel = kernel)
    midx <- which(diff(sign(diff(kd$y))) == -2) + 1
    modes <- kd$x[midx]

    # assign each interval to nearest mode; 0 if farther than max_gap from interval
    lab <- integer(length(s))
    for (i in seq_along(s)) {
      if (!length(modes)) { lab[i] <- unassigned_label; next }
      d <- ifelse(modes >= s[i] & modes <= e[i], 0,
                  pmin(abs(modes - s[i]), abs(modes - e[i])))
      j <- which.min(d)
      lab[i] <- if (length(j)==0 || d[j] > max_gap) unassigned_label else j
    }
    out <- copy(sub)
    out[, cluster_id := as.integer(lab)]

    # split by consecutive gap > max_gap (post)
    out <- .split_by_max_gap(out, "cluster_id", start, end, max_gap)

    out <- .finalize_clusters(out, start, end, chrom_col,
                              warn_if_clusters_closer_than, verbose)
    attr(out, "kde_modes") <- modes
    attr(out, "kde_bw")    <- kd$bw
    out
  }

  res <- .apply_by_chrom(dt, chrom_col, per_chrom, id_scope = id_scope)
  res <- .add_chrom_cluster_id(res,
                               chrom_col = chrom_col,
                               start = start, end = end,
                               id_col = "cluster_id",
                               label_col = "chrom_cluster_id",
                               zero_label = "chrom_0",     # or "NA"
                               pos_stat = "median",
                               pos_round = TRUE)
  res
}


library(data.table)


# Return reads shared between two clusters
# - id_type = "chrom_cluster_id": pass a="Chr1_12345", b="Chr3_67890"
# - id_type = "cluster_id":       pass a=12, b=17; if cluster_id is per-chrom, also set chrom_a / chrom_b
# return = "reads"  -> character vector of read names
#        = "rows"   -> rows from dt for those reads (both clusters)
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
    stop("`chrom_cluster_id` column not found. Build labels before using id_type='chrom_cluster_id'.")
  if (id_type == "cluster_id" && !"cluster_id" %in% names(dt))
    stop("`cluster_id` column not found in dt.")

  # subset rows for cluster A and B
  if (id_type == "chrom_cluster_id") {
    A <- dt[get("chrom_cluster_id") == a]
    B <- dt[get("chrom_cluster_id") == b]
  } else {
    # cluster_id; handle possible per-chrom scopes
    A <- dt[cluster_id == a]
    B <- dt[cluster_id == b]
    # if ambiguous across chroms and user didn't specify, error with hint
    if (!is.null(chrom_a)) A <- A[get(chrom_col) == chrom_a]
    if (!is.null(chrom_b)) B <- B[get(chrom_col) == chrom_b]
    if (!nrow(A)) {
      stop(sprintf("No rows for cluster_id=%s%s.",
                   as.character(a),
                   if (!is.null(chrom_a)) paste0(" (chrom=", chrom_a, ")") else ""))
    }
    if (!nrow(B)) {
      stop(sprintf("No rows for cluster_id=%s%s.",
                   as.character(b),
                   if (!is.null(chrom_b)) paste0(" (chrom=", chrom_b, ")") else ""))
    }
    # if still multiple chroms present and chrom_* not given, warn
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
    # return the rows for those shared reads from BOTH clusters (tag which cluster)
    A2 <- A[get(read_col) %in% shared][, cluster_side := "A"]
    B2 <- B[get(read_col) %in% shared][, cluster_side := "B"]
    return(rbindlist(list(A2, B2), use.names = TRUE, fill = TRUE))
  }
}

# shared <- shared_reads_between(chr_sw,
#                                a = "314",
#                                b = "387",
#                                id_type = "cluster_id",
#                                return = "reads"
# )


# build_cluster_network <- function(dt,
#                                   cluster_col = "chrom_cluster_id",  # or "cluster_id"
#                                   read_col    = "name",
#                                   chrom_col   = "chrom",
#                                   start_col   = "start",
#                                   end_col     = "end",
#                                   id_col      = "cluster_id",
#                                   include_noise = FALSE,
#                                   restrict_within_chrom = FALSE,
#                                   min_shared = 1L) {
#   stopifnot(data.table::is.data.table(dt))
#   x <- copy(dt)
#
#   if (!include_noise && id_col %in% names(x)) x <- x[x[[id_col]] > 0L]
#
#   if (!nrow(x)) {
#     return(list(
#       nodes = data.table(cluster = character(), chrom = character(),
#                          cluster_n = integer(), node_avg_score = numeric(),
#                          node_avg_len = numeric(),
#                          min_start = integer(), max_end = integer(),
#                          edge_count = integer(), total_length = integer()),
#       edges = data.table(cluster1 = character(), cluster2 = character(),
#                          shared_reads = integer(),
#                          shared_avg_score = numeric(), shared_avg_len = numeric()),
#       reads = x[0]
#     ))
#   }
#
#   # Per-read global stats (for edge-level means)
#   read_global <- x[, .(
#     read_avg_score = mean(score, na.rm = TRUE),
#     read_avg_len   = mean(aligned_length, na.rm = TRUE)
#   ), by = read_col]
#   data.table::setnames(read_global, old = read_col, new = "read")
#
#   # Membership table
#   cr <- unique(x[, .(cluster = .SD[[1]], read = .SD[[2]], chrom = .SD[[3]]),
#                  .SDcols = c(cluster_col, read_col, chrom_col)])
#
#   # Per-cluster/per-read detail for node stats
#   cr_detail <- x[, .(
#     read_cluster_score = mean(score, na.rm = TRUE),
#     read_cluster_len   = mean(aligned_length, na.rm = TRUE)
#   ), by = .(cluster = x[[cluster_col]], read = x[[read_col]])]
#
#   # Cluster spans
#   span_dt <- x[, .(
#     min_start = suppressWarnings(min(as.numeric(get(start_col)), na.rm = TRUE)),
#     max_end   = suppressWarnings(max(as.numeric(get(end_col)),   na.rm = TRUE))
#   ), by = .(cluster = get(cluster_col))]
#
#   # Nodes
#   nodes <- cr[, .(cluster_n = .N, chrom = data.table::first(chrom)), by = cluster][
#     cr_detail, on = .(cluster), nomatch = 0][
#       , .(
#         chrom          = data.table::first(chrom),
#         cluster_n      = data.table::first(cluster_n),
#         node_avg_score = mean(read_cluster_score, na.rm = TRUE),
#         node_avg_len   = mean(read_cluster_len,   na.rm = TRUE)
#       ), by = cluster][
#         span_dt, on = .(cluster)]
#
#   # Edges via shared reads
#   data.table::setkey(cr, read)
#   pairs <- cr[cr, allow.cartesian = TRUE, nomatch = 0,
#               .(read, c1 = i.cluster, ch1 = i.chrom,
#                 c2 = x.cluster, ch2 = x.chrom)]
#   pairs <- pairs[c1 != c2]
#   ord <- pairs$c1 > pairs$c2
#   pairs[, `:=`(
#     cluster1 = ifelse(ord, c2, c1),
#     cluster2 = ifelse(ord, c1, c2),
#     chrom1   = ifelse(ord, ch2, ch1),
#     chrom2   = ifelse(ord, ch1, ch2)
#   )][, c("c1","c2","ch1","ch2") := NULL]
#
#   if (restrict_within_chrom) pairs <- pairs[chrom1 == chrom2]
#
#   if (!nrow(pairs)) {
#     return(list(
#       nodes = nodes[0],
#       edges = data.table(cluster1 = character(), cluster2 = character(),
#                          shared_reads = integer(),
#                          shared_avg_score = numeric(), shared_avg_len = numeric()),
#       reads = x[0]
#     ))
#   }
#
#   pairs <- read_global[pairs, on = .(read)]
#   edges <- pairs[, .(
#     shared_reads     = .N,
#     shared_avg_score = mean(read_avg_score, na.rm = TRUE),
#     shared_avg_len   = mean(read_avg_len,   na.rm = TRUE)
#   ), by = .(cluster1, cluster2)]
#
#   if (min_shared > 1L) edges <- edges[shared_reads >= min_shared]
#
#   # Keep only nodes present in edges
#   keep <- unique(c(edges$cluster1, edges$cluster2))
#   nodes <- nodes[cluster %in% keep]
#
#   # Node degree & span length
#   if (nrow(nodes)) {
#     deg_tbl <- rbindlist(list(
#       data.table::data.table(cluster = edges$cluster1),
#       data.table::data.table(cluster = edges$cluster2)
#     ))[, .(edge_count = .N), by = cluster]
#     nodes <- deg_tbl[nodes, on = .(cluster)]
#     nodes[is.na(edge_count), edge_count := 0L]
#     nodes[, total_length := as.integer(max_end - min_start)]
#   }
#
#   # --- READS: only reads present in EXACTLY TWO of the kept clusters ----------
#   # Filter to rows in kept clusters
#   x_keep <- x[get(cluster_col) %chin% keep]
#   if (nrow(x_keep)) {
#     # Count distinct kept clusters per read
#     r2 <- unique(x_keep[, .(read = .SD[[1]], cl = .SD[[2]]),
#                         .SDcols = c(read_col, cluster_col)])[
#                           , .(n_clusters = .N), by = read]
#
#     # Keep reads with exactly 2 kept clusters
#     reads_w2 <- r2[n_clusters == 2L, read]
#     reads_out <- x_keep[get(read_col) %chin% reads_w2]
#   } else {
#     reads_out <- x_keep
#   }
#
#   # Sort outputs
#   nodes <- nodes[order(-cluster_n, cluster)]
#   edges <- edges[order(-shared_reads, cluster1, cluster2)]
#   if (nrow(reads_out)) {
#     reads_out <- reads_out[order(get(chrom_col), get(start_col), get(end_col))]
#   }
#
#   list(nodes = nodes[], edges = edges[], reads = reads_out[])
# }
#


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
  req_cols <- c(cluster_col, read_col, chrom_col, start_col, end_col, "score", "aligned_length")
  miss <- setdiff(req_cols, names(dt))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

  x <- copy(dt)

  # Optional: drop "noise" (cluster_id == 0)
  if (!include_noise && id_col %in% names(x)) {
    x <- x[x[[id_col]] > 0L]
  }

  if (!nrow(x)) {
    return(list(
      nodes = data.table(cluster = character(), chrom = character(),
                         cluster_n = integer(), node_avg_score = numeric(),
                         node_avg_len = numeric(),
                         min_start = integer(), max_end = integer(),
                         edge_count = integer(), total_length = integer()),
      edges = data.table(cluster1 = character(), cluster2 = character(),
                         shared_reads = integer(),
                         shared_avg_score = numeric(), shared_avg_len = numeric()),
      reads = x[0]
    ))
  }

  # Helper: build membership table (unique per cluster-read)
  build_cr <- function(xx) {
    unique(xx[, .(cluster = .SD[[1]], read = .SD[[2]], chrom = .SD[[3]]),
              .SDcols = c(cluster_col, read_col, chrom_col)])
  }

  # Helper: compute edges from cr (+ read-level stats)
  compute_edges <- function(xx, crx) {
    # per-read global stats within xx
    read_global <- xx[, .(
      read_avg_score = mean(score, na.rm = TRUE),
      read_avg_len   = mean(aligned_length, na.rm = TRUE)
    ), by = read_col]
    setnames(read_global, read_col, "read")

    setkey(crx, read)
    pairs <- crx[crx, allow.cartesian = TRUE, nomatch = 0,
                 .(read, c1 = i.cluster, ch1 = i.chrom,
                   c2 = x.cluster, ch2 = x.chrom)]
    pairs <- pairs[c1 != c2]
    if (!nrow(pairs)) return(data.table(cluster1=character(), cluster2=character(),
                                        shared_reads=integer(),
                                        shared_avg_score=numeric(), shared_avg_len=numeric()))
    # canonicalize unordered pair
    ord <- pairs$c1 > pairs$c2
    pairs[, `:=`(
      cluster1 = ifelse(ord, c2, c1),
      cluster2 = ifelse(ord, c1, c2),
      chrom1   = ifelse(ord, ch2, ch1),
      chrom2   = ifelse(ord, ch1, ch2)
    )][, c("c1","c2","ch1","ch2") := NULL]

    if (restrict_within_chrom) pairs <- pairs[chrom1 == chrom2]
    if (!nrow(pairs)) return(data.table(cluster1=character(), cluster2=character(),
                                        shared_reads=integer(),
                                        shared_avg_score=numeric(), shared_avg_len=numeric()))

    pairs <- read_global[pairs, on = .(read)]
    edges <- pairs[, .(
      shared_reads     = .N,
      shared_avg_score = mean(read_avg_score, na.rm = TRUE),
      shared_avg_len   = mean(read_avg_len,   na.rm = TRUE)
    ), by = .(cluster1, cluster2)]

    if (min_shared > 1L) edges <- edges[shared_reads >= min_shared]
    edges[order(-shared_reads, cluster1, cluster2)]
  }

  # Helper: compute nodes from xx constrained to "keep_clusters"
  compute_nodes <- function(xx, keep_clusters) {
    if (!length(keep_clusters)) return(data.table(
      cluster=character(), chrom=character(), cluster_n=integer(),
      node_avg_score=numeric(), node_avg_len=numeric(),
      min_start=integer(), max_end=integer(), edge_count=integer(), total_length=integer()
    ))
    yy <- xx[get(cluster_col) %chin% keep_clusters]

    # per-(cluster, read) detail (for node averages)
    cr_detail <- yy[, .(
      read_cluster_score = mean(score, na.rm = TRUE),
      read_cluster_len   = mean(aligned_length, na.rm = TRUE)
    ), by = .(cluster = yy[[cluster_col]], read = yy[[read_col]])]

    # cluster size = unique reads per cluster (post-filter)
    cl_sizes <- unique(yy[, .(cluster = .SD[[1]], read = .SD[[2]], chrom = .SD[[3]]),
                          .SDcols = c(cluster_col, read_col, chrom_col)])[
                            , .(cluster_n = .N, chrom = first(chrom)), by = cluster]

    # spans
    spans <- yy[, .(
      min_start = suppressWarnings(min(as.numeric(get(start_col)), na.rm = TRUE)),
      max_end   = suppressWarnings(max(as.numeric(get(end_col)),   na.rm = TRUE))
    ), by = .(cluster = get(cluster_col))]

    nodes <- cl_sizes[cr_detail, on = .(cluster), nomatch = 0][
      , .(
        chrom          = first(chrom),
        cluster_n      = first(cluster_n),
        node_avg_score = mean(read_cluster_score, na.rm = TRUE),
        node_avg_len   = mean(read_cluster_len,   na.rm = TRUE)
      ), by = cluster][
        spans, on = .(cluster)]

    nodes
  }

  # -------- Initial pass to get a seed set of clusters ("keep") --------
  cr0    <- build_cr(x)
  edges0 <- compute_edges(x, cr0)
  keep   <- unique(c(edges0$cluster1, edges0$cluster2))

  if (!length(keep)) {
    if (verbose) message("No edges after initial pass; returning empty network.")
    return(list(
      nodes = data.table(cluster=character(), chrom=character(), cluster_n=integer(),
                         node_avg_score=numeric(), node_avg_len=numeric(),
                         min_start=integer(), max_end=integer(),
                         edge_count=integer(), total_length=integer()),
      edges = data.table(cluster1=character(), cluster2=character(),
                         shared_reads=integer(),
                         shared_avg_score=numeric(), shared_avg_len=numeric()),
      reads = x[0]
    ))
  }

  # ---------- Iterate to enforce "reads in exactly 2 kept clusters" ----------
  iter <- 1L
  read_set <- NULL
  repeat {
    # restrict to current keep-set
    x_keep <- x[get(cluster_col) %chin% keep]
    if (!nrow(x_keep)) break

    # membership on keep
    crk <- build_cr(x_keep)

    # reads in exactly 2 of the kept clusters
    r_counts <- crk[, .(n_clusters = .N), by = read]
    reads2   <- r_counts[n_clusters == 2L, read]

    # if no reads satisfy, break
    if (!length(reads2)) {
      if (verbose) message("No reads with exactly 2 clusters among kept; returning empty.")
      return(list(
        nodes = data.table(cluster=character(), chrom=character(), cluster_n=integer(),
                           node_avg_score=numeric(), node_avg_len=numeric(),
                           min_start=integer(), max_end=integer(),
                           edge_count=integer(), total_length=integer()),
        edges = data.table(cluster1=character(), cluster2=character(),
                           shared_reads=integer(),
                           shared_avg_score=numeric(), shared_avg_len=numeric()),
        reads = x[0]
      ))
    }

    # filter rows to those reads
    xf <- x_keep[get(read_col) %chin% reads2]

    # rebuild cr, edges on filtered rows
    crf    <- build_cr(xf)
    edgesf <- compute_edges(xf, crf)
    keep2  <- unique(c(edgesf$cluster1, edgesf$cluster2))

    # Recompute read set: exactly 2 among *new* keep
    crf_keep <- crf[cluster %chin% keep2]
    r2_again <- crf_keep[, .(n_clusters = .N), by = read]
    reads2b  <- r2_again[n_clusters == 2L, read]
    xf2      <- xf[get(read_col) %chin% reads2b]

    # If stable, stop
    if (setequal(keep2, keep) && setequal(reads2b, read_set)) {
      x_final    <- xf2
      edges_final <- edgesf
      keep_final  <- keep2
      break
    }

    # Update and loop
    keep    <- keep2
    read_set <- reads2b
    iter <- iter + 1L
    if (iter > max_iter) {
      if (verbose) message("Reached max_iter; using current filtered network.")
      x_final    <- xf2
      edges_final <- edgesf
      keep_final  <- keep2
      break
    }
  }

  # If edges vanish after filtering, return empties
  if (!nrow(edges_final)) {
    if (verbose) message("No edges after filtering to reads in exactly two clusters.")
    return(list(
      nodes = data.table(cluster=character(), chrom=character(), cluster_n=integer(),
                         node_avg_score=numeric(), node_avg_len=numeric(),
                         min_start=integer(), max_end=integer(),
                         edge_count=integer(), total_length=integer()),
      edges = data.table(cluster1=character(), cluster2=character(),
                         shared_reads=integer(),
                         shared_avg_score=numeric(), shared_avg_len=numeric()),
      reads = x[0]
    ))
  }

  # Recompute nodes on the final filtered rows and final kept clusters
  nodes_final <- compute_nodes(x_final, keep_final)

  # Degree (edge_count) & total_length
  deg_tbl <- rbindlist(list(
    data.table(cluster = edges_final$cluster1),
    data.table(cluster = edges_final$cluster2)
  ))[, .(edge_count = .N), by = cluster]
  nodes_final <- deg_tbl[nodes_final, on = .(cluster)]
  nodes_final[is.na(edge_count), edge_count := 0L]
  nodes_final[, total_length := as.integer(max_end - min_start)]

  # Sort outputs
  nodes_final <- nodes_final[order(-cluster_n, cluster)]
  edges_final <- edges_final[order(-shared_reads, cluster1, cluster2)]
  x_final     <- x_final[order(get(chrom_col), get(start_col), get(end_col))]

  if (verbose) {
    message(sprintf("Network built in %d iteration(s). Kept %d clusters, %d edges, %d reads.",
                    iter, length(keep_final), nrow(edges_final),
                    data.table::uniqueN(x_final[[read_col]])))
  }

  list(nodes = nodes_final[], edges = edges_final[], reads = x_final[])
}


#' Plot a shared-read cluster network (ggraph)
#' - Nodes: size ~ cluster_n (capped), color = chrom, alpha = node_avg_score (default)
#' - Edges: fixed thin black lines (0.25)
#' Filtering rule: keep all nodes with cluster_n >= min_cluster_n, PLUS any node
#' directly connected to one of those (neighbors). Drop other low-only islands.
#'
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
                                        min_cluster_n = 1L,    # threshold
                                        keep_neighbors_of_kept = TRUE) {
  stopifnot(is.list(net), all(c("nodes","edges") %in% names(net)))
  nodes <- data.table::as.data.table(net$nodes)
  edges <- data.table::as.data.table(net$edges)

  # Early exits
  if (!nrow(nodes) || !nrow(edges)) {
    message("No nodes/edges to plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0, y = 0, label = "No edges to plot"))
  }

  # -------- Filter: keep qualifying nodes AND their neighbors --------
  if (!is.null(min_cluster_n) && min_cluster_n > 1L) {
    kept <- nodes[edge_count>min_cluster_n | cluster_n >= min_cluster_n, cluster]

    if (isTRUE(keep_neighbors_of_kept) && length(kept)) {
      nbr1 <- edges[cluster1 %in% kept, cluster2]
      nbr2 <- edges[cluster2 %in% kept, cluster1]
      kept <- unique(c(kept, nbr1, nbr2))
    }

    nodes <- nodes[cluster %in% kept]
    if (!nrow(nodes)) {
      message("No nodes remain after threshold+neighbor filter.")
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::annotate("text", x = 0, y = 0, label = "No nodes after filter"))
    }
    edges <- edges[cluster1 %in% nodes$cluster & cluster2 %in% nodes$cluster]
    if (!nrow(edges)) {
      message("No edges remain after filter; plotting kept nodes without edges.")
    }
  }

  # -------- Alpha mapping (0..1) --------
  if (!is.null(alpha_var) && alpha_var %in% names(nodes)) {
    rng <- range(nodes[[alpha_var]], na.rm = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
      nodes[, node_alpha := (get(alpha_var) - rng[1]) / (rng[2] - rng[1])]
    } else nodes[, node_alpha := 1]
  } else nodes[, node_alpha := 1]

  # -------- Build graph & plot --------
  g <- igraph::graph_from_data_frame(
    d = if (nrow(edges)) edges[, .(from = cluster1, to = cluster2)] else NULL,
    directed = FALSE,
    vertices = nodes[, .(name = if ("name" %in% names(nodes)) name else cluster,
                         cluster = cluster,
                         chrom = chrom,
                         cluster_n = cluster_n,
                         node_alpha = node_alpha)]
  )

  set.seed(seed)
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(color = "black", width = 0.25, alpha = 0.6, show.legend = FALSE) +
    ggraph::geom_node_point(
      ggplot2::aes(size = pmin(cluster_n, cap), color = chrom, alpha = node_alpha),
      stroke = 0.2
    ) +
    ggplot2::scale_size_binned(
      range   = c(0.8, 5),
      breaks  = c(1, 3, 5, 10, 20, cap),
      labels  = c("1","3","5","10","20", paste0(cap, "+")),
      limits  = c(1, cap),
      nice.breaks = FALSE,
      guide   = ggplot2::guide_bins(title = "Cluster size")
    ) +
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::theme_void(base_size = 6) +
    ggplot2::labs(color = "Chromosome")

  if (isTRUE(label)) {
    lab <- if (label_col %in% names(nodes)) label_col else "cluster"
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = .data[[lab]]),
      repel = isTRUE(label_repel),
      size = label_size
    )
  }
  p
}



plot_clusters_raw<-function(cluster_dt){
  cluster_dt<-cluster_dt[cluster_n>1]
  cluster_dt[,ID:=rev(1:nrow(cluster_dt))]
  cluster_dt[,loc:=round(mean((start))), by=chrom_cluster_id]
  ggplot(cluster_dt)+
    geom_segment(aes(x=start, xend=end, y=ID, yend=ID))+
    geom_point(aes(x=start,y=ID), col="red", size=0.1)+
    geom_point(aes(x=end,y=ID), col="red", size=0.1)+
    facet_wrap(~chrom_cluster_id, scales="free")+
    theme(axis.text = element_blank(),
          axis.title = element_blank())
}




plot_clusters_raw<-function(chr_sw){
  chr_sw<-chr_sw[cluster_n>1]
  chr_sw[,ID:=rev(1:nrow(chr_sw))]
  chr_sw[,loc:=round(mean((start))), by=cluster_id]
  ggplot(chr_sw)+
    geom_segment(aes(x=start, xend=end, y=ID, yend=ID))+
    geom_point(aes(x=start,y=ID), col="red", size=0.1)+
    geom_point(aes(x=end,y=ID), col="red", size=0.1)+
    facet_wrap(~chrom_cluster_id, scales="free")+
    theme(axis.text = element_blank(),
          axis.title = element_blank())

}

# nodes_to_fasta:
# - Takes net (with net$nodes) and a genome loaded via seqinr::read.fasta()
# - Extracts [min_start:max_end] for each node's chrom, with optional padding
# - Supports 0-based coords (typical for BED-like ranges) -> convert to 1-based R indexing
# - Returns a named character vector of sequences; can also write a FASTA file.

library(data.table)

nodes_to_fasta <- function(net,
                           genome,
                           clusters = NULL,                 # vector of names to keep (e.g., "Chr2_7354979"); NULL = all
                           chrom_col = "chrom",
                           start_col = "min_start",
                           end_col   = "max_end",
                           name_col  = "cluster",
                           pad_bp    = 0L,                  # add padding on both sides
                           index0    = TRUE,                # TRUE: start/end are 0-based; convert to 1-based for R
                           clamp     = TRUE,                # clamp to chromosome bounds
                           tolower   = FALSE,               # force lowercase output
                           outfile   = NULL,                # path to write FASTA (optional)
                           width     = 80                   # FASTA line width if writing
) {
  stopifnot(is.list(net), "nodes" %in% names(net))
  nodes <- as.data.table(net$nodes)

  # filter clusters if requested
  if (!is.null(clusters)) {
    nodes <- nodes[get(name_col) %in% clusters]
  }
  if (!nrow(nodes)) stop("No nodes to export (after filtering).")

  # sanity checks
  req <- c(chrom_col, start_col, end_col, name_col)
  miss <- setdiff(req, names(nodes))
  if (length(miss)) stop("Missing columns in net$nodes: ", paste(miss, collapse = ", "))

  # build sequences
  seqs <- vector("list", nrow(nodes))
  nms  <- character(nrow(nodes))
  probs <- list()

  for (i in seq_len(nrow(nodes))) {
    chrom <- nodes[[chrom_col]][i]
    if (!chrom %in% names(genome)) {
      probs[[length(probs)+1]] <- sprintf("Chromosome '%s' not found in FASTA; skipping %s",
                                          chrom, nodes[[name_col]][i])
      next
    }

    chr_seq <- genome[[chrom]]        # seqinr::read.fasta: vector of letters
    chr_len <- length(chr_seq)

    # coordinates
    s0 <- suppressWarnings(as.numeric(nodes[[start_col]][i]))
    e0 <- suppressWarnings(as.numeric(nodes[[end_col]][i]))
    if (!is.finite(s0) || !is.finite(e0)) {
      probs[[length(probs)+1]] <- sprintf("Non-finite coords for %s; skipping", nodes[[name_col]][i])
      next
    }
    # convert to 1-based if needed; pad; clamp
    s1 <- as.integer(s0) + if (isTRUE(index0)) 1L else 0L
    e1 <- as.integer(e0) + if (isTRUE(index0)) 1L else 0L
    s  <- s1 - as.integer(pad_bp)
    e  <- e1 + as.integer(pad_bp)

    if (isTRUE(clamp)) {
      s <- max(1L, s)
      e <- min(chr_len, e)
    }

    if (s > e || s < 1L || e > chr_len) {
      probs[[length(probs)+1]] <- sprintf("Out-of-bounds/empty range for %s: %s:%d-%d (len=%d)",
                                          nodes[[name_col]][i], chrom, s, e, chr_len)
      next
    }

    subseq <- chr_seq[s:e]
    if (tolower) subseq <- tolower(subseq)
    seqs[[i]] <- paste0(subseq, collapse = "")
    nms[i]    <- nodes[[name_col]][i]
  }

  # drop failures
  keep <- !vapply(seqs, is.null, logical(1))
  seqs <- unlist(seqs[keep], use.names = FALSE)
  names(seqs) <- nms[keep]

  if (length(seqs) == 0L) stop("No sequences extracted. ", if (length(probs)) paste(unlist(probs), collapse = " | "))

  # optional write to FASTA
  if (!is.null(outfile)) {
    if (!requireNamespace("seqinr", quietly = TRUE)) stop("Please install 'seqinr' to write FASTA.")
    seqinr::write.fasta(sequences = as.list(seqs),
                        names     = names(seqs),
                        nbchar    = as.integer(width),
                        file.out  = outfile)
  }

  # report any issues
  if (length(probs)) message("Warnings: ", paste(unlist(probs), collapse = " | "))

  return(seqs)  # named character vector (FASTA-ready strings)
}


estimate_modes_1d <- function(x, bw = NULL, bw_method = c("SJ","nrd0","ucv","bcv"),
                              kernel = "gaussian", grid_n = 2048) {
  bw_method <- match.arg(bw_method)
  if (is.null(bw)) {
    bw <- switch(bw_method,
                 SJ   = bw.SJ(x, method = "ste"),
                 nrd0 = bw.nrd0(x),
                 ucv  = bw.ucv(x),
                 bcv  = bw.bcv(x))
  }
  kd <- density(x, bw = bw, kernel = kernel, n = grid_n)
  y <- kd$y; idx <- which(diff(sign(diff(y))) == -2) + 1
  list(bw = kd$bw, modes = kd$x[idx], density = kd)
}

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
      SJ   = bw.SJ(x, method = "ste"),
      nrd0 = bw.nrd0(x),
      ucv  = bw.ucv(x),
      bcv  = bw.bcv(x)
    )
  }

  L <- min(x); R <- max(x)

  left_keep  <- x < (L + reflect_multiplier * bw)
  right_keep <- x > (R - reflect_multiplier * bw)

  x_left_ref  <- 2 * L - x[left_keep]
  x_right_ref <- 2 * R - x[right_keep]

  x_ref <- c(x, x_left_ref, x_right_ref)

  kd <- density(
    x_ref,
    bw     = bw,
    kernel = kernel,
    from   = L,
    to     = R,
    n      = grid_n
  )

  # raw (all) modes
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
# assumes `kde_reflect()` is already defined in your session
plot_node_breaks_kde <- function(net,
                                 node,
                                 cluster_col = "chrom_cluster_id",
                                 bw_method   = "SJ",
                                 bins        = 50L,
                                 kernel      = "gaussian") {

  # 1) pull starts/ends for that node
  dt_reads <- as.data.table(net$reads)
  if (!all(c(cluster_col, "start", "end") %in% names(dt_reads))) {
    stop("net$reads must have columns: ", cluster_col, ", start, end")
  }

  node_rows <- dt_reads[get(cluster_col) == node]
  if (!nrow(node_rows)) {
    stop(sprintf("No reads found for node '%s'", node))
  }

  start_ends <- c(node_rows$start, node_rows$end)

  # 2) KDE with reflection
  res_kde <- kde_reflect(start_ends, bw_method = bw_method, kernel = kernel)

  # 3) histogram binwidth from data range / bins
  dt <- data.table(x = start_ends)
  bw_hist <- diff(range(dt$x)) / bins

  # 4) build density df & scale to counts
  dens_df <- data.frame(
    x = res_kde$density$x,
    y = res_kde$density$y
  )
  dens_df$y_scaled <- dens_df$y * nrow(dt) * bw_hist

  # 5) plot
  p <- ggplot(dt, aes(x = x)) +
    geom_histogram(binwidth = bw_hist,
                   fill = "grey90", color = "grey80") +
    geom_line(data = dens_df,
              aes(x = x, y = y_scaled),
              color = "steelblue", linewidth = 1) +
    labs(
      title = paste0("Kernel density for ", node),
      x = "Read breaks",
      y = "Count"
    ) +
    theme_bw(base_size = 6)

  # return both the plot and the data you might want to inspect
  list(
    plot     = p,
    breaks   = start_ends,
    dens_df  = dens_df,
    bw_kde   = res_kde$density$bw,
    bw_hist  = bw_hist
  )
}

#
# estimate_modes_2critical <- function(x,
#                                      k = 2,
#                                      kernel = "gaussian",
#                                      grid_n = 2048,
#                                      h_start = NULL,
#                                      h_max_factor = 1.5,
#                                      step = 1.25) {
#   x <- x[is.finite(x)]
#   if (length(x) < 5L) stop("not enough data")
#
#   # start from something reasonable, but we'll only increase it
#   if (is.null(h_start)) {
#     h_start <- bw.nrd0(x)  # could also use bw.SJ(x)
#   }
#
#   xr <- range(x)
#   h_max <- diff(xr) * h_max_factor  # don't let it grow forever
#
#   count_modes <- function(d) {
#     y <- d$y
#     # interior local maxima of the density
#     sum(diff(sign(diff(y))) == -2)
#   }
#
#   h <- h_start
#   repeat {
#     d <- density(x, bw = h, kernel = kernel, n = grid_n,
#                  from = xr[1], to = xr[2])
#     m <- count_modes(d)
#     if (m <= k || h >= h_max) {
#       return(list(
#         bw = h,
#         modes = d$x[which(diff(sign(diff(d$y))) == -2) + 1],
#         density = d,
#         n_modes = m
#       ))
#     }
#     h <- h * step  # oversmooth a bit more
#   }
# }
