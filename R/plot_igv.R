#' IGV-style read pile-up plot
#'
#' @title Plot BED alignments in IGV style
#' @description Plots reads overlapping a genomic region as horizontal segments
#'   packed into non-overlapping tracks (like IGV). Supports BED12 block
#'   structure and optional colour mapping.
#'
#' @param bed_dt A \code{data.table} from \code{\link{read_bed12}} or any
#'   BED-like table with \code{chrom}, \code{start}, \code{end} columns.
#' @param chrom Character. Target chromosome.
#' @param region_start Integer. Start of the display region (bp).
#' @param region_end Integer. End of the display region (bp).
#' @param max_reads Integer. Maximum number of reads to display (downsampled
#'   randomly if exceeded). Default \code{5000}.
#' @param color_var Column to map to colour (e.g. \code{"strand"} or
#'   \code{"nameN"}). Set to \code{NULL} for no colour mapping. Default
#'   \code{"strand"}.
#' @param show_blocks Logical. If \code{TRUE} and BED12 block columns are
#'   present, draw individual alignment blocks. Default \code{TRUE}.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' bed <- read_bed12(system.file("extdata", "alignVA2.bed12", package = "kolibri"))
#' plot_bed_igv(bed, chrom = "Chr1", region_start = 1e6, region_end = 1.1e6)
#' }
#'
#' @family visualization
#' @importFrom data.table as.data.table
#' @importFrom IRanges IRanges disjointBins
#' @importFrom ggplot2 ggplot aes geom_segment scale_x_continuous scale_y_reverse
#'   theme_classic theme element_blank unit
#' @export
plot_bed_igv <- function(bed_dt,
                         chrom,
                         region_start,
                         region_end,
                         max_reads   = 5000,
                         color_var   = "strand",
                         show_blocks = TRUE) {
  bed_dt <- data.table::as.data.table(bed_dt)

  target_chrom <- chrom

  reads <- bed_dt[
    chrom == target_chrom &
      end   > region_start &
      start < region_end
  ]

  if (!nrow(reads)) stop("No reads found in that region.")

  if (nrow(reads) > max_reads) {
    message("Subsampling ", max_reads, " of ", nrow(reads), " reads for plotting.")
    set.seed(1L)
    reads <- reads[sample(.N, max_reads)]
  }

  ir <- IRanges::IRanges(start = reads$start + 1L, end = reads$end)
  reads[, track := IRanges::disjointBins(ir)]
  reads[, read_id := .I]

  use_color <- FALSE
  if (!is.null(color_var)) {
    if (!is.character(color_var)) color_var <- deparse(substitute(color_var))
    if (color_var %in% names(reads)) {
      reads[, .col := reads[[color_var]]]
      use_color <- TRUE
    } else {
      warning("color_var '", color_var, "' not found; ignoring colour mapping.")
      reads[, .col := NA]
    }
  } else {
    reads[, .col := NA]
  }

  if (isTRUE(show_blocks) &&
      all(c("blockCount","blockSizes","blockStarts") %in% names(reads))) {
    block_dt <- reads[, {
      bc <- as.integer(blockCount)
      if (is.na(bc) || bc <= 1L || is.na(blockSizes) || blockSizes == "") {
        data.table::data.table(block_start = start, block_end = end)
      } else {
        sizes  <- as.integer(strsplit(blockSizes,  ",")[[1]])
        rel_st <- as.integer(strsplit(blockStarts, ",")[[1]])
        sizes  <- sizes[!is.na(sizes)]
        rel_st <- rel_st[!is.na(rel_st)]
        n      <- min(length(sizes), length(rel_st))
        if (n == 0L) {
          data.table::data.table(block_start = start, block_end = end)
        } else {
          data.table::data.table(
            block_start = start + rel_st[seq_len(n)],
            block_end   = start + rel_st[seq_len(n)] + sizes[seq_len(n)]
          )
        }
      }
    }, by = .(read_id, chrom, track, .col)]
  } else {
    block_dt <- reads[, .(
      read_id = read_id, chrom = chrom, track = track, .col = .col,
      block_start = start, block_end = end
    )]
  }

  block_dt[block_start < region_start, block_start := region_start]
  block_dt[block_end   > region_end,   block_end   := region_end]

  if (use_color) {
    p <- ggplot2::ggplot(block_dt) +
      ggplot2::geom_segment(
        ggplot2::aes(x = block_start, xend = block_end,
                     y = track, yend = track, colour = .col),
        linewidth = 0.25, lineend = "round"
      )
  } else {
    p <- ggplot2::ggplot(block_dt) +
      ggplot2::geom_segment(
        ggplot2::aes(x = block_start, xend = block_end,
                     y = track, yend = track),
        linewidth = 0.25, lineend = "round"
      )
  }

  p +
    ggplot2::scale_x_continuous(
      limits = c(region_start, region_end),
      name   = "Position (bp)"
    ) +
    ggplot2::scale_y_reverse(
      name   = "Reads (track)",
      breaks = NULL
    ) +
    ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.25, "cm")
    )
}
