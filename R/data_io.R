#' Read a BED12 file
#'
#' @title Read a BED12 alignment file
#' @description Reads a BED12 file into a \code{data.table}. Multi-block rows
#'   (split reads with \code{blockCount > 1}) are accepted with a message.
#'   Derived columns \code{is_split}, \code{aligned_length}, and \code{nameN}
#'   are added automatically.
#'
#' @param path Character. Path to the BED12 file.
#' @param key_on Character vector of column names to use as the data.table key.
#'   Default is \code{c("chrom","start","end")}. Set to \code{NULL} to skip
#'   keying.
#'
#' @return A keyed \code{data.table} with the 12 standard BED columns plus
#'   \code{is_split} (logical), \code{aligned_length} (integer), and
#'   \code{nameN} (integer count of occurrences of each read name).
#'
#' @examples
#' \dontrun{
#' bed <- read_bed12(system.file("extdata", "alignVA2.bed12", package = "kolibri"))
#' head(bed)
#' }
#'
#' @family data I/O
#' @importFrom data.table fread setkeyv
#' @export
read_bed12 <- function(path, key_on = c("chrom","start","end")) {
  cols <- c("chrom","start","end","name","score","strand",
            "thickStart","thickEnd","itemRgb","blockCount",
            "blockSizes","blockStarts")
  cc <- c("character","integer","integer","character","integer","character",
          "integer","integer","character","integer",
          "character","character")

  dt <- data.table::fread(path, sep = "\t", header = FALSE, fill = TRUE,
                          quote = "", col.names = cols, colClasses = cc,
                          showProgress = FALSE)

  n_multi <- sum(dt$blockCount > 1L, na.rm = TRUE)
  if (n_multi > 0L) {
    message(sprintf("Note: detected %d row(s) with blockCount > 1 (multi-block BED12).",
                    n_multi))
  }

  dt[, is_split       := blockCount > 1L]
  dt[, aligned_length := pmax(0L, end - start)]
  dt[, nameN          := .N, by = name]

  if (!is.null(key_on)) data.table::setkeyv(dt, key_on)
  dt[]
}
