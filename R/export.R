#' Export network node sequences to FASTA
#'
#' @title Export node sequences to FASTA
#' @description Extracts the genomic sequences corresponding to network node
#'   spans from a loaded FASTA genome and returns them as a named character
#'   vector. Optionally writes a FASTA file.
#'
#' @param net A network list from \code{\link{build_cluster_network}}.
#' @param genome A named list of character vectors or \code{seqinr}-style
#'   FASTA sequences (e.g. loaded with \code{seqinr::read.fasta()}).
#' @param clusters Optional character vector of cluster names to export. If
#'   \code{NULL} (default), all nodes are exported.
#' @param chrom_col Name of the chromosome column in \code{net$nodes}. Default
#'   \code{"chrom"}.
#' @param start_col Name of the start column in \code{net$nodes}. Default
#'   \code{"min_start"}.
#' @param end_col Name of the end column in \code{net$nodes}. Default
#'   \code{"max_end"}.
#' @param name_col Name of the node name column. Default \code{"cluster"}.
#' @param pad_bp Integer. Padding in bp to add on each side of the interval.
#'   Default \code{0L}.
#' @param index0 Logical. If \code{TRUE} (default), coordinates are assumed to
#'   be 0-based (BED convention) and converted to 1-based R indexing.
#' @param clamp Logical. Clamp coordinates to chromosome bounds. Default
#'   \code{TRUE}.
#' @param tolower Logical. Convert sequence to lowercase. Default \code{FALSE}.
#' @param outfile Optional path to write a FASTA file. Default \code{NULL}.
#' @param width Integer. FASTA line width when writing. Default \code{80}.
#'
#' @return A named character vector of sequences (one element per node).
#'
#' @examples
#' \dontrun{
#' genome <- seqinr::read.fasta(
#'   system.file("extdata", "Col-CEN_v1.2.fasta.gz", package = "kolibri"))
#' seqs <- nodes_to_fasta(net, genome, outfile = "nodes.fa")
#' }
#'
#' @family export
#' @importFrom data.table as.data.table
#' @export
nodes_to_fasta <- function(net,
                           genome,
                           clusters  = NULL,
                           chrom_col = "chrom",
                           start_col = "min_start",
                           end_col   = "max_end",
                           name_col  = "cluster",
                           pad_bp    = 0L,
                           index0    = TRUE,
                           clamp     = TRUE,
                           tolower   = FALSE,
                           outfile   = NULL,
                           width     = 80) {
  stopifnot(is.list(net), "nodes" %in% names(net))
  nodes <- data.table::as.data.table(net$nodes)

  if (!is.null(clusters)) nodes <- nodes[get(name_col) %in% clusters]
  if (!nrow(nodes)) stop("No nodes to export (after filtering).")

  req  <- c(chrom_col, start_col, end_col, name_col)
  miss <- setdiff(req, names(nodes))
  if (length(miss))
    stop("Missing columns in net$nodes: ", paste(miss, collapse = ", "))

  seqs  <- vector("list", nrow(nodes))
  nms   <- character(nrow(nodes))
  probs <- list()

  for (i in seq_len(nrow(nodes))) {
    chrom   <- nodes[[chrom_col]][i]
    if (!chrom %in% names(genome)) {
      probs[[length(probs)+1]] <- sprintf(
        "Chromosome '%s' not found in FASTA; skipping %s",
        chrom, nodes[[name_col]][i])
      next
    }

    chr_seq <- genome[[chrom]]
    chr_len <- length(chr_seq)

    s0 <- suppressWarnings(as.numeric(nodes[[start_col]][i]))
    e0 <- suppressWarnings(as.numeric(nodes[[end_col]][i]))
    if (!is.finite(s0) || !is.finite(e0)) {
      probs[[length(probs)+1]] <- sprintf("Non-finite coords for %s; skipping",
                                          nodes[[name_col]][i])
      next
    }

    s1 <- as.integer(s0) + if (isTRUE(index0)) 1L else 0L
    e1 <- as.integer(e0) + if (isTRUE(index0)) 1L else 0L
    s  <- s1 - as.integer(pad_bp)
    e  <- e1 + as.integer(pad_bp)

    if (isTRUE(clamp)) { s <- max(1L, s); e <- min(chr_len, e) }

    if (s > e || s < 1L || e > chr_len) {
      probs[[length(probs)+1]] <- sprintf(
        "Out-of-bounds/empty range for %s: %s:%d-%d (len=%d)",
        nodes[[name_col]][i], chrom, s, e, chr_len)
      next
    }

    subseq <- chr_seq[s:e]
    if (tolower) subseq <- tolower(subseq)
    seqs[[i]] <- paste0(subseq, collapse = "")
    nms[i]    <- nodes[[name_col]][i]
  }

  keep <- !vapply(seqs, is.null, logical(1))
  seqs <- unlist(seqs[keep], use.names = FALSE)
  names(seqs) <- nms[keep]

  if (length(seqs) == 0L)
    stop("No sequences extracted. ",
         if (length(probs)) paste(unlist(probs), collapse = " | "))

  if (!is.null(outfile)) {
    if (!requireNamespace("seqinr", quietly = TRUE))
      stop("Please install 'seqinr' to write FASTA.")
    seqinr::write.fasta(sequences = as.list(seqs),
                        names     = names(seqs),
                        nbchar    = as.integer(width),
                        file.out  = outfile)
  }

  if (length(probs))
    message("Warnings: ", paste(unlist(probs), collapse = " | "))

  return(seqs)
}
