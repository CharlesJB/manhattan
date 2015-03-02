#' Produce a manhattan plot.
#'
#' @param gr_snps A \code{GRanges} object.
#' @param title The title of the plot.
#' @param pvalue The name of the column to use as pvalue. Default: NULL (the
#'   first numerical metadata column will be used)
#' @param status Name of the column for the SNP status (Genotype or Imputed)
#' @param r2 The name of the column to use for showing r2. Default: NULL.
#' @param geno_ld The file *.geno.ld produced by `vcftools` or a
#'   \code{data.frame} from the file loaded with read.table.
#'
#' @return \code{} corresponding to the input file table.
#'
#' @examples
#' dataset <- load_dataset()
#'
#' @export
manhattan <- function(gr_snps, title, pvalue = NULL, status = NULL, r2 = NULL,
                      geno_ld = NULL) {
  # TODO: Add parameter validity tests
  # Name the metadata column correctly
  if (is.null(pvalue)) {
    i <- min(which(sapply(GenomicRanges::mcols(gr_snps), is.numeric)))
    pvalue <- colnames(GenomicRanges::mcols(gr_snps))[i]
  }
  gr_snps$pvalue <- GenomicRanges::mcols(gr_snps)[[pvalue]]
  if (!is.null(status)) {
    gr_snps$status <- GenomicRanges::mcols(gr_snps)[[status]]
  }
  if (!is.null(r2)) {
    gr_snps$r2 <- GenomicRanges::mcols(gr_snps)[[r2]]
  }
  if (!is.null(geno_ld)) {
    gr_snps$r2 <- add_r2_from_vcftools(geno_ld = geno_ld, gr_snps = gr_snps,
                    pvalue = pvalue, status = status)
    r2 = "r2"
  }

  ## Produce the manhattan plot
  # Prepare aesthetic
  aesthetic <- get_aes(status = status, r2 = r2)

  # Prepare basic plot
  p <- ggbio::ggplot(gr_snps, aesthetic) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw()

  # Add color gradient (if necessary)
  if (!is.null(r2)) {
    p <- p + ggplot2::scale_color_continuous(low = "#E6E6E6",
                                             high = "steelblue")
  }

  # Add labels
  chr <- unique(GenomicRanges::seqnames(gr_snps))
  stopifnot(length(chr) == 1)
  p <- p + ggplot2::ggtitle(title) +
    ggplot2::xlab(paste0("Position on ", chr)) +
    ggplot2::ylab(paste0("-log(", pvalue, ")"))

  # Return plot
  p
}

get_aes <- function(status, r2) {
  if (is.null(status)) {
    if (is.null(r2)) {
      ggplot2::aes(x = start, y = -log10(pvalue))
    } else {
      ggplot2::aes(x = start, y = -log10(pvalue), color = r2)
    }
  } else {
    if (is.null(r2)) {
      ggplot2::aes(x = start, y = -log10(pvalue), shape = status)
    } else {
      ggplot2::aes(x = start, y = -log10(pvalue), shape = status, color = r2)
    }
  }
}
