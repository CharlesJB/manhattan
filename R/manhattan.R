#' Produce a manhattan plot.
#'
#' @param gr_snps A \code{GRanges} object.
#' @param title The title of the plot.
#' @param pvalue The name of the column to use as pvalue. Default: NULL (the
#'   first numerical metadata column will be used)
#' @param status Name of the column for the SNP status (Genotype or Imputed)
#' @param show_labels SNPs name to add to the graph in a character vector.
#'                    Default: NULL.
#' @param r2 The name of the column to use for showing r2. Default: NULL.
#' @param geno_ld The file *.geno.ld produced by `vcftools` or a
#'   \code{data.frame} from the file loaded with read.table.
#' @param extra_ggplot2 Extra modifications to add to the graph. Must be
#'   valid ggplot2 expressions.
#'
#' @return \code{} corresponding to the input file table.
#'
#' @examples
#' dataset <- load_dataset()
#'
#' @export
manhattan <- function(gr_snps, title, pvalue = NULL, status = NULL,
                      show_labels = NULL, r2 = NULL, geno_ld = NULL,
		      extra_ggplot2 = NULL) {
  # TODO: Add parameter validity tests
  #       Must check if a column
  # Status should be genotype or imputed only
  # Name the metadata column correctly
  metadata <- GenomicRanges::mcols(gr_snps)
  GenomicRanges::mcols(gr_snps) <- NULL
  gr_snps$name <- metadata[["name"]]
  if (is.null(pvalue)) {
    i <- min(which(sapply(metadata, is.numeric)))
    pvalue <- colnames(metadata)[i]
  }
  gr_snps$pvalue <- metadata[[pvalue]]
  if (!is.null(status)) {
    i <- grepl("genotyped", metadata[[status]],
               ignore.case = TRUE)
    gr_snps$status <- 21
    gr_snps$status[i] <- 24
  }
  if (!is.null(show_labels)) {
     i <- metadata[["name"]] %in% show_labels
     gr_snps$labels <- ""
     gr_snps$labels[i] <- as.character(metadata[["name"]])[i]
     gr_snps$labels_color <- "no"
     gr_snps$labels_color[i] <- "yes"
  }
  if (!is.null(r2)) {
    gr_snps$r2 <- metadata[[r2]]
  }
  if (!is.null(geno_ld)) {
    gr_snps$r2 <- add_r2_from_vcftools(geno_ld = geno_ld, gr_snps = gr_snps,
                    pvalue = pvalue, status = status)
    r2 = "r2"
  }

  ## Produce the manhattan plot
  # Prepare aesthetic
  aesthetic <- get_aes(status = status, r2 = r2, show_labels = show_labels)

  # Bind data
  p <- ggbio::ggplot(gr_snps, aesthetic)


  # Add shape (if necessary)
  if (!is.null(status)) {
    p <- p + ggplot2::scale_shape_identity(guide = "legend", breaks = c(21,24),
                                           labels = c("Imputed", "Genotyped"))
  }

  # Add border color (if necessary)
  if (!is.null(show_labels)) {
    p <- p + scale_color_manual(values = c("black", "blue"), guide = FALSE)
  } else {
    p <- p + scale_color_manual(values = "black", guide = FALSE)
  }

  # Add geom_points
  p <- p + get_geom_point(status = status, r2 = r2)

  # Add theme
  p <- p + ggplot2::theme_bw()

  # Add fill gradient (if necessary)
  if (!is.null(r2)) {
    p <- p + ggplot2::scale_fill_continuous(low = "white", high = "red")
  }

  # Show SNPs labels
  if (!is.null(show_labels)) {
    p <- p + geom_text(aes(label = labels), size = 4, hjust = 0.7, vjust = 1.2,
                       angle = 45, color = "blue")
  }

  # Add labels
  chr <- unique(GenomicRanges::seqnames(gr_snps))
  stopifnot(length(chr) == 1)
  p <- p + ggplot2::ggtitle(title) +
    ggplot2::xlab(paste0("Position on ", chr, " (bp)")) +
    ggplot2::ylab(paste0("-log(", pvalue, ")"))

  # Return plot
  if (is.null(extra_ggplot2)) {
    p
  } else {
    p + extra_ggplot2
  }
}

get_aes <- function(status, r2, show_labels) {
  aesthetic <- aes(x = start, y = -log10(pvalue))
  if (!is.null(status)) aesthetic <- c(aesthetic, aes(shape = status))
  if (!is.null(r2)) aesthetic <- c(aesthetic, aes(fill = r2))
  if (!is.null(show_labels)) aesthetic <- c(aesthetic, aes(color = labels_color))
  if (is.list(aesthetic)) class(aesthetic) <- "uneval"
  aesthetic
}

get_geom_point <- function(status, r2) {
  if (is.null(status)) {
    if (is.null(r2)) {
      ggplot2::geom_point(shape = 21, fill = "white", size = 4)
    } else {
      ggplot2::geom_point(shape = 21, size = 4)
    }
  } else {
    if (is.null(r2)) {
      ggplot2::geom_point(fill = "white", size = 4)
    } else {
      ggplot2::geom_point(size = 4)
    }
  }
}
