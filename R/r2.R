#' Add r2 values to GRanges object
#'
#' This function will:
#'   1) Find the best SNPs in \code{gr_snps} using the \code{pvalue} parameter.
#'   2) Will extract relevant lines from \code{geno_ld}.
#'   3) Create \code{r2} vector.
#'
#' @param geno_ld The file *.geno.ld produced by `vcftools` or a
#'   \code{data.frame} from the file loaded with read.table.
#' @param gr_snps A \code{GRanges} object.
#' @param pvalue The name of the column to use as pvalue. Default: NULL (the
#'   first numerical metadata column will be used)
#' @param status Name of the column for the SNP status (Genotype or Imputed).
#'   If value is \code{NULL}, all the SNPs will be used to find the best SNP,
#'   otherwise only SNPs with \code{Genotyped} value will be used.
#'   Default: NULL.
#'
#' @return A \code{numeric} vector corresponding to the pairwise LD of best SNP
#'   versus all the other SNPs.
#'
#' @export
# TODO: Example
add_r2_from_vcftools <- function(geno_ld, gr_snps, pvalue = NULL,
                                 status = NULL) {
  ## 0. Name the metadata column correctly
  if (is.null(pvalue)) {
    i <- min(which(sapply(GenomicRanges::mcols(gr_snps), is.numeric)))
    pvalue <- colnames(GenomicRanges::mcols(gr_snps))[i]
  }
  pvalue <- GenomicRanges::mcols(gr_snps)[[pvalue]]
  if (!is.null(status)) {
    status <- GenomicRanges::mcols(gr_snps)[[status]]
  }

  ## 1. Find the best SNP
  gr_sorted <- gr_snps
  gr_sorted <- gr_sorted[order(pvalue)]
  if (!is.null(status)) {
    gr_sorted <- gr_snps[status == "Genotyped"]
  }
  # TODO: Should use findOverlaps instead of expecting a name metadata column
  best_snp <- as.character(gr_sorted[1]$name)
  i <- which(gr_snps$name == best_snp)

  ## 2. Load the R2 values
  if (is.character(geno_ld)) {
    geno_ld <- read.table(geno_ld, header = TRUE, stringsAsFactors=FALSE)
  }

  ## 3. Subset R2 values
  # The current SNP can be in POS1 or POS2 column
  pos <- GenomicRanges::start(gr_snps)[i]
  subset_R2 <- geno_ld[geno_ld$POS1 == pos | geno_ld$POS2 == pos,]

  # We to make it easier to convert to GRanges, we will make sure all POS1
  # correspond to the best SNP position and all POS2 value correspond to the 2nd
  # SNP
  j <- subset_R2[["POS1"]] != pos
  tmp <- subset_R2[["POS1"]][j]
  subset_R2[["POS1"]][j] <- subset_R2[["POS2"]][j]
  subset_R2[["POS2"]][j] <- tmp

  # Convert to GRanges
  subset_R2 <- GenomicRanges::GRanges(unique(GenomicRanges::seqnames(gr_snps)),
    IRanges::IRanges(subset_R2[["POS2"]], width = 1), R2 = subset_R2[["R.2"]])

  ## 4. Add r2 value to gr_snps
  overlaps <- GenomicRanges::findOverlaps(gr_snps, subset_R2)
  j <- IRanges::queryHits(overlaps)
  k <- IRanges::subjectHits(overlaps)
  r2 <- 1
  r2[j] <- subset_R2$R2[k]

  ## 5. Return result
  r2
}
