#' Exploratory analysis of edgeR data
#'
#' @param disp DGEList object; output of \code{edgeR::estimateDisp}
#' @inheritParams edger_analysis
#' @inheritParams run_dea
#'
#' @importFrom grDevices png dev.off
#' @importFrom edgeR cpm plotBCV

exploratory_analysis_edger <- function(disp, species, metadata) {
  grDevices::png(
    filename = "edgeR_dispersionplot.png",
    width = 8,
    height = 8,
    units = "in",
    res = 500)
  edgeR::plotBCV(disp)
  grDevices::dev.off()

  normalised.counts <- edgeR::cpm(disp)
  normalised.counts_names <- ens2symbol(
    result = normalised.counts,
    species = species,
    columns.of.interest = c("gene", colnames(normalised.counts)),
    colnames = c("gene", colnames(normalised.counts), "symbol")
  )

  make_PCA(normalised.counts, "edgeR", metadata$sampleinfo)
  write.table(
    x = normalised.counts_names,
    file = "edgeR_normalisedcounts_cpm.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    x = normalised.counts_names[, c("gene", "symbol")],
    file = "edgeR_consideredgenes.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
