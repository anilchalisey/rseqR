#' Exploratory analysis of DESeq2 data
#'
#' @param dds DESeqDataSet object created by \code{DESeq2::DESeq}
#' @inheritParams exploratory_analysis_edger
#'
#' @importFrom grDevices png dev.off
#' @importFrom DESeq2 plotDispEsts vst
#' @importFrom SummarizedExperiment assay

exploratory_analysis_deseq2 <- function(dds, species, metadata) {
  grDevices::png(
    filename = "DESeq2_dispersionplot.png",
    width = 8,
    height = 8,
    units = "in",
    res = 500)
  DESeq2::plotDispEsts(dds, genecol = "grey",
                       fitcol = "purple",
                       finalcol = "orange")
  grDevices::dev.off()

  vst <- DESeq2::vst(dds)
  normcounts <- SummarizedExperiment::assay(vst)
  rlddf <- data.frame(normcounts)
  rownames(rlddf) <- names(dds)
  colnames(rlddf) <- dds$sample

  rlddf <- ens2symbol(
    result = rlddf,
    species = species,
    columns.of.interest = c("gene", colnames(rlddf)),
    colnames = c("gene", colnames(rlddf), "symbol")
  )

  make_PCA(normcounts, "DESeq2", metadata$sampleinfo)

  write.table(
    x = rlddf,
    file = "DESeq2_vst_normalisedcounts.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    x = rlddf[, c("gene", "symbol")],
    file = "DESeq2_consideredgenes.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )

}
