#' Plot intersections of list of DE genes
#'
#' @param limma data.frame with the columns gene and padj or FDR; result of Limma-voom analysis
#' @param edger data.frame with the columns gene and padj or FDR; result of edgeR analysis
#' @param deseq2 data.frame with the columns gene and padj or FDR; result of DESeq2 analysis
#'
#' @return An UpSet plot of gene intersections.
#'
#' @importFrom UpSetR upset
#' @importFrom dplyr select
#' @importFrom grDevices png dev.off
#' @export

make_upset <- function(limma, edger, deseq2) {

  dedata <- list(limma, edger, deseq2)

  dedata <- lapply(dedata, function(x) {
    colnames(x) <- gsub("FDR", "padj", colnames(x))
    return(x)
  })

  dedata <- dedata %>% lapply(., function(x) dplyr::select(x, gene, padj))

  # merge all DE genes
  dedata <- Reduce(function(x,y) {
    merge(x, y, by = "gene", all.x = TRUE, all.y = TRUE)
  }, dedata)
  dedata <- unique(dedata)

  # get in format for plot
  rownames(dedata) <- dedata$gene
  dedata <- dedata[2:length(dedata)]
  dedata[dedata < 0.05] <- 1
  dedata[dedata < 1] <- 0
  dedata[is.na(dedata)] <- 0
  colnames(dedata) <- c("Limma-voom", "edgeR", "DESeq2")

  UpSetR::upset(dedata,
                sets = colnames(dedata),
                sets.bar.color = "#ED6925FF",
                main.bar.color = "#781C6DFF",
                order.by = "freq",
                empty.intersections = "on",
                sets.x.label = "Set Size",
                mainbar.y.label = "Intersection Size")
}
