#' Create Volcano Plot
#'
#' @param input result of differential analysis created as part of
#' \code{edger_analysis}, \code{limma_voom_analysis}, or \code{deseq2_analysis}
#' @inheritParams make_MA
#'
#' @import ggplot2

make_volcano <- function(input, proc = c("Limma-voom", "edgeR", "DESeq2")) {

  proc <- match.arg(proc)

  if (proc == "Limma-voom") {
    y <- "padj"
    x <- "logFC"
  }

  if (proc == "edgeR") {
    y = "FDR"
    x = "logFC"
  }

  if (proc == "DESeq2") {
    y = "padj"
    x = "log2FoldChange"
  }

  toptable <- as.data.frame(input)
  toptable$Sig <- "NS"
  toptable$Sig[abs(toptable[, x]) > 2] <- "FC"
  toptable$Sig[toptable[, y] < 0.05] <- "P"
  toptable$Sig[abs(toptable[, x]) > 2 & (toptable[, y] < 0.05)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels=c("NS", "FC", "P", "FC_P"))

  if (min(toptable[, y], na.rm=TRUE) == 0) {
    warning("One or more P values is 0. Converting to minimum possible value...", call. = FALSE)
    toptable[which(toptable[, x] == 0), x] <- .Machine$double.xmin
  }

  if (substr(toptable$gene[1], 1, 4) == 'ENSG') {
    toptable$lab <- toptable$symbol
  } else {
    toptable$lab <- toptable$gene
  }

  toptable$xval <- toptable[, x]
  toptable$yval <- toptable[, y]

  plot <- ggplot2::ggplot(toptable, ggplot2::aes(x = xval, y = -log10(yval))) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(Sig)), alpha = 0.5) +
    ggplot2::scale_colour_manual(values = c(NS = "black", FC = "blue",
                                            P = "green", FC_P = "red"),
                                 labels = c(NS = "NS", FC = "|FC| > 2",
                                            P = "FDR < 0.05", FC_P = "|FC| > 2 & FDR < 0.05"),
                                 name = "") +
    theme_publication() +
    ggplot2::xlab(bquote(~Log[2]~ "fold change")) +
    ggplot2::ylab(bquote(~-Log[10]~italic(P))) +
    ggplot2::xlim(c(min(toptable[,x], na.rm=TRUE),
                    max(toptable[,x], na.rm=TRUE))) +
    ggplot2::ylim(c(0, max(-log10(toptable[,y]), na.rm=TRUE))) +
    ggplot2::geom_vline(xintercept = c(-2, 2),
                        linetype = "longdash",
                        colour = "black",
                        size = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        linetype = "longdash",
                        colour = "black",
                        size = 0.4) +
    ggplot2::geom_text(data = subset(toptable, toptable[, y] < 0.05),
                       ggplot2::aes(label = subset(toptable, toptable[, y] < 0.05)[, "lab"]),
                       size = 2.0, check_overlap = TRUE, vjust = 1.0) +
    ggplot2::theme(legend.position = "top", legend.key.size = ggplot2::unit(0.5, "cm"))

  suppressMessages(ggplot2::ggsave(paste(proc, 'volcano.png', sep="_"), plot, device = "png"))
}

