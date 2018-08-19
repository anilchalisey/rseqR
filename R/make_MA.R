#' Create MA plot
#'
#' Create a scatter plot of log2 fold changes vs the mean expression signal.
#'
#' @param data \code{matrix} or \code{data.frame} with the columns for average
#' gene expression, adjusted P-value, and log fold change; created as part of
#' \code{edger_analysis}, \code{limma_voom_analysis}, or \code{deseq2_analysis}.
#' @param fdr numerical value between 0 and 1 indicating the threshold false
#' discovery rate for discovering differentially expressed genes. [DEFAULT = 0.05]
#' @param label.rectangle logical indicating whether to add rectangles underneath
#' gene names [DEFAULT = FALSE]
#' @param proc character string indicating analysis procedure used.  One of
#' "Limma-voom", "edgeR", "DESeq2". [DEFAULT = "Limma-voom].
#' @param top integer indicating the number of genes to label
#' @param select.method Character string identifying the method of selecting the
#' top genes.  One of \code{padj} or \code{logfc}
#' @param main character string giving the plot main title
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr mutate case_when filter slice
#' @importFrom stats na.omit
#' @importFrom ggrepel geom_label_repel geom_text_repel

make_MA <- function(data, fdr = 0.05, label.rectangle = FALSE,
                    proc = c("Limma-voom", "edgeR", "DESeq2"),
                    top = 0, select.method = c("padj", "logfc"),
                    main = NULL) {

  mean <- lfc <- sig <- gene <- symbol <- padj <- NULL

  proc <- match.arg(proc, c("Limma-voom", "edgeR", "DESeq2"))

  if (proc == "DESeq2")
    res.ma <- data.frame(gene = data$gene, symbol = data$symbol, mean = data$baseMean,
                         lfc = data$log2FoldChange, padj = data$padj) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < fdr & lfc > 0 ~ 1,
                                         padj < fdr & lfc < 0 ~ -1,
                                         TRUE ~ 0)) %>%
    dplyr::mutate(sig = as.factor(sig))

  if (proc == "edgeR")
    res.ma <- data.frame(gene = data$gene, symbol = data$symbol, mean = data$logCPM,
                         lfc = data$logFC, padj = data$FDR) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < fdr & lfc > 0 ~ 1,
                                         padj < fdr & lfc < 0 ~ -1,
                                         TRUE ~ 0)) %>%
    dplyr::mutate(sig = as.factor(sig))

  if (proc == "Limma-voom")
    res.ma <- data.frame(gene = data$gene, symbol = data$symbol, mean = data$AveExpr,
                         lfc = data$logFC, padj = data$padj) %>%
    dplyr::mutate(sig = dplyr::case_when(padj < fdr & lfc > 0 ~ 1,
                                         padj < fdr & lfc < 0 ~ -1,
                                         TRUE ~ 0)) %>%
    dplyr::mutate(sig = as.factor(sig))

  select.method <- match.arg(select.method)
  if (select.method == "padj") {
    res.ma <- res.ma[order(res.ma$padj), ]
  } else {
    res.ma <- res.ma[order(abs(res.ma$lfc), decreasing = TRUE), ]
  }
  labelled <- res.ma %>%
    stats::na.omit() %>%
    dplyr::filter(padj <= 0.05 & gene != "") %>%
    dplyr::slice(1:top)

  set.seed(42)

  if (proc == "DESeq2") p <- ggplot2::ggplot(res.ma, ggplot2::aes(x = log2(mean), y = lfc))
  if (proc == "edgeR") p <- ggplot2::ggplot(res.ma, ggplot2::aes(x = mean, y = lfc))
  if (proc == "Limma-voom") p <- ggplot2::ggplot(res.ma, ggplot2::aes(x = mean, y = lfc))

  p <- p + ggplot2::geom_point(ggplot2::aes(color = sig), alpha = 0.4, size = 1) +
    ggplot2::scale_color_manual(values = c("-1" = "red", "0" = "darkgrey", "1" = "blue"),
                                name = "",
                                breaks = c("1", "0", "-1"),
                                labels = c(paste("Up:", table(res.ma$sig)[3]),
                                           paste("NS:", table(res.ma$sig)[2]),
                                           paste("Down:", table(res.ma$sig)[1]))) +
    ggplot2::labs(x = expression(paste(log[2], " mean expression")),
                  y = expression(paste(log[2], " fold change")),
                  title = main) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 2, colour = "black", size = 0.8) +
    theme_publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", size = 10),
                   axis.text.x = ggplot2::element_text(face = "bold", size = 8),
                   axis.title.y = ggplot2::element_text(face = "bold", size = 10),
                   axis.text.y = ggplot2::element_text(face = "bold", size = 8),
                   legend.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.text = ggplot2::element_text(size = 9))

  if (proc == "DESeq2") p <- p + ggplot2::scale_x_continuous(breaks = seq(0, max(log2(res.ma$mean)), 2))
  if (proc == "edgeR") p <- p + ggplot2::scale_x_continuous(breaks = seq(0, max(res.ma$mean), 2))
  if (proc == "Limma-voom") p <- p + ggplot2::scale_x_continuous(breaks = seq(0, max(res.ma$mean), 2))

  if (top != 0) {
    if (label.rectangle) {
      p <- p + ggrepel::geom_label_repel(data = labelled,
                                         mapping = ggplot2::aes(label = symbol),
                                         box.padding = ggplot2::unit(0.35, "lines"),
                                         point.padding = ggplot2::unit(0.3, "lines"),
                                         force = 1,
                                         size = 3,
                                         color = "black")
    } else {
      p <- p + ggrepel::geom_text_repel(data = labelled,
                                        mapping = ggplot2::aes(label = symbol),
                                        box.padding = ggplot2::unit(0.35, "lines"),
                                        point.padding = ggplot2::unit(0.3, "lines"),
                                        force = 1,
                                        size = 3,
                                        color = "black")
    }
  }
  suppressMessages(ggplot2::ggsave(paste(proc, 'MAplot.png', sep="_"), p, device = "png"))
}
