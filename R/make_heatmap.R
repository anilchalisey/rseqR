#' Create heatmap with clustering of samples and genes for DEGs
#'
#' @inheritParams exploratory_analysis_deseq2
#' @param DEG list of differential genes produced by \code{deseq2_analysis}
#'
#' @importFrom DESeq2 vst
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr mutate select
#' @importFrom magrittr set_rownames
#' @importFrom stats hclust as.dendrogram
#' @importFrom ggdendro segment theme_dendro dendro_data
#' @importFrom viridis scale_fill_viridis
#' @importFrom grid unit.pmax gpar grid.rect grid.draw
#' @importFrom gridExtra arrangeGrob
#' @importFrom gtable gtable_add_cols
#' @import ggplot2

make_heatmap <- function(dds, DEG) {
  deseq.vst <- DESeq2::vst(dds) %>%
    SummarizedExperiment::assay() %>%
    as.data.frame() %>%
    dplyr::mutate(gene = rownames(.))
  deseq.vst <- deseq.vst[deseq.vst$gene %in% DEG$gene, ]

  deseq.vst.df <- reshape2::melt(deseq.vst, id.vars = "gene")

  deseq.vst.matrix <- deseq.vst %>%
    magrittr::set_rownames(.$gene) %>%
    dplyr::select(-gene) %>%
    as.matrix()

  distance.gene <- stats::dist(deseq.vst.matrix)
  distance.sample <- stats::dist(t(deseq.vst.matrix))

  cluster.gene <- stats::hclust(distance.gene, method = "average")
  cluster.sample <- stats::hclust(distance.sample, method = "average")

  sample.model <- stats::as.dendrogram(cluster.sample)
  sample.dendrogram.data <- ggdendro::segment(ggdendro::dendro_data(sample.model, type = "rectangle"))
  sample.dendrogram <- ggplot2::ggplot(sample.dendrogram.data) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggdendro::theme_dendro() +
    ggplot2::scale_x_continuous(expand = c(.0085, .0085)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))

  deseq.vst.df$variable <- factor(deseq.vst.df$variable,
                                  levels = cluster.sample$labels[cluster.sample$order])

  gene.model <- stats::as.dendrogram(cluster.gene)
  gene.dendrogram.data <- ggdendro::segment(ggdendro::dendro_data(gene.model, type = "rectangle"))
  gene.dendrogram <- ggplot2::ggplot(gene.dendrogram.data) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand=c(0, 0)) +
    ggplot2::scale_x_continuous(expand=c(0, 0)) + ggdendro::theme_dendro()

  deseq.vst.df$gene <- factor(deseq.vst.df$gene,
                              levels = cluster.gene$labels[cluster.gene$order])

  heatmap <- ggplot2::ggplot(deseq.vst.df, ggplot2::aes(x = variable, y = gene, fill = value)) +
    ggplot2::geom_raster() +
    viridis::scale_fill_viridis(option = "viridis", trans = "sqrt", alpha = 0.8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::xlab("") + ggplot2::ylab("")


  sample.dendrogram.grob <- ggplot2::ggplotGrob(sample.dendrogram)
  gene.dendrogram.grob <- ggplot2::ggplotGrob(gene.dendrogram +
                                                ggplot2::scale_x_discrete(expand = c(0, 0)))
  heatmap.grob <- ggplot2::ggplotGrob(heatmap)

  sample.dendrogram.grob <- gtable::gtable_add_cols(sample.dendrogram.grob,
                                                    heatmap.grob$widths[7], 6)
  sample.dendrogram.grob <- gtable::gtable_add_cols(sample.dendrogram.grob,
                                                    heatmap.grob$widths[8], 7)

  con.data <- as.data.frame(SummarizedExperiment::colData(dds))
  con.data$sample <- rownames(con.data)
  con.data$sample <- factor(con.data$sample,
                            levels = cluster.sample$labels[cluster.sample$order])
  n <- length(levels(con.data$condition))
  colours <- RColorBrewer::brewer.pal(n, "Set2")
  sample.condition <- ggplot2::ggplot(con.data,
                                      ggplot2::aes(x = sample, y = 1, fill = condition)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(name = "condition", values = colours) +
    ggplot2::theme_void() +
    ggplot2::guides(col = ggplot2::guide_legend(nrow = 2))

  sample.condition.grob <- ggplot2::ggplotGrob(sample.condition)
  maxWidth <- grid::unit.pmax(sample.dendrogram.grob$widths,
                              heatmap.grob$widths,
                              sample.condition.grob$widths)
  sample.dendrogram.grob$widths <- as.list(maxWidth)
  heatmap.grob$widths <- as.list(maxWidth)
  sample.condition.grob$widths <- as.list(maxWidth)

  maxHeight <- grid::unit.pmax(heatmap.grob$heights, gene.dendrogram.grob$heights)
  heatmap.grob$heights <- as.list(maxHeight)
  gene.dendrogram.grob$heights <- as.list(maxHeight)

  # Arrange the grobs into a plot
  blankpanel <- grid::grid.rect(gp = grid::gpar(col = "white"))
  finalGrob <- gridExtra::arrangeGrob(blankpanel, sample.dendrogram.grob,
                                      blankpanel,
                                      sample.condition.grob, gene.dendrogram.grob,
                                      heatmap.grob, ncol = 2, nrow = 3,
                                      widths = c(1, 5), heights = c(2, 0.8, 6))

  grid::grid.draw(finalGrob)

}
