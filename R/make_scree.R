#' Create scree plot from prcomp data
#'
#' @param pca \code{prcomp} object.
#' @inheritParams make_MA
#'
#' @import ggplot2

make_scree <- function(pca, proc) {
  scree_plot <- data.frame(pca$sdev^2/sum(pca$sdev^2))
  scree_plot[1:10,2] <- c(1:10)
  colnames(scree_plot) <- c("variance", "component number")
  scree <- ggplot2::ggplot(scree_plot[1:10,],
                           mapping = ggplot2::aes(x = `component number`, y = variance)) +
    ggplot2::geom_bar(stat = "identity") + theme_publication()
  suppressMessages(ggplot2::ggsave(paste(proc, 'PCA_scree.png', sep="_"), scree, device = "png"))
}
