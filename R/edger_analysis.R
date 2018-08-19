#' Run edgeR analysis on RNA-seq count data
#'
#' @param metadata list created by \code{sanity_check} and with raw
#' count data stored in \code{metadata$count$A} and length-adjusted
#' count data stored in \code{metadata$count$B}.  Created as part of
#' \code{run_dea}
#' @inheritParams run_dea
#'
#' @importFrom edgeR cpm calcNormFactors DGEList scaleOffset
#' @importFrom edgeR estimateDisp glmQLFit topTags
#' @importFrom utils combn
#' @importFrom purrr map
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix
#' @importFrom dplyr as_tibble arrange
#'
#' @export

edger_analysis <- function(metadata, species = c("human", "mouse")) {
  counts <- metadata$counts$A$counts
  group <- metadata$sampleinfo$condition
  minLibSize <- min(colSums(counts))
  minGroupSize <- min(tabulate(group))

  keep <- rowSums(edgeR::cpm(counts) > 10/(minLibSize/1e6)) >= minGroupSize
  counts <- counts[keep, ]

  normmat <- metadata$counts$A$length
  normmat <- normmat[keep, ]
  normmat <- normmat/exp(rowMeans(log(normmat)))
  o <- log(edgeR::calcNormFactors(counts/normmat)) + log(colSums(counts/normmat))
  d <- edgeR::DGEList(counts, group = metadata$sampleinfo$condition)
  d <- edgeR::scaleOffset(d, t(t(log(normmat)) + o))

  design <- stats::model.matrix(metadata$design, data = metadata$sampleinfo)
  disp <- edgeR::estimateDisp(d, design = design, robust = TRUE)
  rownames(design) <- colnames(disp)
  create_dir(file.path(metadata$outdir, "edgeR"))
  out <- file.path(metadata$outdir, "edgeR")
  exploratory_analysis_edger(disp, species = species, metadata = metadata)
  tomove <- list.files(".", pattern = "edgeR_")
  move_file(tomove, out)

  fit <- edgeR::glmQLFit(disp, design = design, robust = TRUE)

  contrast.matrix <-
    colnames(design) %>%
    .[grep("condition", .)] %>%
    factor(., levels = .) %>%
    utils::combn(., 2, simplify = FALSE) %>%
    purrr::map(~paste0(.[2], "-", .[1])) %>%
    limma::makeContrasts(contrasts = ., levels = colnames(design))

  deg <- apply(contrast.matrix, 2, function(x) {
    edgeR::glmQLFTest(fit, contrast = x)
  })

  names(deg) <- gsub("condition", "", names(deg))
  out.sub <- file.path(out, names(deg))
  lapply(out.sub, create_dir)

  DEG <- list()
  for (i in seq_along(deg)) {
    ma.data <- edgeR::topTags(deg[[i]], n = nrow(deg[[i]]$counts))$table
    output <- ens2symbol(
      result = ma.data,
      species = species,
      columns.of.interest = c('gene', 'logFC', 'logCPM', 'PValue', 'FDR'),
      colnames = c("gene", "logFC", "logCPM", "pvalue", "FDR", "symbol"))

    make_MA(output, fdr = 0.05, label.rectangle = TRUE, proc = "edgeR", top = 20, select.method = "logfc")
    DEG[[i]] <- subset(output, FDR < 0.05)

    write.table(
      x = ma.data,
      file = file.path(out.sub[[i]], "edgeR_differential_expression.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE)

    write.table(
      x = as.data.frame(DEG[[i]][, c("gene", "symbol")]),
      file = file.path(out.sub[[i]], "edgeR_DEG.txt"),
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE)

    make_volcano(input = output, proc = "edgeR")
    move_file(tomove = c("edgeR_volcano.png", "edgeR_MAplot.png"), out.sub[[i]])
  }

  names(DEG) <- names(deg)
  lapply(names(DEG), function(x) {
    give_note(paste("\nFound", nrow(subset(DEG[[x]], logFC > 0)), "upregulated genes and",
                      nrow(subset(DEG[[x]], logFC < 0)), "downregulated genes for the contrast", x,
                      "using edgeR-glmQLFTest.\n\n", collapse=" "))
  })

  return(lapply(DEG, function (x) dplyr::as_tibble(x) %>% dplyr::arrange(desc(abs(logFC)), FDR)))
  }

