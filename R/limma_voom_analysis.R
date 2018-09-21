#' Run Limma-voom analysis on RNA-seq count data
#'
#' @inheritParams edger_analysis
#'
#' @importFrom edgeR cpm DGEList calcNormFactors
#' @importFrom limma voom makeContrasts lmFit contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom utils combn write.table
#' @importFrom purrr map
#' @importFrom dplyr as_tibble arrange desc
#'
#' @export

limma_voom_analysis <- function(metadata, species = c("human", "mouse")) {
  
  . <- padj <- logFC <- NULL
  
  design <- stats::model.matrix(metadata$design, data = metadata$sampleinfo)
  counts <- metadata$counts$B$counts
  group <- metadata$sampleinfo$condition
  minLibSize <- min(colSums(counts))
  minGroupSize <- min(tabulate(group))

  keep <- rowSums(edgeR::cpm(counts) > 10/(minLibSize/1e6)) >= minGroupSize
  counts <- counts[keep, ]

  dge <- edgeR::DGEList(counts = counts,
                        group = group)
  dge <- edgeR::calcNormFactors(dge)

  v <- limma::voomWithQualityWeights(dge, design = design, normalize.method = "none")
  normalisedCounts <- ens2symbol(result = v$E,
                                 species = species,
                                 columns.of.interest = c("gene", colnames(v$E)),
                                 colnames = c("gene", colnames(v$E), "symbol"))
  create_dir(file.path(metadata$outdir, "Limma-voom"))
  out <- file.path(metadata$outdir, "Limma-voom")
  utils::write.table(normalisedCounts, file.path(out, "Limma-voom_normalisedcounts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  utils::write.table(normalisedCounts[, c("gene", "symbol")], file.path(out, "Limma-voom_consideredgenes"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  make_PCA(v$E, "Limma-voom", metadata$sampleinfo)
  tomove <- list.files(".", pattern = "Limma-voom_PCA")
  move_file(tomove = tomove, to = out)

  contrast.matrix <-
    colnames(design) %>%
    .[grep("condition", .)] %>%
    factor(., levels = .) %>%
    utils::combn(., 2, simplify = FALSE) %>%
    purrr::map(~paste0(.[2], "-", .[1])) %>%
    limma::makeContrasts(contrasts = ., levels = colnames(design))

  vfit <- limma::lmFit(v)
  vfit <- limma::contrasts.fit(vfit, contrasts = contrast.matrix)
  fit <- limma::eBayes(vfit)

  res.names <- gsub("condition", "", colnames(fit))
  out.sub <- file.path(out, res.names)
  lapply(out.sub, create_dir)

  DEG <- list()
  for (i in seq_len(ncol(fit))) {
    degTable <- limma::topTable(fit, number = Inf, coef = i)
    output <- ens2symbol(result = degTable[order(degTable$adj.P.Val), ],
                         species = species,
                         columns.of.interest = c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val"),
                         colnames = c("gene", "logFC", "AveExpr", "pvalue", "padj", "symbol"))
    make_MA(output, fdr = 0.05, label.rectangle = TRUE, proc = "Limma-voom", top = 20, select.method = "logfc")
    make_volcano(input = output, proc = "Limma-voom")
    move_file(tomove = c("Limma-voom_volcano.png", "Limma-voom_MAplot.png"), out.sub[[i]])
    DEG[[i]] <- subset(output, padj < 0.05)
    utils::write.table(
      x = as.data.frame(output),
      file = file.path(out.sub[[i]], "Limma-voom_differential_expression.txt"),
      sep = "\t",
      row.names = TRUE,
      quote = FALSE)

    utils::write.table(
      x = as.data.frame(DEG[[i]][, c("gene", "symbol")]),
      file = file.path(out.sub[[i]], "Limma-voom_DEG.txt"),
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE)
  }

  names(DEG) <- res.names

  lapply(names(DEG), function(x) {
    return(
      give_note(paste("\nFound", nrow(subset(DEG[[x]], logFC > 0)), "upregulated genes and",
                      nrow(subset(DEG[[x]], logFC < 0)), "downregulated genes for the contrast", x,
                      "using Limma-voom.\n\n", collapse=" "))
    )
  })

  return(lapply(DEG, function (x) dplyr::as_tibble(x) %>% dplyr::arrange(dplyr::desc(abs(logFC)), padj)))
}

