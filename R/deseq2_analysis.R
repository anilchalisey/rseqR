#' Run DESeq2 analysis on RNA-seq count data
#'
#' @inheritParams edger_analysis
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom DESeq2 DESeqDataSetFromTximport counts DESeq results
#' @importFrom edgeR cpm
#' @importFrom utils combn write.table
#' @importFrom grDevices png dev.off
#' @importFrom dplyr as_tibble arrange desc
#'
#' @export

deseq2_analysis <- function(metadata, species) {
  
  log2FoldChange <- padj <- NULL

  threads <- metadata$threads
  doParallel::registerDoParallel(threads)

  deseqdata <- DESeq2::DESeqDataSetFromTximport(
    txi = metadata$counts$A,
    colData = metadata$sampleinfo,
    design = metadata$design
  )

  minLibSize <- min(colSums(DESeq2::counts(deseqdata)))
  minGroupSize <- min(tabulate(deseqdata$condition))
  dds <- deseqdata[rowSums(edgeR::cpm(DESeq2::counts(deseqdata)) > 10/(minLibSize/1e6)) >= minGroupSize, ]

  dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = TRUE)

  create_dir(file.path(metadata$outdir, "DESeq2"))
  out <- file.path(metadata$outdir, "DESeq2")
  exploratory_analysis_deseq2(dds, species = species, metadata = metadata)
  tomove <- list.files(".", pattern = "DESeq2_")
  move_file(tomove, out)

  contrasts <- factor(levels(dds$condition), levels = levels(metadata$sampleinfo$condition))
  con <- utils::combn(contrasts, 2)
  vec <- list()
  for (i in 1:ncol(con)) {
    vec[[i]] <- levels(con[, i, drop = TRUE])
  }
  contrasts <- lapply(vec, function(x) c("condition", rev(x)))
  res.names <- lapply(contrasts, function(x) paste(x[2], x[3], sep = "-"))
  out.sub <- file.path(out, res.names)
  lapply(out.sub, create_dir)

  DEG <- list()
  for (i in seq_along(contrasts)) {
    res <- DESeq2::results(dds, parallel = TRUE,
                           contrast = contrasts[[i]], independentFiltering = FALSE)
    output <- ens2symbol(
      result = res[order(res$padj), ],
      species = "human",
      columns.of.interest = c("gene", "baseMean", "log2FoldChange", "pvalue", "padj"),
      colnames = c("gene", "baseMean", "log2FoldChange", "pvalue", "padj", "symbol")
    )
    make_MA(output, fdr = 0.05, label.rectangle = TRUE,
            proc = "DESeq2", top = 20, select.method = "logfc")
    make_volcano(input = output, proc = "DESeq2")
    move_file(tomove = c("DESeq2_volcano.png", "DESeq2_MAplot.png"), out.sub[[i]])
    DEG[[i]] <- subset(output, padj < 0.05)
    utils::write.table(
      x = as.data.frame(output),
      file = file.path(out.sub[[i]], "DESeq2_differential_expression.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE)

    utils::write.table(
      x = as.data.frame(DEG[[i]][, c("gene", "symbol")]),
      file = file.path(out.sub[[i]], "DESeq2_DEG.txt"),
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE)
  }

  names(DEG) <- res.names

  lapply(names(DEG), function(x) {
    return(
      give_note(paste("\nFound", nrow(subset(DEG[[x]], log2FoldChange > 0)), "upregulated genes and",
                      nrow(subset(DEG[[x]], log2FoldChange < 0)), "downregulated genes for the contrast", x,
                      "using DESeq2.\n\n", collapse=" "))
    )
  })

  for (i in seq_along(DEG)) {
    if (nrow(DEG[[i]]) <= 2) {
      give_warning(paste("No heatmap generated for the contrast", names(DEG)[[i]], "as no DEGs.\n\n"))
    } else {
      grDevices::png(
        filename = file.path(out.sub[[i]], paste0("DESeq2_heatmap_", names(DEG)[[i]], ".png")),
        width = 9,
        height = 9,
        units = "in",
        res = 600)
      make_heatmap(dds, DEG[[i]])
      grDevices::dev.off()
    }
  }
  dev.off()

  return(lapply(DEG, function (x) dplyr::as_tibble(x) %>%
                  dplyr::arrange(dplyr::desc(abs(log2FoldChange)), padj)))
}
