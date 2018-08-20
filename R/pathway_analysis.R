#' Title
#'
#' @param x output from \code{edger_analysis()}, \code{deseq2_analysis},
#' \code{limma_voom_analysis}.
#'
#' @inheritParams run_dea
#'
#' @export
#'
#' @importFrom gage kegg.gsets go.gsets gage
#' @import GO.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom dplyr mutate filter select
#' @importFrom utils write.table


pathway_analysis <- function(x, species = c("human", "mouse")) {
  
  . <- q.val <- kegg.id <- pathway <- GO.id <- NULL
  
  match.arg(species, c("human", "mouse"))
  if (species == "human") {kid <- "hsa"}
  if (species == "mouse") {kid <- "mmu"}

  kg <- gage::kegg.gsets(species = kid)
  kegg.sigmet.gs <- kg$kg.sets[kg$sigmet.idx]
  kegg.dise.gs <- kg$kg.sets[kg$dise.idx]

  go <- gage::go.gsets(species = species)
  go.bp.gs <- go$go.sets[go$go.subs$BP]
  go.mf.gs <- go$go.sets[go$go.subs$MF]
  go.cc.gs <- go$go.sets[go$go.subs$CC]

  res <- as.data.frame(x)
  rownames(res) <- res$gene

  colnames(res) <- gsub("log2FoldChange", "logFC", colnames(res))
  colnames(res) <- gsub("AveExpr", "logFC", colnames(res))

  res <- ens2entrez(
    result = res,
    species = "human",
    columns.of.interest = c("gene", "logFC", "padj"),
    colnames = c("gene", "logFC", "padj", "entrez")
  )

  fc.res <- res$logFC
  names(fc.res) <- res$entrez

  fc.kegg.sigmet.p <- gage::gage(fc.res, gsets = kegg.sigmet.gs)
  fc.kegg.dise.p <- gage::gage(fc.res, gsets = kegg.dise.gs)
  fc.go.bp.p <- gage::gage(fc.res, gsets = go.bp.gs)
  fc.go.mf.p <- gage::gage(fc.res, gsets = go.mf.gs)
  fc.go.cc.p <- gage::gage(fc.res, gsets = go.cc.gs)

  kegg.as.dataframe <- function(x, direction) {
    as.data.frame(x[[direction]]) %>%
    dplyr::mutate(kegg.id = unlist(lapply(strsplit(rownames(.), " "), function(x) x[[1]])),
                  pathway = gsub("hsa.*[[:digit:]] ", "", rownames(.))
    ) %>%
    dplyr::filter(!is.na(q.val)) %>% dplyr::select(kegg.id, pathway, dplyr::everything())
  }

  fc.kegg.sigmet.p.up <- kegg.as.dataframe(fc.kegg.sigmet.p, "greater")
  fc.kegg.sigmet.p.down <- kegg.as.dataframe(fc.kegg.sigmet.p, "less")
  fc.kegg.dise.p.up <- kegg.as.dataframe(fc.kegg.dise.p, "greater")
  fc.kegg.dise.p.up <- kegg.as.dataframe(fc.kegg.dise.p, "less")

  GO.as.dataframe <- function(x, direction) {
    as.data.frame(x[[direction]]) %>%
      dplyr::mutate(GO.id = unlist(lapply(strsplit(rownames(.), " "), function(x) x[[1]])),
                    pathway = gsub("GO:.*[[:digit:]] ", "", rownames(.))
      ) %>%
      dplyr::filter(!is.na(q.val)) %>% dplyr::select(GO.id, pathway, dplyr::everything())
  }

  fc.go.bp.p.up <- GO.as.dataframe(fc.go.bp.p, "greater")
  fc.go.mf.p.up <- GO.as.dataframe(fc.go.mf.p, "greater")
  fc.go.cc.p.up <- GO.as.dataframe(fc.go.cc.p, "greater")
  fc.go.bp.p.down <- GO.as.dataframe(fc.go.bp.p, "less")
  fc.go.mf.p.down <- GO.as.dataframe(fc.go.mf.p, "less")
  fc.go.cc.p.down <- GO.as.dataframe(fc.go.cc.p, "less")

  towrite = ls(pattern = "fc.go.*up|fc.go.*down|fc.kegg.*up|fc.kegg.*down")
  names(towrite) <- towrite
  create_dir("pathway_results")

  sapply(names(towrite), function(x) {
    utils::write.table(towrite[[x]],
                file = file.path("pathway_results", paste(x, "txt", sep = ".")),
                row.names = FALSE, col.names = TRUE, sep = "\t",
                quote = FALSE)
  })
  return(towrite)
}
