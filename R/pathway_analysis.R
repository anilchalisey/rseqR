#' Title
#'
#' @param x data.frame object containing at least the column gene (gene IDs in ENSEMBL format) and 
#' a column of logFC changes (entitled log2FoldChange or logFC)
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
#' @importFrom AnnotationDbi mapIds
#' @import pathview


pathway_analysis <- function(x, species = c("human", "mouse")) {
  
  . <- q.val <- kegg.id <- pathway <- GO.id <- NULL
  
  match.arg(species, c("human", "mouse"))
  if (species == "human") {
    kid <- "hsa"
    org <- org.Hs.eg.db::org.Hs.eg.db
  }
  if (species == "mouse") {
    kid <- "mmu"
    org <- org.Mm.eg.db::org.Mm.eg.db
  }
  
  # KEGG/GO rely on ENTREZ ID, so we convert the ENSEMBL ID
  x$symbol <- AnnotationDbi::mapIds(
    org,
    keys = row.names(x),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  x$entrez <- AnnotationDbi::mapIds(
    org,
    keys = row.names(x),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  x$name <- AnnotationDbi::mapIds(
    org,
    keys = row.names(x),
    column = "GENENAME",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # convert the name of the log2FoldChange column to logFC
  colnames(x) <- gsub("log2FoldChange", "logFC", colnames(res))
  
  # Obtain a vector fold-changes where the name of the values are ENTREZ IDs
  foldchanges <- x$logFC
  names(foldchanges) <- x$entrez
  
  # KEGG pathway analysis first ---------------
  # Generate uptodata KEGG pathway gene sets
  kg <- gage::kegg.gsets(species = kid)
  
  # Subset for signalling and metabolic pathways only
  kegg.sigmet.gs <- kg$kg.sets[kg$sigmet.idx]
  
  # Get the results - here we use same.dir = TRUE to get separate lists for 
  # pathways that are upregulated and downregulated
  keggres <- gage::gage(foldchanges, gsets = kegg.sigmet.gs, same.dir = TRUE)
  
  # if plotting then do so here - we plot the top 3 upregulated and downregulated
  if (plot) {
    keggrespathways.up <- data.frame(id = rownames(keggres$greater), keggres$greater) %>%
      tbl_df() %>% 
      filter(row_number() <= 3) %>% 
      .$id %>% 
      as.character()
    keggrespathways.down <- data.frame(id = rownames(keggres$less), keggres$less) %>%
      tbl_df() %>% 
      filter(row_number() <= 3) %>% 
      .$id %>% 
      as.character()
    keggrespathways <- c(keggrespathways.up, keggrespathways.down)
    keggresids <- substr(keggrespathways, start = 1, stop = 8)
    ## Unfortunately, the code for pathview (specifically the geneannot.map() function)
    ## down't specify AnnotationDbi::select at line 52, and so this conflicts with the 
    ## dplyr namespace, which we need to unload
    if(any(grepl("package:dplyr", search()))) {
      detach("package:dplyr")
      dplyr.unloaded <- TRUE
    }
    tmp <- sapply(keggresids, function(pid) {
      pathview::pathview(gene.data = foldchanges, pathway.id = pid, species = kid)
    })
    ## Re-attach dplyr if it was unloaded
    if (dplyr.unloaded) require("dplyr", quietly = TRUE)
  }
  
  # Gene Ontology analysis; we only do the BP ---------------
  go <- gage::go.gsets(species = species)
  go.bp.gs <- go$go.sets[go$go.subs$BP]
  
  gobres <- gage::gage(foldchanges, gsets = go.bp.gs, same.dir = TRUE)
  
  
  # Write out the results -----------------------------------
  kegg.as.dataframe <- function(x, direction) {
    as.data.frame(x[[direction]]) %>%
    dplyr::mutate(kegg.id = unlist(lapply(strsplit(rownames(.), " "), function(x) x[[1]])),
                  pathway = gsub("hsa.*[[:digit:]] ", "", rownames(.))
    ) %>%
    dplyr::filter(!is.na(q.val)) %>% dplyr::select(kegg.id, pathway, dplyr::everything())
  }

  kegg.up <- kegg.as.dataframe(keggres, "greater")
  kegg.down <- kegg.as.dataframe(keggres, "less")
  
  GO.as.dataframe <- function(x, direction) {
    as.data.frame(x[[direction]]) %>%
      dplyr::mutate(GO.id = unlist(lapply(strsplit(rownames(.), " "), function(x) x[[1]])),
                    pathway = gsub("GO:.*[[:digit:]] ", "", rownames(.))
      ) %>%
      dplyr::filter(!is.na(q.val)) %>% dplyr::select(GO.id, pathway, dplyr::everything())
  }

  GO.up <- GO.as.dataframe(gobres, "greater")
  GO.down <- GO.as.dataframe(gobres, "less")
  
  towrite <- list(kegg.up, kegg.down, GO.up, GO.down)
  names(towrite) <- c("kegg.up", "kegg.down", "GO.up", "GO.down")
  create_dir("pathway_results")

  sapply(names(towrite), function(x) {
    utils::write.table(towrite[[x]],
                file = file.path("pathway_results", paste(x, "txt", sep = ".")),
                row.names = FALSE, col.names = TRUE, sep = "\t",
                quote = FALSE)
  })
  
  tomove1 <- list.files(".", pattern = paste0(kid, ".*\\.png"))
  tomove2 <- list.files(".", pattern = paste0(kid, ".*\\.xml"))
  file.rename(from = tomove1, to = file.path("pathway_results", tomove1))
  file.rename(from = tomove2, to = file.path("pathway_results", tomove2))
  return(towrite)
}
