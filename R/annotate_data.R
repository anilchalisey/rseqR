#' Annotate gene data
#'
#' \code{annotate_data} annotates Ensembl, Entrez or Gene Symbol IDs.
#'
#' @param data matrix or data.frame object whereby rownames are gene IDs
#' @param id.type gene ID type used for rownames.  May be one of \code{ensgene},
#'   \code{symbol}, or \code{entrez}.
#' @param species  name of the species. Only \code{'human'}, \code{'mouse'},
#' and \code{'rat'} are allowed.
#'
#' @export
#'
#' @importFrom dplyr as_tibble select

annotate_data <- function(data, id.type = c("ensgene", "symbol", "entrez"),
                          species = c("human", "mouse")) {

  id.type <- match.arg(id.type)
  species <- match.arg(species)

  if (species == "human") {
    annot_genes <- grch38[match(rownames(data), grch38[[id.type]]), , drop = FALSE]
  }

  if (species == "mouse") {
    annot_genes <- grch38[match(rownames(data), grcm38[[id.type]]), , drop = FALSE]
  }

  enstxp <- NULL

  dplyr::as_tibble(cbind(dplyr::select(annot_genes, id.type), data,
                         dplyr::select(annot_genes, -id.type, -enstxp)))

}
