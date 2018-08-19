#' Wrapper function to run tximport
#'
#' \code{run_tximport} function runs tximport on transcript level
#' abundances from Salmon to summarise to gene level. See Bioconductor
#' package \link[tximport]{tximport} for details.
#'
#' We use Ensembl annotation for both genes and transcripts.
#'
#' @inheritParams run_dea
#' @param salmon.files character vector specifying path to quant.sf files.
#' @param summarise logical, if \code{TRUE} then return gene-level summary, otherwise
#' return transcript level data.
#' @param countsFromAbundance whether to generate counts based on abundance.
#' Available options are: \code{'no'},
#' \code{'scaledTPM'} (abundance based estimated counts scaled up to
#' library size),
#' \code{'lengthScaledTPM'}.
#'
#' @return a list of gene and transcript level estimated counts.
#'
#' @references
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015):
#' Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences. F1000Research.
#' \url{http://dx.doi.org/10.12688/f1000research.7563.1}
#'
#' @importFrom data.table fread
#' @importFrom dplyr select distinct
#' @importFrom tximport tximport summarizeToGene
#'
#' @export

run_tximport <- function(species = c("human", "mouse"),
                         salmon.files, summarise = TRUE,
                         countsFromAbundance=c("no", "scaledTPM", "lengthScaledTPM")) {

  countsFromAbundance <- match.arg(countsFromAbundance,
                                   c("no", "scaledTPM", "lengthScaledTPM"))

  species <- match.arg(species, c("human", "mouse"))

  entrez <- ensgene <- enstxp <- NULL

  if (species == "human") {
    tx2gene <- grch38 %>%
      dplyr::select(-entrez) %>%
      dplyr::distinct() %>%
      dplyr::select(enstxp, ensgene)
    }

  if (species == "mouse") {
    tx2gene <- grcm38 %>%
      dplyr::select(-entrez) %>%
      dplyr::distinct() %>%
      dplyr::select(enstxp, ensgene)
    }

  names(salmon.files) <- dirname(salmon.files) %>% reduce_path()

  give_note("\nGenerating counts table...\n\n")

  tx.t <- tximport::tximport(files = salmon.files, type = "salmon", tx2gene = tx2gene,
                             txOut = TRUE, dropInfReps = FALSE,
                             importer = data.table::fread,
                             countsFromAbundance = countsFromAbundance)

  if (all(apply(is.na(tx.t$counts), 2, any))) {
    txi.t <- tximport::tximport(files = salmon.files, type="salmon", tx2gene = tx2gene,
                                txOut=TRUE, importer = data.table::fread,
                                countsFromAbundance="no")
    countsFromAbundance <- "no"
  } else {
    txi.t <- tx.t
  }

  if (summarise) {
    txi.g <- tximport::summarizeToGene(txi = txi.t, ignoreTxVersion = TRUE,
                                       tx2gene = tx2gene)
    txi.g$countsFromAbundance <- countsFromAbundance
    return(txi.g)
  } else {
    return(txi.t)
  }
}
