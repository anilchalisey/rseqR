#' Import the metadata and check it is correct
#'
#' @inheritParams run_dea
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
#' @return A \code{data.frame}

sanity_check <- function(sample.info, reference = NULL,
                         species = c("human", "mouse"),
                         output.dir, threads = NULL) {
  inputdata <- list()

  # Results directory
  if (!dir.exists(output.dir)) create_dir(output.dir)
  inputdata$outdir <- output.dir

  # Threads
  if (is.null(threads)) threads <- parallel::detectCores() - 1
  inputdata$threads <- threads

  # Annotation
  species <- match.arg(species, c("human", "mouse"))
  inputdata$annotation <- species

  # Check sample info
  sd <- check_sample(sample.info, reference = reference)
  inputdata$sampleinfo <- sd$sample.info
  inputdata$filetype <- sd$file.type

  if (inputdata$filetype == "fastq") inputdata$paired <- sd$paired

  # Make design
  inputdata$design <- make_design(inputdata$sampleinfo)

  return(inputdata)
}
