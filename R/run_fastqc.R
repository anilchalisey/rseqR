#' QC reports for fastq files using FastQC
#'
#' @description Run the FastQC tool to generate HTML reports. FastQC may be installed from
#' \url{https://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
#'
#' @section FastQC path:
#' If the executable is in \code{$PATH}, then the default value for paths
#' (\code{"fastqc"}) will work. If it is not in \code{$PATH}, then the absolute
#' path should be given. If using Windows 10, it is assumed that FASTQC has
#' been installed in WSL, and the same rules apply.
#'
#' @param fastqc.files character vector of paths to fastq files. [DEFAULT = NULL].
#' @param dest.dir character string indicating directory to save results. [DEFAULT = "FASTQC"]
#' @inheritParams run_dea
#'
#' @return \code{character} string specifying the path of the output directory.
#'
#' @examples
#' \dontrun{
#' fqFiles <- list.files(path = "fastqDir", pattern = "*.fastq", full.names = TRUE)
#' run_fastqc(fqFiles)
#' }
#'
#' @importFrom parallel detectCores
#'
#' @export

run_fastqc <- function(fastqc.files = NULL,
                       dest.dir = "FASTQC",
                       threads  = NULL,
                       fastqc   = "fastqc") {

  if (check_cmd(fastqc) == "Not Found") {
    give_error("the fastqc path is invalid or it is not installed correctly.")
  }

  if (is.null(fastqc.files)) {
    give_error("no files specified")
  }

  filesExist <- sapply(fastqc.files, file.exists)
  if (!all(filesExist)) {
    give_error("not all files specified can be found - are all paths correct?")
  }

  fastqc.files <- convert_paths(fastqc.files)

  fq <- paste(fastqc.files, collapse = " ")

  dir.create(path = dest.dir, showWarnings = FALSE)

  if (is.null(threads)) {threads <- parallel::detectCores() - 1}

  cat(paste("Running FastQC... ", Sys.time(), "\n", sep=""))
  fastqc.run <- sprintf('%s %s --outdir=%s --threads=%s',
                        fastqc, fq, dest.dir, threads)
  run_cmd(fastqc.run)
  return(dest.dir)
}
