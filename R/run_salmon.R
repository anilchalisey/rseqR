#' Quantify transcript abundances using Salmon
#'
#' @description Run the abundance quantification tool \code{Salmon} on a set of FASTQ
#' files. Requires \code{Salmon} (\url{https://combine-lab.github.io/salmon/})
#' to be installed and a Salmon transcript index must have been generated prior
#' to using this function. See the Salmon website for installation and basic
#' usage instructions.
#'
#' @return The following items will be returned and saved in the salmon directory:
#' \enumerate{
#' \item quant.sf: plain-text, tab-separated quantification file
#' that contains 5 column: Name,Length,EffectiveLength,TPM, and NumReads.
#' \item quant.sf.bkp: plain-text, tab-separated quantification file
#' that contains 5 column: Name,Length,EffectiveLength,TPM, and NumReads.
#' This is the raw version of the quant.sf file.
#' \item cmd_info.json: A JSON format file that records the main command
#' line parameters with which Salmon was invoked for the run that produced
#' the output in this directory.
#' \item aux_info: This directory will have a number of files (and
#' subfolders) depending on how salmon was invoked.
#' \item meta_info.json: A JSON file that contains meta information about
#' the run, including stats such as the number of observed and mapped
#' fragments, details of the bias modeling etc.
#' \item ambig_info.tsv: This file contains information about the number
#' of uniquely-mapping reads as well as the total number of
#' ambiguously-mapping reads for each transcript.
#' \item lib_format_counts.json: This JSON file reports the number of
#' fragments that had at least one mapping compatible with the designated
#' library format, as well as the number that didn't.
#' \item libParams: The auxiliary directory will contain a text file
#' called flenDist.txt. This file contains an approximation of the
#' observed fragment length distribution.
#' }
#'
#' @references
#' Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and
#' Carl Kingsford (2017): Salmon provides fast and bias-aware
#' quantification of transcript expression. Nature methods, 14(4), 417.
#' \url{https://www.nature.com/articles/nmeth.4197}
#'
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end.
#' @param dest.dir directory where results are to be saved.  If directory does not exist,
#' then it will be created.
#' @param salmon (optional) string giving full command to use to call
#' Salmon, if simply typing "salmon" at the command line does not give the
#' required version of Salmon or does not work. Default is simply "salmon".
#' If used, this argument should give the full path to the desired Salmon
#' binary.
#' @inheritParams run_dea
#'
#' @export

run_salmon <- function(fastq1, fastq2 = NULL, index.dir, dest.dir = "SALMON",
                       salmon = "salmon", threads = NULL, advanced.opts = NULL,
                       bam = FALSE, bootstraps = 0, seqBias = TRUE, gcBias = TRUE,
                       posBias = FALSE, allowOrphans = FALSE) {

  if (is.null(fastq2)) {
    paired = FALSE
  } else {
    if(length(fastq1) != length(fastq2)) {
      stop("The number of forward and reverse reads do not match")
    }
    paired = TRUE
  }

  # Checks
  if (check_cmd(salmon) == "Not Found") {
    stop("the salmon path is invalid or it is not installed correctly.")
  }
  if (!all(file.exists(fastq1))) {
    stop("All of the specified files not found. Is the file path correct?")
  }
  if (paired) {
    if (!all(file.exists(fastq2))) {
      stop("All of the specified files not found. Is the file path correct?")
    }
  }

  # Threads
  if (is.null(threads)) threads <- parallel::detectCores() - 1

  # Labels
  sample.names <- reduce_path(fastq1)

  # Output directory
  if(!dir.exists(dest.dir)) {
    create_dir(dest.dir)
  }
  output.dirs <- file.path(dest.dir, sample.names)

  # Convert paths if necessary
  fastq1 <- convert_paths(fastq1)
  fastq2 <- convert_paths(fastq2)

  ## Generate calls to Salmon
  salmon.args <- paste("quant -i", index.dir, "-l A -o", output.dirs,
                       "--threads", threads,
                       "--numBootstraps", bootstraps)
  names(salmon.args) <- sample.names

  if (seqBias)  salmon.args <- paste0(salmon.args, " --seqBias")
  if (gcBias)   salmon.args <- paste0(salmon.args, " --gcBias")
  if (posBias)  salmon.args <- paste0(salmon.args, " --posBias")

  if (!paired) {
    salmon.args <- paste(salmon.args, "-r", fastq1)
  } else {
    salmon.args <- paste(salmon.args, "-1", fastq1,
                         "-2", fastq2)
  }

  if (allowOrphans) salmon.args <- paste(salmon.args, "--allowOrphans")
  if (!is.null(advanced.opts)) salmon.args <- paste(salmon.args, advanced.opts)

  if (bam) {
    salmon.args <- paste0(salmon, " ", salmon.args, " --writeMappings=",
                          sample.names, "_pseudo.bam")
  } else {
    salmon.args <- paste0(salmon, " ", salmon.args)
  }

  give_note(paste("Analysis started: ", Sys.time(), "\n\n"))
  give_note(paste("Processing", length(sample.names), "samples", "\n\n"))

  for(i in seq_along(salmon.args)) {
    run_cmd(cmd = salmon.args[i], intern = FALSE)
  }

  give_note(paste("Analysis completed: ", Sys.time(), "\n\n"))

  # Clean up the quant.sf files to remove the decimal characters in the transcript names
  # Save the original quant.sf as a backup and then create a clean one
  quant.files <- file.path(output.dirs, "quant.sf")
  lapply(quant.files, function(x) {
    file.rename(x, paste0(x, ".bkp"))
  })
  quant.files.bkp <- paste0(quant.files, ".bkp")
  clean_cmd <- sprintf("cat %s | sed -E 's/\\.[0-9]+//' > %s", quant.files.bkp, quant.files)
  lapply(clean_cmd, function(x) run_cmd(x))
}
