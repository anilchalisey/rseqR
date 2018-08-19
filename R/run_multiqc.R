#' Generate MultiQC report for FastQC
#'
#' \code{run_multiqc} generates a single HTML report from the FastQC results using MultiQC.
#'
#' @section MultiQC path:
#' If the executable is in \code{$PATH}, then the default value for paths
#' (\code{"multiqc"}) will work. If it is not in \code{$PATH}, then the absolute
#' path should be given. If using Windows 10, it is assumed that MultiQC has
#' been installed in WSL, and the same rules apply.
#'
#' @param fastqc.dir directory where all the FastQC files are saved.
#' @param dest.dir directory to save the combined QC report.
#' @inheritParams run_dea
#'
#' @return HTML report.
#'
#' @references
#' Philip Ewels, Mans Magnusson, Sverker Lundin, and Max Kaller (2016):
#' MultiQC: summarize analysis results for multiple tools and samples
#' in a single report. Bioinformatics, 32(19), 3047-3048.
#' \url{https://doi.org/10.1093/bioinformatics/btw354}
#'
#' @examples
#'
#' \dontrun{
#' run_multiqc(fastqc.dir=tempdir(), dest.dir=tempdir())
#' }
#'
#' @export

run_multiqc <- function(fastqc.dir, dest.dir, multiqc = "multiqc") {
  give_note("\nCreating MultiQC report.\n\n")
  file.copy(system.file("extdata/MULTIQC", "multiqc_config.yaml", package = "rseqR"), ".")
  file.copy(system.file("extdata/MULTIQC/Logo.png", package = "rseqR"), ".")
  multiqc.cmd <- sprintf('%s -f %s -c multiqc_config.yaml', multiqc, fastqc.dir)
  run_cmd(multiqc.cmd)
  remove_file("Logo.png")
  remove_file("multiqc_config.yaml")
  dir.create(dest.dir, recursive = TRUE, showWarnings = FALSE)
  file.rename("multiqc_report.html", file.path(dest.dir, "multiqc_fastqc.html"))
  file.rename("multiqc_data", file.path(dest.dir, "multiqc_data"))
}
