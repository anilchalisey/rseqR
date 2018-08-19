#' Run differential expression analysis
#'
#' @param sample.info character string giving the path to a tab-delimited text
#' file with at least the columns <condition> (treatment condition),
#' <sample> (sample name), and <file1> (absolute or relative path to the
#' fastq or salmon quant.sf files).  If fastq files and PE reads, then a column
#' <file2> should also be present.  If a batch effect is to be included in the
#' design, then this should be identified under the column <batch>.
#' @param reference character vector specifying the conditions in order.  For example,
#' c("A", "B", "C", "D") would mean "A" is the reference condition to which "B", "C"
#' and "D" are compared; in addition, "C" and "D" will be compared to "B", and "D"
#' will be compared to "C". If \code{NULL} then the comparisons will be arranged
#' alphabetically.  [DEFAULT = NULL].
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#' @param output.dir character string specifying the directory to which results
#' will be saved.  If the directory does not exist, it will be created.
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#' @param fastqc a character string specifying the path to the fastqc executable.
#' [DEFAULT = "fastqc"].
#' @param multiqc a character string specifying the path to the multiqc executable.
#' [DEFAULT = "multiqc"].
#' @param index.dir directory of the index files needed for read mapping using Salmon.
#' See function \code{'build_index()'}.
#' @param advanced.opts character vector supplying list of advanced option
#' arguments to apply to each Salmon call. For details see Salmon documentation
#' or type \code{salmon quant --help-reads} at the command line.
#' @param bam logical, if \code{TRUE} then create a pseudo-alignment BAM file.
#' [Default = \code{FALSE}]
#' @param bootstraps integer giving the number of bootstrap samples
#' that Salmon should use (default is 0). With bootstrap samples, uncertainty
#' in abundance can be quantified.
#' @param seqBias logical, should Salmon's option be used to model and correct
#' abundances for sequence specific bias? Default is \code{TRUE}.
#' @param gcBias logical, should Salmon's option be used to model and correct
#' abundances for GC content bias? Requires Salmon version 0.7.2 or higher.
#' Default is \code{TRUE}.
#' @param posBias logical, should Salmon's option be used to model and correct
#' abundances for positional biases? Requires Salmon version 0.7.3 or higher.
#' Default is \code{FALSE}.
#' @param allowOrphans logical, if \code{TRUE} then consider orphaned reads as
#' valid hits when performing lightweight-alignment. This option will increase
#' sensitivity (allow more reads to map and more transcripts to be detected), but may
#' decrease specificity as orphaned alignments are more likely to be spurious.
#' For more details see Salmon documentation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' targets.file <- data.frame(
#'   file1 = list.files(system.file("extdata/FASTQ/", package = "rseqR"),
#'                      "*_1.fastq", full.names = TRUE),
#'   file2 = list.files(system.file("extdata/FASTQ/", package = "rseqR"),
#'                      "*_2.fastq", full.names = TRUE),
#'   sample = c("BUFF1", "BUFF2", "OX1", "OX2", "UNT1", "UNT2"),
#'   condition = c(rep("BUFF", 2), rep("OX", 2), rep("UNT", 2)),
#'   batch = rep(1:2, times = 3))
#' write.table(targets.file, "targets.txt", col.names = TRUE,
#'   row.names = FALSE, sep = "\t", quote = FALSE)
#' run_dea(sample.info = "targets.txt", reference = c("UNT", "BUFF", "OX"),
#'         species = "human", output.dir = "results", index.dir = NULL)
#' }

run_dea <- function(
  # Important options
  sample.info, reference = NULL, species = c("human", "mouse"),
  output.dir, threads = NULL,
  # QC options
  fastqc = "fastqc", multiqc = "multiqc",
  # Salmon options
  index.dir = NULL, salmon = "salmon", bam = FALSE,
  bootstraps = 0,  seqBias = TRUE, gcBias = TRUE,
  posBias = FALSE, allowOrphans = FALSE, advanced.opts = NULL) {

  metadata <- sanity_check(sample.info = sample.info, reference = reference,
                           species = species, output.dir = output.dir, threads = threads)

  if (metadata$filetype == "fastq") {
    # Run FastQC
    run_fastqc(fastqc.files = c(metadata$sampleinfo$file1, metadata$sampleinfo$file2),
               dest.dir = file.path(output.dir, "fastqc"), threads = threads, fastqc = fastqc)

    # Run multiQC
    run_multiqc(fastqc.dir = file.path(output.dir, "fastqc"), dest.dir = file.path(output.dir, "multiqc"),
                multiqc = multiqc)

    # Create salmon index
    if (is.null(index.dir)) {
      give_warning(paste("\nIndex directory not specified - creating index for", species,
                         "in", file.path(output.dir, "index, based on Ensembl 93.\n\n")))
      create_dir(file.path(output.dir, "index"))
      build_index(species = species, kmer = 31, ens_release = 93, dest.dir = file.path(output.dir, "index"))
      index.dir <- file.path(output.dir, "index", paste0(species, "transcripts_release93_index"))
    }

    # Run salmon
    index.dir <- convert_paths(index.dir)
    if (metadata$paired) {
      run_salmon(fastq1 = metadata$sampleinfo$file1, fastq2 = metadata$sampleinfo$file2,
                 index.dir = index.dir, dest.dir = file.path(output.dir, "salmon"),
                 salmon = salmon, threads = threads, bam = bam, bootstraps = bootstraps,
                 advanced.opts = advanced.opts, seqBias = seqBias, gcBias = gcBias,
                 posBias = posBias, allowOrphans = allowOrphans)
    } else {
      run_salmon(fastq1 = metadata$sampleinfo$file1, index.dir = index.dir,
                 dest.dir = file.path(output.dir, "salmon"), salmon = salmon, threads = threads,
                 bam = bam, bootstraps = bootstraps, advanced.opts = advanced.opts,
                 seqBias = seqBias, gcBias = gcBias, posBias = posBias,
                 allowOrphans = allowOrphans)
    }


    # Import salmon data

    salmon.files <- list.files(file.path(metadata$outdir, "salmon"), pattern = "quant.sf$",
                               recursive = TRUE, full.names = TRUE)
    metadata$counts$A <- run_tximport(species = species, salmon.files = salmon.files,
                                      summarise = TRUE, countsFromAbundance = "no")
    metadata$counts$B <- run_tximport(species = species, salmon.files = salmon.files,
                                      summarise = TRUE, countsFromAbundance = "lengthScaledTPM")

  }

  if (metadata$filetype == "salmon") {
    salmon.files <- metadata$sampleinfo$file1
    metadata$counts$A <- run_tximport(species = species, salmon.files = salmon.files,
                                      summarise = TRUE, countsFromAbundance = "no")
    metadata$counts$B <- run_tximport(species = species, salmon.files = salmon.files,
                                      summarise = TRUE, countsFromAbundance = "lengthScaledTPM")

  }

  give_note("\nPerforming differential expression analysis with Limma-voom.\n\n")
  DE_limma <- limma_voom_analysis(metadata = metadata, species = species)

  give_note("\nPerforming differential expression analysis with edgeR.\n\n")
  DE_edger <- edger_analysis(metadata = metadata, species = species)

  give_note("\nPerforming differential expression analysis with DESeq2.\n\n")
  DE_deseq2 <- deseq2_analysis(metadata = metadata, species = species)

  for (i in seq_along(DE_limma)) {
    if (sum(nrow(DE_limma[[i]]), nrow(DE_edger[[i]]), nrow(DE_deseq2[[i]])) <= 3) {
      give_warning(paste("No Upset Plot generated for the contrast", names(DE_limma)[[i]], "as no DEGs.\n\n"))
    } else {
      grDevices::png(
        filename = paste0("upset_", names(DE_limma)[[i]], ".png"),
        width = 8,
        height = 8,
        units = "in",
        res = 500)
      make_upset(DE_limma[[i]], DE_edger[[i]], DE_deseq2[[i]])
      grDevices::dev.off()
    }
  }
  move_file(list.files(pattern = "upset_"), output.dir)
  DE_genes <- list(limma = DE_limma, edger = DE_edger, deseq2 = DE_deseq2)
  save(metadata, file = file.path(output.dir, "metadata.rda"))
  save(DE_deseq2, DE_edger, DE_limma, file = file.path(output.dir, "DE_results.rda"))
  return(DE_genes)
}
