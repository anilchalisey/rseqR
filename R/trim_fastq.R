#' An R-based wrapper for Trim Galore!
#'
#' @description Run the Trim Galore! tool
#'
#' @details This script runs the Trim Galore! tool and requires installation
#' of both Cutadapt and Trim Galore!  It is essential that Cutadapt is in the
#' executable path otherwise this tool will not work.
#'
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end.
#' @param adapter1 a character string specifying the adapter sequence to be
#' trimmed. If not specified explicitly, Trim Galore will try to auto-detect
#' whether the Illumina universal, Nextera transposase or Illumina small RNA
#' adapter sequence was used. Also see \code{illumina}, \code{nextera} and
#' \code{small_rna} options. If no adapter can be detected within the first 1
#' million sequences of the first file specified Trim Galore defaults to
#' \code{illumina}.
#' @param adapter2 a character string specifying an optional adapter sequence to
#' be trimmed off read 2 of paired-end files. This option requires paired-end
#' reads.
#' @param illumina a logical specifying that the adapter sequence to be trimmed
#' is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of
#' the default auto-detection of adapter sequence.  Default: \code{FALSE}
#' @param nextera adapter sequence to be trimmed is the first 12bp of the
#' Nextera adapter CTGTCTCTTATA instead of the default auto-detection of adapter
#' sequence.
#' @param small_rna a logical specifying that the adapter sequence to be trimmed
#' is the first 12bp of the Illumina Small RNA 3' Adapter TGGAATTCTCGG instead
#' of the default auto-detection of adapter sequence.  Selecting to trim
#' smallRNA adapters will also lower the \code{length} value to 18bp. If the
#' smallRNA libraries are paired-end then \code{adapter2} will be set to the
#' Illumina small RNA 5' adapter automatically (GATCGTCGGACT) unless
#' \code{adapter2} had been defined explicitly.
#' @param minlength an integer value; reads that become shorter than this length
#' as a result of either quality or adapter trimming are discarded. A value of 0
#' effectively disables this behaviour.  Default: 20 bp.  For paired-end files,
#' both reads of a read-pair need to be longer than bp to be printed out to
#' validated paired-end files. If only one read became too short there is the
#' possibility of keeping such unpaired single-end reads (see
#' \code{retain_unpaired}). Default pair-cutoff: 20 bp.
#' @param minqual an integer value specifying the quality threshold below which
#' to trim low-quality ends from reads in addition to adapter removal. Default
#' Phred score: 20.
#' @param trimN a logical specifying whether to remove Ns from the end of reads.
#' @param retainUnpaired a logical.  If only one of the two paired-end reads
#' become too short, the longer read will be written to either .unpaired_1.fq or
#' .unpaired_2.fq output files. The length cutoff for unpaired single-end reads
#' is governed by the parameters \code{retain1length} and \code{retain2length}.
#' Default: ON.
#' @param retain1length an integer.  Unpaired single-end read length cutoff
#' needed for read 1 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp.
#' @param retain2length an integer.  Unpaired single-end read length cutoff
#' needed for read 2 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp
#' @param clipR1 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 1 (or single-end reads). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clipR2 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 2 (paired-end reads only). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clip3primeR1 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param clip3primeR2 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param robust_check a logical indicating whether to check that the paired
#' files specified are matching and have equal numbers of reads.  Default:
#' \code{FALSE}
#' @param trimgalore a character string specifying the path to the trimgalore executable.
#' On Unix systems, if the executable is in \code{$PATH}, then it may be left as
#' the default. If it is not in \code{$PATH}, then the absolute path should be given.
#` If using the WSL on Windows 10, then the path must be the absolute path in WSL,
#' unless the system has been set up as described in the vignette.
#' @param dest.dir a character string specifying the output directory.  If NULL
#' a directory named "TRIMMED_FASTQC" is created in the current working directory
#' [DEFAULT = NULL].
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#'
#'@param multi.core an integer value indicating the number of cores to be used by cutadapt 
#'multicore functionality requires pigz and python3 installed and accessible.
#'
#'@param cutadapt a character string specifing the path to the cutadapt executable,
#' if the executable is in \code{$PATH} then it may be left as default.
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

trim_fastq <- function(fastq1, fastq2 = NULL, adapter1 = NULL, adapter2 = NULL,
                       illumina = FALSE, nextera = FALSE, small_rna = FALSE,
                       minlength = 20, minqual = 20, trimN = TRUE,
                       retainUnpaired = TRUE, retain1length = 35,
                       retain2length = 35, clipR1 = NULL, clipR2  = NULL,
                       clip3primeR1 = NULL, clip3primeR2 = NULL,
                       robust_check = FALSE, dest.dir = NULL,
                       threads = NULL,
                       multi.core = NULL,
                       cutadapt = NULL,
                       trimgalore = "trim_galore") {

cmd <- paste0(trimgalore)

  if(!is.null(cutadapt)){
    cmd <- paste0(cmd, " --path_to_cutadapt ", cutadapt)
  }
  
if(!is.null(multi.core)){
  cmd <- paste0(cmd, " --cores ", multi.core)
}
  
  #cmd <- paste0(trimgalore, " --path_to_cutadapt /Users/ieo/.local/bin/cutadapt", " --cores 4")



  if (is.null(fastq2)) {
    paired = FALSE
  } else {
    if(length(fastq1) != length(fastq2)) {stop("The number of forward and reverse reads do not match")}
    if (robust_check) {
      fq1lengths <- lapply(fastq1, function (x) {sprintf("gzip -cd %s | wc -l", x)})
      fq2lengths <- lapply(fastq1, function (x) {sprintf("gzip -cd %s | wc -l", x)})
      fq1lengths <- unlist(lapply(fq1lengths, run_cmd, intern = TRUE))
      fq2lengths <- unlist(lapply(fq2lengths, run_cmd, intern = TRUE))
      if (!identical(fq1lengths, fq2lengths)) {
        stop("One or more of the forward and reverse reads pairs have differing number of reads.\n",
             "Are you sure the two lists are in the correct paired order?")
      }
    }
    paired = TRUE
  }

  if (is.null(dest.dir)) dest.dir <- "TRIMMED_FASTQC"
  dir.create(dest.dir, showWarnings = FALSE)

  cmd <- paste(cmd,
               "-q", minqual, "--length", minlength, "-o", dest.dir)
  if(!is.null(adapter1)) cmd <- paste(cmd, "--adapter", adapter1)
  if(!is.null(adapter2)) cmd <- paste(cmd, "--adapter2", adapter2)
  if(illumina) cmd <- paste(cmd, "--illumina")
  if(nextera) cmd <- paste(cmd, "--nextera")
  if(small_rna) cmd <- paste(cmd, "--small_rna")
  if(trimN) cmd <- paste(cmd, "--trim-n")
  if(!is.null(clipR1)) cmd <- paste(cmd, "--clip_R1", clipR1)
  if(!is.null(clipR2)) cmd <- paste(cmd, "--clip_R2", clipR2)
  if(!is.null(clip3primeR1)) cmd <- paste(cmd, "--three_prime_clip_R1", clip3primeR1)
  if(!is.null(clip3primeR2)) cmd <- paste(cmd, "--three_prime_clip_R2", clip3primeR2)
  if(paired) {
    cmd <- paste(cmd, "--paired")
    if(retainUnpaired) {
      cmd <- paste(cmd, "--retain_unpaired",
                   "-r1", retain1length,
                   "-r2", retain2length)
    }
  }

  give_note("\nRemoving adapters and performing quality trimming...\n\n")

  if (is.null(threads)) {threads <- parallel::detectCores() - 1}

  if (threads > 1) {
    cl <- parallel::makeCluster(threads)
    doParallel::registerDoParallel(cl)

    if (paired) {
      foreach (i = seq_along(fastq1)) %dopar% {
        tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
        run_cmd <- function(cmd, intern = FALSE) {
          if (.Platform$OS.type != "windows") {
            system(command = cmd, intern = intern)
          } else {
            shell(cmd = shQuote(cmd), shell = "zsh", intern = intern)
          }
        }
        run_cmd(tgcmd)
      }
    } else {
      foreach (i = seq_along(fastq1)) %dopar% {
        tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
        run_cmd <- function(cmd, intern = FALSE) {
          if (.Platform$OS.type != "windows") {
            system(command = cmd, intern = intern)
          } else {
            shell(cmd = shQuote(cmd), shell = "zsh", intern = intern)
          }
        }
        run_cmd(tgcmd)
      }

    }
    parallel::stopCluster(cl)
  } else {
    if (paired) {
      for (i in seq_along(fastq1)) {
        tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
        run_cmd(tgcmd)
      }
    } else {
      for (i in seq_along(fastq1)) {
        tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
        run_cmd(tgcmd)
      }
    }
  }

  if (paired) {
    trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
    lapply(trimmed.files, function(x) {
     file.rename(from = x, to = gsub("_val_[0-9]", "_trimmed", x))
    })
    unpaired.files <- list.files(path = dest.dir, pattern = "*unpaired", full.names = TRUE)
    lapply(unpaired.files, function(x) {
      file.rename(from = x, to = gsub("_unpaired_[0-9]", "_unpaired", x))
    })
  } else {
    trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
    lapply(trimmed.files, function(x) {
      file.rename(from = x, to = gsub("_val_[0-9]", "_trimmed", x))
    })
  }
}
