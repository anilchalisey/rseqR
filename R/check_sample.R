#' Check that targets file contains the necessary info in the right format
#'
#' @inheritParams run_dea
#'
#' @importFrom utils read.table

check_sample <- function(sample.info, reference = NULL) {
  if (!file.exists(sample.info)) {
    give_error("SAMPLE.INFO ERROR:\nsample information file does not exist or path is incorrect.")
  }
  sample.info <- utils::read.table(sample.info, header = TRUE, stringsAsFactors = FALSE)
  names(sample.info) <- tolower(names(sample.info))

  for (required in c("condition", "sample", "file1")) {
    if (!required %in% names(sample.info)) {
      give_error(paste0("SAMPLE.INFO ERROR:\nCould not find required field '",
                         required, "' in sample information file."))
    }
  }

  if (min(table(sample.info$condition)) < 2) {
    give_error("SAMPLE.INFO ERROR:\nLess than 2 replicates in smallest group from <condition>.")
  }

  if (length(unique(sample.info$condition)) < 2) {
    give_error("SAMPLE.INFO ERROR:\nfield <condition> needs at least two different values/groups.")
  }

  for (unique_values_required_field in c("sample", "file1")) {
    if (anyDuplicated(sample.info[, unique_values_required_field])) {
      give_error(paste("SAMPLE.INFO ERROR:\nValues in field",
                       unique_values_required_field,
                       "in sample info file are not unique."))}
  }

  for (path in sample.info$file1) {
    if (!file.exists(path)) {
      give_error(paste("\nIncorrect path to", path, "\nFile not found.\n"))
    }
  }

  if (all(grepl("fastq|fq", sample.info$file1))) {
    give_note("Files are of type fastq")
    file.type <- "fastq"
  } else {
    if (any(grepl("fastq|fq", sample.info$file1))) {
      give_error("\nMix of file types - unable to proceed\n")
    } else {
      if (all(grepl(".sf", sample.info$file1))) {
        give_note("\nFiles are of type salmon quantification\n")
        file.type <- "salmon"
      } else {
        if (any(grepl(".sf", sample.info$file1))) {
          give_error("Mix of file types - unable to proceed\n")
        } else {
          give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq> or <.sf>\n")
        }
      }
    }
  }

  if (file.type == "fastq") {
    if ("file2" %in% names(sample.info)) {
      if (!all(grepl("fastq|fq", sample.info$file2))) {
        give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq>\n")
      }
      for (path in sample.info$file2) {
        if (!file.exists(path)) {
          give_error(paste0("\nIncorrect path to ", path, ". File not found.\n"))
        }
      }
      give_note("\nEach sample has 2 files - assuming PE reads\n")
      paired <- TRUE
    } else {
      give_note("\nOnly one file per sample - assuming SE reads\n")
      paired <- FALSE
    }
  }

  if (is.null(reference)) {
    sample.info$condition <- factor(sample.info$condition)
  } else {
      if (length(reference) == 1) {
      sample.info$condition <- relevel(factor(sample.info$condition), ref = reference)
    } else {
      if (length(reference) > 1) {
        if (length(setdiff(sample.info$condition, reference)) > 0) {
          give_error("ERROR:\n The factors specified in <reference> are those in the condition column do not match\n")
        }
        if (length(setdiff(reference, sample.info$condition)) > 0) {
          give_error("ERROR:\n The factors specified in <reference> are those in the condition column do not match\n")
        }
        sample.info$condition <- factor(sample.info$condition, levels = reference)
      }
    }
  }

  if ("batch" %in% names(sample.info)) sample.info$batch <- as.factor(sample.info$batch)

  if (file.type == "fastq") {
    return(list(sample.info = sample.info, file.type = file.type, paired = paired))
  } else {
    return(list(sample.info = sample.info, file.type = file.type))
  }

}
