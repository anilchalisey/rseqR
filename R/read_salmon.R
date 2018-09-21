#' Read the run info from Salmon output
#' 
#' Function to import metadata from Salmon into R
#'
#' @author Anil Chalisey
#'
#' @param data.dir Path to directories containing the Salmon results
#'
#' @importFrom rjson fromJSON
#' @importFrom dplyr bind_rows
#' @importFrom lubridate mdy
#' @importFrom magrittr %>%

read_salmon <- function(data.dir = NULL) {
  directories <- data.dir
  ## extract run info
  json.file <- lapply(directories, function(x) file.path(x, "aux_info", "meta_info.json"))
  jnames <- c("salmon_version", "frag_dist_length", "seq_bias_correct",
              "gc_bias_correct", "num_bias_bins", "mapping_type",
              "num_targets", "num_processed", "num_mapped",
              "percent_mapped", "start_time")
  if (!all(sapply(json.file, file.exists))) {
    stop("all meta_info.json files not found or do not exist.")
  }
  
  tmp.run <- lapply(json.file, function(x) {
    tmp <- rjson::fromJSON(file = x)
    tmp <- tmp[names(tmp) %in% jnames]
    return(tmp)
  })
  
  run.info <- dplyr::bind_rows(tmp.run) %>%
    data.frame(., stringsAsFactors = FALSE)
  run.info$date <- lubridate::mdy(gsub("[[:digit:]]+\\:*[[:digit:]]+\\:[[:digit:]]*", "",
                                       run.info$start_time))
  run.info$start_time <- unlist(regmatches(run.info$start_time,
                                           gregexpr("[[:digit:]]+\\:*[[:digit:]]+\\:[[:digit:]]*",
                                                    run.info$start_time)))
  
  ## extract library format info
  libformat.json.file <- file.path(directories, "lib_format_counts.json")
  lnames <- c("read_files", "expected_format", "compatible_fragment_ratio",
              "num_compatible_fragments", "num_assigned_fragments",
              "num_consistent_mappings", "num_inconsistent_mappings",
              "strand_mapping_bias")
  if (!all(sapply(libformat.json.file, file.exists))) {
    stop("all lib_format_counts.json files not found or do not exist.")
  }
  tmp.lib <- lapply(libformat.json.file, function(x) {
    tmp <- rjson::fromJSON(file = x)
    tmp <- tmp[names(tmp) %in% lnames]
    return(tmp)
  })
  
  lib.info <- dplyr::bind_rows(tmp.lib) %>%
    data.frame(., stringsAsFactors = FALSE)
  
  ## extract command info
  cmd.json.file <- file.path(directories, "cmd_info.json")
  cnames <- c("index", "libType", "mates1", "mates2", "output", "auxDir")
  if (!all(sapply(cmd.json.file, file.exists))) {
    stop("all cmd_info.json files not found or do not exist.")
  }
  
  tmp.cmd <- lapply(cmd.json.file, function(x) {
    tmp <- rjson::fromJSON(file = x)
    tmp <- tmp[names(tmp) %in% cnames]
    return(tmp)
  })
  
  cmd.info <- dplyr::bind_rows(tmp.cmd) %>%
    data.frame(., stringsAsFactors = FALSE)
  
  output <- list(run.info = run.info, libformat.info = lib.info, cmd.info = cmd.info)
  output <- lapply(output, function(x) {
    rownames(x) <- targets$sample.names
    return(x)
  })
  
  # return the dataframes
  return(output)
}
