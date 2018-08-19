# Internal functions not for exporting

#--------- run_cmd: run commands via bash or wsl ---------#
run_cmd <- function(cmd, intern = FALSE) {
  if (.Platform$OS.type != "windows") {
    system(command = cmd, intern = intern)
  } else {
    shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
  }
}

#--------- check_cmd: check whether a command exists in bash/wsl ---------#
check_cmd <- function(cmd) {
  test <- sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", cmd)
  run_cmd(test, intern = TRUE)
}

#--------- extract_numbers: extract numbers from a list of strings ---------#
extract_numbers <- function(string) {
  x <- unlist(regmatches(string, gregexpr('\\(?[0-9,.]+', string)))
  x <- as.numeric(gsub('\\(', '', gsub(',', '', x)))
  return(x)
}

#--------- reduce_path: remove file path and extension from file name ---------#
#' @importFrom tools file_path_sans_ext
reduce_path <- function(filename) {
  pattern1  <- "\\.(fastq|fq|fasta|fa|sam|bam|bed|bedgraph|bg|bigwig|bw)"
  pattern2  <- "\\.(txt|tar|bz2|bz|csv|zip|tsv).*"
  pattern3  <- "_fastqc$"
  reduced   <- tools::file_path_sans_ext(basename(filename))
  reduced   <- sub(pattern1, "", reduced)
  reduced   <- sub(pattern2, "", reduced)
  (reduced  <- sub(pattern3, "", reduced))
}

#--------- move_file: convenience function to move files ---------#
move_file <- function(tomove, to) {
  if (length(tomove) > 1) {
    lapply(tomove, function(x) {
      file.rename(from = x, to = file.path(to, x))
    })
  } else {
    file.rename(from = tomove, to = file.path(to, tomove))
  }
}

#--------- remove_file: convenience function to remove a file or directory ---------#
remove_file <- function(file) {
  if (file.exists(file)) {
    file.remove(file)
  } else {
    return("file not found")
  }
}

#--------- create_dir: convenience function to create a dir ---------#
create_dir <- function(dir, recurs = TRUE) {
  if (!dir.exists(dir)) {
    dir.create(path = dir, showWarnings = TRUE, recursive = recurs)
    return("directory created")
  } else {
    return("directory already exists")
  }
}

#--------- error messages: custom messages ----------#
#' @importFrom crayon red bold yellow underline green italic

give_error <- function(message) {stop(paste0("\n", crayon::red(crayon::bold(message), "\n")))}
give_warning <- function(message) {cat(crayon::yellow(crayon::underline(message)))}
give_note <- function(message) {cat(crayon::green(crayon::italic(message)))}


#-------- create design matrix ----------#
#' @importFrom stats formula

make_design <- function(sample.info) {
  if ("batch" %in% colnames(sample.info)) batch <- TRUE
  design <- stats::formula(paste0("~0+", "condition",
                                  ifelse(batch, paste0("+batch"), "")))
  design
}

#--------- convert ENSEMBL to SYMBOL --------#
#' @importFrom dplyr select distinct left_join
#' @importFrom magrittr set_colnames

ens2symbol <- function(result, columns.of.interest, colnames, species) {
  if (species == "human") ann <- grch38 %>% dplyr::select(gene = "ensgene", symbol) %>% dplyr::distinct()
  if (species == "mouse") ann <- grcm38 %>% dplyr::select(gene = "ensgene", symbol) %>% dplyr::distinct()
  output <- cbind(gene = row.names(result), as.data.frame(result), stringsAsFactors = FALSE)[, columns.of.interest] %>%
    dplyr::left_join(ann, by = "gene") %>% magrittr::set_colnames(colnames)
  output$symbol[which(output$symbol == "")] <- NA
  colnames(output) <- colnames
  return(output)
}

#--------- convert ENSEMBL to ENTREZ --------#
#' @importFrom dplyr select distinct left_join
#' @importFrom magrittr set_colnames

ens2entrez <- function(result, columns.of.interest, colnames, species) {
  if (species == "human") ann <- grch38 %>% dplyr::select(gene = "ensgene", "entrez") %>% dplyr::distinct()
  if (species == "mouse") ann <- grcm38 %>% dplyr::select(gene = "ensgene", "entrez") %>% dplyr::distinct()
  output <- cbind(gene = row.names(result), as.data.frame(result), stringsAsFactors = FALSE)[, columns.of.interest] %>%
    dplyr::left_join(ann, by = "gene") %>% magrittr::set_colnames(colnames)
  output$entrez[which(output$symbol == "")] <- NA
  colnames(output) <- colnames
  return(output)
}

#--------- common theme for plots -------------#
#' @import ggplot2

theme_publication <- function() {
  (ggplot2::theme_bw(base_size   = 14) +
     ggplot2::theme(
       plot.title  = ggplot2::element_text(face   = "bold",
                                           size   = ggplot2::rel(1.2),
                                           hjust  = 0.5),
       text                 = ggplot2::element_text(size   = 12),
       panel.background     = ggplot2::element_rect(colour = NA),
       plot.background      = ggplot2::element_rect(colour = NA),
       panel.border         = ggplot2::element_rect(colour = NA),
       axis.title           = ggplot2::element_text(face   = "bold",
                                                    size   = ggplot2::rel(1)),
       axis.title.y         = ggplot2::element_text(angle  = 90,
                                                    vjust  = 2),
       axis.title.x         = ggplot2::element_text(vjust  = -0.2),
       axis.text            = ggplot2::element_text(),
       axis.line            = ggplot2::element_line(colour = "black"),
       axis.ticks           = ggplot2::element_line(),
       panel.grid.major     = ggplot2::element_line(colour = "#f0f0f0"),
       panel.grid.minor     = ggplot2::element_blank(),
       legend.key           = ggplot2::element_rect(colour = NA),
       legend.key.size      = ggplot2::unit(0.7, "lines"),
       legend.spacing       = ggplot2::unit(0, "cm"),
       legend.title         = ggplot2::element_text(face   = "italic"),
       plot.margin          = ggplot2::unit(c(10, 5, 5, 5), "mm"),
       strip.background     = ggplot2::element_rect(colour = "#f0f0f0",
                                                    fill   = "#f0f0f0"),
       strip.text           = ggplot2::element_text(face   = "bold"))
  )
}

#----------------- W2U paths: Convert WIndows paths to linux for WSL ---------------#

convert_paths <- function(x) {
  if (.Platform$OS.type == "windows") {
    paths <- lapply(x, function(y) {
      cmd <- sprintf("wslpath %s", y)
      run_cmd(cmd, intern = TRUE)
    }) %>% unlist()
    return(paths)
  } else {
    return(x)
  }
}


