% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_analysis.R
\name{pathway_analysis}
\alias{pathway_analysis}
\title{Title}
\usage{
pathway_analysis(x, species = c("human", "mouse"))
}
\arguments{
\item{x}{data.frame object containing at least the column gene (gene IDs in ENSEMBL format) and 
a column of logFC changes (entitled log2FoldChange or logFC)}

\item{species}{character string specifying the name of the species. Only
\code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].}
}
\description{
Title
}
