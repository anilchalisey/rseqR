% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_MA.R
\name{make_MA}
\alias{make_MA}
\title{Create MA plot}
\usage{
make_MA(data, fdr = 0.05, label.rectangle = FALSE,
  proc = c("Limma-voom", "edgeR", "DESeq2"), top = 0,
  select.method = c("padj", "logfc"), main = NULL)
}
\arguments{
\item{data}{\code{matrix} or \code{data.frame} with the columns for average
gene expression, adjusted P-value, and log fold change; created as part of
\code{edger_analysis}, \code{limma_voom_analysis}, or \code{deseq2_analysis}.}

\item{fdr}{numerical value between 0 and 1 indicating the threshold false
discovery rate for discovering differentially expressed genes. [DEFAULT = 0.05]}

\item{label.rectangle}{logical indicating whether to add rectangles underneath
gene names [DEFAULT = FALSE]}

\item{proc}{character string indicating analysis procedure used.  One of
"Limma-voom", "edgeR", "DESeq2". [DEFAULT = "Limma-voom].}

\item{top}{integer indicating the number of genes to label}

\item{select.method}{Character string identifying the method of selecting the
top genes.  One of \code{padj} or \code{logfc}}

\item{main}{character string giving the plot main title}
}
\value{
A ggplot2 object
}
\description{
Create a scatter plot of log2 fold changes vs the mean expression signal.
}
