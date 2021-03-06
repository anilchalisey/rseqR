% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_upset.R
\name{make_upset}
\alias{make_upset}
\title{Plot intersections of list of DE genes}
\usage{
make_upset(limma, edger, deseq2)
}
\arguments{
\item{limma}{data.frame with the columns gene and padj or FDR; result of Limma-voom analysis}

\item{edger}{data.frame with the columns gene and padj or FDR; result of edgeR analysis}

\item{deseq2}{data.frame with the columns gene and padj or FDR; result of DESeq2 analysis}
}
\value{
An UpSet plot of gene intersections.
}
\description{
Plot intersections of list of DE genes
}
