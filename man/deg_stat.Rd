% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.R
\name{deg_stat}
\alias{deg_stat}
\title{Compute the DEGs}
\usage{
deg_stat(
  expr_mat,
  groups,
  pCutoff = 0.05,
  log2FCCutoff = 1,
  method = "t.test",
  min.count = 0
)
}
\arguments{
\item{expr_mat}{a data frame which columns are genes, rows are samples. the values are counts, tpm or fpkm.}

\item{groups}{a factor which contains group information.}

\item{pCutoff}{a cutoff for p.value.}

\item{log2FCCutoff}{a cutoff for log2FC.}

\item{method}{one of the t.test, wilcox.test, limma.}

\item{min.count}{numeric. Minimum count required for at least some samples (when counts).}
}
\value{
a data frame which contains gene, P.Value, statistic, logFC and AveExpr columns.
}
\description{
Compute the DEGs
}
