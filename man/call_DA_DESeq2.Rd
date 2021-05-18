% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_DA.R
\name{call_DA_DESeq2}
\alias{call_DA_DESeq2}
\title{Evaluate differential abundance with DESeq2 (via Seurat)}
\usage{
call_DA_DESeq2(data, groups)
}
\arguments{
\item{data}{simulated data set (samples x features)}

\item{groups}{group (cohort) labels}
}
\value{
unadjusted p-values for all features
}
\description{
Evaluate differential abundance with DESeq2 (via Seurat)
}