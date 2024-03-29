% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_DA.R
\name{call_DA_DESeq2}
\alias{call_DA_DESeq2}
\title{Evaluate differential abundance with DESeq2 (via Seurat)}
\usage{
call_DA_DESeq2(data, groups, control_indices = NULL)
}
\arguments{
\item{data}{simulated data set (samples x features)}

\item{groups}{group (cohort) labels}

\item{control_indices}{if specified, these features are used as the stable 
features against which to normalize observed abundances (DESeq2-only)}
}
\value{
unadjusted p-values and "effect sizes" for DA for all features
}
\description{
Evaluate differential abundance with DESeq2 (via Seurat)
}
