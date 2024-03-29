% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_DA.R
\name{DA_wrapper}
\alias{DA_wrapper}
\title{Wrapper function for differential abundance calling}
\usage{
DA_wrapper(
  ref_data,
  data,
  groups,
  method = "NBGLM",
  oracle_calls = NULL,
  control_indices = NULL
)
}
\arguments{
\item{ref_data}{absolute abundance data set (samples x features)}

\item{data}{relative abundance data set (samples x features)}

\item{groups}{group (cohort) labels}

\item{method}{differential abundance calling method (e.g. "DESeq2")}

\item{oracle_calls}{optional baseline (true) differential abundance calls}

\item{control_indices}{if specified, these features are used as the stable 
features against which to normalize observed abundances (DESeq2-only)}
}
\value{
named list of baseline (true) differential abundance calls and
calls made on relative abundances (the observed data)
}
\description{
Wrapper function for differential abundance calling
}
