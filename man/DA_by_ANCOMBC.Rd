% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_DA.R
\name{DA_by_ANCOMBC}
\alias{DA_by_ANCOMBC}
\title{Wrapper for call_DA_ANCOMBC generating reference and target calls}
\usage{
DA_by_ANCOMBC(ref_data, data, groups, oracle_calls = NULL)
}
\arguments{
\item{ref_data}{absolute abundance data set (samples x features)}

\item{data}{relative abundance data set (samples x features)}

\item{groups}{group (cohort) labels}

\item{oracle_calls}{optional baseline (true) differential abundance calls}
}
\value{
baseline (true) differential abundance calls and calls made on 
relative abundances
}
\description{
Wrapper for call_DA_ANCOMBC generating reference and target calls
}
