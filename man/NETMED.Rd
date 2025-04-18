% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NETMED.R
\name{NETMED}
\alias{NETMED}
\title{NETMED: Network-Based Mediation Analysis}
\usage{
NETMED(A, M, Y, E, C = NULL, C_names = NULL)
}
\arguments{
\item{A}{Numeric vector of exposures.}

\item{M}{Numeric vector of mediators.}

\item{Y}{Numeric vector of outcomes.}

\item{E}{Symmetric adjacency matrix (entries between 0 and 1; diagonal entries must be 1).}

\item{C}{Optional numeric matrix of covariates.}

\item{C_names}{Optional character vector of covariate names.}
}
\value{
List containing the results of the network mediation analysis.
}
\description{
This function performs network-based mediation analysis using adjacency matrices,
mediators, exposures, and outcomes.
}
