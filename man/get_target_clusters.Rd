% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_target_clusters}
\alias{get_target_clusters}
\title{Extract target gene clusters for given penalization parameters}
\usage{
get_target_clusters(fit, penalization = NULL)
}
\arguments{
\item{fit}{An object of class \code{scregclust}}

\item{penalization}{A numeric vector of penalization parameters.
The penalization parameters specificed here must have
been used used during fitting of the \code{fit} object.}
}
\value{
A list of lists of final target clusters. One list for each
parameter in \code{penalization}. The lists contain the clustering of
target genes for each final configuration.
}
\description{
Extract target gene clusters for given penalization parameters
}
\concept{utilities}
