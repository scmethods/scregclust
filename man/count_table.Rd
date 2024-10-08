% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{count_table}
\alias{count_table}
\title{Format count table nicely}
\usage{
count_table(counts, title, row_names, col_width = 5)
}
\arguments{
\item{counts}{a list of count vectors with \code{1 + n_cl} entries each.
\code{NA} values are replaced with \code{-}}

\item{title}{title above the table}

\item{row_names}{a vector of row names, one for each count vector}

\item{col_width}{minimum width for columns}
}
\value{
A string formatted as a table
}
\description{
Format count table nicely
}
\keyword{internal}
