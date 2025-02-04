% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDistMat.R
\name{calculateDistMatOneModel}
\alias{calculateDistMatOneModel}
\title{Calculate Distance Matrix Using a Single Model}
\usage{
calculateDistMatOneModel(
  matrix,
  metadata,
  verbose = TRUE,
  method = "voom",
  additional.variate = NULL
)
}
\arguments{
\item{matrix}{A count matrix with rows representing genes or features and columns representing samples or cells.}

\item{metadata}{A data frame containing metadata corresponding to the samples or cells. 
Each row should match a column in \code{matrix}.}

\item{verbose}{Logical. If \code{TRUE}, the function displays progress messages 
and a progress bar. The default is \code{TRUE}.}

\item{method}{A character string specifying the method for differential 
expression analysis. Options are "voom" or "trend", with "trend" as the default.}

\item{additional.variate}{A character vector of additional variates to 
include in the linear model for regression.}
}
\value{
A similarity matrix.
}
\description{
This function computes a similarity matrix by utilising 
a single linear model for differential expression analysis.
}
\seealso{
\code{\link{getDistMat}}
}
