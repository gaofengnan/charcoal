% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithms.R
\name{orthogonalProjection}
\alias{orthogonalProjection}
\title{Computing the orthonormal matrix spanning the orthogonal complement of the column span of X}
\usage{
orthogonalProjection(X)
}
\arguments{
\item{X}{an n x p matrix with n > p}
}
\value{
matrix A of size n x (n-p) with orthonormal columns
}
\description{
it seems that QR decomposition in R is slow when the dimension is large e.g. n and p > 5000; when working with large matrices, we recommend use either the python or matlab version of this package.
}
