% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_val.R
\name{p_val}
\alias{p_val}
\title{This function gives us test statistic and corresponding p-value for testing if a scalar predictor is correlated to a functional response.}
\usage{
p_val(data, kernel = 0, use.sig = 1)
}
\arguments{
\item{data}{has n data points comprising X, Y, Z, and t}

\item{X}{is a matrix. Each column of X corresponds to q many covariates of an individual}

\item{Z}{is a matrix. Each column of Z corresponds to r many nuisance covariates of an individual.}

\item{t}{is the list of n time points at which the data is available}

\item{Y}{is a list of n responses at time points corresponding to t}

\item{ker}{= 1 corresponds to linear kernel for X
= 2 corresponds to quadratic kernel for X
= 0 will correspond to gaussian kernel in the future}

\item{sig}{is indicator for whether we are using sigma in the test statistic}
}
\value{
p.val is the p-value for the given data.

test.stat is the test statistic obtained.
}
\description{
This function gives us test statistic and corresponding p-value for testing if a scalar predictor is correlated to a functional response.
}
\details{
A functional response variable Y is regressed on scalar variables X and Z to determine if X is a significant contributor to Y.
}
\author{
Mityl Biswas
}
