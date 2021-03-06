% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/optFunctions.R
\name{getBIC}
\alias{getBIC}
\title{Get BIC of a given network configuration}
\usage{
getBIC(S, n, Omega, id)
}
\arguments{
\item{S}{sample covariance matrix}

\item{n}{number of samples}

\item{Omega}{estimated precision}

\item{id}{vector grouping elements of S, Omega}
}
\value{
Bayesian information criterion for model implied by Omega
}
\description{
In a zero-mean Gaussian graphical model, calculate the Bayesian information
criterion (BIC) of a model with the proposed precision matrix Omega.
}

