% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/optFunctions.R
\name{getAIC}
\alias{getAIC}
\title{Get AIC of a given network configuration}
\usage{
getAIC(S, n, Omega, id)
}
\arguments{
\item{S}{sample covariance matrix}

\item{n}{number of samples}

\item{Omega}{estimated precision}

\item{id}{vector grouping elements of S, Omega}
}
\value{
Akaike information criterion for model implied by Omega
}
\description{
In a zero-mean Gaussian graphical model, calculate the Akaike information
criterion (AIC) of a model with the proposed precision matrix Omega.
}

