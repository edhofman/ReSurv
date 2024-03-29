% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotReSurvFit.R
\name{plot.ReSurvFit}
\alias{plot.ReSurvFit}
\title{Plot for machine learning models}
\usage{
\method{plot}{ReSurvFit}(object, nsamples = NULL)
}
\arguments{
\item{object}{\code{ReSurvFit} object.}

\item{nsamples}{\code{integer}, number of observations to sample for neural networks features importance plot.}
}
\value{
predictions
}
\description{
This function plots the mean absolute SHAP values for the ReSurv fits of machine learning models.
}
