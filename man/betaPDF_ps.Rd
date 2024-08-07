% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaPDF.R
\name{betaPDF_ps}
\alias{betaPDF_ps}
\title{Performs betaPDF calculation on Parameter Space row}
\usage{
betaPDF_ps(psv, k = 1, interval = 0.001)
}
\arguments{
\item{psv}{(row) vector, or single row data frame from parameter space}

\item{k}{Scale parameter, usually left at 1 by default}

\item{interval}{Interval length for steps between a and b. By default 0.001}
}
\value{
Data frame containing X and Y values for the desired intervals between a and b
}
\description{
Performs betaPDF calculation on Parameter Space row
}
\examples{
a <- 1; b <- 2
alpha = gamma <- 1:2
par_space <- set_parameter_space(a, b, alpha, gamma)
betaPDF_ps(psv = par_space[1, ])
}
