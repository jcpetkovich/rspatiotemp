% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{updateLifeTab}
\alias{updateLifeTab}
\title{Update Life Table (updateLifeTab)}
\usage{
updateLifeTab(lifeTab, accDeg)
}
\arguments{
\item{lifeTab}{The life table created from createLifeTab or updateLifeTab.}

\item{accDeg}{The output of the 'accDegrad' function. A vector containing the Accumulated Degredation values of the time series data.}
}
\value{
A Life Table used for RUL computation.
}
\description{
Updates Life Table produced from 'createLifeTab' function
}

