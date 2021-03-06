% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RULHMM.R
\name{createModel}
\alias{createModel}
\title{Create HMM (createModel)}
\usage{
createModel(timeSeries, SAXalphabetSize, SAXgroupSize, exp2, hidAlphabetSize)
}
\arguments{
\item{timeSeries}{The time series to be converted.}

\item{SAXalphabetSize}{The alphabet size used for the SAX conversion.}

\item{SAXgroupSize}{The size of each group to be converted to a single SAX symbol.}

\item{exp2}{Two to the power of exp2 will be the size.}

\item{hidAlphabetSize}{The alphabet size of the hidden sequence.}
}
\value{
A HMM generated from the time series.
}
\description{
Create a single hidden markov model from a given time series.
}

