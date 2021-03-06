% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RULHMM.R
\name{trainHMM}
\alias{trainHMM}
\title{Train Hidden Markov Model (trainHMM)}
\usage{
trainHMM(ObsSeq, obsAlphabetSize, hidAlphabetSize)
}
\arguments{
\item{ObsSeq}{A discrete observed time series. Output from toWPDSAX.}

\item{obsAlphabetSize}{The alphabet size of the observed sequence.}

\item{hidAlphabetSize}{The alphabet size of the hidden sequence.}
}
\value{
The probability matrices, mean, standard deviation vectors and the final failure state of the hidden markov model.
}
\description{
Train Hidden Markov Models using the Baum Welch Algorithm
}

