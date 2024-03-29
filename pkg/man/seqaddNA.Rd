% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missings_generation.R
\name{seqaddNA}
\alias{seqaddNA}
\title{Generation of missing data under the form of gaps, which
is the typical form of missing data with longitudinal data.
A missing completely at random (MAR) mechanism is used.}
\usage{
seqaddNA(
  data,
  states.high = NULL,
  propdata = 1,
  pstart.high = 0.1,
  pstart.low = 0.005,
  maxgap = 3
)
}
\arguments{
\item{data}{either a data frame containing sequences of a multinomial variable with missing data (coded as \code{NA}) or
a state sequence object built with the TraMineR package}

\item{states.high}{list of states that will have a larger probability to trigger a
subsequent gap of missing data}

\item{propdata}{proportion of the observations on which missing data will be simulated}

\item{pstart.high}{probability to start a missing data for the specified states}

\item{pstart.low}{probability to start a missing data for the other states}

\item{maxgap}{maximum length of a gap of missing data}
}
\value{
Returns either a data frame or a state sequence object, depending
the type of data that was provided to the function
}
\description{
Generation of missing data under the form of gaps, which
is the typical form of missing data with longitudinal data.
A missing completely at random (MAR) mechanism is used.
}
