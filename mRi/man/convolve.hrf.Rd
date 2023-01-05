\name{convolve.hrf}
\alias{convolve.hrf}
\title{Convolve With Canonical HRF (from AFNI)}
\description{}
\usage{
  convolve.hrf(onset, duration=NULL, amplitude=NULL, tr, n.tr, st.betas, scale=FALSE, precision=2)
}
\arguments{
  \item{onset}{
  a vector of stimulus n onset, nearest 3 digits If duration are not provided, onset are assumed to be finite impulses.
  }
  \item{duration}{
  default = NULL, can be a scalar if stimulus duration is the same at every onset or a vector of values to specify a duration at each onset timepoint.
  }
  \item{amplitude}{
  default = NULL, can be a scalar if amplitude modulation is the same at every onset or a vector of values to specify the amplitude at each onset timepoint.
  }
  \item{tr}{repetition time, s}
  \item{n.tr}{number of volumes in timeseries}
  \item{st.betas}{logical indicating whether or not to convolve each onset individually, Default = TRUE}
  \item{scale}{Logical indicating whether or not to scale the design matrix. Default=FALSE}
  \item{precision}{}
}
\value{}
\author{
  Timothy R. Koscik <timkoscik+mripkg@gmail.com>
}
\examples{}