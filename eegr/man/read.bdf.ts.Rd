\name{read.bdf.ts}
\alias{read.bdf.ts}
\title{Read BDF Timeseries}
\description{Read BDF Timeseries}
\usage{read.bdf.ts(bdf.file, start = NULL, stop = NULL, channels = "all")}
\arguments{
  \item{bdf.file}{Full file listing to BDF file}
  \item{start}{Start time of section to read, default 0.}
  \item{stop}{Stop time of section to read, default full file}
  \item{channels}{Numeric indicating which channel(s) to read. Default = 'all' to read all channels}
}
\value{A matrix with # of channels acros the columns and number of samples down the rows.
}
\references{
<http://www.biosemi.com/faq/file_format.htm>
}
\author{
Timothy R. Koscik <timkoscik+eegr@gmail.com>
}