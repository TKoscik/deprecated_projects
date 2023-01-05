\name{read.bdf.events}
\alias{read.bdf.events}
\title{Read Events from BDF File}
\description{Read Events from BDF File}
\usage{
read.bdf.events(bdf.file, event.ts = FALSE)
}
\arguments{
  \item{bdf.file}{Full file listing to BDF file}
  \item{event.ts}{Logical indicating whetehr or not to return the Events as a timeseries in addition to normal output}
}
\value{
  \item{sample}{A list containing sample numbers during which each marked event occur. One list per event type }
  \item{time}{A list containing the time in s when each event occured}
  \item{duration}{The length of between an event and the following event}
  \item{ts}{A time series of event markers}
}
\references{
<http://www.biosemi.com/faq/file_format.htm>
}
\author{
Timothy R. Koscik <timkoscik+eegr@gmail.com>
}