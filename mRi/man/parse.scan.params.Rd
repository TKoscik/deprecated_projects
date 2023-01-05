\name{parse.scan.params}
\alias{parse.scan.params}
\title{
Parse Scan Parameters Input
}
\description{
Parse Scan Parameters Input
}
\usage{
parse.scan.params(scan.params, sjx, blx)
}
\arguments{
  \item{scan.params}{
Full directory listing to scan parameters csv file
}
  \item{sjx}{
number of subjects
}
  \item{blx}{
number of scans per subject
}
}
\value{
A data frame for generating convolved CSV paradigm files.
}
\author{
Timothy R. Koscik <timkoscik+mripkg@gmail.com>
}
\note{
mRi internal function
}