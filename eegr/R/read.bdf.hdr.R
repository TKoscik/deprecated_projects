read.bdf.hdr <-
function(bdf.file) {

# Check inputs -------------------------------------------------------------------------
  stopifnot(!missing(bdf.file), file.exists(bdf.file))

# Initialize BDF file ------------------------------------------------------------------
  fid <- file(bdf.file, 'rb')
  endian <- .Platform$endian
  
# Read Header --------------------------------------------------------------------------
  hdr <- list()
  hdr$id.ascii <- strtoi(readBin(fid, raw(), n=1, size=1, endian=endian), 16L)
  hdr$id.code <- rawToChar(readBin(fid, raw(), n=7, size=1, endian=endian))
  hdr$subject.id <- rawToChar(readBin(fid, raw(), n=80, size=1, endian=endian))
  hdr$record.id <- rawToChar(readBin(fid, raw(), n=80, size=1, endian=endian))
  hdr$start.date <- rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian))
  hdr$start.time <- rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian))
  hdr$hdr.bytes <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  hdr$data.format <- rawToChar(readBin(fid, raw(), n=44, size=1, endian=endian))
  hdr$data.format <- gsub("[[:space:]]*$", "", hdr$data.format)
  hdr$n.records <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  hdr$rec.dur.s <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  hdr$n.channels <- as.numeric(rawToChar(readBin(fid, raw(), n=4, size=1, endian=endian)))
  hdr$channel.labels <- character(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$channel.labels[i] <- rawToChar(readBin(fid, raw(), n=16, size=1, endian=endian))
    hdr$channel.labels[i] <- gsub("[[:space:]]*$", "", hdr$channel.labels[i])
  }
  hdr$transducer <- character(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$transducer[i] <- rawToChar(readBin(fid, raw(), n=80, size=1, endian=endian))
    hdr$transducer[i] <- gsub("[[:space:]]*$", "", hdr$transducer[i])
  }
  hdr$channel.dim <- character(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$channel.dim[i] <- rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian))
    hdr$channel.dim[i] <- gsub("[[:space:]]*$", "", hdr$channel.dim[i])
  }
  hdr$phys.min<- numeric(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$phys.min[i] <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  }
  hdr$phys.max <- numeric(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$phys.max[i] <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  }
  hdr$digi.min <- numeric(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$digi.min[i] <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  }
  hdr$digi.max <- numeric(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$digi.max[i] <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  }
  hdr$prefilter <- character(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$prefilter[i] <- rawToChar(readBin(fid, raw(), n=80, size=1, endian=endian))
    hdr$prefilter[i] <- gsub("[[:space:]]*$", "", hdr$prefilter[i])
    hdr$prefilter[i] <- gsub(": ", ":", hdr$prefilter[i])
  }
  hdr$n.samples <- numeric(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$n.samples[i] <- as.numeric(rawToChar(readBin(fid, raw(), n=8, size=1, endian=endian)))
  }
  hdr$reserved <- character(hdr$n.channels)
  for (i in 1:hdr$n.channels) {
    hdr$reserved[i] <- rawToChar(readBin(fid, raw(), n=32, size=1, endian=endian))
    hdr$reserved[i] <- gsub("[[:space:]]*$", "", hdr$reserved[i])
  }
  close(fid)
  return(hdr)
}
