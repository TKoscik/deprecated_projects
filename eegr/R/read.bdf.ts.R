read.bdf.ts <-
function(bdf.file, start=NULL, stop=NULL, channels="all") {
# Check inputs -------------------------------------------------------------------------
  stopifnot(!missing(bdf.file), file.exists(bdf.file))
  
# Retreive BDF Header ------------------------------------------------------------------
  hdr <- read.bdf.hdr(bdf.file)
  
# Parse Channel Input ------------------------------------------------------------------
  if (is.character(channels)) {
    chan.temp <- channels
    channels <- numeric()
    for (i in 1:length(chan.temp)) {
      if (chan.temp[i] %in% c("all","eeg", "misc", "trig")) {
        channels <- switch(chan.temp[i],
          `all` = c(channels, 1:hdr$n.channels),
          `eeg` = c(channels, which(hdr$reserved=="MON")),
          `misc` = c(channels, which(hdr$reserved!="MON" & hdr$reserved!="TRI")),
          `trig` = c(channels, which(hdr$reserved=="TRI")))
      } else if (chan.temp[i] %in% hdr$channel.labels) {
        channels <- c(channels, which(chan.temp[i] == hdr$channel.labels))
      } else {
        stop("Unknown channel specified")
      }
    }
  } else if (is.numeric(channels)) {
    stopifnot(all(channels %in% 1:hdr$n.channels))
  } else {
    stop("Cannot parse channel input")
  }
  
# Parse start and stop times -----------------------------------------------------------
  if (is.null(start)) {
    start <- 0
  } else {
    stopifnot(start < hdr$n.records - 1)
  }
  
  if (is.null(stop)) {
    stop <- hdr$n.records - 1
  } else {
    stopifnot(stop <= hdr$n.records - 1)
  }
  
# Initialize BDF file ------------------------------------------------------------------
  fid <- file(bdf.file, 'rb')
  endian <- .Platform$endian
  
# Read data ----------------------------------------------------------------------------
  sample.rate <- unique(hdr$n.samples[channels])
  if (length(sample.rate) > 1) {
    stop("Some of the specified channels have different sample rates,
         please read these in separately")
  }
  
  recs.read <- floor(start):floor(stop)
  n.recs <- length(recs.read)
  
  ts <- matrix(0, nrow=sample.rate*n.recs, ncol=length(channels))
  colnames(ts) <- hdr$channel.labels[channels]
  
  for (i in 1:n.recs) {
    offset <- hdr$hdr.bytes + sum(hdr$n.samples)*recs.read[i]*3
    seek(fid, where=offset, origin="start", size=1)
    temp <- strtoi(readBin(fid, raw(), n=3*sum(hdr$n.samples),
                           size=1, endian=endian), 16L)
    temp <- t(matrix(temp, nrow=3))
    temp <- c(temp[1:((hdr$n.channels-1)*sample.rate),1] + 
                temp[1:((hdr$n.channels-1)*sample.rate),2] * 256 + 
                temp[1:((hdr$n.channels-1)*sample.rate),3] * 256 * 256,
              temp[((hdr$n.channels-1)*sample.rate+1):nrow(temp),1] +
                temp[((hdr$n.channels-1)*sample.rate+1):nrow(temp),2] * 256)
    temp <- matrix(temp, ncol=hdr$n.channels)
    # temp[ ,80] <- temp[ ,80] - 153*256*256
    temp <- temp[ ,channels]
    ineg <- temp >= 256*256*128
    temp[ineg] <- temp[ineg] - 256*256*256
    gain <- (hdr$phys.max-hdr$phys.min)/(hdr$digi.max-hdr$digi.min)/1000000
    gain[hdr$n.channels] <- 1
    gain <- gain[channels]
    temp = sweep(temp, MARGIN=2, gain, '*')
    # temp <- temp[ ,channels]
    ts[((i-1)*sample.rate+1):(i*sample.rate), ] <- temp
  }
  
# Trim TS to user precision
  trim.start <- round((start - floor(start)) * sample.rate)
  trim.stop <- sample.rate - round((stop - floor(stop)) * sample.rate)
  ts <- ts[trim.start:(nrow(ts)-trim.stop), ]
  
  close(fid)
  return(ts)
}
