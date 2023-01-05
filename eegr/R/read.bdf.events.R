read.bdf.events <-
function(bdf.file, event.ts=FALSE) {

# Check inputs -------------------------------------------------------------------------
  stopifnot(!missing(bdf.file), file.exists(bdf.file))

# Retreive BDF Header ------------------------------------------------------------------
  hdr <- read.bdf.hdr(bdf.file)
  
# Initialize BDF file ------------------------------------------------------------------
  fid <- file(bdf.file, 'rb')
  endian <- .Platform$endian
  
# Read Events --------------------------------------------------------------------------
  ts <- numeric(hdr$n.records * hdr$n.samples[hdr$n.channels])
  for (i in 1:hdr$n.records) {
    offset <- hdr$hdr.bytes +
      sum(hdr$n.samples[-hdr$n.channels])*3 +
      (sum(hdr$n.samples))*(i-1)*3
    seek(fid, where=offset, origin="start", size=1)
    temp <- strtoi(readBin(fid,
                           raw(),
                           n=3*hdr$n.samples[hdr$n.channels],
                           size=1,
                           endian=endian),
                   16L)
    temp <- t(matrix(temp, nrow=3))
    temp <- temp[ ,-3]
    temp <- temp[ ,1] + temp[ ,2] * 256
    ts[((i-1)*hdr$n.samples[hdr$n.channels]+1):(i*hdr$n.samples[hdr$n.channels])] <-temp 
  }
  close(fid)
  
# Make event list ----------------------------------------------------------------------
  event.ls <- unique(ts)
  event.ls <- sort(event.ls[-which(event.ls==0)])
  events <- vector("list", length=length(event.ls))
  for (i in 1:length(event.ls)) {
    events[[i]] <- list()
    events[[i]]$marker <- event.ls[i]
    ev.temp <- which(ts==event.ls[i])
    ev.diff <- c(ev.temp[1], diff(ev.temp))
    ev.locs <- which(ev.diff > 1)
    events[[i]]$sample <- ev.temp[ev.locs]
    events[[i]]$time <- (1/hdr$n.samples[hdr$n.channels]) * (events[[i]]$sample - 1)
    events[[i]]$duration <- c(diff(ev.locs), length(ev.temp) - ev.locs[length(ev.locs)])
  }
  
# Optionally append timeseries to output -----------------------------------------------
  if (event.ts) {
    events$ts <- ts
  }
  
  return(events)
}
