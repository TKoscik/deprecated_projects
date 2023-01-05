convolve.hrf <- function (onset, duration=NULL, amplitude=NULL,
                          tr, n.tr, st.betas, scale=FALSE, precision=2) {
  #-----------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------

# Check onset -------------------------------------------------------------------------
  onset <- unlist(round(onset, digits=precision)) + 1
  onset0 <- which(onset < 0)
  onset.neg <- onset * (onset < 0)
  onset[onset0] <- 0

# Setup Basis Function --- (currently only SPMG1 is implemented) -----------------------
  params <- list()
  params$precision <- 10^precision
  params$kernel.length <- 32
  params$a1 <- 0.0083333333
  params$a2 <- 1.274527e-13
  params$p1 <- 5
  params$p2 <- 15
  params$t <- seq(from=0, to=params$kernel.length, by=1/params$precision)
  params$hrf <- exp(-params$t)*(params$a1*params$t^params$p1-params$a2*params$t^params$p2)
  params$hrf <- params$hrf[3:length(params$hrf)] # should this be a function of precision?
  params$hrf <- params$hrf / max(params$hrf);
  hrfs <- list()

# Convolve on square wave to account for stimulus duration ----------------------------
  if (!is.null(duration)) {
    if (length(duration) == 1) {
      duration <- round(duration * params$precision)
      sq.wave <- numeric(duration)
      if (duration == 0 ) {
        sq.wave[1] <- 1
      } else {
        sq.wave[1:duration] <- 1
      }
      hrfs[[1]] <- convolve(sq.wave, rev(params$hrf), type="open")
      hrfs[[1]] <- hrfs[[1]] / max(hrfs[[1]])
    } else if (length(duration) == length(onset)) {
      for (i in 1:length(duration)) {
        dur.temp <- round(duration[i] * params$precision)
        sq.wave <- numeric(dur.temp)
        if (duration[i] == 0 ) {
          sq.wave[1] <- 1
        } else {
          sq.wave[1:dur.temp] <- 1
        }
        hrfs[[i]] <- convolve(sq.wave, rev(params$hrf), type="open")
        hrfs[[i]] <- hrfs[[i]] / max(hrfs[[i]])
      }
    } else
      stop("duration must be a scalar value or a vector of the same length as onset")
  } else {
    hrfs[[1]] <- params$hrf
  }

# Amplitude Modulation -----------------------------------------------------------------
  if (is.null(amplitude)) {
    amplitude <- rep(1,length(onset))
  }

# Convolve HRF -------------------------------------------------------------------------
  n.ts <- length(onset)
  ts <- matrix(0, nrow=(n.tr*tr+params$kernel.length)*params$precision, ncol=n.ts)
  for (i in 1:n.ts) {
    temp.ts <- numeric(nrow(ts))
    temp.ts[onset[i]*params$precision+1] <- amplitude[i]
    if (length(hrfs) == 1) {
      temp.hrf <- convolve(temp.ts, rev(hrfs[[1]]), type="open")
    } else {
      temp.hrf <- convolve(temp.ts, rev(hrfs[[i]]), type="open")
    }
    if (i %in% onset0) {
      temp.hrf <- c(temp.hrf[-(1:(abs(onset.neg[onset0])*params$precision))],
                    rep(0, abs(onset.neg[onset0])*params$precision))
    }

    ts[ ,i] <- temp.hrf[1:nrow(ts)]
  }

  if (!st.betas) {
    ts <- rowSums(ts)
  }

  if (is.vector(ts)) {
    ts <- ts[(params$precision*tr):length(ts)]
    ts <- ts[seq(from=1, to=(n.tr*tr*params$precision), by=tr*params$precision)]
  } else if (is.matrix(ts)) {
    ts <- ts[(params$precision*tr):nrow(ts), ]
    ts <- ts[seq(from=1, to=(n.tr*tr*params$precision), by=tr*params$precision), ]
  }

  output <- list()
  output$params <- params
  output$hrfs <- hrfs
  output$ts <- ts
  return(output)
}
