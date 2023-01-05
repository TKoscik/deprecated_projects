convolve.hrf.2 <- function(onset,
                           duration = NULL,
                           amplitude = NULL,
                           tr, n.tr,
                           conv.fxn = "spm",
                           stbs = TRUE,
                           precision = 2,
                           num.cores = detectCores()/2) {
  
  if (missing(onset)) { stop("you must provide onset times")}
  if (missing(tr)) { stop("you must provide tr")}
  if (missing(n.tr)) { stop("you must provide n.tr")}
  
  if (conv.fxn == "spm") {
    kernel.length <- 32
    a1 <- 0.0083333333
    a2 <- 1.274527e-13
    p1 <- 5
    p2 <- 15
    t <- seq(from = 0, to = kernel.length, by = 1/(10^precision))
    conv.fxn <- exp(-t) * (a1*t^p1 - a2*t^p2)
    conv.fxn <- conv.fxn[3:length(conv.fxn)] # should this be a function of precision?
    conv.fxn <- conv.fxn / max(conv.fxn)
  }
  
  # precision onsets
  onset <- round(onset*10^precision) + 1
  
  # Amplitude Modulation
  if (is.null(amplitude)) {
    amplitude <- rep(1, length(onset))
  }
  
  # Duration Modulation
  if (is.null(duration)) {
    duration <- rep(1,length(onset))
  } else {
    duration <- round(duration * 10^precision)
  }
  
  # make timeseries
  n.samples <- tr * n.tr * 10^precision
  ts <- matrix(0, nrow = n.samples, ncol = length(onset))
  for (i in 1:length(onset)) {
    ts[(onset[i]+1):(onset[i]+duration[i]), i] <- rep(amplitude[i], duration[i])
  }
  ts.all <- rowSums(ts)
  
  # NO STBetas
  if (!stbs) {
    hrf <- convolve(ts.all, rev(conv.fxn), type="open")#[1:n.samples]
    hrf <- hrf[(10^precision*tr):length(hrf)]
    hrf <- hrf[seq(1, n.samples, tr*10^precision)]
  }
  
  # YES STBetas
  if (stbs) {
    my.fxn <- function(X, ...) {
      output <- convolve(ts[ ,X], rev(conv.fxn), type="open")#[1:n.samples]
    }
    hrf <- mcmapply(X = 1:length(onset),
                    FUN = my.fxn,
                    mc.preschedule=FALSE,
                    mc.cores=num.cores)
    hrf <- hrf[(10^precision*tr):nrow(hrf), ]
    hrf <- hrf[seq(1, n.samples, tr*10^precision), ]
  }
  
  return(hrf)
}