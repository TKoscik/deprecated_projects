convolve.paradigm <- function(data.onsets, scan.params, conv.type, save.dir, hrf.precision, ...) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # check save.dir ---------------------------------------------------------------------
  if (!dir.exists(save.dir)) {dir.create(save.dir)}
  
  # Parse onset files ------------------------------------------------------------------
  onset.fls <- parse.onsets(data.onsets)
  
  # Load onset files -------------------------------------------------------------------
  data.onsets <- vector("list", length(onset.fls))
  for (i in 1:length(onset.fls)) {
    data.onsets[[i]] <- read.csv(onset.fls[i])
  }
  tf <- unique(data.onsets[[1]][ ,1:2])
  sjx <- as.character(tf$subject)
  blx <- as.character(tf$block)
  
  # Parse scan parameters --------------------------------------------------------------
  scan.params <- parse.scan.params(scan.params, sjx, blx)
  
  # Initialize paradigm.df -------------------------------------------------------------
  subject <- character()
  block <- character()
  tr <- numeric()
  tr2 <- numeric()
  
  for (i in 1:nrow(scan.params)) {
    subject <- c(subject, rep(as.character(scan.params$subject[i]), scan.params$n.tr[i]))
    block <- c(block, rep(as.character(scan.params$block[i]), scan.params$n.tr[i]))
    tr <- c(tr, 1:scan.params$n.tr[i])
  }
  tr2 <- tr^2
  
  pdm.df <- data.frame(subject, block, tr, tr2)
  
  # Create convolved variables ---------------------------------------------------------
  if (conv.type=="hrf") {
    if (missing(hrf.precision)) {hrf.precision = 2L}
    for (i in 1:length(data.onsets)) {
      for (l in 1:(ncol(data.onsets[[i]])-4)) {
        conv.ts <- numeric()
        count <- 0
        for (j in 1:length(unique(sjx))) {
          which.sjx <- unique(sjx)[j]
          for (k in 1:length(unique(blx[which(sjx==which.sjx)]))) {
          count <- count +1
            which.blx <- blx[which(sjx==which.sjx)][k]
            tf <- data.onsets[[i]][which(data.onsets[[i]]$subject==which.sjx & data.onsets[[i]]$block==which.blx), ]
            conv.temp <- convolve.hrf(onset=tf$onset,
                                      duration=tf$duration,
                                      amplitude=tf[ ,(l+4)],
                                      tr=scan.params$tr[count],
                                      n.tr=scan.params$n.tr[count],
                                      st.betas=FALSE,
                                      scale=FALSE,
                                      precision=hrf.precision)
            conv.ts <- c(conv.ts, conv.temp$ts)
          }
        }
        pdm.df[colnames(data.onsets[[i]])[l+4]] <- conv.ts
      }
    }
  } else if (conv.type=="fir"){
    stop("FIR convolution not implemented yet")
  }
  
  write.table(pdm.df,
              file=paste(save.dir, "conv.paradigm.csv", sep="/"),
              row.names=FALSE,
              col.names=TRUE,
              quote=FALSE,
              sep=",")
}
