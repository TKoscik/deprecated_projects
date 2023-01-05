parse.scan.params <- function(scan.params, sjx, blx) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(file.exists(scan.params))
  scan.temp <- read.csv(scan.params, as.is=TRUE)
  scan.params <- data.frame(subject=sjx, block=blx, tr=numeric(length(sjx)), n.tr=numeric(length(sjx)))
  for (i in 1:length(sjx)) {
    if (scan.temp$subject[1]=="all" & scan.temp$block[1]=="all") {
      scan.params$tr[i] <- scan.temp$tr[1]
      scan.params$n.tr[i] <- scan.temp$n.tr[1]
    } else if (scan.temp$subject[1]=="all") {
      scan.params$tr[i] <- scan.temp$tr[which(scan.temp$block==blx[i])]
      scan.params$n.tr[i] <- scan.temp$n.tr[which(scan.temp$block==blx[i])]
    } else if (scan.temp$block[1]=="all") {
      scan.params$tr[i] <- scan.temp$tr[which(scan.temp$subject==sjx[i])]
      scan.params$n.tr[i] <- scan.temp$n.tr[which(scan.temp$subject==sjx[i])]
    } else {
      scan.params$tr[i] <- scan.temp$tr[which(scan.temp$subject==sjx[i] & scan.temp$block==blx[i])]
      scan.params$n.tr[i] <- scan.temp$n.tr[which(scan.temp$subject==sjx[i] & scan.temp$block==blx[i])]
    }
  }
  scan.params <- unique(scan.params)
  return(scan.params)
}