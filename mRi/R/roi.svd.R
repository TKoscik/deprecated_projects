roi.svd <- function(data.4d, data.mask, data.paradigm, scan.start.vol=1, var.name="roi.svd") {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  data.4d <- parse.4d(data.4d)
  n.4d <- length(data.4d[[1]])
  
  n.mask <- length(data.mask)
  
  if (length(roi.svd) < n.mask) {
    for (i in length(roi.svd):n.mask) {
      roi.svd <- c(roi.svd, paste0("roi.svd.", i))
    }
  }
  
  paradigm.df <- parse.paradigm(data.paradigm)
  
  for (i in 1:n.mask) {
    mask <- read.nii.volume(data.mask[i], 1)
    s <- numeric()
    for (j in 1:n.4d) {
      if (length(scan.start.vol) == 1) {
        n.scans <- 1
        scan.stop.vol <- nii.dims(data.4d[[1]][j])[4]
      } else {
        n.scans <- length(scan.start.vol)
        scan.stop.vol <- c(scan.start.vol[2:n.scans]-1, nii.dims(data.4d[[1]][j])[4])
      }
      for (k in 1:n.scans) {
        tf <- data.frame(matrix(0, nrow=sum(mask), ncol=(scan.stop.vol-scan.start.vol)+1))
        count <- 0
        for (l in scan.start.vol[k]:scan.stop.vol[k]) {
          count <- count + 1
          uf <- read.nii.volume(data.4d[[1]][j], l)
          tf[ ,count] <- uf[mask==1]
        }
        
        s.temp <- svd(tf)$v[ ,1]
        m <- colMeans(tf)
        r <- cor(s.temp,m)
        if (sign(r)==-1) {
          s.temp <- s.temp * -1
        }
        
        s <- c(s, s.temp)
        # write.table(s, file=paste0(save.dir, "/", fname[i]), sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
      }
    }
    paradigm.df <- cbind(paradigm.df, s)
    colnames(paradigm.df)[ncol(paradigm.df)] <- var.name[i]
    write.table(paradigm.df, file=data.paradigm, sep=",", row.names=F, col.names=T)
  }
}
