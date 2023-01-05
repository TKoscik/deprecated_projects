deconvolve <- function(data.4d,
                       data.onsets,
                       scan.params,
                       data.mask,
                       save.dir=getwd(),
                       num.cores,
                       incl.trend = c(TRUE, FALSE)) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  data.4d <- parse.4d(data.4d)
  n.scans <- length(data.4d[[1]])
  
  data.onsets <- parse.onsets(data.onsets)
  if (length(data.onsets) == 1) {
    data.onsets <- rep(data.onsets, n.scans)
  } 
  stopifnot(length(data.onsets) == n.scans) 
  
  scan.params <- parse.scan.params(scan.params, 1:n.scans, 1L)
  
  stopifnot(length(data.mask)==1)
  which.voxels <- load.mask(data.mask)
  n.masks <- length(which.voxels)
  
# Deconvolution function ---------------------------------------------------------------
  deconv.fun <- function(X, ...) {
    print(sprintf("%0.0f of %0.0f voxels", X, nrow(which.voxels[[1]])))
    coords <<- c(which.voxels[[1]][X,1],
                 which.voxels[[1]][X,2],
                 which.voxels[[1]][X,3])
    
    fmri <- read.nii.voxel(data.4d[[1]][i], coords=c(coords, Inf))
    
    if (sum(abs(fmri)) != 0) {
      df$fmri <- fmri
      for (j in 1:n.trials) {
        if (all(incl.trend)) {
          mdl <- lm(df$fmri ~ df[ ,j] + df$trend.1 + df$trend.2)
        } else if (any(incl.trend)) {
          mdl <- lm(df$fmri ~ df[ ,j] + df$trend)
        } else {
          md  <- lm(df$fmri ~ df[ ,j])
        }
        beta <- summary(mdl)$coef[2,1]
        write.nii.voxel(nii.file=paste0(save.file, ".nii"), coords=c(coords, j), value=beta)
      }
    }
  }
  
# Loop through scans -------------------------------------------------------------------
  for (i in 1:n.scans) {
    design <- read.csv(data.onsets[i])
    hrf <- convolve.hrf(onset=design$onset,
                        duration=design$duration,
                        amplitude=design$amplitude,
                        tr=scan.params$tr[i],
                        n.tr=scan.params$n.tr[i],
                        st.betas=TRUE, scale=TRUE, precision=2L)
    
    img.dims <- nii.dims(data.4d[[1]][i])
    n.voxels <- prod(img.dims[1:3])
    n.trials <- nrow(design)
    
    trend.regs <- numeric()
    if (incl.trend[1]) {
      trend.regs <- cbind(trend.regs, 1:scan.params$n.tr[i])
    }
    if (incl.trend[2]) {
      trend.regs <- cbind(trend.regs, (1:scan.params$n.tr[i])^2)
    }
    if (any(incl.trend)) {
      df <- data.frame(hrf=hrf$ts, trend=trend.regs)
    } else {
      df <- data.frame(hrf=hrf$ts)
    }
    
    # initialize output
    fname <- strsplit(data.4d[[1]][i], "[/]")[[1]]
    fname <- fname[length(fname)]
    fname <- strsplit(fname, "[.]")[[1]]
    fname <- fname[-length(fname)]
    fname <- paste(fname, collapse=".")
    save.file <- paste0(save.dir, "/", fname, ".nii")
    init.nii(file.name = save.file,
             dims = c(img.dims[1:3], n.trials),
             pixdim = unlist(nii.hdr(data.4d[[1]][i], "pixdim")),
             orient = nii.orient(data.4d[[1]][i]))
    
    mclapply(X=seq_len(nrow(which.voxels[[1]])),
             FUN=deconv.fun,
             mc.preschedule=FALSE,
             mc.cores = num.cores)
  }
}