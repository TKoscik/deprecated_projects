run.voxel.doPar <- function(data.4d,
                      data.paradigm,
                      data.mask,
                      my.function,
                      my.command,
                      save.dir,
                      file.prefix,
                      num.cores,
                      log.name,
                      verbose, ...) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(!missing(data.4d), !missing(data.paradigm), !missing(my.function))
  if (missing(save.dir)) { save.dir <- getwd() }
  if (missing(save.dir)) { file.prefix <- "model" }
  if (missing(log.name)) { log.name <- paste0(save.dir, "/log.nii")}
  
  # Parse 4D data ------------------------------------------------------------------------
  data.4d <- parse.4d(data.4d)
  n.4d <- length(data.4d)
  
  # Parse paradigm data ------------------------------------------------------------------
  paradigm.df <- parse.paradigm(data.paradigm)
  
  # Check save directory -----------------------------------------------------------------
  # check.dir(save.dir, paste0(file.prefix, "*.nii$"))
  
  # Check mask ---------------------------------------------------------------------------
  which.voxels <- load.mask(data.mask)
  n.masks <- length(which.voxels)
  
  # Check num.cores ----------------------------------------------------------------------
  if (missing(num.cores)) { num.cores <- detectCores() - 1 }
  # Check for Intel Math Kit and disable multithreading, MKL interferes with mclapply
  if (exists("getMKLthreads") & num.cores != 1) { setMKLthreads(1L) }
  
  # Check model for output boilerplate ---------------------------------------------------
  if (missing(my.command) | is.null(my.command)) {
    my.command <- "my.function(df, coords, img.dims, save.dir, prefix)"
  }
  
  # Initialize Log -----------------------------------------------------------------------
  if (!file.exists(log.name)) {
    img.dims <<- nii.dims(data.4d[[1]][1])
    pixdim <- unlist(nii.hdr(data.4d[[1]][1], field="pixdim"))
    orient <- nii.orient(data.4d[[1]][1])
    init.nii(log.name, dims=img.dims[1:3], pixdim=pixdim, orient=orient)
  }
  
  # Run model on voxels ------------------------------------------------------------------
  model.fxn <- function(X, ...) {
    
    # Load brain data
    if (verbose) { print(sprintf("%0.0f of %0.0f voxels",
                                 X, nrow(which.voxels[[mask.count]]))) }
    coords <<- c(which.voxels[[mask.count]][X,1],
                 which.voxels[[mask.count]][X,2],
                 which.voxels[[mask.count]][X,3])
    
    status <- as.logical(read.nii.voxel(log.name, coords))
    
    if (!status | is.na(status)) {
      fmri <- load.voxel(data.4d, coords)
      
      if (!any(colSums(abs(fmri))==0)) {
        df <- data.frame(fmri, paradigm.df)
        prefix <<- paste0(file.prefix, ".mask", mask.count)
        img.dims <<- nii.dims(data.4d[[1]][1])[-4]
        eval(parse(text=my.command))
      }
      
      write.nii.voxel(log.name, coords, 1)
    }
  }
  
  # Run models ---------------------------------------------------------------------------
  oldw <- getOption("warn")
  options(warn = -1)
  print(num.cores)
  if (num.cores > 1) {
    for (mask.count in 1:n.masks) {
      # mclapply(X=seq_len(nrow(which.voxels[[mask.count]])),
      #          FUN=model.fxn,
      #          mc.preschedule=FALSE,
      #          mc.cores=num.cores)
      if (num.cores > 1) {
        registerDoParallel(num.cores)
        for (mask.count in 1:n.masks) {
          foreach(X=seq_len(nrow(which.voxels[[mask.count]]))) %dopar% model.fxn(X)
        }
        stopImplicitCluster()
      }
      
    }
  } else {
    for (mask.count in 1:n.masks) {
      lapply(X=seq_len(nrow(which.voxels[[mask.count]])),
             FUN=model.fxn)
    }
  }
  options(warn=oldw)
}
