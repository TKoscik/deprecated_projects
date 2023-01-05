run.cluster <- function(data.4d,
                        data.paradigm,
                        data.cluster, which.cluster,
                        my.command,
                        my.function,
                        save.dir,
                        file.prefix,
                        num.cores,
                        save.mean.ts=FALSE,
                        verbose=FALSE, ...) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
# Parse 4D data ------------------------------------------------------------------------
  # data.4d <- parse.4d(data.4d)
  # n.4d <- length(data.4d)
  
# Parse paradigm data ------------------------------------------------------------------
  paradigm.df <- parse.paradigm(data.paradigm)
  
# Load mean time series ----------------------------------------------------------------
  if (file.exists(paste0(save.dir, "/", file.prefix, ".cluster.mean.ts.RData"))) {
    load(paste0(save.dir, "/", file.prefix, ".cluster.mean.ts.RData"))
    fmri <- df[ ,which.cluster]
    rm(list="df")
  } else {
    fmri <- load.cluster(data.4d, data.cluster, which.cluster, verbose=verbose)
  }
  n.cluster <- ncol(fmri)
  
# Save mean time series if desired -----------------------------------------------------
  if (save.mean.ts) {
    df <- data.frame(fmri, paradigm.df)
    save("df", file=paste0(save.dir, "/", file.prefix, ".cluster.mean.ts.RData"))
  }

# Check num.cores ----------------------------------------------------------------------
  if (missing(num.cores)) { num.cores <- detectCores() - 1 }
  # Check for Intel Math Kit and disable multithreading, MKL interferes with mclapply
  if (exists("getMKLthreads") & num.cores != 1) { setMKLthreads(1L) }
  
# Check model for output boilerplate ---------------------------------------------------
  if (missing(my.command) | is.null(my.command)) {
    my.command <- "my.function(df, save.dir, prefix, cluster.name, cluster.vol)"
  }
  
# Run model on clusters ------------------------------------------------------------------
  model.fxn <- function(X, ...) {
    if (verbose) { print(sprintf("%0.0f of %0.0f clusters", X, n.cluster)) }
    
    df <- data.frame(fmri=fmri[ ,X], paradigm.df)
    cluster.name <- data.cluster
    cluster.vol <- X
    prefix <<- paste0(file.prefix, ".cluster", X)
    eval(parse(text=my.command))
  }
  
  # Run models ---------------------------------------------------------------------------
  oldw <- getOption("warn")
  options(warn = -1)
  print(num.cores)
  if (num.cores > 1) {
    # for (mask.count in 1:n.masks) {
      mclapply(X=seq_len(n.cluster),
               FUN=model.fxn,
               mc.preschedule=FALSE,
               mc.cores=num.cores)
    # }
  } else {
    # for (mask.count in 1:n.masks) {
      lapply(X=seq_len(n.cluster),
             FUN=model.fxn)
    # }
  }
  options(warn=oldw)
}