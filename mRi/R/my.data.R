my.data <- function(data.4d = NULL,
                           data.paradigm,
                           data.mask = NULL,
                           vxl.order = "rand",
                           save.dir,
                           log.name = "log",
                           prefix,
                           num.cores = detectCores()-1,
                           verbose = TRUE) {
  #-----------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-----------------------------------------------------------------------------

  my.data <- list()

  # Get scan list --------------------------------------------------------------
  if (!is.null(data.4d)) { my.data$scans <- parse.4d(data.4d) }

  # Get paradigm data ----------------------------------------------------------
  my.data$paradigm <- parse.paradigm(data.paradigm)

  # Get voxels in mask ---------------------------------------------------------
  if (!is.null(data.mask)) {
    my.data$vxls <- load.mask(data.mask)
    my.data$n.masks <- length(my.data$vxls)
  } else if (!is.null(data.4d)) {
    my.data$vxls <- which(array(TRUE, dim=nii.dims(my.data$scans[[1]][1])),
                            arr.ind = TRUE)
    my.data$n.masks <- 1
  }

  my.data$vxl.order <- vxl.order

  # Setup save directory -------------------------------------------------------
  count <- 0
  while (dir.exists(save.dir)) {
    count <- count + 1
    save.dir <- paste0(save.dir, count)
  }
  dir.create(save.dir)
  my.data$save.dir <- save.dir

  # set log name ---------------------------------------------------------------
  my.data$log <- paste(save.dir, log.name, sep="/")

  # store other varaiables -----------------------------------------------------
  my.data$prefix <- prefix
  my.data$num.cores <- num.cores
  my.data$verbose <- verbose

  # Return study setup ---------------------------------------------------------
  return(my.data)
}
