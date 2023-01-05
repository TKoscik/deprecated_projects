mri.make <- function(data.4d, data.mask, save.dir, file.name) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # Parse Input ------------------------------------------------------------------------
  data.4d <- parse.4d(data.4d)
  fname <- paste0(save.dir, "/", file.name, ".mri")
  if (file.exists(fname)) {
    stop("File exists, please choos another filename and/or location")
  }
  
  # Gather information before writing to file ------------------------------------------
  img.dims <- unlist(nii.hdr(data.4d[[1]][1], "dim"))
  pixdim <- unlist(nii.hdr(data.4d[[1]][1], "pixdim"))
  orient <- nii.orient(data.4d[[1]][1])
  
  n.ts <- length(data.4d)
  ts.length <- length(load.voxel(data.4d, c(1,1,1)))
  
  which.voxels <- load.mask(data.mask)
  n.masks <- length(which.voxels)
  n.voxels <- 0
  for (i in 1:n.masks) {
    n.voxels <- n.voxels + nrow(which.voxels[[i]])
  }
  
  # Initialize File for Writing --------------------------------------------------------
  fid <- file(fname, "w+b")
  
  writeBin(168L, fid, size = 4) # Header length
  suppressWarnings(writeChar("mRi Study File", fid, nchars = 16, eos = NULL)) # File Description
  writeBin(as.double(n.voxels), fid, size = 8)  # Number of voxels
  writeBin(as.double(n.ts), fid, size = 8)  # Number of timeseries
  writeBin(as.double(ts.length), fid, size = 8)# Length of concatenated time series
  writeBin(as.integer(img.dims), fid, size = 2) # NII dimensions for saving NII output
  writeBin(as.double(pixdim), fid, size = 4)    # NII pixdim for saving NII Output
  writeBin(as.integer(orient$qform_code), fid, size = 2) # NII orient for daving NII output
  writeBin(as.integer(orient$sform_code), fid, size = 2)
  writeBin(as.double(orient$quatern_b), fid, size = 4)
  writeBin(as.double(orient$quatern_c), fid, size = 4)
  writeBin(as.double(orient$quatern_d), fid, size = 4)
  writeBin(as.double(orient$qoffset_x), fid, size = 4)
  writeBin(as.double(orient$qoffset_y), fid, size = 4)
  writeBin(as.double(orient$qoffset_z), fid, size = 4)
  writeBin(as.double(orient$srow_x), fid, size = 4)
  writeBin(as.double(orient$srow_y), fid, size = 4)
  writeBin(as.double(orient$srow_z), fid, size = 4)
  
  for (i in 1:n.masks) {
    n.voxels <- nrow(which.voxels[[i]])
    for (j in 1:n.voxels) {
      fmri <- load.voxel(data.4d, which.voxels[[i]][j, ])
      data <- c(which.voxels[[i]][j, ], fmri)
      writeBin(as.double(data), fid, size=8)
      print(paste(j, "of", n.voxels, sep=" "))
    }
  }
  
  close(fid)
}