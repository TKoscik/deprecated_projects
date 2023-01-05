mri.load <- function(mri.file, vxl.num) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # Check inputs -----------------------------------------------------------------------
  stopifnot(!missing(mri.file), file.exists(mri.file), !missing(vxl.num))
  
  # Get necessary NII file info --------------------------------------------------------
  hdr <- mri.hdr(mri.file, field = "all")
  
  # Check if Voxel number is in range --------------------------------------------------
  stopifnot(vxl.num <= hdr$n.voxels)
  
  # Initialize NII file for binary reading ---------------------------------------------
  fid <- file(mri.file, "rb")
  endian <- .Platform$endian
  
  # Get data ---------------------------------------------------------------------------
  seek(fid, where=(hdr$hdr.size + (vxl.num - 1) * hdr$ts.length), origin = "start")
  data.temp <- readBin(fid, double(), hdr$ts.length+3, 8, endian=endian)
  
  close(fid)
  
  # Format data
  data <- list()
  data$coords <- data.temp[1:3]
  data$img.dims <- hdr$img.dims
  data$pixdim <- hdr$pixdim
  data$orient <- list()
  data$orient$qform_code <- hdr$qform_code
  data$orient$sform_code <- hdr$sform_code
  data$orient$quatern_b <- hdr$quatern_b
  data$orient$quatern_c <- hdr$quatern_c
  data$orient$quatern_d <- hdr$quatern_d
  data$orient$qoffset_x <- hdr$qoffset_x
  data$orient$qoffset_y <- hdr$qoffset_y
  data$orient$qoffset_z <- hdr$qoffset_z
  data$orient$srow_x <- hdr$srow_x
  data$orient$srow_y <- hdr$srow_y
  data$orient$srow_z <- hdr$srow_z
  
  data$fmri <- matrix(0, nrow=hdr$ts.length, ncol=hdr$n.ts)
  if (hdr$n.ts == 1) {
    colnames(data$fmri) <- "fmri"
  } else {
    colnames(data$fmri) <- paste0("fmri", 1:hdr$n.ts)
  }
  for (i in 1:hdr$n.ts) {
    data$fmri[ ,i] <- data.temp[(4+(i-1)*hdr$ts.length):(3+i*hdr$ts.length)]
  }
  
  return(data)
}