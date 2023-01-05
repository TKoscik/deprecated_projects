mri.hdr <- function(mri.file, field = "all") {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(file.exists(mri.file))
  
  fid <- file(mri.file, "rb")
  endian <- .Platform$endian
  
  if (field[1] == "all") {
    field <- c("hdr.size", "descrip",
               "n.voxels", "n.ts", "ts.length",
               "img.dims", "pixdim",
               "qform_code", "sform_code",
               "quatern_b", "quatern_c", "quatern_d", "qoffset_x", "qoffset_y",
               "qoffset_z", "srow_x", "srow_y", "srow_z", "intent_name", "magic")
  }
  
  hdr <- list()
  
  if ("hdr.size" %in% field) {
    invisible(seek(fid, 0, "start", "rb"))
    hdr$hdr.size <- readBin(fid, integer(), size=4, endian=endian)
  }
  
  if ("descrip" %in% field) {
    invisible(seek(fid, 4, "start", "rb"))
    txt <- readBin(fid, "raw", 16)
    hdr$descrip <- iconv(rawToChar(txt[txt!=as.raw(0)]), to="UTF-8")
  }
  
  if ("n.voxels" %in% field) {
    invisible(seek(fid, 20, "start", "rb"))
    hdr$n.voxels <- readBin(fid, double(), size = 8, endian = endian)
  }
  
  if ("n.ts" %in% field) {
    invisible(seek(fid, 28, "start", "rb"))
    hdr$n.ts <- readBin(fid, double(), size = 8, endian = endian)
  }
  
  if ("ts.length" %in% field) {
    invisible(seek(fid, 36, "start", "rb"))
    hdr$ts.length <- readBin(fid, double(), size = 8, endian = endian)
  }
  
  if ("img.dims" %in% field) {
    invisible(seek(fid, 44, "start", "rb"))
    hdr$img.dims <- readBin(fid, integer(), 8, size = 2, endian = endian)
  }
  
  if ("pixdim" %in% field) {
    invisible(seek(fid, 60, "start", "rb"))
    hdr$pixdim <- readBin(fid, numeric(), 8, size = 4, endian = endian)
  }
  bad_pixdim = !is.finite(hdr$pixdim)
  if (any(bad_pixdim)) { hdr$pixdim[bad_pixdim] = 1 }
  
  if ("qform_code" %in% field) {
    invisible(seek(fid, 92, "start", "rb"))
    hdr$qform_code <- readBin(fid, integer(), size = 2, endian = endian)
  }
  
  if ("sform_code" %in% field) {
    invisible(seek(fid, 94, "start", "rb"))
    hdr$sform_code <- readBin(fid, integer(), size = 2, endian = endian)
  }
  
  if ("quatern_b" %in% field) {
    invisible(seek(fid, 96, "start", "rb"))
    hdr$quatern_b <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("quatern_c" %in% field) {
    invisible(seek(fid, 100, "start", "rb"))
    hdr$quatern_c <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("quatern_d" %in% field) {
    invisible(seek(fid, 104, "start", "rb"))
    hdr$quatern_d <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("qoffset_x" %in% field) {
    invisible(seek(fid, 108, "start", "rb"))
    hdr$qoffset_x <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("qoffset_y" %in% field) {
    invisible(seek(fid, 112, "start", "rb"))
    hdr$qoffset_y <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("qoffset_z" %in% field) {
    invisible(seek(fid, 116, "start", "rb"))
    hdr$qoffset_z <- readBin(fid, numeric(), size = 4, endian = endian)
  }
  
  if ("srow_x" %in% field) {
    invisible(seek(fid, 120, "start", "rb"))
    hdr$srow_x <- readBin(fid, numeric(), 4, size = 4, endian = endian)
  }
  
  if ("srow_y" %in% field) {
    invisible(seek(fid, 136, "start", "rb"))
    hdr$srow_y <- readBin(fid, numeric(), 4, size = 4, endian = endian)
  }
  
  if ("srow_z" %in% field) {
    invisible(seek(fid, 152, "start", "rb"))
    hdr$srow_z <- readBin(fid, numeric(), 4, size = 4, endian = endian)
  }
  
  close(fid)
  return(hdr)
}