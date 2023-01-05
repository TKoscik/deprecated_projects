cluster.bilateral <- function(nii.file, nii.vol, radius) {

  # DEBUG ----
  nii.file <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/supplementals/fdr.maps/badEV.nii"
  nii.vol <- 1
  radius <- 8
  spatial.scalar <- 1
  value.scalar <- 1
  # ----

  in.vol <- read.nii.volume(nii.file, nii.vol)
  out.vol <- array(0, dim=dim(in.vol))

  nii.info <- nii.hdr(nii.file, c("dim", "pixdim"))
  dir.mx <- get.voxel.sphere(sphereRadius = radius,
                             voxelSz = nii.info$pixdim[2:4],
                             center = c(0,0,0),
                             volumeDims = nii.info$dim[2:4])
  x.max <- max(dir.mx[ ,1])
  y.max <- max(dir.mx[ ,2])
  z.max <- max(dir.mx[ ,3])

  pad.vol <- array(0, dim=nii.info$dim[2:4]+c(x.max,y.max,z.max)*2)
  ctr.x <- seq(from = 1+x.max, by = 1, length.out = nii.info$dim[2])
  ctr.y <- seq(from = 1+y.max, by = 1, length.out = nii.info$dim[3])
  ctr.z <- seq(from = 1+z.max, by = 1, length.out = nii.info$dim[4])
  pad.vol[ctr.x, ctr.y, ctr.z] <- in.vol

  for (i in 1:nrow(dir.mx)) {
    tvol <- array(0, dim=nii.info$dim[2:4]+c(x.max,y.max,z.max)*2)
    idx.x <- seq(from = 1+x.max, by = 1, length.out = nii.info$dim[2]) + dir.mx[i,1]
    idx.y <- seq(from = 1+y.max, by = 1, length.out = nii.info$dim[3]) + dir.mx[i,2]
    idx.z <- seq(from = 1+z.max, by = 1, length.out = nii.info$dim[4]) + dir.mx[i,3]

    tvol[idx.x, idx.y, idx.z] <- in.vol
    if (sum(dir.mx[i, ]) != 0) {
    svol <- tvol / (sqrt(sum(dir.mx[i, ]^2)) * spatial.scalar)
    }
    vvol <- abs(tvol - pad.vol) * value.scalar

    tvol <- (svol * vvol) / nrow(dir.mx)
    tvol <- tvol[ctr.x, ctr.y, ctr.z]
    if (any(is.nan(tvol))) { stop("crap")}
    out.vol <- out.vol + tvol
  }
    # out.vol <- out.vol / nrow(dir.mx)
    out.vol <- 1 - out.vol
    init.nii(file.name = "C:/Users/tim.koscik/Documents/Data/bi.filt.test.nii", dims = c(91,109,91))
    write.nii.volume(nii.file = "C:/Users/tim.koscik/Documents/Data/bi.filt.test.nii",
                     vol.num = 1, values = out.vol)
}
