cluster.lmax <- function(tmap, tvol,
                         pmap, pvol,
                         mask, mvol,
                         alpha = 0.05,
                         radius = 5,
                         dist.weight = 1,
                         val.weight = 1,
                         save.dir,
                         file.name = NULL) {

  luf <- tmap

  # load maps
  tmap <- read.nii.volume(tmap, tvol)
  pmap <- read.nii.volume(pmap, pvol)
  mask <- read.nii.volume(mask, mvol)

  tmap[is.na(tmap)] <- 0

  # get nifti info
  img.dims <- dim(tmap)
  vxl.size <- unlist(nii.hdr(luf, "pixdim"))
  orient <- nii.orient(luf)

  # Create masks
  pmask <- (pmap < alpha) * mask
  pos.mask <- (tmap > 0) * mask
  neg.mask <- (tmap < 0) * mask

  # create maps
  pos.map <- tmap * pos.mask * pmask
  neg.map <- tmap * neg.mask * pmask * -1

  # get file names
  fname <- unlist(strsplit(luf, "[/]"))
  fname <- fname[(length(fname))]
  fname <- unlist(strsplit(fname, "[.]"))
  fname <- paste(fname[-length(fname)], collapse=".")
  if (is.null(file.name)) {
    fname <- paste0(save.dir, "/", fname, ".vol", tvol,
                    ".a", alpha*100, ".r", radius,
                    ".d", dist.weight, ".v", val.weight)
  } else {
    fname <- paste0(save.dir, "/", file.name)
  }
  fname <- c(paste0(fname, ".pos.clusters.nii"), paste0(fname, ".pos.peaks.nii"),
             paste0(fname, ".neg.clusters.nii"), paste0(fname, ".neg.peaks.nii"))


  if (max(range(pos.map)) != 0) {
    which.vxls <- which(pos.map > 0, arr.ind = TRUE)
    cl.out <- array(0, dim=dim(pos.map))
    pk.out <- array(0, dim=dim(pos.map))
    count <- 0
    for (i in 1:nrow(which.vxls)) {
      roi <- get.voxel.sphere(sphereRadius = radius, voxelSz = vxl.size[2:4], center = which.vxls[i, ], volumeDims = img.dims)
      local.max <- max(pos.map[roi])
      if (pos.map[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] == local.max) {
        count <- count + 1
        pk.out[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] <- count
      }
    }
    lmax.ls <- which(pk.out > 0, arr.ind=T)
    for (i in 1:nrow(which.vxls)) {
      dist.vec <- sqrt(rowSums(cbind((which.vxls[i,1] - lmax.ls[ ,1])^2,
                                     (which.vxls[i,2] - lmax.ls[ ,2])^2,
                                     (which.vxls[i,3] - lmax.ls[ ,3])^2)))
      val.vec <- abs(pos.map[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] - pos.map[lmax.ls])
      # val.vec = 1
      w.factor <- sqrt(sum(dist.weight * dist.vec^2 + val.weight * val.vec^2))
      cl.out[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] <- which(w.factor == min(w.factor))[1]
    }

    init.nii(fname[1], dims = img.dims, pixdim = vxl.size, orient = orient)
    init.nii(fname[2], dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(fname[1], 1, cl.out)
    write.nii.volume(fname[2], 1, pk.out)
  }

  if (max(range(neg.map)) != 0) {
    which.vxls <- which(neg.map > 0, arr.ind = TRUE)
    cl.out <- array(0, dim=dim(neg.map))
    pk.out <- array(0, dim=dim(neg.map))
    count <- 0
    for (i in 1:nrow(which.vxls)) {
      roi <- get.voxel.sphere(sphereRadius = radius, voxelSz = vxl.size[2:4], center = which.vxls[i, ], volumeDims = img.dims)
      local.max <- max(neg.map[roi])
      if (neg.map[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] == local.max) {
        count <- count + 1
        pk.out[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] <- count
      }
    }
    lmax.ls <- which(pk.out > 0, arr.ind=T)
    for (i in 1:nrow(which.vxls)) {
      dist.vec <- sqrt(rowSums(cbind((which.vxls[i,1] - lmax.ls[ ,1])^2,
                                     (which.vxls[i,2] - lmax.ls[ ,2])^2,
                                     (which.vxls[i,3] - lmax.ls[ ,3])^2)))
      val.vec <- abs(neg.map[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] - neg.map[lmax.ls])
      # val.vec = 1
      w.factor <- sqrt(sum(dist.weight * dist.vec^2 + val.weight * val.vec^2))
      cl.out[which.vxls[i,1], which.vxls[i,2], which.vxls[i,3]] <- which(w.factor == min(w.factor))[1]
    }

    init.nii(fname[3], dims = img.dims, pixdim = vxl.size, orient = orient)
    init.nii(fname[4], dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(fname[3], 1, cl.out)
    write.nii.volume(fname[4], 1, pk.out)
  }
}
