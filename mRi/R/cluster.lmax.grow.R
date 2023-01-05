cluster.lmax.grow <- function(tmap, tvol,
                         pmap, pvol,
                         mask, mvol,
                         alpha = 0.05,
                         lmax.radius = NULL,
                         min.size = 10,
                         min.size.action = "merge",
                         save.dir,
                         file.name = NULL,
                         verbose = FALSE) {

  luf <- tmap   # use tmap to look up header info

  # load maps
  tmap <- read.nii.volume(tmap, tvol)
  pmap <- read.nii.volume(pmap, pvol)
  mask <- read.nii.volume(mask, mvol)

  tmap[is.na(tmap)] <- 0
  pmap[is.na(pmap)] <- 1
  mask <- (mask > 0) * 1

  # get nifti info
  img.dims <- dim(tmap)
  vxl.size <- unlist(nii.hdr(luf, "pixdim"))
  orient <- nii.orient(luf)
  if (is.null(lmax.radius)) {
    lmax.radius <- mean(vxl.size[2:4])
  }

  # Create masks
  pmask <- (pmap < alpha) * mask
  pos.mask <- (tmap > 0) * mask
  neg.mask <- (tmap < 0) * mask

  # create maps
  pos.map <- tmap * pos.mask * pmask
  neg.map <- tmap * neg.mask * pmask * (-1)

  #
  neighborhood <- t(matrix(c(-1,-1,-1,-1,-1,0,-1,-1,1,-1,0,-1,-1,0,0,-1,0,1,-1,
                             1,-1,-1,1,0,-1,1,1,0,-1,-1,0,-1,0,0,-1,1,0,0,-1,0,
                             0,1,0,1,-1,0,1,0,0,1,1,1,-1,-1,1,-1,0,1,-1,1,1,0,
                             -1,1,0,0,1,0,1,1,1,-1,1,1,0,1,1,1), nrow=3))

  # get file names
  if (is.null(file.name)) {
    fname <- unlist(strsplit(luf, "[/]"))
    fname <- fname[(length(fname))]
    fname <- unlist(strsplit(fname, "[.]"))
    fname <- paste(fname[-length(fname)], collapse=".")
    fname <- paste0(save.dir, "/", fname, ".vol", tvol,
                    ".a", alpha*100, ".r", lmax.radius)
  } else {
    fname <- paste0(save.dir, "/", file.name)
  }

# FIND POSITIVE CLUSTERS
  if (max(range(pos.map, na.rm=TRUE)) != 0) {
    pos.exist <- TRUE
    which.vxls <- which(pos.map > 0, arr.ind = TRUE)
    pos.pk <- array(0, dim=dim(pos.map))

    # Find local maxima
    count <- 0
    for (i in 1:nrow(which.vxls)) {
      roi <- get.voxel.sphere(sphereRadius = lmax.radius,
                              voxelSz = vxl.size[2:4],
                              center = which.vxls[i, ],
                              volumeDims = img.dims)
      local.max <- max(pos.map[roi])
      cur.vxl <- matrix(which.vxls[i, ], ncol=3)
      if (pos.map[cur.vxl] == local.max) {
        count <- count + 1
        pos.pk[cur.vxl] <- count
      }
    }
    pos.lmax <- which(pos.pk > 0, arr.ind=T)
    pos.cl <- pos.pk

    no.grow <- (pos.map > 0) * 1
    no.grow[no.grow == 0] <- NA

    grow.ls <- pos.lmax
    while (nrow(grow.ls) > 0) {
      cur.vxl <- matrix(grow.ls[1, ], ncol=3)
      cur.cl <- pos.cl[cur.vxl]

      roi <- rbind(cur.vxl, cbind(neighborhood[ ,1] + cur.vxl[ ,1],
                                  neighborhood[ ,2] + cur.vxl[ ,2],
                                  neighborhood[ ,3] + cur.vxl[ ,3]))

      grow.mask <- no.grow[roi]

      if (any(is.na(grow.mask))) {
        roi <- matrix(roi[-which(is.na(grow.mask)), ], ncol=3)
      }
      which.cur <- 1

      roi.vals <- pos.map[roi]
      roi.vals <- (roi.vals - roi.vals[which.cur])
      grow.vals <- which(roi.vals <= 0)
      stop.vals <- which(roi.vals > 0)

      pos.cl[matrix(roi[grow.vals, ], ncol=3)] <- cur.cl

      no.grow[cur.vxl] <- NA
      if (length(stop.vals) > 0) {
        no.grow[matrix(roi[stop.vals, ], ncol=3)] <- NA
      }

      grow.ls <- matrix(grow.ls[-1, ], ncol=3)
      grow.vals <- grow.vals[-which.cur]
      grow.vals <- grow.vals[as.logical(no.grow[matrix(roi[grow.vals, ], ncol=3)])]
      if (length(grow.vals) > 0) {
        add.vxls <- matrix(roi[grow.vals, ], ncol=3)
        grow.ls <- unique(matrix(rbind(grow.ls, add.vxls), ncol=3))
      }
      if (verbose) {
        print(sprintf("grow.ls: %0.0f; remaining voxels: %0.0f", nrow(grow.ls), sum(no.grow, na.rm=T)))
      }
    }
    pos.cl <- pos.cl * mask
  } else { pos.exist <- FALSE }

# FIND NEGATIVE CLUSTERS
  if (max(range(neg.map, na.rm=TRUE)) != 0) {
    neg.exist <- TRUE
    which.vxls <- which(neg.map > 0, arr.ind = TRUE)
    neg.pk <- array(0, dim=dim(neg.map))

    # Find local maxima
    count <- 0
    for (i in 1:nrow(which.vxls)) {
      roi <- get.voxel.sphere(sphereRadius = lmax.radius,
                              voxelSz = vxl.size[2:4],
                              center = which.vxls[i, ],
                              volumeDims = img.dims)
      local.max <- max(neg.map[roi])
      cur.vxl <- matrix(which.vxls[i, ], ncol=3)
      if (neg.map[cur.vxl] == local.max) {
        count <- count + 1
        neg.pk[cur.vxl] <- count
      }
    }
    neg.lmax <- which(neg.pk > 0, arr.ind=T)
    neg.cl <- neg.pk

    no.grow <- (neg.map > 0) * 1
    no.grow[no.grow == 0] <- NA

    grow.ls <- neg.lmax
    while (nrow(grow.ls) > 0) {
      cur.vxl <- matrix(grow.ls[1, ], ncol=3)
      cur.cl <- neg.cl[cur.vxl]

      roi <- rbind(cur.vxl, cbind(neighborhood[ ,1] + cur.vxl[ ,1],
                                  neighborhood[ ,2] + cur.vxl[ ,2],
                                  neighborhood[ ,3] + cur.vxl[ ,3]))

      grow.mask <- no.grow[roi]

      if (any(is.na(grow.mask))) {
        roi <- matrix(roi[-which(is.na(grow.mask)), ], ncol=3)
      }
      which.cur <- 1

      roi.vals <- neg.map[roi]
      roi.vals <- (roi.vals - roi.vals[which.cur])
      grow.vals <- which(roi.vals <= 0)
      stop.vals <- which(roi.vals > 0)

      neg.cl[matrix(roi[grow.vals, ], ncol=3)] <- cur.cl

      no.grow[cur.vxl] <- NA
      if (length(stop.vals) > 0) {
        no.grow[matrix(roi[stop.vals, ], ncol=3)] <- NA
      }

      grow.ls <- matrix(grow.ls[-1, ], ncol=3)
      grow.vals <- grow.vals[-which.cur]
      grow.vals <- grow.vals[as.logical(no.grow[matrix(roi[grow.vals, ], ncol=3)])]
      if (length(grow.vals) > 0) {
        add.vxls <- matrix(roi[grow.vals, ], ncol=3)
        grow.ls <- unique(matrix(rbind(grow.ls, add.vxls), ncol=3))
      }
      if (verbose) {
        print(sprintf("grow.ls: %0.0f; remaining voxels: %0.0f", nrow(grow.ls), sum(no.grow, na.rm=T)))
      }
    }
    neg.cl <- neg.cl * mask
  } else { neg.exist <- FALSE}

  # save raw clusters
  if (pos.exist) {
    tname <- paste0(fname, ".pos.cluster.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, pos.cl)

    tname <- paste0(fname, ".pos.peak.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, pos.pk)
  }
  if (neg.exist) {
    tname <- paste0(fname, ".neg.cluster.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, neg.cl)

    tname <- paste0(fname, ".neg.peak.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, neg.pk)
  }

  # Summarize clusters
  if (pos.exist) {
    n.cl <- max(range(pos.pk))
    pos.df <- data.frame(cluster = 1:n.cl,
                        n.voxels = numeric(n.cl),
                        volume = numeric(n.cl),
                        peak.x = c(pos.lmax[ ,1]),
                        peak.y = c(pos.lmax[ ,2]),
                        peak.z = c(pos.lmax[ ,3]),
                        peak.x.mni = numeric(n.cl),
                        peak.y.mni = numeric(n.cl),
                        peak.z.mni = numeric(n.cl),
                        peak.value = c(pos.map[pos.lmax]),
                        cog.x = numeric(n.cl),
                        cog.y = numeric(n.cl),
                        cog.z = numeric(n.cl),
                        cog.x.mni = numeric(n.cl),
                        cog.y.mni = numeric(n.cl),
                        cog.z.mni = numeric(n.cl),
                        cog.value = numeric(n.cl))

    for (i in 1:max(range(pos.pk))) {
      pos.df$n.voxels[i] <- sum(pos.cl == i, na.rm=TRUE)
      pos.df$volume[i] <- pos.df$n.voxels[i] * prod(vxl.size[2:4])
    }
  }
  if (neg.exist) {
    n.cl <- max(range(neg.pk))
    neg.df <- data.frame(cluster = 1:n.cl,
                         n.voxels = numeric(n.cl),
                         volume = numeric(n.cl),
                         peak.x = c(neg.lmax[ ,1]),
                         peak.y = c(neg.lmax[ ,2]),
                         peak.z = c(neg.lmax[ ,3]),
                         peak.x.mni = numeric(n.cl),
                         peak.y.mni = numeric(n.cl),
                         peak.z.mni = numeric(n.cl),
                         peak.value = c(neg.map[neg.lmax]),
                         cog.x = numeric(n.cl),
                         cog.y = numeric(n.cl),
                         cog.z = numeric(n.cl),
                         cog.x.mni = numeric(n.cl),
                         cog.y.mni = numeric(n.cl),
                         cog.z.mni = numeric(n.cl),
                         cog.value = numeric(n.cl))

    for (i in 1:max(range(neg.pk))) {
      neg.df$n.voxels[i] <- sum(neg.cl == i, na.rm=TRUE)
      neg.df$volume[i] <- neg.df$n.voxels[i] * prod(vxl.size[2:4])
    }
  }

  # Calculate COG
  if (pos.exist) {
    for (i in 1:max(range(pos.pk))) {
      cl.temp <- matrix(which(pos.cl == i, arr.ind = TRUE), ncol=3)
      cl.vals <- pos.map[cl.temp]
      pos.df$cog.x[i] <- round(weighted.mean(cl.temp[ ,1], cl.vals))
      pos.df$cog.y[i] <- round(weighted.mean(cl.temp[ ,2], cl.vals))
      pos.df$cog.z[i] <- round(weighted.mean(cl.temp[ ,3], cl.vals))
    }
  }
  if (neg.exist) {
    for (i in 1:max(range(neg.pk))) {
      cl.temp <- matrix(which(neg.cl == i, arr.ind = TRUE), ncol=3)
      cl.vals <- neg.map[cl.temp]
      neg.df$cog.x[i] <- round(weighted.mean(cl.temp[ ,1], cl.vals))
      neg.df$cog.y[i] <- round(weighted.mean(cl.temp[ ,2], cl.vals))
      neg.df$cog.z[i] <- round(weighted.mean(cl.temp[ ,3], cl.vals))
    }
  }

  # Clean up small clusters
  if (min.size.action == "merge") {
    if (pos.exist) {
      zero.cl <- which(pos.df$n.voxels == 0)
      pos.cl[pos.cl %in% zero.cl] <- 0
      pos.df <- pos.df[-zero.cl, ]

      small <- TRUE
      while (small) {
        smallest <- which(pos.df$n.voxels == min(pos.df$n.voxels))[1]
        cog.dist <- sqrt((pos.df$cog.x[smallest] - pos.df$cog.x)^2 +
                           (pos.df$cog.y[smallest] - pos.df$cog.y)^2 +
                           (pos.df$cog.z[smallest] - pos.df$cog.z)^2)
        cog.dist[smallest] <- NA
        nearest <- which(cog.dist == min(cog.dist, na.rm=TRUE))[1]
        pos.cl[pos.cl == pos.df$cluster[smallest]] <- pos.df$cluster[nearest]
        pos.df$n.voxels[nearest] <- pos.df$n.voxels[nearest] + pos.df$n.voxels[smallest]
        pos.df$volume[nearest] <- pos.df$volume[nearest] + pos.df$volume[smallest]
        if (pos.df$peak.value[nearest] < pos.df$peak.value[smallest]) {
          pos.df$peak.x[nearest] <- pos.df$peak.x[smallest]
          pos.df$peak.y[nearest] <- pos.df$peak.y[smallest]
          pos.df$peak.z[nearest] <- pos.df$peak.z[smallest]
        }
        pos.df <- pos.df[-smallest, ]
        if (min(pos.df$n.voxels) >= min.size) {
          small <- FALSE
        }
      }
    }
    if (neg.exist) {
      zero.cl <- which(neg.df$n.voxels == 0)
      neg.cl[neg.cl %in% zero.cl] <- 0
      neg.df <- neg.df[-zero.cl, ]

      small <- TRUE
      while (small) {
        smallest <- which(neg.df$n.voxels == min(neg.df$n.voxels))[1]
        cog.dist <- sqrt((neg.df$cog.x[smallest] - neg.df$cog.x)^2 +
                           (neg.df$cog.y[smallest] - neg.df$cog.y)^2 +
                           (neg.df$cog.z[smallest] - neg.df$cog.z)^2)
        cog.dist[smallest] <- NA
        nearest <- which(cog.dist == min(cog.dist, na.rm=TRUE))[1]
        neg.cl[neg.cl == neg.df$cluster[smallest]] <- neg.df$cluster[nearest]
        neg.df$n.voxels[nearest] <- neg.df$n.voxels[nearest] + neg.df$n.voxels[smallest]
        neg.df$volume[nearest] <- neg.df$volume[nearest] + neg.df$volume[smallest]
        if (neg.df$peak.value[nearest] < neg.df$peak.value[smallest]) {
          neg.df$peak.x[nearest] <- neg.df$peak.x[smallest]
          neg.df$peak.y[nearest] <- neg.df$peak.y[smallest]
          neg.df$peak.z[nearest] <- neg.df$peak.z[smallest]
        }
        neg.df <- neg.df[-smallest, ]
        if (min(neg.df$n.voxels) >= min.size) {
          small <- FALSE
        }
      }
    }
  } else if (min.size.action == "delete") {
    if (pos.exist) {
      small.cl <- which(pos.df$n.voxels < min.size)
      if (length(small.cl) > 0) {
        pos.cl[pos.cl %in% small.cl] <- 0
        pos.df <- pos.df[-small.cl]
      }
    }
    if (neg.exist) {
      small.cl <- which(neg.df$n.voxels < min.size)
      if (length(small.cl) > 0) {
        neg.cl[neg.cl %in% small.cl] <- 0
        neg.df <- neg.df[-small.cl]
      }
    }
  }

  # Renumber Clusters and sort by and size
  if (pos.exist) {
    pos.df <- pos.df[order(pos.df$n.voxels, decreasing = TRUE), ]
    pos.cl <- pos.cl * (-1)
    for (i in 1:nrow(pos.df)) {
      pos.cl[pos.cl == (pos.df$cluster[i] * (-1))] <- i
      pos.df$cluster[i] <- i
    }
  }
  if (neg.exist) {
    neg.df <- neg.df[order(neg.df$n.voxels, decreasing = TRUE), ]
    neg.cl <- neg.cl * (-1)
    for (i in 1:nrow(neg.df)) {
      neg.cl[neg.cl == (neg.df$cluster[i] * (-1))] <- i
      neg.df$cluster[i] <- i
    }
  }

  # Re-Calculate COG
  if (pos.exist) {
    for (i in 1:max(range(pos.cl))) {
      cl.temp <- matrix(which(pos.cl == i, arr.ind = TRUE), ncol=3)
      cl.vals <- pos.map[cl.temp]
      pos.df$cog.x[i] <- round(weighted.mean(cl.temp[ ,1], cl.vals))
      pos.df$cog.y[i] <- round(weighted.mean(cl.temp[ ,2], cl.vals))
      pos.df$cog.z[i] <- round(weighted.mean(cl.temp[ ,3], cl.vals))
    }
  }
  if (neg.exist) {
    for (i in 1:max(range(neg.cl))) {
      cl.temp <- matrix(which(neg.cl == i, arr.ind = TRUE), ncol=3)
      cl.vals <- neg.map[cl.temp]
      neg.df$cog.x[i] <- round(weighted.mean(cl.temp[ ,1], cl.vals))
      neg.df$cog.y[i] <- round(weighted.mean(cl.temp[ ,2], cl.vals))
      neg.df$cog.z[i] <- round(weighted.mean(cl.temp[ ,3], cl.vals))
    }
    neg.df$cog.value <- neg.map[matrix(unlist(neg.df[ ,11:13]), ncol=3)]
  }

  # Convert Coordinates to MNI
  tform <- nii.hdr(nii.file = luf,
                   field = c("srow_x", "srow_y", "srow_z"))
  tform <- rbind(tform$srow_x, tform$srow_y, tform$srow_z, c(0,0,0,1))
  if (pos.exist) {
    pos.df[ ,7:9]<- t(tform %*% rbind(pos.df$peak.x-1,
                                      pos.df$peak.y-1,
                                      pos.df$peak.z-1,
                                      rep(1, nrow(pos.df))))[ ,-4]
    pos.df[ ,14:16]<- t(tform %*% rbind(pos.df$cog.x-1,
                                        pos.df$cog.y-1,
                                        pos.df$cog.z-1,
                                        rep(1, nrow(pos.df))))[ ,-4]
  }
  if (neg.exist) {
    neg.df[ ,7:9]<- t(tform %*% rbind(neg.df$peak.x-1,
                                      neg.df$peak.y-1,
                                      neg.df$peak.z-1,
                                      rep(1, nrow(neg.df))))[ ,-4]
    neg.df[ ,14:16]<- t(tform %*% rbind(neg.df$cog.x-1,
                                        neg.df$cog.y-1,
                                        neg.df$cog.z-1,
                                        rep(1, nrow(neg.df))))[ ,-4]
  }

  # Save output
  if (pos.exist) {
    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".pos.cluster.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, pos.cl)

    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".pos.peak.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, pos.pk)

    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".pos.cluster.csv")
    write.table(pos.df, file = tname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep=",")
  }
  if (neg.exist) {
    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".neg.cluster.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, neg.cl)

    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".neg.peak.nii")
    init.nii(tname, dims = img.dims, pixdim = vxl.size, orient = orient)
    write.nii.volume(tname, 1, neg.pk)

    tname <- paste0(fname, ".sz", min.size, ".", min.size.action, ".neg.cluster.csv")
    write.table(neg.df, file = tname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep=",")
  }
}
