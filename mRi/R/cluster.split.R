cluster.split <- function(cluster.nii,
                          cluster.vol,
                          value.nii,
                          value.vol,
                          peaks = TRUE,
                          min.peak.distance = 5,
                          split.if.bigger = 1000,
                          save.dir) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------

  # Debug ----
  # rm(list=ls())
  # gc()
  # library(mRi)
  # cluster.nii <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/supplementals/clusters/goodOutcome.neg.nii"
  # cluster.vol <- "all"
  # value.nii <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/supplementals/fdr.maps/goodOutcome.nii"
  # value.vol <- 1
  # peaks <- FALSE
  # min.peak.distance = 10
  # save.dir = "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/analyses/supplementals/clusters.split"
  # ----

  # Check cluster volumes ----
  n.clusters <- nii.dims(cluster.nii)[4]
  if (is.numeric(cluster.vol)) {
    stopifnot(all(cluster.vol %in% 1:n.clusters))
    n.clusters <- length(cluster.vol)
  } else {
    cluster.vol <- 1:n.clusters
  }

  cl.ls <- load.mask(cluster.nii)
  nii.info <- nii.hdr(cluster.nii, c("dim", "pixdim"))
  orient <- nii.orient(cluster.nii)

  value.nii <- read.nii.volume(value.nii, value.vol)

  new.clusters <- numeric(n.clusters)
  init.nii(file.name = paste0(save.dir, "/temp.nii"),
           dims = nii.info$dim[2:5],
           pixdim = nii.info$pixdim,
           orient = orient)

  for (i in 1:n.clusters) {
    c.map <- array(0, dim = nii.info$dim[2:4])
    c.map[cl.ls[[i]]] <- 1
    if (nrow(cl.ls[[i]]) >= split.if.bigger) {
      # find local maxima
      for (j in 1:nrow(cl.ls[[i]])) {
        print(sprintf("%0.0f of %0.0f", j, nrow(cl.ls[[i]])))

        l.map <- array(0, dim = nii.info$dim[2:4])

        # c.map[cl.ls[[i]][j,1], cl.ls[[i]][j,2], cl.ls[[i]][j,3]] <- 1
        locs <- get.voxel.sphere(sphereRadius = min.peak.distance,
                                 voxelSz = c(2,2,2), # nii.info$pixdim[2:4]
                                 center = cl.ls[[i]][j, ],
                                 volumeDims = nii.info$dim[2:4])
        chk.locs <- c(which(locs[ ,1] > nii.info$dim[2],1),
                      which(locs[ ,2] > nii.info$dim[3],2),
                      which(locs[ ,3] > nii.info$dim[4],3))
        if (length(chk.locs) > 0) {
          locs <- locs[-chk.locs]
        }

        l.map[locs] <- 1
        l.map <- l.map * c.map
        in.cl <- which(l.map == 1, arr.ind = TRUE)

        # in.cl <- numeric(nrow(locs))
        # for (k in 1:nrow(locs)) {
        #   in.cl[k] <- sum((locs[k,1] == cl.ls[[i]][ ,1]) *
        #       (locs[k,2] == cl.ls[[i]][ ,2]) *
        #       (locs[k,3] == cl.ls[[i]][ ,3])) != 0
        # }

        # vals <- value.nii[locs[which(in.cl==1),]]
        vals <- value.nii[in.cl]
        ctr <- value.nii[cl.ls[[i]][j,1], cl.ls[[i]][j,2], cl.ls[[i]][j,3]]

        if (peaks) {
          if (ctr != max(vals)) {
            c.map[cl.ls[[i]][j,1], cl.ls[[i]][j,2], cl.ls[[i]][j,3]] <- 0
          }
        } else {
          if (ctr != min(vals)) {
            c.map[cl.ls[[i]][j,1], cl.ls[[i]][j,2], cl.ls[[i]][j,3]] <- 0
          }
        }
      }

      # assign to nearest maxima
      pks <- which(c.map == 1, arr.ind = TRUE)
      new.clusters[i] <- nrow(pks)
      for (j in 1:nrow(cl.ls[[i]])) {
        print(sprintf("%0.0f of %0.0f", j, nrow(cl.ls[[i]])))
        d2peak <- unname(as.matrix(dist(x = rbind(cl.ls[[i]][j, ], pks)))[-1,1])
        c.map[cl.ls[[i]][j,1], cl.ls[[i]][j,2], cl.ls[[i]][j,3]] <-
          which(d2peak == min(d2peak))[1]
      }
    }
    write.nii.volume(paste0(save.dir, "/temp.nii"), i, c.map)
  }
  fname <- unlist(strsplit(cluster.nii, split = "/"))
  fname <- unlist(strsplit(fname[length(fname)], split = "[.]"))
  fname <- paste(c(fname[-length(fname)], "split.nii"), collapse =".")

  init.nii(file.name = paste0(save.dir, "/", fname),
           dims = c(nii.info$dim[2:4],sum(new.clusters)),
           pixdim = nii.info$pixdim,
           orient = orient)
  count <- 0
  for (i in 1:n.clusters) {
    tvol <- read.nii.volume(paste0(save.dir, "/temp.nii"), i)
    for (j in 1:new.clusters[i]) {
      count <- count + 1
      new.vol <- (tvol == j)*1
      write.nii.volume(paste0(save.dir, "/", fname), count, new.vol)
    }
  }
  file.remove(paste0(save.dir, "/temp.nii"))
}
