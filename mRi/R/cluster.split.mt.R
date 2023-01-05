cluster.split.mt <- function(data.cluster,
                             which.cluster = "all",
                             data.value,
                             which.value,
                             mt.min = 0.25,
                             mt.max = 0.95,
                             cluster.sz = 10,
                             thresh.dir = "gt",
                             connectivity = 26L,
                             save.dir,
                             file.name = NULL) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2018 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  # Debug
  rm(list=ls())
  gc()
  library("mRi")
  library("tk.r.misc")
  library("raster")
  library("gridExtra")
  library("fastcluster")
  library("dynamicTreeCut")
  library("ggdendro")
  home <- "/Shared/nopoulos/resting_state/varianceFiles"
  tdir <- paste0(home, "/analysis.20180529")
  i=2
  data.cluster <- paste0(tdir, "/cluster/kidsHD.var.mask1.mdl3.coef.FDR.vol", i, ".cl26.p0.05.sz25.neg.nii")
  which.cluster = "all"
  data.value = paste0(tdir, "/kidsHD.var.mask1.mdl3.coef.tvalue.nii")
  which.value = i
  mt.min = 0.25
  cluster.sz = 10
  thresh.dir = "lt"
  connectivity = 26L
  save.dir = paste0(tdir, "/cluster")
  file.name = NULL
  #

  stopifnot(file.exists(data.cluster), file.exists(data.value))
  img.dims <- nii.dims(data.cluster)[1:3]
  n.clusters <- nii.dims(data.cluster)[4]

  if(is.character(which.cluster) & which.cluster == "all") {
    cl.split <- rep(1, n.clusters)
  } else {
    cl.split <- rep(0, n.clusters)
    cl.split[which.cluster] <- 1
  }

  # get file names
  if (is.null(file.name)) {
    fname <- unlist(strsplit(data.cluster, "[/]"))
    fname <- fname[(length(fname))]
    fname <- unlist(strsplit(fname, "[.]"))
    fname <- paste(fname[-length(fname)], collapse=".")
    fname <- paste0(save.dir, "/", fname,
                    ".mtmin", mt.min,
                    ".", thresh.dir,
                    ".sz", cluster.sz,
                    ".mtsplit.nii")
  } else {
    fname <- paste0(save.dir, "/", file.name)
  }

  # Initialize connectivity neighbourhood ----
  if (is.null(connectivity)) {connectivity <- 26}
  n <- 3
  dimorder <- 1:3
  cdim <- cumprod(c(1, img.dims[dimorder][-n]))
  neighborhood <- switch(
    as.character(connectivity),
    `6` = t(matrix(c(1,2,2,2,1,2,2,2,1,3,2,2,2,3,2,2,2,3), nrow=3)),
    `18` = t(matrix(
      c(1,1,2,1,2,1,1,2,2,1,2,3,1,3,2,2,1,1,2,1,2,2,1,3,2,2,1,2,2,2,
        2,2,3,2,3,1,2,3,2,2,3,3,3,1,2,3,2,1,3,2,2,3,2,3,3,3,2), nrow=3)),
    `26` = t(matrix(
      c(1,1,1,1,1,2,1,1,3,1,2,1,1,2,2,1,2,3,1,3,1,1,3,2,1,3,3,2,1,1,
        2,1,2,2,1,3,2,2,1,2,2,3,2,3,1,2,3,2,2,3,3,3,1,1,3,1,2,3,1,3,
        3,2,1,3,2,2,3,2,3,3,3,1,3,3,2,3,3,3), nrow=3)),
    stop("Unrecognized connectivity"))
  center.pt <- t(matrix(c(2,2,2), nrow=3))
  offsets <- as.integer(colSums(t(neighborhood[ ,1:3, drop=FALSE]-1)*cdim) + 1L) -
    as.integer(colSums(t(center.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)

  neighbor6 <- t(matrix(c(1,2,2,2,1,2,2,2,1,3,2,2,2,3,2,2,2,3), nrow=3))
  offsets6 <- as.integer(colSums(t(neighbor6[ ,1:3, drop=FALSE]-1)*cdim) + 1L) -
    as.integer(colSums(t(center.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)

  new.clusters <- vector("list", 0)
  cl.count <- 0
  cl.vxls <- load.mask(data.cluster)
  cl.value <- read.nii.volume(data.value, which.value)

  thresh <- seq(mt.min, mt.max, 0.01)

  # dbscan eps
  eps.value <- sqrt(sum(unlist(nii.hdr(data.cluster, "pixdim"))[2:4]^2))

  cl.start.val <- 0
  init.nii(fname, c(img.dims),
           pixdim = unlist(nii.hdr(data.cluster, "pixdim")),
           orient = nii.orient(data.cluster))
  for (i in 1:n.clusters) {
    cl.orig <- read.nii.volume(data.cluster, i)

    value.map <- cl.value * cl.orig
    value.ls <- cl.value[cl.vxls[[i]]]
    value.q <- quantile(value.ls, thresh)
    if (thresh.dir == "lt") {
      value.q <- rev(value.q)
    }

    cl.num <- numeric(length(value.q))
    mt.connected <- array(0, dim=c(img.dims, length(value.q)))
    for (j in 1:length(thresh)) {
      volume <- switch(thresh.dir,
                       `gt`=array((value.map > value.q[j])*1, dim=img.dims[1:3]),
                       `lt`=array((value.map < value.q[j])*1, dim=img.dims[1:3]),
                       error("Cannot parse threshold direction. either gt or lt"))
      volume[is.na(volume)] <- 0
      connected <- array(0, dim=img.dims[1:3])

      num.clusters <- 0
      idx <- numeric()
      for (x in 1:img.dims[1]) {
        for (y in 1:img.dims[2]) {
          for (z in 1:img.dims[3]) {
            if (volume[x,y,z] == 1) {
              num.clusters <- num.clusters + 1
              current.pt <- t(as.matrix(c(x,y,z)))
              idx <- as.integer(colSums(t(current.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
              connected[idx] <- num.clusters
              while (length(idx)!=0) {
                volume[idx] <- 0
                neighbors = as.vector(apply(as.matrix(idx), 1, '+', offsets))
                neighbors = unique(neighbors[which(neighbors > 0)])
                idx = neighbors[which(volume[neighbors]!=0)]
                connected[idx] <- num.clusters
              }
            }
          }
        }
      }
      mt.connected[ , , ,j] <-  connected
      cl.num[j] <- num.clusters
    }

    # get matrix of cluster assignments per voxel
    cl.assign <- matrix(0, nrow=nrow(cl.vxls[[i]]), ncol=length(thresh))
    for (j in 1:nrow(cl.vxls[[i]])) {
      cl.assign[j, ] <- mt.connected[cl.vxls[[i]][j,1], cl.vxls[[i]][j,2], cl.vxls[[i]][j,3], ]
    }
    # cl.assign[cl.assign ==  0] <- NA
    cl.na <- which(rowSums((cl.assign == 0)) == ncol(cl.assign))
    cl.cl <- which(rowSums((cl.assign == 0)) != ncol(cl.assign))


    ctrs <- matrix(0, nrow=0, ncol=3)
    ###
    if (length(cl.cl) > 10000) {
      sub.sample <- sample(cl.cl, 10000, replace=FALSE)
    } else {
      sub.sample <- cl.cl
    }

    cl.assign.sub <- cl.assign[sub.sample, ]
    cl.dist <- matrix(0, nrow=nrow(cl.assign.sub), ncol=nrow(cl.assign.sub))
    for (j in 1:nrow(cl.assign.sub)) {
      for (k in 1:(j-1)) {
        cl.dist[j,k] <- (1 - sum(cl.assign.sub[j, ] == cl.assign.sub[k, ], na.rm=TRUE) / ncol(cl.assign.sub))^2
        cl.dist[k,j] <- cl.dist[j,k]
      }
      print(c(j,k))
    }
    cl.dist.obj <- as.dist(cl.dist)
    cl.hclust <- hclust(cl.dist.obj)
    cl.cut <- unname(cutreeHybrid(cl.hclust, dist=as.matrix(cl.dist.obj), minClusterSize=2, deepSplit=0)$labels)

    # get centroids to initialize k-means
    for (j in 1:max(cl.cut)) {
      tidx <- cl.vxls[[i]][sub.sample]
      ctrs <- rbind(ctrs, matrix)
    }

    ###

    # cl.dist <- matrix(0, nrow=nrow(cl.assign), ncol=nrow(cl.assign))
    # for (j in 1:nrow(cl.assign)) {
    #   for (k in 1:(j-1)) {
    #       cl.dist[j,k] <- (1 - sum(cl.assign[j, ] == cl.assign[k, ], na.rm=TRUE) / ncol(cl.assign))^2
    #       cl.dist[k,j] <- cl.dist[j,k]
    #   }
    # }
    # cl.dist[upper.tri(cl.dist)] <- t(cl.dist)[upper.tri(cl.dist)]
    # cl.dist <- as.dist(cl.dist)

    sub.sample <- sample(1:nrow(cl.assign), 10000, replace=FALSE)
    cl.dist.temp <- as.dist(cl.dist[sub.sample, sub.sample])
    cl.hclust <- hclust(cl.dist.temp)
    cl.cut <- unname(cutreeHybrid(cl.hclust, dist=as.matrix(cl.dist.temp), minClusterSize=2, deepSplit=0)$labels)



    # cl.hclust <- hclust(cl.dist)
    # cl.cut <- unname(cutreeHybrid(cl.hclust, dist=as.matrix(cl.dist), minClusterSize=2, deepSplit=0)$labels)
    cl.vxls[[i]] <- cbind(cl.vxls[[i]], cl.cut)
    cl.count <- cl.count + max(cl.cut)
  }

  init.nii(fname, c(img.dims),
           pixdim = unlist(nii.hdr(data.cluster, "pixdim")),
           orient = nii.orient(data.cluster))

  unlabelled <- array(0, dim=img.dims)
  loop.count <- 0
  out.nii <- array(0, dim=img.dims)
  for (i in 1:length(cl.vxls)) {
    for (j in 1:max(cl.vxls[[i]][ ,4])) {
      new.cluster <- array(0, dim=img.dims)
      list.o.vxls <- cl.vxls[[i]][cl.vxls[[i]][ ,4]==j, ]
      list.o.vxls <- list.o.vxls[ ,-4]
      new.cluster[matrix(list.o.vxls, ncol=3)] <- 1

      connected <- array(0, dim=img.dims[1:3])
      num.clusters <- 0
      idx <- numeric()
      for (x in 1:img.dims[1]) {
        for (y in 1:img.dims[2]) {
          for (z in 1:img.dims[3]) {
            if (new.cluster[x,y,z] == 1) {
              num.clusters <- num.clusters + 1
              current.pt <- t(as.matrix(c(x,y,z)))
              idx <- as.integer(colSums(t(current.pt[ ,1:3, drop=FALSE]-1)*cdim) + 1L)
              connected[idx] <- num.clusters
              while (length(idx)!=0) {
                new.cluster[idx] <- 0
                neighbors = as.vector(apply(as.matrix(idx), 1, '+', offsets6))
                neighbors = unique(neighbors[which(neighbors > 0)])
                idx = neighbors[which(new.cluster[neighbors]!=0)]
                connected[idx] <- num.clusters
              }
            }
          }
        }
      }
      connected.size <- as.data.frame(table(connected))
      for (k in 1:nrow(connected.size)) {
        if (connected.size$Freq[k] >= cluster.sz & connected.size$connected[k] != 0) {
          loop.count <- loop.count + 1
          out.nii <- out.nii + (connected == connected.size$connected[k]) * loop.count
        }
      }
    }
  }

  # Sort unlabelled to nearest label
  # unlabelled <- array(0, dim=img.dims)
  # for (i in 1:n.clusters) {
  #   unlabelled <- unlabelled + read.nii.volume(data.cluster, i)
  # }
  # unlabelled <- unlabelled * (out.nii == 0)
  # unlabel.ls <- which(unlabelled != 0, arr.ind=TRUE)
  # label.ls <- which(out.nii != 0, arr.ind=TRUE)
  #
  # for (i in 1:nrow(unlabel.ls)) {
  #   dist.temp <- sqrt(
  #     (unlabel.ls[i,1] - label.ls[ ,1])^2 +
  #     (unlabel.ls[i,2] - label.ls[ ,2])^2 +
  #     (unlabel.ls[i,3] - label.ls[ ,3])^2)
  #   nearest <- which(dist.temp == min(dist.temp))[1]
  #   out.nii[matrix(unlabel.ls[i, ], ncol=3)] <- out.nii[matrix(label.ls[nearest, ], ncol=3)]
  # }

  write.nii.volume(fname, 1, out.nii)
}
