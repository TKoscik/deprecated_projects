cluster.split.anatomy <- function(data.cluster, data.anatomy, which.cluster="all", cluster.sz) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(file.exists(data.cluster), file.exists(data.anatomy))
  img.dims <- nii.dims(data.cluster)[1:3]
  n.clusters <- nii.dims(data.cluster)[4]
  data.anatomy <- read.nii.volume(data.anatomy, 1)
  
  new.cluster <- vector("list", n.clusters)
  n.new <- numeric(n.clusters)
  for (i in 1:n.clusters) {
    cluster.temp <- read.nii.volume(data.cluster, i)
    cluster.temp <- data.anatomy * cluster.temp
    roi <- unique(as.vector(cluster.temp))
    roi <- roi[-which(roi==0)]
    if (length(roi) != 0) {
      new.cluster[[i]] <- array(0, dim=c(img.dims, length(roi)))
      sub.cluster.sz <- numeric(length(roi))
      for (j in 1:length(roi)) {
        sub.cluster.sz[j] <- sum(cluster.temp==roi[j])
        if (sub.cluster.sz[j] >= cluster.sz) {
          n.new[i] <- n.new[i] + 1
          new.cluster[[i]][ , , ,n.new[i]] <- cluster.temp == roi[j]
        }
      }
    } else {
      n.new[i] <- n.new[i] + 1
      new.cluster[[i]] <- array(read.nii.volume(data.cluster, i), dim=c(img.dims, 1))
    }
  }
  comb.cluster <- array(0, dim=c(img.dims, sum(n.new)))
  clust.count <- 0
  for (i in 1:n.clusters) {
    if (n.new[i] !=0) {
      for (j in 1:n.new[i]) {
        clust.count <- clust.count +1
        comb.cluster[ , , ,clust.count] <- new.cluster[[i]][ , , ,j]
      }
    }
  }
  
  
  #n <- 3
  #dimorder <- 1:3
  #cdim <- cumprod(c(1, img.dims[dimorder][-n]))
  #neighborhood <- t(matrix(
  #      c(1,1,1,1,1,2,1,1,3,1,2,1,1,2,2,1,2,3,1,3,1,1,3,2,1,3,3,2,1,1,
  #        2,1,2,2,1,3,2,2,1,2,2,3,2,3,1,2,3,2,2,3,3,3,1,1,3,1,2,3,1,3,
  #        3,2,1,3,2,2,3,2,3,3,3,1,3,3,2,3,3,3), nrow=3))
  #center.pt <- t(matrix(c(2, 2, 2), nrow = 3))
  #offsets <- as.integer(colSums(t(neighborhood[, 1:3, drop = FALSE] - 1) * cdim) + 1L) -
  #  as.integer(colSums(t(center.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
  #neighborhood <- neighborhood - 2
  #
  #num.cl <- dim(comb.cluster)[4]
  #dup.cl <- numeric()
  #for (i in 1:num.cl) {
  #  temp.cl <- comb.cluster[ , , ,i]
  #  connected <- array(0, dim = img.dims[1:3])
  #  num.clusters <- 0
  #  idx <- numeric()
  #  for (x in 1:img.dims[1]) {
  #    for (y in 1:img.dims[2]) {
  #      for (z in 1:img.dims[3]) {
  #        if (temp.cl[x, y, z] == 1) {
  #          num.clusters <- num.clusters + 1
  #          current.pt <- t(as.matrix(c(x, y, z)))
  #          idx <- as.integer(colSums(t(current.pt[, 1:3, drop = FALSE] - 1) * cdim) + 1L)
  #          connected[idx] <- num.clusters
  #          while (length(idx) != 0) {
  #            temp.cl[idx] <- 0
  #            neighbors = as.vector(apply(as.matrix(idx), 1, "+", offsets))
  #            neighbors = unique(neighbors[which(neighbors > 0)])
  #            idx = neighbors[which(temp.cl[neighbors] != 0)]
  #            connected[idx] <- num.clusters
  #          }
  #        }
  #      }
  #    }
  #  }
  #  cl.ls <- as.data.frame(table(connected))
  #  cl.ls$connected <- as.numeric(levels(cl.ls$connected))
  #  cl.ls <- cl.ls[-1, ]
  #  remove.ls <- as.numeric(init.cl.ls$connected[which(cl.ls$Freq < cluster.sz)])
  #  connected[connected %in% remove.ls] <- 0
  #  cl.ls <- cl.ls[!(cl.ls$connected %in% remove.ls), ]
  #  if (nrow(cl.ls) > 1) { 
  #    for (j in 2:nrow(cl.ls)) {
  #      new.array <- array(0,dim=c(dim(comb.cluster)[1:3], dim(comb.cluster)[4]+1))
  #      new.array[ , , ,1:dim(comb.cluster)[4]] <- comb.cluster
  #      new.array[ , , ,dim(comb.cluster)[4]+1] <- (connected == j) * 1
  #      comb.cluster <- new.array
  #    }
  #    dup.cl <- c(dup.cl, i)
  #  }
  #}
  #comb.cluster <- comb.cluster[ , , ,-dup.cl]
 # 
  #cl.size <- numeric(dim(comb.cluster)[4])
  #for (i in 1:dim(comb.cluster)[4]) {
   # cl.size[i] <- sum(comb.cluster[ , , , i])
  #}
  #comb.cluster <- comb.cluster[ , , , order(cl.size, decreasing=TRUE)]
  
  init.nii(file.name = paste0(substr(data.cluster,1,(nchar(data.cluster)-3)), "anat.nii"),
           dims = dim(comb.cluster),
           pixdim = unlist(nii.hdr(data.cluster, "pixdim")),
           orient = nii.orient(data.cluster))
  for (i in 1:dim(comb.cluster)[4]) {
    write.nii.volume(nii.file = paste0(substr(data.cluster,1,(nchar(data.cluster)-3)), "anat.nii"),
                     vol.num = i,
                     values = comb.cluster[ , , ,i])
  }
}
