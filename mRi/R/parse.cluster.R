parse.cluster <- function(data.cluster, which.cluster="all") {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  cluster.ls <- character()
  if (is.character(data.cluster)) {
    for (i in 1:length(data.cluster)) {
      if (dir.exists(data.cluster[i])) {
        cluster.ls <- c(cluster.ls,
          list.files(data.cluster[i], pattern="*.nii$", full.names=TRUE))
      } else if (file.access(data.cluster[i])) {
        cluster.ls <- c(cluster.ls, data.cluster[i])
      } else {
        stop("Cluster files and/or directories not found.")
      }
    }
  }
  n.maps <- length(cluster.ls)
  
  temp.cluster <- which.cluster
  which.cluster <- vector("list", n.maps)
  if (is.list(temp.cluster)) {
    if (length(temp.cluster) == n.maps) {
      for (i in 1:n.maps) {
        if (is.character(temp.cluster[[i]])) {
          if (temp.cluster[[i]]=="all") {
            which.cluster[[i]] <- seq(1,nii.dims(cluster.ls[i])[4])
          } else {
            stop("Cannot parse cluster inputs.")
          }
        } else if (is.numeric(temp.cluster[[i]])) {
          poss.clusters <- seq(1,nii.dims(cluster.ls[i])[4])
          if (all(temp.cluster %in% poss.clusters)) {
            which.cluster[[i]] <- temp.cluster[[i]]
          } else {
            stop("Specified cluster out of bounds of clusters in maps.")
          }
        }
      }
    } else {
      stop("Cluster maps and volume lists must match.")
    }
  } else if (is.character(temp.cluster) & length(temp.cluster)==1 & temp.cluster=="all") {
    for (i in 1:n.maps) {
      which.cluster[[i]] <- seq(1,nii.dims(cluster.ls[i])[4])
    }
  } else if (is.numeric(temp.cluster)) {
    for (i in 1:n.maps) {
      which.cluster[[i]] <- seq(1,nii.dims(cluster.ls[i])[4])
      if (all(temp.cluster %in% poss.clusters)) {
        which.cluster[[i]] <- temp.cluster
      } else {
        stop("Specified cluster out of bounds of clusters in maps.")
      }
    }
  }
  
  cluster <- vector("list", n.maps)
  for (i in 1:n.maps) {
    cluster[[i]]$file <- cluster.ls[i]
    cluster[[i]]$volumes <- which.cluster[[i]]
  }
  
  return(cluster)
}