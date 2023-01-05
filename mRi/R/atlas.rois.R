atlas.rois <- function(cluster.nii,
                       cluster.vol="all",
                       atlas.nii,
                       atlas.csv,
                       return.table = TRUE,
                       save.table = FALSE,
                       save.dir = NULL,
                       prefix = NULL) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------

  if (cluster.vol=="all") { cluster.vol <- 1:(nii.dims(cluster.nii)[4]) }

  atlas.nii <- read.nii.volume(atlas.nii, 1)
  atlas.csv <- read.csv(atlas.csv, header=TRUE, stringsAsFactors=FALSE)

  cluster.id <- numeric(0L)
  cluster.csv <- atlas.csv[0, ]
  for (i in cluster.vol) {
    temp.nii <- read.nii.volume(cluster.nii, i)
    atlas.nums <- unique(atlas.nii[temp.nii==1])

    for (j in 1:length(atlas.nums)) {
      cluster.id <- c(cluster.id, i)
      if (atlas.nums[j] != 0) {
        cluster.csv <- rbind(cluster.csv, atlas.csv[atlas.nums[j], ])
      } else {
        cluster.csv[nrow(cluster.csv)+1, ] <- rep(NA, ncol(atlas.csv))
      }
    }
  }
  cluster.csv <- data.frame(Cluster.ID=cluster.id, cluster.csv)

  if (save.table) {
    write.table(cluster.csv, file=paste0(save.dir, "/", prefix, ".csv"),
                row.names=FALSE, sep=",", quote=FALSE)
  }

  if (return.table) {
    return(cluster.csv)
  }
}
