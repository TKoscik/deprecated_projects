merge.masks <- function(data.dir, save.dir, prefix = NULL) {

  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------

  fls <- list.files(path = data.dir, pattern = "*nii")
  fls.sub <- unlist(strsplit(x = fls, split = "mask"))
  if (is.null(prefix)) {
    prefix <- substr(fls.sub[1], start = 1, stop = nchar(fls.sub) - 1)
  }

  to.merge <- substr(fls.sub[seq(2,length(fls.sub),2)], start = 3, stop = 1000)
  unique.fls <- unique(to.merge)

  if (length(unique.fls) == length(to.merge)) {
    warning('Did not find files to merge.')
  } else {
    for (i in 1:length(unique.fls)) {
      n.masks <- sum(to.merge == unique.fls[i])
      if (n.masks > 1) {
        img.dims <- nii.dims(paste0(data.dir, "/", prefix, ".mask1.", unique.fls[i]))
        if (length(img.dims == 4)) {
          n.vols <- img.dims[4]
        } else {
          n.vols <- 1
        }

        init.nii(file.name = paste0(save.dir, "/", prefix, ".", unique.fls[i]),
                 dims = img.dims,
                 pixdim = unlist(nii.hdr(nii.file = paste0(data.dir, "/", prefix, ".mask1.", unique.fls[i]),
                                         field = "pixdim")),
                 orient = nii.orient(nii.file = paste0(data.dir, "/", prefix, ".mask1.", unique.fls[i])))

        for (j in 1:n.vols) {
          b <- array(0, dim = img.dims[1:3])
          for (k in 1:n.masks) {
            a <- read.nii.volume(nii.file = paste0(data.dir, "/", prefix, ".mask", k, ".", unique.fls[i]),
                                      vol.num = j)
            b[a != 0 & !is.na(a)] <- a[a != 0 & !is.na(a)]
          }
          write.nii.volume(nii.file = paste0(save.dir, "/", prefix, ".", unique.fls[i]),
                           vol.num = j,
                           values = b)
        }
      } else {
        file.copy(from = paste0(data.dir, "/", prefix, ".mask1.", unique.fls[i]),
                  to = paste0(data.dir, "/", prefix, ".", unique.fls[i]))
      }
    }
  }
}
