parse.4d <- function(data.4d) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(!missing(data.4d))
  if (is.character(data.4d)) {
    temp.4d <- data.4d
    data.4d <- list()
    if (length(temp.4d) == 1) { # Get filelist from specified directory
      if(!dir.exists(temp.4d)) {
        data.4d[[1]] <-temp.4d
      } else {
        data.4d[[1]] <- list.files(path=temp.4d, pattern="*.nii$", full.names=TRUE)
      }
    } else {
      data.4d[[1]] <- temp.4d
    }
    for (i in 1:length(data.4d[[1]])) { # check if 4d files exist
      stopifnot(file.exists(data.4d[[1]][i]))
    }
  } else if (is.list(data.4d)) {
    n.4d <- length(data.4d)
    temp.4d <- data.4d
    data.4d <- vector("list", n.4d)
    for (i in 1:n.4d) {
      if (length(temp.4d[[i]]) == 1) { # Get filelist from specified directory
        if(!dir.exists(temp.4d[[1]])) {
          data.4d[[1]] <- temp.4d[[i]]
        } else {
          data.4d[[1]] <- list.files(path=temp.4d[[i]], pattern="*.nii$", full.names=TRUE)
        }
        # data.4d[[i]] <- list.files(path = temp.4d[[i]], pattern="*.nii$", full.names=TRUE)
      } else {
        data.4d[[i]] <- temp.4d[[i]]
      }
      for (j in 1:length(data.4d[[i]])) { # check if 4d files exist
        stopifnot(file.exists(data.4d[[i]][j]))
      }
    }
  }
  return(data.4d)
}