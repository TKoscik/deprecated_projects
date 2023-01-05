parse.onsets <- function(data.onsets) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(!missing(data.onsets))
  
  temp.onsets <- character()
  if (is.list(data.onsets)) {
    for (i in 1:length(data.onsets)) {
      for (j in 1:length(data.onsets[[i]])) {
        if (dir.exists(data.onsets[[i]][j])) {
          fls <- list.files(path=data.onsets[[i]][j], pattern=".csv", full.names=TRUE)
          for (k in 1:length(fls)) {
            if (!grepl("scan.params", fls[k])) {
              temp.onsets <- c(temp.onsets, fls[k])
            }
          }
        } else if (file.exists(data.onsets[[i]][j])) {
          temp.onsets <- c(temp.onsets, data.onsets[[i]][j])
        } else {
          stop("Cannot parse onsets input.")
        }
      }
    }
  } else if (is.character(data.onsets)) {
    for (i in 1:length(data.onsets)) {
      if (dir.exists(data.onsets[i])) {
        fls <- list.files(path=data.onsets[i], pattern=".csv", full.names=TRUE)
        for (k in 1:length(fls)) {
          if (!grepl("scan.params", fls[k])) {
            temp.onsets <- c(temp.onsets, fls[k])
          }
        }
      } else if (file.exists(data.onsets[i])) {
        temp.onsets <- c(temp.onsets, data.onsets[i])
      } else {
        stop("Cannot parse onsets input.")
      }
    }
  }
  data.onsets <- temp.onsets
  return(data.onsets)
}