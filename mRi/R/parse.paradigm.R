parse.paradigm <- function(data.paradigm) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(!missing(data.paradigm))
  
  if (is.character(data.paradigm)) {
    if (file.exists(data.paradigm)) {
      paradigm.df <- read.csv(data.paradigm, header=TRUE)
    } else if (dir.exists(data.paradigm)) {
      fls <- list.files(path=data.paradigm, pattern=c(".csv|.txt"), full.names=TRUE)
      paradigm.df <- data.frame()
      for (i in 1:length(fls)){
        temp <- read.csv(data.paradigm[i], header=TRUE)
        paradigm.df <- cbind(paradigm.df, temp)
      }
    } else {
      stop("Cannot parse paradigm input")
    }
  } else if (is.data.frame(data.paradigm)) {
    paradigm.df <- data.paradigm
  } else if (is.matrix(data.paradigm)) {
    paradigm.df <- as.data.frame(data.paradigm)
  } else {
    stop("Cannot parse paradigm input")
  }
  
  # var.names <- colnames(paradigm.df)
  return(paradigm.df)
}