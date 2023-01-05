faster.cor <- function(inputMx, method="pearson") {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
if (!is.matrix(inputMx)) {
  inputMx <- as.matrix(inputMx)
}

if (method != "pearson") {
  inputMx <- apply(inputMx, 2L, rank, na.last = "keep")
}

inputMx <- t(inputMx)
inputMx <- inputMx - rowMeans(inputMx);
inputMx <- inputMx / sqrt(rowSums(inputMx^2))
outputMx <- tcrossprod(inputMx)
gc()
return(outputMx)

}