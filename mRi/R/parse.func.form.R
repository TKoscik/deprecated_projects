parse.func.form <- function(FUNC, FORM) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  stopifnot(!missing(FUNC), is.character(FUNC),
            !missing(FORM), is.character(FORM))
  n.func <- length(FUNC)
  n.form <- length(FORM)
  if (n.func != n.form) {
    if (n.func < n.form) {
      FUNC <- rep(FUNC, n.form)
    } else {
      FORM <- rep(FORM, n.func)
    }
  }
  model <- paste0(FUNC, "(", FORM, ", data=df)")
}