parse.model <- function(FUNC, FORM,
                        out.rm.subject, out.rm.trim=matrix(c(3,3), ncol=2),
                        which.effects="all") {
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
  
  temp.model <- paste0(FUNC, "(", FORM, ", data=df)")
  n.model <- length(temp.model)
  
  model <- vector("list", n.model)
  for (i in 1:n.model) {
    model[[i]]$model <- temp.model[i]
  }
  
  if (length(out.rm.subject) != n.model) {
    if (length(out.rm.subject) == 1) {
      out.rm.subject <- rep(out.rm.subject, n.model)
    } else {
      stop("Cannot parse outlier removal.")
    }
  }
  for (i in 1:n.model) {
    model[[i]]$out.rm.subject <- out.rm.subject[i]
  }
  
  if (nrow(out.rm.trim) != n.model) {
    if (nrow(out.rm.trim) == 1) {
      out.rm.trim <- matrix(rep(out.rm.trim, each=n.model), nrow=n.model)
    } else {
      stop("Cannot parse outlier removal.")
    }
  }
  for (i in 1:n.model) {
    model[[i]]$out.rm.trim <- as.vector(out.rm.trim[i, ])
  }
  
  if (is.list(which.effects)) {
    if (length(which.effects) == n.model) {
      for (i in 1:n.model) {
        temp.fx <- colnames(attr(terms(formula(FORM)), "factors"))
        temp.fx <- temp.fx[-grep("[|]", temp.fx)]
        if (which.effects[[i]] == "all") {
          model[[i]]$effects <- temp.fx
        } else {
          if (all(which.effects[[i]] %in% temp.fx)) {
            model[[i]]$effects <- which.effects[[i]]
          } else {
            stop("Specified effects not in model")
          }
        }
      }
    } else {
      stop("List of effects must match the specified models")
    }
  } else if (is.character(which.effects) & length(which.effects)==1 & which.effects=="all") {
    for (i in 1:n.model) {
      temp.fx <- colnames(attr(terms(formula(FORM)), "factors"))
      model[[i]]$effects <- temp.fx[-grep("[|]", temp.fx)]
    }
  } else {
    stop("Cannot parse effects.")
  }
}