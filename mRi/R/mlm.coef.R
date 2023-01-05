mlm.coef <- function(mdl,
                     df = c("glm-df", "conservative")) {

  mdl.coef <- as.data.frame(summary(mdl, ddf = "lme4")$coef)

  n.efs <- nrow(mdl.coef)
  n.obs <- nrow(mdl@pp$V)
  t.coef <- coef(mdl)
  n.lvls <- numeric(length(t.coef))
  for (i in 1:length(t.coef)) {
    n.lvls[i] <- nrow(t.coef[[i]])
  }
  n.unique <- numeric(n.efs)
  for (i in 1:n.efs) {
    n.unique[i] <- length(unique(mdl@pp$V[ ,i]))
  }

  if ("glm-df" %in% df) {
    mdl.coef$glm.df <- c(n.obs-n.efs,
                         rep(n.obs - n.efs, n.efs-1))
    mdl.coef$glm.p <- 2 * (1 - pt(abs(mdl.coef$`t value`), df=mdl.coef$glm.df))
  }

  if ("conservative" %in% df) {
    df.temp <- numeric(n.efs)
    df.temp[1] <- (n.obs-n.efs)/mean(n.lvls)
    for (i in 2:n.efs) {
      if (grepl(":", rownames(mdl.coef)[i])) {
        df.temp[i] <- (n.obs - n.efs - mean(n.lvls))
      } else {
        df.temp[i] <- (n.obs - n.efs - mean(n.lvls)) * (n.unique[i]/n.obs)
      }
    }

    mdl.coef$con.df <- df.temp
    mdl.coef$con.p <- 2 * (1 - pt(abs(mdl.coef$`t value`), df=mdl.coef$con.df))
  }

  return(mdl.coef)
}
