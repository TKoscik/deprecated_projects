cluster.results <- function(cluster.nii,
                            cluster.vol = "all",
                            save.dir,
                            file.name,
                            opts = list(ts.analysis = FALSE,
                                        peak.find = FALSE)) {
  #-----------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-----------------------------------------------------------------------------

  # # Debug ----
  # rm(list=ls())
  # gc()
  # library("parallel")
  # library("lmerTest")
  #
  # cluster.nii <- "C:/Users/tim.koscik/Dropbox/Analyses/duplex.mni/decision/evANDdecision_ROIs.nii"
  # cluster.vol <- "all"
  # save.dir <- "C:/Users/tim.koscik/Documents/Analyses/duplex.cluster.tables"
  # file.name <- "cluster.results.test.RData"
  # opts <- list()
  # opts$ts.analysis <- TRUE
  # opts$cluster.ts <- "C:/Users/tim.koscik/Dropbox/Analyses/duplex.mni/decision/evANDdecision.clusterTS.RData"
  #
  # opts$peak.find <- TRUE
  # opts$peak.analysis <- TRUE
  # opts$value.nii <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/models/20170602/duplex.mdl3.coef.zvalue.nii"
  # opts$value.vol <- 4
  # opts$data.4d <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/stbs_RT"
  # opts$data.paradigm <- "C:/Users/tim.koscik/Documents/Data/Duplex.Decision/paradigm.onset.csv"
  # opts$data.paradigm <- parse.paradigm(opts$data.paradigm)
  # opts$data.paradigm$run <- opts$data.paradigm$block
  # opts$data.paradigm$subject.run <- paste0(opts$data.paradigm$subject, opts$data.paradigm$run)
  # opts$data.paradigm$response10 <- (opts$data.paradigm$response == 1)*1
  #
  # opts$my.command <- "my.function(df)"
  # opts$my.function <- function(df) {
  #   df <- df[df$task == "display", ]
  #   ex1 <- outliers.by.subject(df, "fmri", "subject.run", 3)
  #   df <- ex1$data
  #
  #   output <- list()
  #   output$f1 <- formula("fmri ~ ev.good.c + ev.bad.c + (1|subject/run)")
  #   output$f2 <- formula("response10 ~ (fmri.grand.c + fmri.subject.c + fmri.subject.run.c) + (ev.good.c + ev.bad.c) + (1|subject/run)")
  #
  #   tf <- quantile.scale(df, "fmri", center=F, prob=0.9)
  #   tf <- decompose.var(tf, "fmri", grouping=c("grand", "subject", "subject.run"))
  #
  #   mdl1 <- lmer(output$f1, df)
  #   mdl2 <- glmer(output$f2, tf, family=binomial,
  #                 nAGQ=1, control=glmerControl(calc.derivs=FALSE, optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
  #
  #   ex2.1 <- outliers.by.std.resid(model=mdl1, data=df, trim=3)
  #   ex2.2 <- outliers.by.std.resid(model=mdl2, data=tf, trim=3)
  #
  #   tf1 <- ex2.1$data
  #   tf2 <- ex2.2$data[ ,-seq(ncol(ex2.2$data)-2, ncol(ex2.2$data), 1)]
  #
  #   tf2 <- decompose.var(tf2, "fmri", grouping=c("grand", "subject", "subject.run"))
  #
  #   mdl1 <- lmer(output$f1, tf1)
  #   mdl2 <- glmer(output$f2, tf2, family=binomial,
  #                 nAGQ=1, control=glmerControl(calc.derivs=FALSE, optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
  #
  #   output$coef1 <- as.data.frame(summary(mdl1)$coef)
  #   output$coef2 <- as.data.frame(summary(mdl2)$coef)
  #
  #   # Save residual standard deviation for each model
  #   output$ressd <- matrix(c(sd(resid(mdl1)), sd(resid(mdl2))), ncol=1)
  #   rownames(output$ressd) <- c("mdl1", "mdl2")
  #   colnames(output$ressd) <- "sigma"
  #
  #   output$out.rm <- matrix(c(ex1$n.removed, ex2.1$n.removed, ex2.2$n.removed,
  #                      ex1$percent.removed, ex2.1$percent.removed, ex2.2$percent.removed), ncol=2)
  #   rownames(output$out.rm) <- c("By Subject", "By Std Resid 1", "By Std Resid 2")
  #   colnames(output$out.rm) <- c("N", "PCT")
  #
  #   return(output)
  # }
  # opts$num.cores <- 1
  # # ----

  # Parse Options ----
  opts.ls <- names(opts)

  if (opts$peak.find) {
    if (!("value.nii" %in% opts.ls)) { stop("must provide a value map (opts$value.nii) to locate peak voxels") }
    if (!("value.vol" %in% opts.ls)) { stop("must provide a volume index for value map (opts$value.vol) to locate peak voxels") }
    if (opts$peak.analysis) {
      if (!("data.4d" %in% opts.ls)) { stop("must provide location of 4D data for modelling peak voxel") }
      if (!("data.paradigm" %in% opts.ls)) { stop("must provide paradigm data for modelling peak voxel") }
      if (!("my.command" %in% opts.ls)) { opts$my.command <- "my.function(df)" }
      if (!("my.function" %in% opts.ls)) { stop("must provide modelling function") }
    }
  } else {
    opts$value.nii <- NULL
    opts$value.vol <- NULL
  }

  if (opts$ts.analysis) {
    if (!"cluster.ts" %in% opts.ls) { stop("must provide cluster time series") }
    if (!("my.command" %in% opts.ls)) { opts$my.command <- "my.function(df)" }
    if (!("my.function" %in% opts.ls)) { stop("must provide modelling function") }
  }

  if (("my.function" %in% opts.ls)) { my.function <- opts$my.function }

  # Retrieve Cluster Coordinates ----
  coords <- cluster.coord(cluster.nii, cluster.vol, opts$value.nii, opts$value.vol)
  n.clusters <- length(coords$n.voxels)

  # Calculate Model at Peak cluster ----
  if (opts$peak.analysis) {
    data.4d <- parse.4d(opts$data.4d)
    pf <- parse.paradigm(opts$data.paradigm)
    coord.ls <- as.matrix(coords$peak[ ,1:3])
    temp.fxn <- function(X, ...) {
      fmri <- load.voxel(data.4d, coord.ls[X, ])
      df <- data.frame(fmri, pf)
      output <- eval(parse(text=opts$my.command))
    }
    peak.results <- mclapply(X=seq_len(n.clusters),
                             FUN=temp.fxn,
                             mc.preschedule=FALSE,
                             mc.cores=opts$num.cores)
  }

  # Calculate model for cluster time series
  if (opts$ts.analysis) {
    tname <- load(opts$cluster.ts)
    fmri.ts <- get(tname)
    rm("tname")
    pf <- parse.paradigm(opts$data.paradigm)
    
    temp.fxn <- function(X, ...) {
      fmri <- fmri.ts[[X]]
      df <- data.frame(fmri, pf)
      output <- eval(parse(text=opts$my.command))
    }
    ts.results <- mclapply(X=seq_len(n.clusters),
                             FUN=temp.fxn,
                             mc.preschedule=FALSE,
                             mc.cores=opts$num.cores)
  }

  save(list=c("coords", "peak.results", "ts.results"),
       file = paste0(save.dir, "/", file.name))
}
