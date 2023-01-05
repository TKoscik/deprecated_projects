draw.summary <- function(saved.model,
                         overlay.opts=list(),
                         effect.opts=list(),
                         table.opts=list(),
                         save.dir) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # debug
  # rm(list=ls())
  # saved.model <- "C:/Users/tim.koscik/Documents/Data/duplex.decision.20170227.nAGQ0/clusterAnalysis/clusterTest.cluster1.mdl1.RData"
  # overlay.opts=list()
  # effect.opts=list()
  # table.opts=list()
  # overlay.opts <- list(anat.nii="C:/Users/tim.koscik/Documents/Data/brains/MNI152_T1_2mm_brain.nii")
  # save.dir <- "C:/Users/tim.koscik/Documents/Data/duplex.decision.20170227.nAGQ0/summaries"
  #
  require(rmarkdown)
  
  load(saved.model)
  var.names <- ls()
  model.name <- var.names[!(var.names %in% c("cluster.name", "cluster.vol", "saved.model"))]
  
  # Initialize Rmd File --------------------------------------------------------
  out.name <- unlist(strsplit(saved.model, "[/]"))
  rmd.title <- unlist(strsplit(out.name[length(out.name)], split="[.]"))
  rmd.title <- paste0(rmd.title[-length(rmd.title)], collapse=".")
  out.name <- paste0(save.dir, "/", out.name[length(out.name)])
  out.name <- unlist(strsplit(out.name, split="[.]"))
  out.name <- paste0(paste(out.name[-length(out.name)], collapse="."), ".Rmd")
  fid <- file(out.name, "wt")
  cat('---\n', file=fid)
  cat(paste0('title: "', rmd.title, '"\n'), file=fid)
  cat('author: "mRi"\n', file=fid)
  cat('date: "', format(Sys.Date(), "%B %e, %Y"), '"\n', file=fid)
  cat('output: html_document\n', file=fid)
  cat('---\n', file=fid)
  cat('\n```{r setup, include=FALSE}\n', file=fid)
  cat('knitr::opts_chunk$set(echo = FALSE)\n', file=fid)
  cat('library(mRi)\n', file=fid)
  cat('library(ggplot2)\n', file=fid)
  cat('library(reshape2)\n', file=fid)
  cat('library(plotly)\n', file=fid)
  cat('library(gridExtra)\n', file=fid)
  cat('library(raster)\n', file=fid)
  cat('library(rgeos)\n', file=fid)
  cat('library(effects)\n', file=fid)
  cat('library(DT)\n', file=fid)
  cat('```\n', file=fid)
  
  # Parse Overlay Options ------------------------------------------------------
  if (!("over.nii" %in% names(overlay.opts))) {
    path.temp <- unlist(strsplit(cluster.name, split="[/]"))
    fname <- path.temp[length(path.temp)]
    path.temp <- paste0(path.temp[-length(path.temp)], collapse="/")
    fname <- unlist(strsplit(fname, split=".vol"))
    overlay.opts$over.nii <- paste0(path.temp, "/", fname[1], ".nii")
    overlay.opts$over.vol <- as.numeric(substr(fname[2],1,1))
  }
  if (!("over.vol" %in% names(overlay.opts))) {
    overlay.opts$over.vol <- nii.dims(overlay.opts$over.nii)[4]
  }
  if (!("over.color" %in% names(overlay.opts))) {
    overlay.opts$over.color <- c("#640000", "#ff0000", "#ffff00", "#ffff64")
  }
  if (!("mask.nii" %in% names(overlay.opts))) {
    overlay.opts$mask.nii <- cluster.name
  }
  if (!("mask.vol" %in% names(overlay.opts))) {
    overlay.opts$mask.vol <- cluster.vol
  }
  if (!("roi.nii" %in% names(overlay.opts))) {
    overlay.opts$roi.nii <- NULL
  }
  if (!("roi.val" %in% names(overlay.opts))) {
    overlay.opts$roi.val <- NULL
  }
  if (!("roi.color" %in% names(overlay.opts))) {
    overlay.opts$roi.color <- NULL
  }
  if (!("orientation" %in% names(overlay.opts))) {
    overlay.opts$orientation <- c("coronal", "axial", "sagittal")
  }
  if ("save.plot" %in% names(overlay.opts)) {
    if (overlay.opts$save.plot == TRUE) {
      if (!("save.dir" %in% names(overlay.opts))) {
        overlay.opts$save.dir <- save.dir
      }
      if (!("file.name" %in% names(overlay.opts))) {
        overlay.opts$file.name <- paste0(model.name, ".overlay.")
      }
      if (!("img.format" %in% names(overlay.opts))) {
        overlay.opts$img.format <- "png"
      }
      if (!("img.w" %in% names(overlay.opts))) {
        overlay.opts$img.w <- 10
      }
      if (!("img.unit" %in% names(overlay.opts))) {
        overlay.opts$img.unit <- "cm"
      }
      if (!("img.dpi" %in% names(overlay.opts))) {
        overlay.opts$img.dpi <- 600
      }
    } else {
      overlay.opts$save.plot <- FALSE
      overlay.opts$save.dir <- FALSE
      overlay.opts$file.name <- FALSE
      overlay.opts$img.format <- FALSE
      overlay.opts$img.w <- FALSE
      overlay.opts$img.unit <- FALSE
      overlay.opts$img.dpi <- FALSE
    }
  } else {
    overlay.opts$save.plot <- FALSE
    overlay.opts$save.dir <- FALSE
    overlay.opts$file.name <- FALSE
    overlay.opts$img.format <- FALSE
    overlay.opts$img.w <- FALSE
    overlay.opts$img.unit <- FALSE
    overlay.opts$img.dpi <- FALSE
  }
  
  # Plot Overlays --------------------------------------------------------------
  if ("coronal" %in% overlay.opts$orientation) {
    cat('\n```{r coronal overlay}\n', file=fid)
    cat(paste0('the.plot <- draw.overlay(anat.nii = "', overlay.opts$anat.nii, '",\n'), file=fid)
    cat(paste0('  over.nii = "', overlay.opts$over.nii, '",\n'), file=fid)
    cat(paste0('  over.vol = "', overlay.opts$over.vol, '",\n'), file=fid)
    cat(paste0('  over.color = list(c(', paste(sprintf("'%s'", overlay.opts$over.color), collapse=", "), ')),\n'), file=fid)
    cat(paste0('  mask.nii = "', overlay.opts$mask.nii, '",\n'), file=fid)
    cat(paste0('  mask.vol = "', overlay.opts$mask.vol, '",\n'), file=fid)
    cat(paste0('  roi.nii = "', overlay.opts$roi.nii, '",\n'), file=fid)
    cat(paste0('  roi.val = "', overlay.opts$roi.val, '",\n'), file=fid)
    cat(paste0('  roi.color = "', overlay.opts$roi.color, '",\n'), file=fid)
    cat(paste0('  orientation = "coronal",\n'), file=fid)
    cat(paste0('  idx.slice = TRUE,\n'), file=fid)
    cat(paste0('  cbars = TRUE,\n'), file=fid)
    cat(paste0('  save.dir = "', overlay.opts$save.dir, '",\n'), file=fid)
    cat(paste0('  file.name = "', paste0(overlay.opts$filename, "coronal"), '",\n'), file=fid)
    cat(paste0('  img.format = "', overlay.opts$img.format, '",\n'), file=fid)
    cat(paste0('  img.w = "', overlay.opts$img.w, '",\n'), file=fid)
    cat(paste0('  img.unit = "', overlay.opts$img.unit, '",\n'), file=fid)
    cat(paste0('  img.dpi = "', overlay.opts$img.dpi, '",\n'), file=fid)
    cat(paste0('  save.plot = "', overlay.opts$save.plot, '",\n'), file=fid)
    cat(paste0('  return.plot = TRUE)\n'), file=fid)
    cat('grid.newpage()\n', file=fid)
    cat('grid.draw(the.plot)\n', file=fid)
    cat('```\n', file=fid)
  }
  if ("axial" %in% overlay.opts$orientation) {
    cat('\n```{r axial overlay}\n', file=fid)
    cat(paste0('the.plot <- draw.overlay(anat.nii = "', overlay.opts$anat.nii, '",\n'), file=fid)
    cat(paste0('  over.nii = "', overlay.opts$over.nii, '",\n'), file=fid)
    cat(paste0('  over.vol = "', overlay.opts$over.vol, '",\n'), file=fid)
    cat(paste0('  over.color = list(c(', paste(sprintf("'%s'", overlay.opts$over.color), collapse=", "), ')),\n'), file=fid)
    cat(paste0('  mask.nii = "', overlay.opts$mask.nii, '",\n'), file=fid)
    cat(paste0('  mask.vol = "', overlay.opts$mask.vol, '",\n'), file=fid)
    cat(paste0('  roi.nii = "', overlay.opts$roi.nii, '",\n'), file=fid)
    cat(paste0('  roi.val = "', overlay.opts$roi.val, '",\n'), file=fid)
    cat(paste0('  roi.color = "', overlay.opts$roi.color, '",\n'), file=fid)
    cat(paste0('  orientation = "axial",\n'), file=fid)
    cat(paste0('  idx.slice = TRUE,\n'), file=fid)
    cat(paste0('  cbars = TRUE,\n'), file=fid)
    cat(paste0('  save.dir = "', overlay.opts$save.dir, '",\n'), file=fid)
    cat(paste0('  file.name = "', paste0(overlay.opts$filename, "axial"), '",\n'), file=fid)
    cat(paste0('  img.format = "', overlay.opts$img.format, '",\n'), file=fid)
    cat(paste0('  img.w = "', overlay.opts$img.w, '",\n'), file=fid)
    cat(paste0('  img.unit = "', overlay.opts$img.unit, '",\n'), file=fid)
    cat(paste0('  img.dpi = "', overlay.opts$img.dpi, '",\n'), file=fid)
    cat(paste0('  save.plot = "', overlay.opts$save.plot, '",\n'), file=fid)
    cat(paste0('  return.plot = TRUE)\n'), file=fid)
    cat('grid.newpage()\n', file=fid)
    cat('grid.draw(the.plot)\n', file=fid)
    cat('```\n', file=fid)
  }
  if ("sagittal" %in% overlay.opts$orientation) {
    cat('\n```{r sagittal overlay}\n', file=fid)
    cat(paste0('the.plot <- draw.overlay(anat.nii = "', overlay.opts$anat.nii, '",\n'), file=fid)
    cat(paste0('  over.nii = "', overlay.opts$over.nii, '",\n'), file=fid)
    cat(paste0('  over.vol = "', overlay.opts$over.vol, '",\n'), file=fid)
    cat(paste0('  over.color = list(c(', paste(sprintf("'%s'", overlay.opts$over.color), collapse=", "), ')),\n'), file=fid)
    cat(paste0('  mask.nii = "', overlay.opts$mask.nii, '",\n'), file=fid)
    cat(paste0('  mask.vol = "', overlay.opts$mask.vol, '",\n'), file=fid)
    cat(paste0('  roi.nii = "', overlay.opts$roi.nii, '",\n'), file=fid)
    cat(paste0('  roi.val = "', overlay.opts$roi.val, '",\n'), file=fid)
    cat(paste0('  roi.color = "', overlay.opts$roi.color, '",\n'), file=fid)
    cat(paste0('  orientation = "sagittal",\n'), file=fid)
    cat(paste0('  idx.slice = TRUE,\n'), file=fid)
    cat(paste0('  cbars = TRUE,\n'), file=fid)
    cat(paste0('  save.dir = "', overlay.opts$save.dir, '",\n'), file=fid)
    cat(paste0('  file.name = "', paste0(overlay.opts$filename, "sagittal"), '",\n'), file=fid)
    cat(paste0('  img.format = "', overlay.opts$img.format, '",\n'), file=fid)
    cat(paste0('  img.w = "', overlay.opts$img.w, '",\n'), file=fid)
    cat(paste0('  img.unit = "', overlay.opts$img.unit, '",\n'), file=fid)
    cat(paste0('  img.dpi = "', overlay.opts$img.dpi, '",\n'), file=fid)
    cat(paste0('  save.plot = "', overlay.opts$save.plot, '",\n'), file=fid)
    cat(paste0('  return.plot = TRUE)\n'), file=fid)
    cat('grid.newpage()\n', file=fid)
    cat('grid.draw(the.plot)\n', file=fid)
    cat('```\n', file=fid)
  }
  
  # Parse Effects Options ------------------------------------------------------
  if (!("effect.name" %in% effect.opts)) {
    effect.opts$effect.name <- labels(terms(eval(parse(text=model.name))))
  }
  n.terms <- length(effect.opts$effect.name)
  if (!("axes" %in% effect.opts)) {
    effect.opts$axes <- list()
    for (i in 1:n.terms) {
      ef.vars <- unlist(strsplit(effect.opts$effect.name[i], split=":"))
      ef.num <- length(ef.vars)
      effect.opts$axes[[i]] <- switch(as.character(ef.num),
        `1` = c(ef.vars[1], "fit"),
        `2` = c(ef.vars[1], "fit", ef.vars[2]),
        `3` = c(ef.vars[1], "fit", ef.vars[2], ef.vars[3]),
        `4` = c(ef.vars[1], "fit", ef.vars[2], ef.vars[3], ef.vars[4]),
        stop("Cannot plot interactions larger than 4-way"))
    }
  }
  if (!("labels" %in% effect.opts)) {
    effect.opts$labels <- list()
    for (i in 1:n.terms) {
      s.temp <- unlist(strsplit(effect.opts$axes[[i]], " "))
      effect.opts$labels[[i]] <- paste(toupper(substring(s.temp,1,1)), substring(s.temp,2), sep="")
    }
  }
  if (!("plot.colors" %in% effect.opts)) {
    effect.opts$plot.colors=list()
    for (i in 1:n.terms) {
      effect.opts$plot.colors[[i]] <- get.color()
    }
  }
  if ("save.plot" %in% names(effect.opts)) {
    if (effect.opts$save.plot == TRUE) {
      if (!("save.dir" %in% names(overlay.opts))) {
        effect.opts$save.dir <- save.dir
      }
      if (!("file.name" %in% names(overlay.opts))) {
        effect.opts$file.name <- paste0(model.name, ".overlay.")
      }
      if (!("img.format" %in% names(overlay.opts))) {
        effect.opts$img.format <- "png"
      }
      if (!("img.w" %in% names(overlay.opts))) {
        effect.opts$img.w <- 10
      }
      if (!("img.unit" %in% names(overlay.opts))) {
        effect.opts$img.unit <- "cm"
      }
      if (!("img.dpi" %in% names(overlay.opts))) {
        effect.opts$img.dpi <- 600
      }
    }
  }
    
  # Plot Effects ---------------------------------------------------------------
  for (i in 1:n.terms) {
    cat(paste0('\n```{r effect', i, '}\n'), file=fid)
    cat(paste0('draw.effect(model = ', model.name, ',\n'), file=fid)
    cat(paste0("  effect.name = '", effect.opts$effect.name[i], "',\n"), file=fid)
    cat(paste0('  axes = c(', paste(sprintf("'%s'", effect.opts$axes[[i]]), collapse=", "), '),\n'), file=fid)
    cat(paste0('  labels = c(', paste(sprintf("'%s'", effect.opts$labels[[i]]), collapse=", "), '),\n'), file=fid)
    cat(paste0('  plot.colors = c(', paste(sprintf("'%s'", effect.opts$plot.colors[[i]]), collapse=", "), '),\n'), file=fid)
    cat(paste0('  save.dir = "', effect.opts$save.dir, '",\n'), file=fid)
    cat(paste0('  file.name = "', paste0(effect.opts$filename, "coronal"), '",\n'), file=fid)
    cat(paste0('  img.format = "', effect.opts$img.format, '",\n'), file=fid)
    cat(paste0('  img.w = "', effect.opts$img.w, '",\n'), file=fid)
    cat(paste0('  img.unit = "', effect.opts$img.unit, '",\n'), file=fid)
    cat(paste0('  img.dpi = "', effect.opts$img.dpi, '",\n'), file=fid)
    cat(paste0('  save.plot = "', effect.opts$save.plot, '",\n'), file=fid)
    cat(paste0('  return.plot = TRUE)\n'), file=fid)
    cat('grid.newpage()\n', file=fid)
    cat('grid.draw(the.plot)\n', file=fid)
    cat('```\n', file=fid)
  }
  
  # Print Tables ---------------------------------------------------------------
  if (!("which.tables" %in% table.opts)) {
    table.opts$which.tables <- c("coef", "aov") # other options? lsm (lsmeans), dlsm (difflsmeans)
  }
  if ("coef" %in% table.opts$which.table) {
    cat(paste0('\n```{r fixed effects table}\n'), file=fid)
    cat(paste0('table.temp <- as.data.frame(summary(', model.name, ')$coef)\n'), file=fid)
    cat('datatable(table.temp, row.names=FALSE, caption = "Fixed Effects", filter="bottom")\n', file=fid)
    cat('```\n', file=fid)
  }
  if ("aov" %in% table.opts$which.table) {
    cat(paste0('\n```{r anova table}\n'), file=fid)
    cat(paste0('table.temp <- as.data.frame(Anova(', model.name, ', type=3L))\n'), file=fid)
    cat('datatable(table.temp, row.names=FALSE, caption = "ANOVA (car)", filter="bottom")\n', file=fid)
    cat('```\n', file=fid)
  }
  if ("lsm" %in% table.opts$which.table) {
    cat(paste0('\n```{r lsm table}\n'), file=fid)
    cat(paste0('table.temp <- lsmeans(', model.name, ')$lsmeans.table)\n'), file=fid)
    cat('datatable(table.temp, row.names=FALSE, caption = "Least-Squares Means", filter="bottom")\n', file=fid)
    cat('```\n', file=fid)
  }
  if ("dlsm" %in% table.opts$which.table) {
    cat(paste0('\n```{r dlsm table}\n'), file=fid)
    cat(paste0('table.temp <- difflsmeans(', model.name, ')$diffs.lsmeans.table)\n'), file=fid)
    cat('datatable(table.temp, row.names=FALSE, caption = "Difference of Least-Squares Means", filter="bottom")\n', file=fid)
    cat('```\n', file=fid)
  }
  close(fid)
  render(out.name)
}
