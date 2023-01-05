easy.mlm.cluster <- function(
  FUNC,
  FORM,
  df,
  model.name = "mdl",
  outliers.subject = c("fmri", "subject.run", 3),
  outliers.model = 3,
  scale.vars = c("fmri", FALSE, 0.9),
  decompose.vars = c("fmri", "grand", "subject", "subject.run"),
  save.dir,
  prefix,
  cluster.name, cluster.vol, ...) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  # Remove outliers by subject -------------------------------------------------
  if (!is.null(outliers.subject[1])) {
    ex1 <- lme.outliers.by.subject(
      df,
      outliers.subject[1],
      outliers.subject[2],
      as.numeric(outliers.subject[3]))
    tf1 <- ex1$data
  } else {
    tf1 <- df
  }
  
  # Scale Variables ------------------------------------------------------------
  if (is.character(scale.vars)) {
    if (!is.null(scale.vars[1])) {
      tf1 <- quantile.scale(
        tf1,
        scale.vars[1],
        center=as.logical(scale.vars[2]),
        prob=as.numeric(scale.vars[3]))
    }
  } else if (is.list(scale.vars)) {
    for (i in 1:length(scale.vars)) {
      tf1 <- quantile.scale(
        tf1,
        scale.vars[[i]][1],
        center=as.logical(scale.vars[[i]][2]),
        prob=as.numeric(scale.vars[[i]][3]))
    }
  } else {
    error("Cannot parse scale.vars input.")
  }
  
  # Decompose Variables --------------------------------------------------------
  if (is.character(decompose.vars)) {
    if (!is.null(decompose.vars[1])) {
      tf1 <- decompose.var(tf1,
                           decompose.vars[1],
                           grouping=decompose.vars[2:length(decompose.vars)])
    }
  } else if (is.list(desompose.vars)) {
    for (i in 1:length(decompose.vars)) {
      tf1 <- decompose.var(
        tf1,
        decompose.vars[[i]][1],
        grouping=decompose.vars[[i]][2:length(decompose.vars)])
    }
  } else {
    error("Cannot parse decompose.vars input.")
  }
  
  # Run Model ------------------------------------------------------------------
  model.str <- switch(FUNC,
                      `lmer` = paste0("lmer(", formula(FORM), ", data=tf1)"),
                      `glmer` = paste0("glmer(",
                                       FORM, ", ",
                                       "data=tf1, ",
                                       "family=binomial, ",
                                       "nAGQ=0, ",
                                       "control=glmerControl(",
                                       "calc.derivs=FALSE, ",
                                       "optimizer='bobyqa', ",
                                       "optCtrl=list(maxfun=1000000))) "))
  assign(model.name, eval(parse(text=model.str)))
  
  # Outlier Removal by Model ---------------------------------------------------
  if (!is.null(outliers.model[1])) {
    if (exists(ex1)) {
      ex2 <- lme.outliers.by.std.resid(model=eval(parse(text=model.name)),
                                       data=ex1$data,
                                       trim=outliers.model[1])
    } else {
      ex2 <- lme.outliers.by.std.resid(model=eval(parse(text=model.name)),
                                       data=df,
                                       trim=outliers.model[1])
    }
    tf2 <- ex2$data
    
    ## Re-Scale Variables ------------------------------------------------------
    if (is.character(scale.vars)) {
      if (!is.null(scale.vars[1])) {
        tf2 <- quantile.scale(
          tf2,
          scale.vars[1],
          center=as.logical(scale.vars[2]),
          prob=as.numeric(scale.vars[3]))
      }
    } else if (is.list(scale.vars)) {
      for (i in 1:length(scale.vars)) {
        tf2 <- quantile.scale(
          tf2,
          scale.vars[[i]][1],
          center=as.logical(scale.vars[[i]][2]),
          prob=as.numeric(scale.vars[[i]][3]))
      }
    } else {
      error("Cannot parse scale.vars input.")
    }
    
    ## Re-Decompose Variables --------------------------------------------------
    if (is.character(decompose.vars)) {
      if (!is.null(decompose.vars[1])) {
        tf2 <- decompose.var(tf2,
                             decompose.vars[1],
                             grouping=decompose.vars[2:length(decompose.vars)])
      }
    } else if (is.list(desompose.vars)) {
      for (i in 1:length(decompose.vars)) {
        tf2 <- decompose.var(
          tf2,
          decompose.vars[[i]][1],
          grouping=decompose.vars[[i]][2:length(decompose.vars)])
      }
    } else {
      error("Cannot parse decompose.vars input.")
    }
    ## Re-Run Model ------------------------------------------------------------
    model.str <- switch(FUNC,
                        `lmer` = paste0("lmer(", formula(FORM), ", data=tf2)"),
                        `glmer` = paste0("glmer(",
                                         FORM, ", ",
                                         "data=tf2, ",
                                         "family=binomial, ",
                                         "nAGQ=0, ",
                                         "control=glmerControl(",
                                         "calc.derivs=FALSE, ",
                                         "optimizer='bobyqa', ",
                                         "optCtrl=list(maxfun=1000000))) "))
    assign(model.name, eval(parse(text=model.str)))
  }
  
  # Save Output ----------------------------------------------------------------
  save(list = c("model.name", "cluster.name", "cluster.vol"),
       file = paste0(save.dir, "/", prefix, ".", model.name, ".RData"))
}