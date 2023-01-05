# ====================================================================================
# File:   			setup.design.R
#	Type:					R function
#	Description:	Create design matrices for deconvolution.
# NOTES:  2016-06-22  INPUT CHECKS AND ALL FUNCTIONALITY HAVE BEEN IMPLEMENTED YET - TK
# Usage ----
#               trials      [required]  Either a # x 3 matrix or the full path to a .csv
#                                       file contianing the same information. The format
#                                       is identical to FSL's 3-column format:
#                                       onset (s), duration (s), amplitude.  Note,
#                                       single trial betas will be based on this input.
#               n.tr        [required]  The number of volumes in the dataset.
#               tr          [required]  The repetition time of the collected data.
#               method      [optional]  The method for creating the design matrix,
#                                       options: "lsa" (equivalent AFNI 3dDeconvolve),
#                                       "lss" (equivalent to Mumford, et al. (2012)),
#                                       "lsr" [Default]
#               nuisance    [optional]  LSR method only.
#                                       A list of # x 3 matrices or full file paths to
#                                       .csv files.  Data is in the same 3 column
#                                       format.  Note: this input will generate
#                                       regressors that have been zeroed where the
#                                       current trial is present, which makes this
#                                       innappropriate to regress out noise (e.g., due
#                                       to movement)
#               noise       [optional]  A list of # x 3 matrices or full file paths to
#                                       .csv files.  Data is in the same 3 column
#                                       format.  Note: these regressors will not be
#                                       zeroed for the current trial when using LSR.
#               num.cores   [optional]  Default = num.cores - 1
#               save.design [optional]  Whether or not to save design
#               save.name   [optional]  full file path to save design data
# NOTE: INPUT CHECKS ARE NOT IMPLEMENTED
# ----
#	Project:			mRi
#	Author:				TR Koscik
#	Created:			2016-06-22 -- 15:37
#	Revisions:
# Citations:    LSS Method:
#               Mumford JA, Turner BO, Ashby FG, & Poldrack RA (2012) Deconvolving BOLD
#               activation in event-related designs for multivoxel pattern
#               classification analyses. NeuroImage, 59, 2636-3643.
# ======================================================================================

setup.design <- function (trials, n.tr, tr,
                          method = "lsr",
                          nuisance = NULL,
                          noise = NULL,
                          num.cores = NULL,
                          save.design = FALSE,
                          save.name = NULL) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  

# Parse Inputs -------------------------------------------------------------------------
  if (is.character(trials)) {
    trials <- read.csv(trials, header=FALSE, as.is=TRUE)
  }
  if (!is.null(nuisance)) {
    if (is.list(nuisance)) {
      for (i in 1:length(nuisance)) {
        if (is.character(nuisance[[i]])) {
          nuisance[[i]] <- read.csv(nuisance[[i]], header=FALSE, as.is=TRUE)
        } else if (is.matrix(nuisance[[i]])) {
        } else {
          stop("cannot parse nuisance regressors")
        }
      }
    } else if (is.character(nuisance)) {
      nuisance <- list(read.csv(nuisance, header=FALSE, as.is=TRUE))
    } else if (is.matrix(nuisance)) {
      nuisance <- list(nuisance)
    } else {
      stop("cannot parse nuisance regressors")
    }
  }
  if (!is.null(noise)) {
    if (is.list(noise)) {
      for (i in 1:length(nuisance)) {
        if (is.character(noise[[i]])) {
          noise[[i]] <- read.csv(noise[[i]], header=FALSE, as.is=TRUE)
        } else if (is.matrix(noise[[i]])) {
        } else {
          stop("cannot parse noise regressors")
        }
      }
    } else if (is.character(noise)) {
      noise <- list(read.csv(noise, header=FALSE, as.is=TRUE))
    } else if (is.matrix(noise)) {
      noise <- list(noise)
    } else {
      stop("cannot parse noise regressors")
    }
  }

# Setup Design Matrices ----------------------------------------------------------------
  hrf <- convolve.hrf(onset=trials[ ,1],
                      duration=trials[ ,2],
                      amplitude=trials[ ,3],
                      tr=tr, n.tr=n.tr,
                      st.betas=TRUE, scale=FALSE)

  if (!is.null(noise)) {
    noise.ts <- matrix(0, nrow=n.tr, ncol=length(noise))
    for (i in 1:length(noise)) {
      noise.temp <- convolve.hrf(onset=noise[[i]][ ,1],
                                 duration=noise[[i]][ ,2],
                                 amplitude=noise[[i]][ ,3],
                                 tr=tr, n.tr=n.tr,
                                 st.betas=TRUE, scale=FALSE)
      noise.ts[ ,i] <- rowSums(noise.temp$ts)
    }
  }

  if (method == "lsa") {
    design <- hrf$ts
    if (!is.null(noise)) {
      design <- cbind(design, noise.ts)
    }
    success <- TRUE
  } else if (method == "lss") {
    design <- list()
    for (i in 1:nrow(trials)) {
      design[[i]] <- cbind(hrf$ts[ ,i], rowSums(hrf$ts[ ,-i]))
      if (!is.null(noise)) {
        design[[i]] <- cbind(design[[i]], noise.ts)
      }
    }
    success <- TRUE
  } else if (method == "lsr") {
    if (!is.null(nuisance)) {
      # nuisance.ts <- matrix(0, nrow=n.tr, ncol=length(nuisance))
      nuisance.ts <- vector("list", length(nuisance))
      for (i in 1:length(nuisance)) {
        nuisance.temp <- convolve.hrf(onset=nuisance[[i]][ ,1],
                                      duration=nuisance[[i]][ ,2],
                                      amplitude=nuisance[[i]][ ,3],
                                      tr=tr, n.tr=n.tr,
                                      st.betas=TRUE, scale=FALSE)
        nuisance.ts[[i]] <- nuisance.temp$ts
      }
    }

    design <- list()
    for (i in 1:nrow(trials)) {
      design[[i]] <- cbind(hrf$ts[ ,i], rowSums(hrf$ts[ ,-i]))
      if (!is.null(nuisance)) {
        for (j in 1:length(nuisance)) {
          design[[i]] <- cbind(design[[i]], rowSums(nuisance.ts[[j]][ ,-i]))
        }
      }
      if (!is.null(noise)) {
        design[[i]] <- cbind(design[[i]], noise.ts)
      }
    }
    success <- TRUE
  } else {
    success <- FALSE
  }

  if (success) {
    if (save.design) {
      save("design", file=save.name)
    }
    return(design)
  }
}