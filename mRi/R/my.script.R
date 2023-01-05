my.script <- function(documentation = list(user.name = "Anonymous", description = "No Description"),
                      my.data,
                      script.name,
                      save.dir,
                      model.function = "lmer",
                      model.formula,
                      model.opts = NULL,
                      data.subset = NULL,
                      data.functions = NULL,
                      outliers.extreme = TRUE,
                      outliers.group = TRUE,
                      outliers.model = TRUE) {
  #-----------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-----------------------------------------------------------------------------

  # Debug ----
  # rm(list=ls())
  # gc()
  # script.name <- "scriptTest.R"
  # save.dir <- "C:/Users/tim.koscik/Dropbox (Personal)/Programs/script.write.test"

  # ----

  # fname <- paste(save.dir, script.name, sep="/") # deal with file extension
  # count <- 0
  # while (file.exists(fname)) {
#     count <- count + 1
#     fname <- paste0(fname, count)
#   }
#   fid <- file(fname, "wb")
#   writeLines("#===============================================================================", fid)
#   writeLines(sprintf("# This script was generated on %s", date()), fid)
#   writeLines(sprintf("# This script is provided without warranty and without licence. All rights are reserved."), fid)
#   writeLines("#===============================================================================", fid)
#   close(fid)

}
