\name{run.voxel}
\alias{run.voxel}
\title{Run mRi Models Voxelwise}
\description{
A function that handles mRi input (4D, paradigm, masks, etc.) and applies the user defined models (my.function) over all voxels specified by voxels contained within masks.
}
\usage{
  run.voxel(data.4d, data.paradigm, data.mask, my.function, my.command, save.dir, file.prefix, num.cores, log.name, verbose, ...)
}
\arguments{
  \item{data.4d}{
  A character vector specifying the full directory listing to the folder containing 4D *.nii files for processing. Input may also be a character vector with the full directory listing for all 4D files (to allow specifying files in different directories). Input may also be list of directories and file paths.
  }
  \item{data.paradigm}{
  Full directory listing to a *.csv or *.txt file containing information for the paradigm (i.e., all other variables that are desired to be modelled).  If a directory is given, then all *.csv and *.txt files will be read and assumed to be represent data for individual subjects/runs in the same order as 4D files.
  }
  \item{data.mask}{}
  \item{my.function}{
  User specified models to applied to each voxel. Output from statistical models can be output using the table.to.nii function.
  }
  \item{my.command}{
  A character string indicating how to call my.function. Default my.command <- "my.function(df, coords, img.dims, save.dir, prefix)"
  }
  \item{save.dir}{
  Full directory listing of the folder in which to save files
  }
  \item{file.prefix}{
  A character string to append to filenames. filenames will be <file.prefix>.mask<#>.<table.name>.<var.name>.nii
  }
  \item{num.cores}{
  An integer indicating the number of cores to use for parallelization, Default value or NULL uses detectCores() - 1 to determine the number of cores to use. NOT IMPLEMENTED IN WINDOWS
  }
  \item{log.name}{A character string (without file extension) to name the output log, including directory listing.  Default is <save.dir>/log}
  \item{verbose}{}
  \item{...}{}
}
\details{
NIFTI files need to be unzipped before using this function or any other portions of the mRi package. This is necessary given the inconsistent way in which gzipped files are indexed (Some information on this is given in the documentation for the readBin function).

Scripts that call run.model must preload all necessary libraries for whatever functions they need.
}
\value{
Values are output directly to NIFTI files.
}
\author{
  Timothy R. Koscik <timkoscik+mripkg@gmail.com>
}
\note{
BLAS libraries may interfere with parallelization, currently the function detects if the Intel Math Kit (MKL) is installed and will set the number of threads to 1 for the MKL. It is currently not known if other BLAS libraries like OpenBLAS will interefere.
}
\examples{}