\name{nii.qc.snr}
\alias{nii.qc.snr}
\title{Regional Signal-to-Noise Ratio}
\description{
Calculate the signal to noise ratio within the masked region.
\deqn{\displaystyle SNR = \frac{\mu_ROI}{\sigma_ROI\sqrt{\frac{n}{\left(n-1\right)}}}}{SNR = \mu_ROI/[\sigma_ROI\sqrt(n/(n-1))]}
}
\usage{
nii.qc.snr(img.nii, img.vol=1L, mask.nii, mask.vol=1L, mask.dir="gt", mask.thresh=0)
}
\arguments{
  \item{img.nii}{Full directory listing to desired nifti image (unzipped)}
  \item{img.vol=1L}{For multi volume images; an integer indicating the volume to load (1..N)}
  \item{mask.nii}{Full directory listing to desired nifti image to use for masking (unzipped). Must be identical pixel dimensions to base image.}
  \item{mask.vol=1L}{For multi volume images; an integer indicating the volume to load (1..N)}
  \item{mask.dir="gt"}{Direction to threshold the mask image to generate the desired ROI. Options include ("gt") [>] Default, "ge" [>=], "lt" [<], "le" [<=], "eq" [==]}
  \item{mask.thresh=0}{Value at which to threshold the mask image to generate an ROI. E.g., for binary masks, where the value 1 indicates the ROI, default of >0, will generate the appropriate mask}
}
\value{
a numeric value indicating signal-to-noise-ratio
}
\author{
  Timothy R. Koscik <timkoscik+niiqc@gmail.com>
}
\examples{}
