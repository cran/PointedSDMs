% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockedCV.R
\name{blockedCV}
\alias{blockedCV}
\title{\emph{blockedCV}: run spatial blocked cross-validation on the integrated model.}
\usage{
blockedCV(data, options = list())
}
\arguments{
\item{data}{An object produced by either \code{\link{startISDM}} of \code{\link{startSpecies}}. Requires the slot function, \code{.$spatialBlock} to be run first in order to specify how the data in the model is blocked.}

\item{options}{A list of \pkg{INLA} or \pkg{inlabru} options to be used in the model. Defaults to \code{list()}.}
}
\value{
An object of class \code{blockedCV}, which is essentially a list of DIC values obtained from each iteration of the model.
}
\description{
This function is used to perform spatial blocked cross-validation with regards to model selection for the integrated model. It does so by leaving out a block of data in the full model, running a model with the remaining data, and then calculating the deviance information criteria (DIC) as a score of model fit.
}
\examples{

\dontrun{
 if(requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                            responsePA = 'Present',
                            Projection = proj)
 
 #Set up spatial block
 organizedData$spatialBlock(k = 2, rows = 2, cols = 1)
 
 #Run spatial block cross-validation
 blocked <- blockedCV(organizedData)
 
 #Print summary
 blocked
   
 }
}

}
