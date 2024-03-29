% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel.R
\name{runModel}
\alias{runModel}
\title{\emph{runModel}: function used to run the integrated model. Note that this function is depreciated, and will be removed in a later version of the package.}
\usage{
runModel(data, options = list())
}
\arguments{
\item{data}{A intModel object to be used in the integrated model.}

\item{options}{A list of INLA options used in the model. Defaults to \code{list()}.}
}
\value{
An inlabru model with additional lists containing some more metadata attached.
}
\description{
This function takes a \code{intModel} object and produces an \code{inlabru} model object with additional lists and meta-data added.
}
\examples{

\dontrun{
 
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
                             Projection = proj, responsePA = 'Present')
 
  ##Run the model
  modelRun <- runModel(organizedData, 
  options = list(control.inla = list(int.strategy = 'eb')))
   
  #Print summary of model
  modelRun
   
 }
}

}
