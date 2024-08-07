% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasetOut.R
\name{datasetOut}
\alias{datasetOut}
\title{\emph{datasetOut}: function that removes a dataset out of the main model, and calculates some cross-validation score.}
\usage{
datasetOut(model, dataset, predictions = TRUE)
}
\arguments{
\item{model}{Model of class modISDM run with multiple datasets.}

\item{dataset}{Names of the datasets to leave out. If missing, will run for all datasets used in the full model.}

\item{predictions}{Will new models be used for predictions. If \code{TRUE} returns marginals and bru_info in model. Defaults to \code{TRUE}.}
}
\value{
A list of inlabru models with the specified dataset left out. If predictions is \code{FALSE}, these objects will be missing their \code{bru_info} and \code{call} lists.
}
\description{
This function calculates the difference in covariate values between a full integrated model and a model with one dataset left out, as well as some cross-validation score, which is used to obtain a score of the relative importance of the dataset in the full model. The score is calculated as follows:
\enumerate{

  \item Running a new model with one less dataset (from the main model) -- resulting in a reduced model,
  \item predicting the intensity function at the locations of the left-out dataset with the reduced model,
  \item using the predicted values as an offset in a new model,
  \item finding the difference between the marginal-likelihood of the main model (ie the model with all the datasets considered) and the marginal-likelihood of the offset model.

}
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
 organizedData <- startISDM(data, Mesh = mesh,
                            Projection = proj, 
                            responsePA = 'Present')
 
  ##Run the model
  modelRun <- fitISDM(organizedData,
              options = list(control.inla = list(int.strategy = 'eb')))
   
  #Choose dataset to leave out
  eBirdOut <- datasetOut(modelRun, dataset = 'eBird')
  
  #Print datasetOut summary
  eBirdOut
  
  
   
 }
}
}
