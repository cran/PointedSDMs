## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## ----setup, warning = FALSE, message = FALSE----------------------------------
#  
#  library(PointedSDMs)
#  library(terra)
#  library(ggpolypath)
#  library(INLA)
#  library(ggplot2)
#  

## ----safe, include = FALSE----------------------------------------------------
#  
#  bru_options_set(inla.mode = "experimental")
#  

## ----load data----------------------------------------------------------------
#  
#  data('SolitaryTinamou')
#  projection <- "+proj=longlat +ellps=WGS84"
#  
#  covariates <- terra::rast(system.file('extdata/SolitaryTinamouCovariates.tif',
#                                        package = "PointedSDMs"))
#  
#  datasets <- SolitaryTinamou$datasets
#  region <- SolitaryTinamou$region
#  mesh <- SolitaryTinamou$mesh
#  

## ----look at data-------------------------------------------------------------
#  
#  str(datasets)
#  class(region)
#  

## ----covariates, fig.width=8, fig.height=5------------------------------------
#  
#  covariates <- scale(covariates)
#  crs(covariates) <- projection
#  plot(covariates)
#  

## ----mesh, fig.width=8, fig.height=5------------------------------------------
#  
#  ggplot() + gg(mesh)
#  

## ----set up base model, warning = FALSE, message = FALSE----------------------
#  
#  base <- intModel(datasets, spatialCovariates = covariates, Coordinates = c('X', 'Y'),
#                   Projection = projection, responsePA = 'Present', Offset = 'area',
#                   Mesh = mesh, pointsSpatial = NULL)
#  

## ----data, fig.width=8, fig.height=5------------------------------------------
#  
#  base$plot(Boundary = FALSE) +
#    geom_sf(data = st_boundary(region)) +
#    ggtitle('Plot of the species locations by dataset')
#  

## ----priorsFixed--------------------------------------------------------------
#  
#  base$priorsFixed(Effect = 'Forest', mean.linear = 0.5, prec.linear = 0.01)
#  

## ----run base model, warning = FALSE, message = FALSE-------------------------
#  
#  baseModel <- fitISDM(data = base)
#  baseModel
#  

## ----set up model with fields, warning = FALSE, message = FALSE---------------
#  
#  fields <- intModel(datasets, Coordinates = c('X', 'Y'),
#                     Projection = projection, Mesh = mesh, responsePA = 'Present',
#                     pointsIntercept = FALSE)
#  

## ----specifySpatial-----------------------------------------------------------
#  
#  fields$specifySpatial(sharedSpatial = TRUE, prior.range = c(50,0.01),
#                        prior.sigma = c(0.1, 0.01))
#  

## ----addBias------------------------------------------------------------------
#  
#  fields$addBias('eBird')
#  

## ----run fields model, warning = FALSE, message = FALSE-----------------------
#  
#  fieldsModel <- fitISDM(fields, options = list(control.inla = list(int.strategy = 'eb')))
#  fieldsModel
#  

## ----copy model---------------------------------------------------------------
#  
#  copy <- intModel(datasets, Coordinates = c('X', 'Y'),
#                   Projection = projection, Mesh = mesh, responsePA = 'Present',
#                   pointsSpatial = 'copy')
#  
#  copy$specifySpatial(sharedSpatial = TRUE, prior.range = c(20,0.01),
#                        prior.sigma = c(0.1, 0.01))
#  
#  copy$changeComponents()
#  

## ----run copy model-----------------------------------------------------------
#  
#  copyModel <- fitISDM(copy, options = list(control.inla = list(int.strategy = 'eb')))
#  copyModel
#  

## ----predict spatial, warning = FALSE, message = FALSE------------------------
#  
#  spatial_predictions <- predict(fieldsModel, mesh = mesh,
#                         mask = region,
#                         spatial = TRUE,
#                         fun = 'linear')
#  

## ----spatial, fig.width=8, fig.height=5---------------------------------------
#  
#  plot(spatial_predictions)
#  

## ----predict bias, warning = FALSE, message = FALSE---------------------------
#  
#  bias_predictions <- predict(fieldsModel,
#                      mesh = mesh,
#                      mask = region,
#                      biasfield = TRUE,
#                      fun = 'linear')
#  

## ----bias, fig.width=8, fig.height=5------------------------------------------
#  
#  plot(bias_predictions)
#  

## ----datasetOut, warning = FALSE, message = FALSE-----------------------------
#  
#  eBird_out <- datasetOut(model = fieldsModel, dataset = 'eBird')
#  

## ----print datasetOut---------------------------------------------------------
#  
#  eBird_out
#  

