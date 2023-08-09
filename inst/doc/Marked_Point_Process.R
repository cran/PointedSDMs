## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  warning = FALSE,
  message = FALSE
)


## ---- setup-------------------------------------------------------------------
#  library(spatstat)
#  library(PointedSDMs)
#  library(sf)
#  library(sp)
#  library(ggplot2)
#  library(inlabru)
#  library(INLA)

## ----load_data----------------------------------------------------------------
#  
#  data(Koala)
#  eucTrees <- Koala$eucTrees
#  boundary <- Koala$boundary
#  

## ---- clean_data, echo = FALSE,fig.width=7, fig.height=5----------------------
#  
#  proj <- "+init=epsg:27700"
#  
#  boundary <- as(boundary, 'sf')
#  st_crs(boundary) <- proj
#  
#  euc <- st_as_sf(x = eucTrees,
#                  coords = c('E', 'N'),
#                  crs = proj)
#  
#  euc$food <- euc$FOOD/1000
#  euc <- euc[unlist(st_intersects(boundary, euc)),]
#  
#  class(trees)
#  
#  mesh = inla.mesh.2d(boundary = boundary, max.edge = 20)
#  mesh$crs <- st_crs(proj)
#  
#  ggplot() +
#    geom_sf(data = st_boundary(boundary)) +
#    geom_sf(data = euc, aes(color = koala)) +
#    ggtitle('Plot showing number of koalas at each site')
#  
#  ggplot() +
#    geom_sf(data = st_boundary(boundary)) +
#    geom_sf(data = euc, aes(color = food)) +
#    ggtitle('Plot showing the food value index at each site')
#  
#  

## ---- analysis_of_data, eval = FALSE, echo= FALSE-----------------------------
#  
#  data(euc) ##will add this in the future when data is on archive
#  

## ---- points_only,fig.width=7, fig.height=5-----------------------------------
#  
#  points <- intModel(euc, Coordinates = c('x', 'y'),
#                    Projection = proj, Mesh = mesh)
#  
#  pointsModel <- fitISDM(points, options = list(control.inla = list(int.strategy = 'eb')))
#  
#  pointsPredictions <- predict(pointsModel, mask = boundary,
#                               mesh = mesh, predictor = TRUE)
#  
#  plot(pointsPredictions)

## ---- include_marks,fig.width=7, fig.height=5---------------------------------
#  
#  marks <- intModel(euc, Coordinates = c('x', 'y'), Projection = proj,
#                    markNames = c('food', 'koala'), markFamily = c('gamma', 'poisson'),
#                    Mesh = mesh)
#  
#  marksModel <- fitISDM(marks, options = list(control.inla = list(int.strategy = 'eb'),
#                                               safe = TRUE))
#  
#  foodPredictions <- predict(marksModel, mask = boundary,
#                             mesh = mesh, marks = 'food', spatial = TRUE)
#  
#  koalaPredictions <- predict(marksModel, mask = boundary,
#                              mesh = mesh, marks = 'koala', spatial = TRUE)
#  
#  plot(foodPredictions)
#  plot(koalaPredictions)

## ---- marks_add_scaling,fig.width=7, fig.height=5-----------------------------
#  
#  marks2 <- intModel(euc, Coordinates = c('x', 'y'), Projection = proj,
#                    markNames = 'food', markFamily = 'gaussian',
#                    Mesh = mesh, pointsSpatial = 'individual')
#  
#  marks2$updateFormula(markName = 'food',
#        newFormula = ~ exp(food_intercept + (euc_spatial + 1e-6)*scaling + food_spatial))
#  
#  marks2$changeComponents(addComponent = 'scaling')
#  
#  marks2$specifySpatial(datasetName = 'euc',
#                        prior.sigma = c(0.1, 0.01),
#                        prior.range = c(10, 0.01))
#  
#  marks2$specifySpatial(Mark = 'food',
#                        prior.sigma = c(0.1, 0.01),
#                        prior.range = c(10, 0.01))
#  
#  marksModel2 <- fitISDM(marks2, options = list(control.inla = list(int.strategy = 'eb'),
#                                                bru_max_iter = 2, safe = TRUE))
#  
#  
#  predsMarks2 <- predict(marksModel2, mask = boundary, mesh = mesh,
#      formula =  ~ (food_intercept + (euc_spatial + 1e-6)*scaling + food_spatial))
#  
#  plot(predsMarks2)

