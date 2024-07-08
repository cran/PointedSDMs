## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## ----Install PointedSDMs, warning = FALSE, message = FALSE, eval = TRUE-------

##Install if need be
library(PointedSDMs)


## ----startISDM----------------------------------------------------------------
#  
#  args(startISDM)
#  

## ----startSpecies-------------------------------------------------------------
#  
#  args(startSpecies)
#  

## ----fitISDM------------------------------------------------------------------
#  
#  args(fitISDM)
#  

## ----args for blockedCV-------------------------------------------------------
#  
#  args(blockedCV)
#  

## ----datasetOut---------------------------------------------------------------
#  
#  args(datasetOut)
#  

## ----Load packages, message=FALSE, warning=FALSE------------------------------
#  
#  library(INLA)
#  library(inlabru)
#  library(USAboundaries)
#  library(sf)
#  library(blockCV)
#  library(ggmap)
#  library(sn)
#  library(terra)
#  library(RColorBrewer)
#  library(cowplot)
#  library(knitr)
#  library(kableExtra)
#  library(dplyr)
#  library(spocc)
#  

## ----Map of PA----------------------------------------------------------------
#  
#  proj <- "+proj=utm +zone=17 +datum=WGS84 +units=km"
#  
#  PA <- USAboundaries::us_states(states = "Pennsylvania")
#  
#  PA <- st_transform(PA, proj)
#  

## ----get_eBird----------------------------------------------------------------
#  
#  species <- c('caerulescens', 'fusca', 'magnolia')
#  
#  dataSets <- list()
#  for (bird in species) {
#  
#    raw_data <- spocc::occ(
#                query = paste('Setophaga', bird),
#                from = "gbif",
#                date = c("2005-01-01", "2005-12-31"),
#                geometry = st_bbox(st_transform(PA,
#                  '+proj=longlat +datum=WGS84 +no_defs')))$gbif
#  
#    rows <- grep("EBIRD", raw_data$data[[paste0('Setophaga_', bird)]]$collectionCode)
#  
#    raw_data <- data.frame(raw_data$data[[1]][rows, ])
#    raw_data$Species_name <- rep(bird, nrow(raw_data))
#  
#    data_sp <- st_as_sf(
#      x = raw_data[, names(raw_data) %in% c("longitude", "latitude", 'Species_name')],
#      coords = c('longitude', 'latitude'),
#      crs = '+proj=longlat +datum=WGS84 +no_defs')
#    data_sp <- st_transform(data_sp, proj)
#  
#    dataSets[[paste0('eBird_', bird)]] <- data_sp[unlist(st_intersects(PA, data_sp)),]
#  
#    }
#  

## ----Load points--------------------------------------------------------------
#  
#  data('SetophagaData')
#  dataSets[['BBA']] <- SetophagaData$BBA
#  dataSets[['BBS']] <- SetophagaData$BBS
#  

## ----Covariate data, message = FALSE, warning = FALSE-------------------------
#  
#  covariates <- scale(terra::rast(system.file('extdata/SetophagaCovariates.tif',
#                                        package = "PointedSDMs")))
#  names(covariates) <- c('elevation', 'canopy')
#  
#  plot(covariates)
#  

## ----Mesh, warning = FALSE, message = FALSE, fig.width=8, fig.height=5--------
#  
#  mesh <- inla.mesh.2d(boundary = inla.sp2segment(PA),
#                       cutoff = 0.2 * 5,
#                       max.edge = c(0.1, 0.24) * 40, #120
#                       offset = c(0.1, 0.4) * 100,
#                       crs = st_crs(proj))
#  
#  mesh_plot <- ggplot() +
#               gg(mesh) +
#               ggtitle('Plot of mesh') +
#               theme_bw() +
#               theme(plot.title = element_text(hjust = 0.5))
#  mesh_plot
#  

## ----modelOptions-------------------------------------------------------------
#  
#  modelOptions <- list(control.inla =
#                         list(int.strategy = 'eb',
#                              diagonal = 0.1),
#                              verbose = TRUE,
#                              safe = TRUE)
#  

## ----Model prep, warning = FALSE, message = FALSE-----------------------------
#  
#  caerulescensData <- dataSets[c(1,4,5)]
#  
#  caerulescensModel <- startISDM(caerulescensData, Boundary = PA,
#                            Projection = proj, Mesh = mesh,
#                            responsePA = 'NPres', responseCounts = 'Counts',
#                            spatialCovariates = covariates,
#                            Formulas =
#                            list(
#            covariateFormula = ~ elevation + I(elevation^2) + canopy + I(canopy^2))
#                               )
#  

## ----help, eval = FALSE-------------------------------------------------------
#  
#  caerulescensModel$help()
#  

## ----dataset plot, fig.width=8, fig.height=5----------------------------------
#  
#  caerulescensModel$plot() +
#    theme_bw() +
#    ggtitle('Plot of the datasets') +
#    theme(plot.title = element_text(hjust = 0.5))
#  

## ----specifySpatial-----------------------------------------------------------
#  
#  caerulescensModel$specifySpatial(sharedSpatial = TRUE,
#                                   prior.sigma = c(1, 0.1),
#                                   prior.range = c(15, 0.1))
#  

## ----bias fields, eval = FALSE------------------------------------------------
#  
#  caerulescensModel$addBias(datasetNames = 'eBird_caerulescens')
#  
#  caerulescensModel$specifySpatial(Bias = TRUE,
#                                   prior.sigma = c(1, 0.1),
#                                   prior.range = c(15, 0.1))
#  

## ----priorsFixed--------------------------------------------------------------
#  
#  caerulescensModel$priorsFixed(Effect = 'Intercept',
#                                mean.linear = 0,
#                                prec.linear = 0.1)
#  

## ----changeComponents---------------------------------------------------------
#  
#  caerulescensModel$changeComponents()
#  

## ----specifyRandom------------------------------------------------------------
#  
#  caerulescensModel$specifyRandom(copyModel = list(beta = list(fixed = TRUE)))
#  
#  caerulescensModel$changeComponents()
#  

## ----fitISDM run--------------------------------------------------------------
#  
#  caerulescensEst <- fitISDM(data = caerulescensModel,
#                     options = modelOptions)
#  

## ----predict and plot---------------------------------------------------------
#  
#  caerulescensPredictions <- predict(caerulescensEst,
#                                     data = fm_pixels(mesh = mesh,
#                                                      mask = PA),
#                                     spatial = TRUE,
#                                     n.samples = 1000)
#  
#  plot(caerulescensPredictions, variable = c('mean', 'sd'))
#  
#  

## ----startSpeciesStart--------------------------------------------------------
#  
#  speciesModel <- startSpecies(dataSets, Boundary = PA, pointsSpatial = NULL,
#                               Projection = proj, Mesh = mesh,
#                               responsePA = 'NPres', responseCounts = 'Counts',
#                               spatialCovariates = covariates,
#                               speciesName = 'Species_name')
#  

## ----species help, eval = FALSE-----------------------------------------------
#  
#  speciesModel$help()
#  

## ----specifySpecies-----------------------------------------------------------
#  
#  speciesModel$specifySpatial(Species  = TRUE,
#                              prior.sigma = c(1, 0.1),
#                              prior.range = c(15, 0.1))
#  
#  speciesModel$priorsFixed(Effect = 'Intercept',
#                           mean.linear = 0,
#                           prec.linear = 0.1)
#  
#  speciesModel$specifyRandom(speciesGroup = list(model = "iid",
#                                                 hyper = list(prec = list(prior = "pc.prec",
#                                                 param = c(0.1, 0.1)))),
#                             speciesIntercepts = list(prior = 'pc.prec',
#                                                      param = c(0.1, 0.1)))
#  
#  

## ----fitSpecies---------------------------------------------------------------
#  
#  speciesEst <- fitISDM(data = speciesModel,
#                        options = modelOptions)
#  
#  summary(speciesEst)
#  

## ----predictionsSpecies-------------------------------------------------------
#  
#  speciesPredictions <- predict(speciesEst,
#                                     data = fm_pixels(mesh = mesh,
#                                                      mask = PA),
#                                     spatial = TRUE,
#                                     n.samples = 1000)
#  
#  plot(speciesPredictions)
#  

## ----spatialBlock, warning = FALSE, message = FALSE,  fig.width=8, fig.height=5----
#  
#  caerulescensModel$spatialBlock(k = 2, rows_cols = c(2, 2), plot = TRUE) + theme_bw()
#  

## ----blockedCV, warning = FALSE, eval = FALSE---------------------------------
#  
#  spatialBlocked <- blockedCV(data = caerulescensModel,
#                              options = modelOptions)
#  

## ----print spatialBlocked-----------------------------------------------------
#  
#  spatialBlocked
#  

## ----No fields model, message = FALSE, warning = FALSE------------------------
#  
#  no_fields <- startISDM(dataSets,
#                        pointsSpatial = NULL,
#                        Projection = proj, Mesh = mesh,
#                        responsePA = 'NPres', responseCounts = 'Counts',
#                        spatialCovariates = covariates)
#  
#  no_fields$spatialBlock(k = 2, rows_cols = c(2, 2), plot = TRUE) + theme_bw()
#  

## ----spatialBlocked_no_fields, eval = FALSE-----------------------------------
#  
#  spatialBlocked_no_fields <- blockedCV(data = no_fields,
#                                        options = modelOptions)
#  

## ----print spatialBlocked_no_fields-------------------------------------------
#  
#  spatialBlocked_no_fields
#  

## ----Summary of model, message = FALSE, warning = FALSE, echo = TRUE,fig.width=7, fig.height=5----
#  
#  results_plot <- joint_model$summary.fixed %>%
#                  mutate(species = gsub('_.*$','',
#                                        row.names(joint_model$summary.fixed))) %>%
#                  mutate(coefficient = row.names(joint_model$summary.fixed))
#  
#  
#  coefficient_plot <- ggplot(results_plot, aes(x = coefficient, y = mean)) +
#                      geom_hline(yintercept = 0, colour = grey(0.25), lty = 2) +
#                      geom_point(aes(x = coefficient,
#                                     y = mean)) +
#                      geom_linerange(aes(x = coefficient,
#                                         ymin = `0.025quant`,
#                                         ymax = `0.975quant`,
#                                         col = species),
#                                         lwd = 1) +
#                                         theme_bw() +
#                      scale_colour_manual(values = c('#003f5c', '#bc5090','#ffa600')) +
#                      theme(legend.position="bottom",
#                      plot.title = element_text(hjust = 0.5)) +
#                      ggtitle("95% credibility intervals of the fixed effects\n
#                              for the three studied species") +
#                      labs(x = 'Variable', y = 'Coefficient value') +
#                      coord_flip()
#  
#  coefficient_plot
#  

## ----Leave one out, message = FALSE, warning = FALSE, eval = FALSE------------
#  
#  dataset_out <- datasetOut(model = caerulescensEst,
#                            dataset = "BBA",
#                            predictions = TRUE)
#  
#  dataset_out
#  

