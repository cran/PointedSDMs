## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## ----Install PointedSDMs, warning = FALSE, message = FALSE, eval = TRUE-------

##Install if need be
library(PointedSDMs)


## ----intModel-----------------------------------------------------------------
#  
#  args(intModel)
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
#  library(sp)
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
#  proj <- '+proj=longlat +datum=WGS84 +no_defs'
#  PA <- USAboundaries::us_states(states = "Pennsylvania")
#  #PA <- PA$geometry
#  PA <- as(PA, "sf")
#  st_crs(PA) <- proj
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
#                geometry = PA)$gbif
#  
#    rows <- grep("EBIRD", raw_data$data[[paste0('Setophaga_', bird)]]$collectionCode)
#  
#    raw_data <- data.frame(raw_data$data[[1]][rows, ])
#    raw_data$Species_name <- rep(bird, nrow(raw_data))
#  
#    data_sp <- st_as_sf(
#      x = raw_data[, names(raw_data) %in% c("longitude", "latitude", 'Species_name')],
#      coords = c('longitude', 'latitude'),
#      crs = proj)
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
#  crs(covariates) <- proj
#  

## ----Mesh, warning = FALSE, message = FALSE, fig.width=8, fig.height=5--------
#  
#  mesh <- inla.mesh.2d(boundary = inla.sp2segment(PA),
#                       cutoff = 0.2,
#                       max.edge = c(0.1, 0.24),
#                       offset = c(0.1, 0.4),
#                       crs = st_crs(proj))
#  
#  mesh_plot <- ggplot() +
#               gg(mesh) +
#               ggtitle('Plot of mesh') +
#               theme_bw() +
#               theme(plot.title = element_text(hjust = 0.5))
#  mesh_plot
#  

## ----Model prep, warning = FALSE, message = FALSE-----------------------------
#  
#  spatial_data <- intModel(dataSets,
#                          Coordinates = c('X', 'Y'),
#                          Projection = proj, Mesh = mesh,
#                          responsePA = 'NPres', responseCounts = 'Counts',
#                          spatialCovariates = covariates, speciesName = 'Species_name')
#  

## ----dataset plot, fig.width=8, fig.height=5----------------------------------
#  
#  spatial_data$plot(Boundary = FALSE) +
#    geom_sf(data = PA, fill = 'black', alpha = 0.15, lwd = 1.2) +
#    theme_bw() +
#    ggtitle('Plot of the datasets') +
#    theme(plot.title = element_text(hjust = 0.5))
#  

## ----species plot, fig.width=8, fig.height=5----------------------------------
#  
#  spatial_data$plot(Species = TRUE, Boundary = FALSE) +
#    geom_sf(data = PA, fill = 'black', alpha = 0.15, lwd = 1.2) +
#    theme_bw() +
#    ggtitle('Plot of the species') +
#    theme(plot.title = element_text(hjust = 0.5))
#  

## ----specifySpatial-----------------------------------------------------------
#  
#  spatial_data$specifySpatial(sharedSpatial = TRUE,
#                              prior.sigma = c(5, 0.01),
#                              prior.range = c(1, 0.01))
#  

## ----bias fields, eval = FALSE------------------------------------------------
#  
#  spatial_data$addBias('eBird_caerulescens')
#  spatial_data$addBias('eBird_fusca')
#  spatial_data$addBias('eBird_magnolia')
#  

## ----priorsFixed--------------------------------------------------------------
#  
#  spatial_data$priorsFixed(Effect = 'elevation', Species = 'fusca',
#                           mean.linear = 2, prec.linear = 0.05)
#  

## ----changeComponents---------------------------------------------------------
#  
#  spatial_data$changeComponents()
#  

## ----spatialBlock, warning = FALSE, message = FALSE,  fig.width=8, fig.height=5----
#  
#  spatial_data$spatialBlock(k = 4, rows_cols = c(2, 2), plot = TRUE) + theme_bw()
#  

## ----blockedCV, warning = FALSE, eval = FALSE---------------------------------
#  
#  spatialBlocked <- blockedCV(data = spatial_data, options = list(control.inla = list(int.strategy = 'eb')))
#  

## ----print spatialBlocked-----------------------------------------------------
#  
#  spatialBlocked
#  

## ----No fields model, message = FALSE, warning = FALSE------------------------
#  
#  no_fields <- intModel(dataSets,
#                        Coordinates = c('X', 'Y'),
#                        pointsSpatial = NULL,
#                        Projection = proj, Mesh = mesh,
#                        responsePA = 'NPres', responseCounts = 'Counts',
#                        spatialCovariates = covariates, speciesName = 'Species_name')
#  
#  no_fields$spatialBlock(k = 4, rows = 2, cols = 2)
#  

## ----spatialBlocked_no_fields, eval = FALSE-----------------------------------
#  
#  spatialBlocked_no_fields <- blockedCV(data = no_fields, options = list(control.inla = list(int.strategy = 'eb')))
#  

## ----print spatialBlocked_no_fields-------------------------------------------
#  
#  spatialBlocked_no_fields
#  

## ----Running model, message=FALSE, warning=FALSE, eval = FALSE----------------
#  
#  joint_model <- fitISDM(data = spatial_data,
#                         options = list(control.inla = list(int.strategy = 'eb')))
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
#  dataset_out <- datasetOut(model = joint_model,
#                            dataset = "BBA",
#                            predictions = TRUE)
#  
#  dataset_out
#  

## ----Projections, message = FALSE, warning = FALSE, eval = FALSE--------------
#  
#  projections <- predict(joint_model, mesh = mesh, mask = PA,
#                         spatial = TRUE,
#                         fun = 'linear', n.samples = 1000)
#  

## ----Plots, fig.width=8, fig.height=5, message = FALSE, warning = FALSE-------
#  
#  plot(projections, whattoplot = 'mean',
#       colourLow = 'orange', colourHigh = 'dark red')
#  

