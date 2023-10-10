## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## ----load_packages------------------------------------------------------------
#  
#  library(PointedSDMs)
#  library(inlabru)
#  library(ggplot2)
#  library(spocc)
#  library(INLA)
#  library(dplyr)
#  library(sp)
#  library(sf)
#  

## ----Alabama_map--------------------------------------------------------------
#  
#  proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#  
#  AL <- USAboundaries::us_states(states = "Alabama", resolution = 'high')
#  AL <- as(AL, "sf")
#  st_crs(AL) <- proj
#  
#  mesh <- inla.mesh.2d(boundary = inla.sp2segment(AL[1]),
#                       cutoff = 0.1,
#                       max.edge = c(0.2, 0.8),
#                       offset = c(0.1, 0.2),
#                       crs = st_crs(proj))
#  
#  
#  

## ----get_routes_data----------------------------------------------------------
#  
#  data("BBSColinusVirginianus")
#  

## ----get_eBird_data-----------------------------------------------------------
#  
#  eBird2015 <- spocc::occ(
#    query = 'Colinus virginianus',
#    from = 'gbif',
#    date = c("2015-01-01", "2015-12-31"),
#    geometry = AL
#  )$gbif
#  
#  eBird2016 <- spocc::occ(
#    query = 'Colinus virginianus',
#    from = 'gbif',
#    date = c("2016-01-01", "2016-12-31"),
#    geometry = AL
#  )$gbif
#  
#  eBird2017 <- spocc::occ(
#    query = 'Colinus virginianus',
#    from = 'gbif',
#    date = c("2017-01-01", "2017-12-31"),
#    geometry = AL
#  )$gbif
#  
#  eBird <- data.frame(eBird2015$data[[1]]) %>%
#    bind_rows(data.frame(eBird2016$data[[1]])) %>%
#    bind_rows(data.frame(eBird2017$data[[1]]))
#  
#  
#  eBird <- st_as_sf(x = eBird,
#                    coords = c('longitude', 'latitude'),
#                    crs = proj)
#  
#  eBird$Year <- eBird$year
#  
#  eBird <- eBird[unlist(st_intersects(AL, eBird)),]
#  

## ----setup_model--------------------------------------------------------------
#  
#  hyperParams <- list(rho = list(prior = "pc.prec", param = c(0.01, 0.01)))
#  
#  modelSetup <- intModel(eBird, BBSColinusVirginianus,
#                         Coordinates = c('Longitude', 'Latitude'),                                                    temporalName = 'Year',
#                         Projection =  proj, Mesh = mesh,
#                         responsePA =  'NPres', trialsPA = 'Ntrials',
#                         temporalModel = list(model = 'ar1', hyper = hyperParams))
#  
#  modelSetup$specifySpatial(sharedSpatial = TRUE, prior.sigma = c(0.2, 0.01),
#                            prior.range = c(0.4, 0.01))
#  

## ----data_plot,fig.width=8, fig.height=5--------------------------------------
#  
#  modelSetup$plot()
#  

## ----model_components---------------------------------------------------------
#  
#  modelSetup$changeComponents()
#  

## ----run_model----------------------------------------------------------------
#  
#  mod <- fitISDM(modelSetup,
#                  options = list(control.inla = list(int.strategy = 'eb')))
#  

## ----predictions, fig.width=8, fig.height=5, message = FALSE------------------
#  
#  preds <- predict(mod, mask = AL, mesh = mesh, temporal = TRUE, fun = '')
#  
#  plot_preds <- plot(preds, whattoplot = 'median', plot = FALSE)
#  
#  plot_preds +
#    geom_sf(data = st_boundary(AL), lwd = 1.2) +
#    scico::scale_fill_scico(palette = "lajolla") +
#    theme_minimal()
#  

