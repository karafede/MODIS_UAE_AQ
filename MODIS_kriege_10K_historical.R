
library(RCurl)
library(stringr)
library(plyr)
library(dplyr)
library(threadr)
library(gdalUtils)
library(rgdal)
library(raster)
library(RNetCDF)
library(readr)
library(gstat)
library(curl)
library(leaflet)
library(webshot)
library(htmlwidgets)

### data ectraction from HDF files -----------------------------------------------


# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2017_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2017_Terra_Aqua")

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2016_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2016_Terra_Aqua")

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2015_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2015_Terra_Aqua")

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2014_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2014_Terra_Aqua")

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2013_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2013_Terra_Aqua")

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2012_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2012_Terra_Aqua")

setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2011_Terra_Aqua")
wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2011_Terra_Aqua")

# list data for each directory
DAYS <- str_sub(list.dirs(), start = 3, end = -1)
DAYS <- DAYS[-1]

# make a look for each directory that corresponds to each day

#-----START of the LOOP for all files------------------------
# file <- filenames[1]
# i <- 1

 for (i in 1:length(DAYS)) {
  date <- DAYS[i]
  setwd(paste0(wd,"/",DAYS[i]))
  filenames <- list.files(pattern = "\\.hdf$") 

  
  ###############################################################################
  # collate the tiles together
  
  filenames <- list.files(pattern = "\\.hdf$")
  
  dir <- "/disk3/fkaragulian/MODIS_AOD/UAE_boundary"
  ### shapefile for UAE
  shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")
  
  # ----- Transform to EPSG 4326 - WGS84 (required)
  shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
  # names(shp)
  
  shp_UAE@data$name <- 1:nrow(shp_UAE)
  # plot(shp_UAE)
  
  
  # # make a function to generate all .csv files with Lon, lat and value ##----
  ######-----------------------------------------------------------------------
  
  
  extract_HDF <- function (file) {      ## this is the filenames 
    
    # get list of field names
    nome <- str_sub(file, start = 1, end = -9)
   
    # working with FWtools library--------
    system(paste0('gdal_translate_FWT ', get_subdatasets(file)[64],' AOD.tif'))
    
    # AOD 10km 
    r <- raster("AOD.tif")
    values <- rasterToPoints(r)
    colnames(values) <- c("x", "y", "values")
    values <- as.data.frame (values) 
    values$values <- (values$values)/1000
    
    system(paste0('gdal_translate_FWT ', get_subdatasets(file)[71],' lon.tif'))
    
    # lon 10km
    lon <- raster(paste0("lon.tif"))
    # data values for longitude
    longitude <- rasterToPoints(lon)  
    colnames(longitude) <- c("x", "y", "lon")
    longitude <- as.data.frame (longitude)
    
    
   #  lat 10km
    system(paste0('gdal_translate_FWT ', get_subdatasets(file)[72],' lat.tif'))
    
    lat <- raster(paste0("lat.tif"))
    # data values for longitude
    latitude <- rasterToPoints(lat)  
    colnames(latitude) <- c("x", "y", "lat")
    latitude <- as.data.frame (latitude)
    
    # Join  lat, lon 
    Lat_Lon <- latitude %>% 
      inner_join(longitude, c("x", "y"))
    
    Lat_Lon_Values <- Lat_Lon %>% 
      inner_join(values, c("x", "y"))  
    
    MODIS_data <- Lat_Lon_Values %>%
      dplyr:: select(lon, lat, values)
    MODIS_data <- na.omit(MODIS_data)
    
    write.csv(MODIS_data, file = paste(nome,".csv", sep = ""), row.names=FALSE)
    
  }  
  
  
  BBB <- lapply(filenames, extract_HDF)
  
  # delete HDF files
  if (file.exists(filenames)) file.remove(filenames)  
  
  ######################################################################################
  ######################################################################################
  
  # collate the tiles together ####-------------------------------------
  
  filenames_tiles <- list.files(pattern = "\\.csv$")
  
  LAT = NULL
  LON = NULL
  aod = NULL
  
  
  ## Bind all data together 
  for (i in 1:length(filenames_tiles)) {
    lon <- read_csv(filenames_tiles[i])[,1]
    lat <- read_csv(filenames_tiles[i])[,2]
    AOD <- read_csv(filenames_tiles[i])[,3]
    LON = rbind(LON, data.frame(lon))
    LAT = rbind(LAT, data.frame(lat))
    aod = rbind(aod, data.frame(AOD))
  }
  
  MODIS04_data <- cbind(LON, LAT, aod)
  MODIS04_data <- subset(MODIS04_data, !is.na(values) & !lat == -999 & !lon == -999)
  
  write.csv(MODIS04_data, paste0("AOD_MOD04_10km_UAE","_",date,".csv"))
  
  head(MODIS04_data)
  
  
  ################################################################################
  # subset data for a selected region
  # UAE
  
  DATA_EMIRATES <- subset(MODIS04_data, lon <= 80 & lon >= 10 & lat >= 10 & lat <= 55)
  
  head(DATA_EMIRATES)
  
  write.csv(DATA_EMIRATES, paste0("AOD_MOD04_10km_UAE","_",date,".csv"))
  DATA_EMIRATES <- read.csv(paste0("AOD_MOD04_10km_UAE","_",date,".csv"))[-1]
  MODIS04_data <- read.csv(paste0("AOD_MOD04_10km_UAE","_",date,".csv"))[2:4]
  
  # converta data into a spatialpolygon dataframe------------------------------------
  
  DATA_EMIRATES$x <- DATA_EMIRATES$lon
  DATA_EMIRATES$y <- DATA_EMIRATES$lat
  
  coordinates(DATA_EMIRATES) = ~x + y  ## Set spatial coordinates to create a Spatial object:
  
  
  
  # make a variogram----------------------------------------------------------------
  
  vargram_AOD <- variogram(values ~ 1, DATA_EMIRATES) # calculates sample variogram values 
  
  # fite the variogram
  vargram_AOD_fit  <- fit.variogram(vargram_AOD, fit.ranges = FALSE, fit.sills = FALSE,
                                    vgm(1, c("Sph", "Exp", "Mat")), fit.kappa = TRUE)
  
 # plot(vargram_AOD, vargram_AOD_fit) # plot the sample values, along with the fit model
  
  
  # make a regular empty grid
  x.range <- as.numeric(c(10, 85))  # min/max longitude of the interpolation area
  y.range <- as.numeric(c(-1.7, 60))  # min/max latitude of the interpolation area
  
  ## grid at 10km resolution
  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
                     y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE

  f.1 <- as.formula(Precip_in ~ X + Y)
  # perform kriging
  dat.krg <- gstat::krige(values ~ 1, DATA_EMIRATES, grd, vargram_AOD_fit, nmax = 50)
  r <- raster(dat.krg)
  projection(r) <- CRS("+proj=longlat +datum=WGS84")
  r <- crop(r, extent(shp_UAE))
  r <- mask(r, shp_UAE)
  
  
  writeRaster(r, paste0("AOD_MOD04_10km_UAE","_",date,".tif"), overwrite = TRUE)
  Data_Emirates_10km_tif <- raster(paste0("AOD_MOD04_10km_UAE","_",date,".tif"))
  
  #########-------------------------------------------------------------------------
  # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
  
  PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
  projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
  
  ### crop raster over the UAE shp file  ###############################
  PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
  PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
  
  file.tiff_GWR_PM25_1km_2015 <- 'GWR_PM25_1km_2015.tif'
  GWR_PM25_1km_2015_tiff <- writeRaster(PM25_2015_GRW, filename = file.tiff_GWR_PM25_1km_2015, format = 'GTiff', overwrite = T)

  
  ###################
  
  ## make resolution of MODIS-data as the one of GWR-------------------------------
  Data_Emirates_tif_1km = projectRaster(Data_Emirates_10km_tif, PM25_2015_GRW)
  Data_Emirates_tif_1km <- crop(Data_Emirates_tif_1km, extent(shp_UAE))
  Data_Emirates_tif_1km <- mask(Data_Emirates_tif_1km, shp_UAE)

  file.Emirates_tif_1km <- paste0("AOD_MOD04_1km_UAE","_",date,".tif")
  Emirates_tif_1km_tiff <- writeRaster(Data_Emirates_tif_1km, filename = file.Emirates_tif_1km, format = 'GTiff', overwrite = T)

  ### Extract points from raster tiff ############################################
  
  Emirates_tif_1km_pts <- rasterToPoints(Emirates_tif_1km_tiff)
  colnames(Emirates_tif_1km_pts) <- c("Lon", "Lat", "AOD_1km")
  Emirates_tif_1km_pts <- as.data.frame(Emirates_tif_1km_pts) 
  Emirates_tif_1km_pts <- subset(Emirates_tif_1km_pts, !is.na(AOD_1km) & AOD_1km>0)
  
  write_csv(Emirates_tif_1km_pts , paste0("Emirates_tif_1km.csv","_",date,".csv"))
  
  
  ###########################################################################
  #### make Leafleft map with current-AOD data 1km ###-----------------------
  
  library(leaflet)
  
 
  # define color palette
  rast_pal_EMIRATES <- colorNumeric(c("#9999FF", "#ffd699", "#FFFF00", "#ffbf00", "#ffc700", "#FF0000", "#994c00"),
                                    getValues(Emirates_tif_1km_tiff),na.color = "transparent")
  
  
  
  map <- leaflet() %>% 
    # setView(45, 25, 5) %>%
    addTiles(group = "OSM (default)") %>%
    addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
    addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
    addRasterImage(Emirates_tif_1km_tiff, 
                   colors = rast_pal_EMIRATES, 
                   opacity = 0.6,
                   group = "AOD_EMIRATES") %>%
    addLegend("bottomright", pal = rast_pal_EMIRATES, values = values(Emirates_tif_1km_tiff),# values = c(MIN_data, MAX_data),
              title = paste("<br><strong>AOD (MODIS-5km)",date,"</strong>"),
              labFormat = labelFormat(prefix = ""),
              opacity = 0.6) %>%
    addPolygons(stroke = TRUE, data = shp_UAE,
                weight = 1.5, color = "black",
                fillOpacity = 0,
                group = "shape_UAE") %>%
    addLayersControl(
      baseGroups = c("Road map", "Topographical", "Satellite", "Toner Lite"),
      overlayGroups = c("AOD_EMIRATES"),
      options = layersControlOptions(collapsed = TRUE)) 
  # map
  
  # save map
  saveWidget(map, paste0("AOD_Emirates","_",date,".html"), selfcontained = FALSE)

}  # END of the LOOP!
  
  
  
  
  
  
  
  