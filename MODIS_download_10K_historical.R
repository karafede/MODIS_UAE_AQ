
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

# setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2016_Terra_Aqua")
# wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2016_Terra_Aqua")

setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2015_Terra_Aqua")
wd <- setwd("/disk3/fkaragulian/MODIS_AOD/MODIS_LAADS_NASA/2015_Terra_Aqua")


# list data for each directory
DAYS <- str_sub(list.dirs(), start = 3, end = -1)
DAYS <- DAYS[-1]

# make a look for each directory that corresponds to each day

#-----START of the LOOP for all files------------------------

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
  
  MODIS04_data <- read.csv(paste0("AOD_MOD04_10km_UAE","_",date,".csv"))[2:4]
  DATA_EMIRATES <- subset(MODIS04_data, lon <= 80 & lon >= 10 & lat >= 10 & lat <= 55)
  
  head(DATA_EMIRATES)
  
  write.csv(DATA_EMIRATES, paste0("AOD_MOD04_10km_UAE","_",date,".csv"))
  DATA_EMIRATES <- read.csv(paste0("AOD_MOD04_10km_UAE","_",date,".csv"))[-1]
  
  ############  Make a data interpolation on a regular grid ########################
  
  DATA_EMIRATES$x <- DATA_EMIRATES$lon
  DATA_EMIRATES$y <- DATA_EMIRATES$lat
  
  coordinates(DATA_EMIRATES) = ~x + y  ## Set spatial coordinates to create a Spatial object:
  
  x.range <- as.numeric(c(20, 85))  # min/max longitude of the interpolation area
  y.range <- as.numeric(c(-1.7, 60))  # min/max latitude of the interpolation area
  
  ## grid at 5km resolution
  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.05),
                     y = seq(from = y.range[1], to = y.range[2], by = 0.05))  # expand points to grid
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE
  
  # plot(grd, cex = 1.5, col = "grey")
  # points(DATA_EMIRATES, pch = 1, col = "red", cex = 1)
  
  idw <- idw(formula = values ~ 1, locations = DATA_EMIRATES, 
             newdata = grd)  # apply idw model for the data (interpolation)
  
  idw.output = as.data.frame(idw)  # output is defined as a data table
  names(idw.output)[1:3] <- c("Lon", "Lat", "values")  # give names to the modelled variables
  
 # write.csv(idw.output, file = paste0("AOD_MOD04_5km_UAE","_",date,"_interp.csv"), row.names=FALSE)
  
  ##### make raster with interpolated data----------------------------------------
  
#  DATA_EMIRATES_interp <- read.csv(paste0("AOD_MOD04_5km_UAE","_",date,"_interp.csv"))
  
  DATA_EMIRATES_interp <- subset(idw.output, Lon <= 57 & Lon >= 50 & Lat >= 21 & Lat <= 56)
  write.csv(DATA_EMIRATES_interp, file = paste0("AOD_MOD04_5km_UAE","_",date,"_interp.csv"), row.names=FALSE)
  
  coordinates(DATA_EMIRATES_interp) <- ~ Lon + Lat
  # coerce to SpatialPixelsDataFrame
  gridded(DATA_EMIRATES_interp) <- TRUE
  raster_DATA_EMIRATES_interp <- raster(DATA_EMIRATES_interp)
  projection(raster_DATA_EMIRATES_interp) <- CRS("+proj=longlat +datum=WGS84")
  # plot(raster_DATA_EMIRATES_interp)
  
  raster_DATA_EMIRATES_interp <- crop(raster_DATA_EMIRATES_interp, extent(shp_UAE))
  raster_DATA_EMIRATES_interp <- mask(raster_DATA_EMIRATES_interp, shp_UAE)
  
  writeRaster(raster_DATA_EMIRATES_interp, paste0("AOD_MOD04_5km_UAE","_",date,".tif"), overwrite = TRUE)
  
  Data_Emirates_tif <- raster(paste0("AOD_MOD04_5km_UAE","_",date,".tif"))
  # plot(Data_Emirates_tif)
  
  ##########################################################################
  ##########################################################################
  
  # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
  
  PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
  projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
  
  ### crop raster over the UAE shp file  ###############################
  PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
  PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
  
  file.tiff_GWR_PM25_1km_2015 <- 'GWR_PM25_1km_2015.tif'
  GWR_PM25_1km_2015_tiff <- writeRaster(PM25_2015_GRW, filename = file.tiff_GWR_PM25_1km_2015, format = 'GTiff', overwrite = T)

  # check resoltuion of Donkelaar data and of ECMWF data
  res(GWR_PM25_1km_2015_tiff) # 0.01
  res(Data_Emirates_tif) # 0.1
  
  ###################
  
  ## make resolution of MODIS-data as the one of GWR-------------------------------
  Data_Emirates_tif_1km = projectRaster(Data_Emirates_tif, PM25_2015_GRW)
  Data_Emirates_tif_1km <- crop(Data_Emirates_tif_1km, extent(shp_UAE))
  Data_Emirates_tif_1km <- mask(Data_Emirates_tif_1km, shp_UAE)

  file.Emirates_tif_1km <- paste0("AOD_MOD04_1km_UAE","_",date,".tif")
  Emirates_tif_1km_tiff <- writeRaster(Data_Emirates_tif_1km, filename = file.Emirates_tif_1km, format = 'GTiff', overwrite = T)

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

}  # END of the LOOP
  
  
  
  
  
  
  
  