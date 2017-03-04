
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
library(sp)
library(tools)

options(warn=-1)

setwd("/disk3/fkaragulian/MODIS_AOD/")
# setwd("/disk3/fkaragulian/MODIS_AOD/2017/001/")
# folder_day <- as.character("001")

wd <- getwd()
Sys.time()
current_date <- str_sub(Sys.time(), start = 1, end = -10)
str(current_date)

## 3km
# url_list_data = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/"
## 10km MODIS Terra
url_list_data_Terra = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2017/" 
## 10km MODIS Aqua
url_list_data_Aqua = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MYD04_L2/2017/" 


# Terra
list_data_Terra = getURL(url_list_data_Terra, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
list_data_Terra = paste(url_list_data_Terra, strsplit(list_data_Terra, "\r*\n")[[1]], sep = "") 
# last element (date) of the list in order of time 
list_data_Terra <- sort(list_data_Terra)
n <- as.numeric(length(list_data_Terra))
list_data_Terra[n]
folder <- unlist(list_data_Terra[n])

# Acqua
list_data_Aqua = getURL(url_list_data_Aqua, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
list_data_Aqua = paste(url_list_data_Aqua, strsplit(list_data_Aqua, "\r*\n")[[1]], sep = "") 
# last element (date) of the list in order of time 
list_data_Aqua <- sort(list_data_Aqua)
n <- as.numeric(length(list_data_Aqua))
list_data_Aqua[n]
folder <- unlist(list_data_Aqua[n])


# isolate folder name
folder_year <- str_sub(folder, start = 74, end = -5)
folder_day <- str_sub(folder, start = 79, end = -1) 
# folder_day <- "001"

# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/",folder_day, "/", sep = "")
# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/","349", "/", sep = "")

# Terra
url_Terra = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2017/",folder_day, "/", sep = "")
# Aqua
url_Aqua = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MYD04_L2/2017/",folder_day, "/", sep = "")
# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/","350", "/", sep = "")

# Terra
filenames_Terra = getURL(url_Terra, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
filenames_Terra = paste(url_Terra, strsplit(filenames_Terra, "\r*\n")[[1]], sep = "") 

# Aqua
filenames_Aqua = getURL(url_Aqua, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
filenames_Aqua = paste(url_Aqua, strsplit(filenames_Aqua, "\r*\n")[[1]], sep = "") 


# select only files in the time range as in the Overpass time

Overpass_times_Terra <- read.csv(paste0("/disk3/fkaragulian/MODIS_AOD/Overpass_times_Terra_",
                                        current_date,".csv"))
Overpass_times_Aqua <- read.csv(paste0("/disk3/fkaragulian/MODIS_AOD/Overpass_times_Aqua_",
                                       current_date,".csv"))

matches_Terra <- unique (grep(paste(Overpass_times_Terra$x,collapse="|"), 
                        filenames_Terra, value=TRUE))

matches_Aqua <- unique (grep(paste(Overpass_times_Aqua$x,collapse="|"), 
                              filenames_Aqua, value=TRUE))


# create a new working directory for each download date
dir.create(folder_year)
# setwd(folder)

# download data
setwd(paste0(wd,"/",folder_year))
dir.create(folder_day)
setwd(paste0(wd,"/",folder_year,"/", folder_day))

filenames_MODIS_10k_Terra <- unlist(str_extract_all(matches_Terra, ".+(.hdf$)"))
filenames_MODIS_10k_Terra <- sort(filenames_MODIS_10k_Terra)

filenames_MODIS_10k_Aqua <- unlist(str_extract_all(matches_Aqua, ".+(.hdf$)"))
filenames_MODIS_10k_Aqua <- sort(filenames_MODIS_10k_Aqua)

# start downloading data in the main directory -----------------------
mapply(download.file, filenames_MODIS_10k_Terra,basename(filenames_MODIS_10k_Terra)) 
mapply(download.file, filenames_MODIS_10k_Aqua,basename(filenames_MODIS_10k_Aqua)) 


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
##  info <- gdalinfo(file) 
 #  sds <- get_subdatasets(file)

#  AOD <- sds[64]  ## field value AOD AOD_550_Dark_Target_Deep_Blue_Combined # 10km
#  AOD <- sds[12]  ## # AOT Image_Optical_Depth_Land_And_Ocean #3km
  
  # filename <- rasterTmpFile()
  # extension(filename) <- 'tif'
  
 #  gdal_translate(sds[12], dst_dataset = filename) #3km
 #  gdal_translate(sds[64], dst_dataset = filename) #10km
 #  gdal_translate(sds[64], dst_dataset = "/home/fkaragulian/MODIS_AOD/blank.tif") #10km
   
   #  system(paste0('gdal_translate ', get_subdatasets("/disk3/fkaragulian/MODIS_AOD/2016/363/MOD04_L2.A2016363.0730.006.NRT.hdf")[64],' myfile_', band,'.tif'))
   # not working
   # system(paste0('gdal_translate ', get_subdatasets(filenames[1])[64],' myfile_', 64,'.tif'))
  
    # working with FWtools library--------
   system(paste0('gdal_translate_FWT ', get_subdatasets(file)[64],' AOD.tif'))
   
   r <- raster("AOD.tif")
   values <- rasterToPoints(r)
   colnames(values) <- c("x", "y", "values")
   values <- as.data.frame (values) 
   values$values <- (values$values)/1000
  
  
#  lon <- sds[71]  #longitude # 10km
#  lon <- sds[51]  #longitude #3km
#  filename <- rasterTmpFile()
#  extension(filename) <- 'tif'
# gdal_translate(sds[51], dst_dataset = filename) #3km
# gdal_translate(sds[71], dst_dataset = filename) #10km

 system(paste0('gdal_translate_FWT ', get_subdatasets(file)[71],' lon.tif'))

   lon <- raster(paste0("lon.tif"))
  # data values for longitude
  longitude <- rasterToPoints(lon)  
  colnames(longitude) <- c("x", "y", "lon")
  longitude <- as.data.frame (longitude)
  
  
#  lat <- sds[72] # latitude  #10km
#  lat <- sds[52] # latitude  #3km
#  filename <- rasterTmpFile()
#  extension(filename) <- 'tif'
# gdal_translate(sds[52], dst_dataset = filename) #3km
#  gdal_translate(sds[72], dst_dataset = filename) #10km

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
  #  Lat_Lon_Values <- merge(Lat_Lon, values, all=TRUE)
  
  MODIS_data <- Lat_Lon_Values %>%
    dplyr:: select(lon, lat, values)
  MODIS_data <- na.omit(MODIS_data)
  
  write.csv(MODIS_data, file = paste(nome,".csv", sep = ""), row.names=FALSE)
  
}  


BBB <- lapply(filenames, extract_HDF)


######################################################################################
######################################################################################

# collate the tiles together ####-------------------------------------

filenames_tiles <- list.files(pattern = "\\.csv$")

LAT = NULL
LON = NULL
aod = NULL

# filenames_tiles <- filenames_tiles[1:5]

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

write.csv(MODIS04_data, paste0("AOD_MOD04_10km_UAE","_",folder_day,".csv"))

head(MODIS04_data)


################################################################################
# subset data for a selected region
# UAE

DATA_EMIRATES <- subset(MODIS04_data, lon <= 80 & lon >= 10 & lat >= 10 & lat <= 55)

head(DATA_EMIRATES)

write.csv(DATA_EMIRATES, paste0("AOD_MOD04_10km_UAE","_",folder_day,".csv"))
DATA_EMIRATES <- read.csv(paste0("AOD_MOD04_10km_UAE","_",folder_day,".csv"))[-1]
MODIS04_data <- read.csv(paste0("AOD_MOD04_10km_UAE","_",folder_day,".csv"))[2:4]

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


writeRaster(r, paste0("AOD_MOD04_10km_UAE","_",folder_day,".tif"), overwrite = TRUE)
Data_Emirates_10km_tif <- raster(paste0("AOD_MOD04_10km_UAE","_",folder_day,".tif"))

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

file.Emirates_tif_1km <- paste0("AOD_MOD04_1km_UAE","_",folder_day,".tif")
Emirates_tif_1km_tiff <- writeRaster(Data_Emirates_tif_1km, filename = file.Emirates_tif_1km, format = 'GTiff', overwrite = T)

### Extract points from raster tiff ############################################

Emirates_tif_1km_pts <- rasterToPoints(Emirates_tif_1km_tiff)
colnames(Emirates_tif_1km_pts) <- c("Lon", "Lat", "AOD_1km")
Emirates_tif_1km_pts <- as.data.frame(Emirates_tif_1km_pts) 
Emirates_tif_1km_pts <- subset(Emirates_tif_1km_pts, !is.na(AOD_1km) & AOD_1km>0)

write_csv(Emirates_tif_1km_pts , paste0("Emirates_tif_1km.csv","_",folder_day,".csv"))


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
            title = paste("<br><strong>AOD (MODIS-5km)",folder_day,"</strong>"),
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
saveWidget(map, paste0("AOD_Emirates","_",folder_day,".html"), selfcontained = FALSE)




# 
# 
# 
# 
# 
# ################################################################################
# # subset data for a selected region
# # United Arab Emirates
# 
# MODIS04_data <- read.csv(paste0("AOD_MOD04_5km_UAE","_",folder_day,".csv"))[2:4]
# DATA_EMIRATES <- subset(MODIS04_data, lon <= 80 & lon >= 10 & lat >= 10 & lat <= 55)
# 
# head(DATA_EMIRATES)
# 
# write.csv(DATA_EMIRATES, paste0("AOD_MOD04_5km_UAE","_",folder_day,".csv"))
# DATA_EMIRATES <- read.csv(paste0("AOD_MOD04_5km_UAE","_",folder_day,".csv"))[-1]
# 
# ############  Make a data interpolation on a regular grid ########################
# 
# # min(DATA_EMIRATES$lat)
# # max(DATA_EMIRATES$lat)
# # min(DATA_EMIRATES$lon)
# # max(DATA_EMIRATES$lon)
# 
# DATA_EMIRATES$x <- DATA_EMIRATES$lon
# DATA_EMIRATES$y <- DATA_EMIRATES$lat
# 
# coordinates(DATA_EMIRATES) = ~x + y  ## Set spatial coordinates to create a Spatial object:
# 
# x.range <- as.numeric(c(20, 85))  # min/max longitude of the interpolation area
# y.range <- as.numeric(c(-1.7, 60))  # min/max latitude of the interpolation area
# 
# ## grid at 5km resolution
# grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.05),
#                    y = seq(from = y.range[1], to = y.range[2], by = 0.05))  # expand points to grid
# coordinates(grd) <- ~x + y
# gridded(grd) <- TRUE
# 
# # plot(grd, cex = 1.5, col = "grey")
# # points(DATA_EMIRATES, pch = 1, col = "red", cex = 1)
# 
# idw <- idw(formula = values ~ 1, locations = DATA_EMIRATES, 
#            newdata = grd)  # apply idw model for the data (interpolation)
# 
# idw.output = as.data.frame(idw)  # output is defined as a data table
# names(idw.output)[1:3] <- c("Lon", "Lat", "values")  # give names to the modelled variables
# 
# # write.csv(idw.output, file = paste0("AOD_MOD04_5km_UAE","_",folder_day,"_interp.csv"), row.names=FALSE)
# 
# ##### make raster with interpolated data----------------------------------------
# 
# # DATA_EMIRATES_interp <- read.csv(paste0("AOD_MOD04_5km_UAE","_",folder_day,"_interp.csv"))
# 
# DATA_EMIRATES_interp <- subset(idw.output, Lon <= 57 & Lon >= 50 & Lat >= 21 & Lat <= 56)
# write.csv(DATA_EMIRATES_interp, file = paste0("AOD_MOD04_5km_UAE","_",folder_day,"_interp.csv"), row.names=FALSE)
# 
# coordinates(DATA_EMIRATES_interp) <- ~ Lon + Lat
# # coerce to SpatialPixelsDataFrame
# gridded(DATA_EMIRATES_interp) <- TRUE
# raster_DATA_EMIRATES_interp <- raster(DATA_EMIRATES_interp)
# projection(raster_DATA_EMIRATES_interp) <- CRS("+proj=longlat +datum=WGS84")
# # plot(raster_DATA_EMIRATES_interp)
# 
# raster_DATA_EMIRATES_interp <- crop(raster_DATA_EMIRATES_interp, extent(shp_UAE))
# raster_DATA_EMIRATES_interp <- mask(raster_DATA_EMIRATES_interp, shp_UAE)
# 
# writeRaster(raster_DATA_EMIRATES_interp, paste0("AOD_MOD04_5km_UAE","_",folder_day,".tif"), overwrite = TRUE)
# 
# Data_Emirates_tif <- raster(paste0("AOD_MOD04_5km_UAE","_",folder_day,".tif"))
# # plot(Data_Emirates_tif)
# 
# 
# # Data_Emirates_10km_nc <- writeRaster(Data_Emirates_tif,
# #                                         filename="Data_Emirates_10km.nc",
# #                                         format="CDF", overwrite=TRUE) 
# 
# ##########################################################################
# ##########################################################################
# 
# # load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)
# 
# PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
# projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")
# 
# ### crop raster over the UAE shp file  ###############################
# PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
# PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)
# 
# file.tiff_GWR_PM25_1km_2015 <- 'GWR_PM25_1km_2015.tif'
# GWR_PM25_1km_2015_tiff <- writeRaster(PM25_2015_GRW, filename = file.tiff_GWR_PM25_1km_2015, format = 'GTiff', overwrite = T)
# 
# # check resoltuion of Donkelaar data and of ECMWF data
# res(GWR_PM25_1km_2015_tiff) # 0.01
# res(Data_Emirates_tif) # 0.1
# 
# ###################
# 
# ## make resolution of MODIS-data as the one of GWR-------------------------------
# Data_Emirates_tif_1km = projectRaster(Data_Emirates_tif, PM25_2015_GRW)
# Data_Emirates_tif_1km <- crop(Data_Emirates_tif_1km, extent(shp_UAE))
# Data_Emirates_tif_1km <- mask(Data_Emirates_tif_1km, shp_UAE)
# 
# file.Emirates_tif_1km <- paste0("AOD_MOD04_1km_UAE","_",folder_day,".tif")
# Emirates_tif_1km_tiff <- writeRaster(Data_Emirates_tif_1km, filename = file.Emirates_tif_1km, format = 'GTiff', overwrite = T)
# 
# 
# #### make Leafleft map with current-AOD data 1km ###-------------------------------
# 
# library(leaflet)
# 
# # Emirates_tif_1km_tif <- raster("C:/MI Drive/MODIS_AOD/2016/355/Emirates_tif_1km.tif")
# # plot(Emirates_tif_1km_tif)
# 
# # MIN_data <- min(minValue(Emirates_tif_1km_tif))
# # MAX_data <- max(maxValue(Emirates_tif_1km_tif))   
# # 
# # 
# # ## create an unique legend and colorbar for Data Emirates
# # rast_pal_EMIRATES <- colorNumeric(c("#ffffff", "#b7b700", "#e50000"), 
# #                                   c(MIN_data,MAX_data),
# #                                   na.color = "transparent")
# 
# # define color palette
# rast_pal_EMIRATES <- colorNumeric(c("#9999FF", "#ffd699", "#FFFF00", "#ffbf00", "#ffc700", "#FF0000", "#994c00"),
#                                         getValues(Emirates_tif_1km_tiff),na.color = "transparent")
# 
# 
# 
# map <- leaflet() %>% 
#   # setView(45, 25, 5) %>%
#   addTiles(group = "OSM (default)") %>%
#   addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
#   addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
#   addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
#   addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
#   addRasterImage(Emirates_tif_1km_tiff, 
#                  colors = rast_pal_EMIRATES, 
#                  opacity = 0.6,
#                  group = "AOD_EMIRATES") %>%
#   addLegend("bottomright", pal = rast_pal_EMIRATES, values = values(Emirates_tif_1km_tiff),# values = c(MIN_data, MAX_data),
#             title = paste("<br><strong>AOD (MODIS-5km)",current_folder_day,"</strong>"),
#             labFormat = labelFormat(prefix = ""),
#             opacity = 0.6) %>%
#   addPolygons(stroke = TRUE, data = shp_UAE,
#               weight = 1.5, color = "black",
#               fillOpacity = 0,
#               group = "shape_UAE") %>%
#   addLayersControl(
#     baseGroups = c("Road map", "Topographical", "Satellite", "Toner Lite"),
#     overlayGroups = c("AOD_EMIRATES"),
#     options = layersControlOptions(collapsed = TRUE)) 
# # map
# 
# # save map
# saveWidget(map, paste0("AOD_Emirates","_",current_folder_day,".html"), selfcontained = FALSE)
# # webshot(paste0("AOD_Emirates","_",current_date,".html"), file=paste0("AOD_Emirates","_",current_date,".png"),
# #         vwidth = 1200, vheight = 800, 
# #         cliprect = 'viewport')
# 
# 
# 





