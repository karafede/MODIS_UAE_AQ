
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
###############################################################################


# collate the tiles together

################
### TERRA ######
################

filenames <- list.files(pattern = c("MOD04" ,"\\.hdf$"))

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

# delete HDF files
if (file.exists(filenames)) file.remove(filenames)  

######################################################################################
######################################################################################

# collate the tiles together ####-------------------------------------

filenames_tiles <- list.files(pattern = c("MOD04", "\\.csv$"))

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

write.csv(MODIS04_data, paste0("AOD_TERRA_10km_UAE","_",folder_day,".csv"))

head(MODIS04_data)

###################################################################################
###################################################################################
###################################################################################


################
### AQUA ######
################


setwd("E:/MASDAR_FK/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2017_Terra_Aqua/027")
filenames <- list.files(pattern = c("MYD04" ,"\\.hdf$"))

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
  
  r <- raster("AOD.tif")
  values <- rasterToPoints(r)
  colnames(values) <- c("x", "y", "values")
  values <- as.data.frame (values) 
  values$values <- (values$values)/1000
  
  system(paste0('gdal_translate_FWT ', get_subdatasets(file)[71],' lon.tif'))
  
  lon <- raster(paste0("lon.tif"))
  # data values for longitude
  longitude <- rasterToPoints(lon)  
  colnames(longitude) <- c("x", "y", "lon")
  longitude <- as.data.frame (longitude)
  
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

# delete HDF files
if (file.exists(filenames)) file.remove(filenames)  

######################################################################################
######################################################################################

# collate the tiles together ####-------------------------------------

filenames_tiles <- list.files(pattern = c("MYD04", "\\.csv$"))

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

write.csv(MODIS04_data, paste0("AOD_AQUA_10km_UAE","_",folder_day,".csv"))

head(MODIS04_data)

###################################################################################
###################################################################################
###################################################################################

dir <- "D:/MODIS_AOD/UAE_boundary"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

setwd("D:/361")
folder_day <- "361"

dir <- "D:/Dust_Event_UAE_2015/WRFChem_domain"

# larger WRF domain
shp_WRF <- readOGR(dsn = dir, layer = "domain_d02_4km_WRFChem_small")
plot(shp_WRF)


#####################################################################################################################
#####################################################################################################################
############## RASTERIZE MODIS data #################################################################################
#####################################################################################################################


############
# TERRA ####
############

  
  filenames_TERRA <- list.files(pattern = "AOD_TERRA_10km_")
  
  ##### make a function to create a raster from an irregular data frame (lat , lon , AOD)

    date <- str_sub(filenames_TERRA, start = 20, end = -5)
    AOD_TERRA <- read_csv(filenames_TERRA)[-1]  # load the one without NA values
    # subset data for a selected region
    # UAE
    
    AOD_TERRA <- subset(AOD_TERRA, lon <= 57 & lon >= 47 & lat >= 21 & lat <= 27)
    colnames(AOD_TERRA) <- c('x', 'y', 'z')
    
    
    xmn = min(AOD_TERRA$x)
    xmx = max(AOD_TERRA$x)
    ymn = min(AOD_TERRA$y)
    ymx = max(AOD_TERRA$y)
    
    
    x.range <- as.numeric(c(xmn,xmx ))
    y.range <- as.numeric(c(ymn,ymx))
    
    
    grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.20),
                       y = seq(from = y.range[1], to = y.range[2], by = 0.20))  # expand points to grid 10km (0.1)
    plot(grd)
    
    grd_1 <- dplyr::filter(grd, grd$x == xmn)
    nrow(grd_1)
    grd_2<- dplyr::filter(grd, grd$y == ymn)
    nrow(grd_2)
    
    r <- raster(xmn=min(AOD_TERRA$x), xmx=max(AOD_TERRA$x), ymn=min(AOD_TERRA$y),
                ymx=max(AOD_TERRA$y), ncol=nrow(grd_2), nrow= nrow(grd_1))
    
    
     raster_TERRA <- rasterize(AOD_TERRA[, 1:2], r, AOD_TERRA[,3], fun=mean)
     
    res(raster_TERRA)
    plot(raster_TERRA)
    projection(raster_TERRA) <- CRS("+proj=longlat +datum=WGS84")
    plot(shp_WRF, add=TRUE, lwd=1)
    
    # crop over the UAE area
    raster_TERRA <- crop(raster_TERRA, extent(shp_WRF))
    raster_TERRA <- mask(raster_TERRA, shp_WRF)
    plot(raster_TERRA)
    
    
    writeRaster(raster_TERRA, paste0("AOD_TERRA_10km_UAE_", date,".tif"), overwrite = TRUE)   

#######################################################################################
#######################################################################################

    ############
    # AQUA #####
    ############
    
    
    filenames_AQUA <- list.files(pattern = "AOD_AQUA_10km_")
    
    ##### make a function to create a raster from an irregular data frame (lat , lon , AOD)
    
    date <- str_sub(filenames_AQUA, start = 20, end = -5)
    AOD_AQUA <- read_csv(filenames_AQUA)[-1]  # load the one without NA values
    # subset data for a selected region
    # UAE
    
    AOD_AQUA <- subset(AOD_AQUA, lon <= 57 & lon >= 47 & lat >= 21 & lat <= 27)
    colnames(AOD_AQUA) <- c('x', 'y', 'z')
    
    
    xmn = min(AOD_AQUA$x)
    xmx = max(AOD_AQUA$x)
    ymn = min(AOD_AQUA$y)
    ymx = max(AOD_AQUA$y)
    
    
    x.range <- as.numeric(c(xmn,xmx ))
    y.range <- as.numeric(c(ymn,ymx))
    
    
    grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.20),
                       y = seq(from = y.range[1], to = y.range[2], by = 0.20))  # expand points to grid 10km (0.1)
    plot(grd)
    
    grd_1 <- dplyr::filter(grd, grd$x == xmn)
    nrow(grd_1)
    grd_2<- dplyr::filter(grd, grd$y == ymn)
    nrow(grd_2)
    
    r <- raster(xmn=min(AOD_AQUA$x), xmx=max(AOD_AQUA$x), ymn=min(AOD_AQUA$y),
                ymx=max(AOD_AQUA$y), ncol=nrow(grd_2), nrow= nrow(grd_1))
    
    
    raster_AQUA <- rasterize(AOD_AQUA[, 1:2], r, AOD_AQUA[,3], fun=mean)
    
    res(raster_AQUA)
    plot(raster_AQUA)
    projection(raster_AQUA) <- CRS("+proj=longlat +datum=WGS84")
    plot(shp_WRF, add=TRUE, lwd=1)
    
    # crop over the UAE area
    raster_AQUA <- crop(raster_AQUA, extent(shp_WRF))
    raster_AQUA <- mask(raster_AQUA, shp_WRF)
    plot(raster_AQUA)
    
    
    writeRaster(raster_AQUA, paste0("AOD_AQUA_10km_UAE_", date,".tif"), overwrite = TRUE)   
    
    
    
# ############## AVERAGE MODIS Terra and MODIS Aqua ################################################
# ##############--------------------------------------- ############################################
 
    
reference <- raster_AQUA
# reproject each raster with the same extent and resolution of the reference raster above
raster_TERRA = projectRaster(raster_TERRA, reference)

stack_MODIS <- stack(raster_TERRA, raster_AQUA)
# calculate the mean of all the rasters (2) in the stack
mean_MODIS <- mean(stack_MODIS, na.rm=TRUE)
plot(mean_MODIS)



#########-------------------------------------------------------------------------
# load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)

# PM25_2015_GRW <- raster("/disk3/fkaragulian/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
PM25_2015_GRW <- raster("D:/MODIS_AOD/GlobalGWR_PM25_GL_201501_201512-RH35.nc")


projection(PM25_2015_GRW) <- CRS("+proj=longlat +datum=WGS84")

### crop raster over the UAE shp file  ###############################
PM25_2015_GRW <- crop(PM25_2015_GRW, extent(shp_UAE))
PM25_2015_GRW <- mask(PM25_2015_GRW, shp_UAE)

file.tiff_GWR_PM25_1km_2015 <- 'GWR_PM25_1km_2015.tif'
GWR_PM25_1km_2015_tiff <- writeRaster(PM25_2015_GRW, filename = file.tiff_GWR_PM25_1km_2015, format = 'GTiff', overwrite = T)


###################

## make resolution of MODIS-data as the one of GWR-------------------------------
Data_Emirates_tif_1km = projectRaster(mean_MODIS, PM25_2015_GRW)

# Data_Emirates_tif_1km = projectRaster(raster_AQUA, PM25_2015_GRW)

Data_Emirates_tif_1km <- crop(Data_Emirates_tif_1km, extent(shp_UAE))
Data_Emirates_tif_1km <- mask(Data_Emirates_tif_1km, shp_UAE)

plot(Data_Emirates_tif_1km)

file.Emirates_tif_1km <- paste0("AOD_MODIS_1km_UAE","_",folder_day,".tif")
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
  #         title = "<br><strong>PM<sub>10</sub> (<font face=symbol>m</font>g/m<sup>3</sup>) MAIAC: </strong>",
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


