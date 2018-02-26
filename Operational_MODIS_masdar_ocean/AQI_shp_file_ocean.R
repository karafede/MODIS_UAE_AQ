
library(readr)
library(sp)
library(raster)
library(gstat)
library(rgdal)
library(RNetCDF)
library(ncdf4)
library(stringr)
library(rgeos)
library(leaflet)
library(htmlwidgets)
library(dplyr)
library(spatialEco)
library(RCurl)


# load function to calculte Air Quality Indexes
source("/home/mariners/MODIS_AOD/aqi_fun_UAE.R")

Sys.time()
current_date <- str_sub(Sys.time(), start = 1, end = -10)
str(current_date)
folder_year <- str_sub(current_date, start = 1, end = -7)

## 3km
# url_list_data = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/61/MOD04_3K/2016/"
## 10km MODIS Terra
url_list_data_Terra = paste0("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/61/MOD04_L2/",folder_year,"/") 
## 10km MODIS Aqua
url_list_data_Aqua = paste0("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/61/MYD04_L2/",folder_year,"/")


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
folder_year <- str_sub(folder, start = 75, end = -5)
folder_day <- str_sub(folder, start = 80, end = -1) 


#########################
### shapefile UAE #######
#########################

dir <- "/home/mariners/MODIS_AOD/UAE_boundary/"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
# names(shp)

shp_UAE@data$name <- 1:nrow(shp_UAE)
# plot(shp_rect)
# plot(shp_UAE, add=TRUE, lwd=1)


##################################################
# make a point for a location in the UAE #########
##################################################

require(sf)

# load coordinates of all UAE Airports
coordinates_UAE <- read.table(text="
                              longitude    latitude ID
                              54.65    24.433  AbuDhabi
                              55.333   25.25  Dubai 
                              55.609   24.262  AlAin
                              55.517   25.329  Sharjah
                              56.20   25.16  Fujairah
                              53.633  23.633  MedinaZayed
                              54.548  24.248  AlDhafra
                              54.458  24.428  Bateen
                              55.172  24.886  AlMaktoum
                              55.939  25.613  Rak
                              53.383  23.617 Buhasa
                              53.65   23.0355 LiwaOasis
                              55.71   25.52 UmmAlQuwain
                              52.72   24.104 Ruwais",
                              header=TRUE)



coordinates_UAE <- st_as_sf(x = coordinates_UAE, 
                        coords = c("longitude", "latitude"),
                        crs = "+proj=longlat +datum=WGS84")

# convert to sp object if needed
coordinates_UAE_point <- as(coordinates_UAE, "Spatial")
setwd("/home/mariners/MODIS_AOD/shapes/")
shapefile(coordinates_UAE_point, "points.shp", overwrite=TRUE)

dir <- "/home/mariners/MODIS_AOD/shapes/"
points <- readOGR(dsn = dir, layer = "points")
# plot(points)
# add a buffer around each point
shp_buff <- gBuffer(points, width=0.075, byid=TRUE)
shp_buff <- spTransform(shp_buff, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(shp_buff)
# plot(shp_UAE, add=TRUE, lwd=1)

# save shp file for the circular buffer
setwd("/home/mariners/MODIS_AOD/shapes/")
shapefile(shp_buff, "circular_buffer.shp", overwrite=TRUE)
dir <- "/home/mariners/MODIS_AOD/shapes/"
# reload and plot domain
shp_buff <- readOGR(dsn = dir, layer = "circular_buffer")


#######################################
# load processed MODIS point data #####
#######################################
# folder_year <- "2018"
# folder_day <- "056"

setwd(paste0("/home/mariners/MODIS_AOD/", folder_year, "/", folder_day))

values <- read.csv(paste0("Emirates_tif_1km.csv","_",folder_day,".csv"))
str(values)

# get points from the raster (lat, lon, points)

head(values)

crs <- projection(shp_buff) ### get projections from shp file
# make a spatial object 
values <- SpatialPointsDataFrame(values[,1:2], values, 
                                 proj4string=CRS(crs))

# get NUMBER of POINTS that fall into the circular buffer zone
# pts_in_buffer <- sp::over(values, shp_buff, fun = NULL)
pts_in_buffer_ID <- over(values, shp_buff[, "ID"])
pts_in_buffer_ID <- na.omit(pts_in_buffer_ID)


library(spatialEco)
# find points inside the buffer
# pts_in_buffer <- point.in.poly(values, shp_buff)
pts_in_buffer <- values[shp_buff,] 
head(pts_in_buffer@data)

data_points <- pts_in_buffer@data 
names(data_points)
# Join ID
data_points <- cbind(data_points, pts_in_buffer_ID)


# Aggregate by zone/polygon
# data_points <- pts_in_buffer@data 
names(data_points)

data_points <- data_points %>% 
  dplyr::group_by(ID) %>% 
  dplyr::summarise(AVERAGE_PM25 = mean(AOD_1km)) %>% 
  dplyr::ungroup()


##############################################################################
# calculate Air Quality Index for the only PM2.5 from satellite MODIS data ###
##############################################################################


PM25_data <- as.vector(data_points$AVERAGE_PM25)
PM25_data <- round(PM25_data, digits = 0)
PM25_data <- as.numeric(PM25_data)


# calculate Air Quality index for PM2.5
aqi_PM25 <- lapply(PM25_data, aqi_PM25_fun)

aqi_PM25 <- as.numeric(aqi_PM25)
aqi_PM25 <- round(aqi_PM25)
aqi_PM25 <- as.data.frame(aqi_PM25)

# join AQi to original dataset
data_points <- cbind(data_points, aqi_PM25)
str(data_points)


# Join aggregation to polygons 
shp_buff@data <- shp_buff@data %>% 
  left_join(data_points, "ID")

# Filter out polygons with no data
data_points <- subset(data_points, !is.na(AVERAGE_PM25))


# Export shp file
head(shp_buff@data)


#### Write GeoJSON for Leaflet application ############################

# ----- Write data to GeoJSON
dir <- paste0("/home/mariners/MODIS_AOD/", folder_year, "/", folder_day)
leafdat <- paste(dir, "/",  "AQI", "_", folder_day, ".geojson", sep="") 
# leafdat
writeOGR(shp_buff, leafdat, layer="", driver="GeoJSON")  ## erase existing .geojson file when re-runing code 

#### Plot Jeojson with leaflet
# AQI <- readOGR(paste0(dir,"/", "AQI", "_", folder_day,".geojson"), "OGRGeoJSON")
AQI <- readOGR(paste0(dir,"/", "AQI", "_", folder_day, ".geojson"))


#### colors for map 
# qpal_AQI <- colorQuantile("Reds", AQI$aqi_PM25, n = 7)  (quantiles)

breaks <- c(0, 50, 100, 200, 300, 500)
colori <- c("#00CD00","#ffff00", "#e59400","#ff0000", "#800080", "#800000")
pal_AQI <- colorBin(colori, bins=breaks)
pal_AQI(shp_buff$aqi_PM25)


###########################################################################################################
# add line for the AQI level and the Category
# geom_hline(yintercept =50, col="#00CD00", lty=1, size = 2) +
#   #  geom_text(aes(x = min  , y = 55, label = "Good"), size = 6) +
#   
#   geom_hline(yintercept = 100, col="#ffff00", lty=1, size=2) +
#   #  geom_text(aes(x = min , y = 105, label = "Moderate"), size = 6) +
#   
#   geom_hline(yintercept = 150, col="#e59400", lty=1, size=2)  +
#   #  geom_text(aes(x = min , y = 155, label = "Unhealthy for Sensitive Groups"), size = 6) +
#   
#   geom_hline(yintercept = 200, col="#ff0000", lty=1, size=2) +
#   #  geom_text(aes(x = 0.7 , y = 220, label = "Unhealthy"), size = 3) +
#   
#   geom_hline(yintercept = 300, col="#800080", lty=1, size=2) +
#   #  geom_text(aes(x = 0.7 , y = 220, label = "Very Unhealthy"), size = 3) +
#   
#   geom_hline(yintercept = 500, col="#800000", lty=1, size=2) +
#   #  geom_text(aes(x = 0.7 , y = 220, label = "Hazardous), size = 3) +
###########################################################################################################

popup_AQI <- paste0("<p>", 
                    AQI$ID, 
                    "<br><strong>AQI PM<sub>2.5</sub>: </strong>", 
                    AQI$aqi_PM25)


####### CREATE map #######################################################

map = leaflet(AQI) %>% 
  addTiles(group = "OSM (default)") %>%
  addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
  addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
  addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
  
  
  # PM2.5 satellite data 
  addPolygons(stroke = TRUE, smoothFactor = 0.2, 
              fillOpacity = 0.5, 
#              color = ~ qpal_AQI(aqi_PM25),
              color = ~ pal_AQI(aqi_PM25),
              weight = 2,
              popup = popup_AQI,
              group = "AQI") %>%


  # Layers control
  addLayersControl(
    baseGroups = c("Toner Lite", "Road map", "Topographical", "Satellite"),
    overlayGroups = c("AQI"),
    options = layersControlOptions(collapsed = FALSE))

# map


#### to export into html file, use the button in the Wiever window: 
#### "save as Web Page"

saveWidget(map,
           file = paste0(dir,"/", "AQI", "_", folder_day, ".html"),
           selfcontained = FALSE)

