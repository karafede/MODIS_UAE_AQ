
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


setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_AOD_csv_interp_5km")
wd <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA"


dir <- paste0(wd,"/","UAE_boundary")
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
# names(shp)

shp_UAE@data$name <- 1:nrow(shp_UAE)
plot(shp_UAE)

filenames <- list.files(pattern = "\\.csv$")



for (i in 1:length(filenames)) {
file <- read_csv(filenames[i])
date <- str_sub(filenames[i], start = 1, end = -5)
coordinates(file) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(file) <- TRUE
raster_DATA_EMIRATES_interp <- raster(file)
projection(raster_DATA_EMIRATES_interp) <- CRS("+proj=longlat +datum=WGS84")

raster_DATA_EMIRATES_interp <- crop(raster_DATA_EMIRATES_interp, extent(shp_UAE))
raster_DATA_EMIRATES_interp <- mask(raster_DATA_EMIRATES_interp, shp_UAE)

writeRaster(raster_DATA_EMIRATES_interp, paste0("AOD_MOD04_5km_UAE","_",date,".tif"), overwrite = TRUE)

Data_Emirates_tif <- raster(paste0("AOD_MOD04_5km_UAE","_",date,".tif"))



# load Donkelaar Satellite-Derived PM2.5, 2015, at 35% RH [ug/m3] (with Geographical regression adjustment)

PM25_2015_GRW <- raster("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/GlobalGWR_PM25_GL_201501_201512-RH35.nc")
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

file.Emirates_tif_1km <- paste0(date,".tif")
Emirates_tif_1km_tiff <- writeRaster(Data_Emirates_tif_1km, filename = file.Emirates_tif_1km, format = 'GTiff', overwrite = T)

###########################################################################
#### make Leafleft map with current-AOD data 1km ###-----------------------

# # define color palette
# rast_pal_EMIRATES <- colorNumeric(c("#9999FF", "#ffd699", "#FFFF00", "#ffbf00", "#ffc700", "#FF0000", "#994c00"),
#                                   getValues(Emirates_tif_1km_tiff),na.color = "transparent")
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
#             title = paste("<br><strong>AOD (MODIS-5km)",date,"</strong>"),
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
# saveWidget(map, paste0("AOD_Emirates","_",date,".html"), selfcontained = FALSE)


}  
  
  
  