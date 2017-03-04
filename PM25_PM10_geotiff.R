
library(raster) 
library(stringr)

# setwd("D:/website_MODIS/AOD_web/2015_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/2015_AOD_tiff_1km")
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/2013_AOD_tiff_1km")


filenames <- list.files(pattern = "\\.tif$")

### Function to extract points from raster AOD ########################

PM25_PM10_tiff <- function (file) {      ## this is the filenames 
  
 # file <- filenames[6]
  
  name <- str_sub(file, start = 1, end = -5)
  

AOD_pts <- rasterToPoints(raster(file))

colnames(AOD_pts) <- c("Lon", "Lat", "AOD_1km")
AOD_pts <- as.data.frame (AOD_pts)
AOD_pts <- subset(AOD_pts, !is.na(AOD_1km) & AOD_1km>0)

AOD_pts$AOD_PM25 <- (AOD_pts$AOD_1km)*93
AOD_pts$AOD_PM10 <- (AOD_pts$AOD_1km)*279

AOD_PM25 <- AOD_pts %>%
  dplyr::select(Lon,
         Lat,
         AOD_PM25)


coordinates(AOD_PM25) <- ~ Lon + Lat
gridded(AOD_PM25) <- TRUE
raster_AOD_PM25 <- raster(AOD_PM25)
projection(raster_AOD_PM25) <- CRS("+proj=longlat +datum=WGS84")
# plot(raster_AOD_PM25)
writeRaster(raster_AOD_PM25, paste0("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/2013_AOD_PM25_tiff_1km/",name,".tif"), overwrite = TRUE)


AOD_PM10 <- AOD_pts %>%
  dplyr::select(Lon,
                Lat,
                AOD_PM10)

coordinates(AOD_PM10) <- ~ Lon + Lat
gridded(AOD_PM10) <- TRUE
raster_AOD_PM10 <- raster(AOD_PM10)
projection(raster_AOD_PM10) <- CRS("+proj=longlat +datum=WGS84")
# plot(AOD_PM10)
writeRaster(raster_AOD_PM10, paste0("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/2013_AOD_PM10_tiff_1km/",name,".tif"), overwrite = TRUE)

}

  
# extract_points_to_csv(filenames[1])
# apply the function to exctract point from raster from all the list of geotif files
BBB <- lapply(filenames, PM25_PM10_tiff)






