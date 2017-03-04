
library(rgdal)
library(ggplot2)
library(gstat)
library(sp)
library(raster)   
library(rgeos)
library(plyr)
library(dplyr)
library(leaflet)
library(htmltools)
library(ncdf4)
library(RNetCDF)
library(fields)
library(readr)
library(threadr)
library(htmlwidgets)


# load shp files for UAE

dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/website_MODIS/UAE_boundary"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
# names(shp)

shp_UAE@data$name <- 1:nrow(shp_UAE)
# plot(shp_UAE)


## list tiff data------------------------------------------

# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/2013_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/2014_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/2015_AOD_tiff_1km")
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/2016_AOD_tiff_1km")

filenames <- list.files(pattern = "\\.tif$")


## make resolution of raster of MODIS data at 1km as the one of ECMWF (40km)--------------------------------
# load a sample ECMWF raster to get its resolution

ECMWF_sample <- raster("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/Jun_2015.nc")
plot(ECMWF_sample)

ECMWF_sample <- crop(ECMWF_sample, extent(shp_UAE))
ECMWF_sample <- mask(ECMWF_sample, shp_UAE)
plot(ECMWF_sample)

# create directory-- 2013
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/csv_40km")
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/csv_40km"


# create directory-- 2014
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/csv_40km")
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/csv_40km"


# create directory--2015
# dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/csv_40km")
# new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/csv_40km"


# create directory-- 2016
dir.create("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/csv_40km")
new_dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/csv_40km"


# function to regrid all tiff file to 40km
# file <- filenames[1]
regrid_40km <- function (file) {
  date <- str_sub(file, start = 1, end = -5)
  
  AOD_1km <- raster(file)
    # change resolution as the one of the ECMWF
  AOD_40km = projectRaster(AOD_1km, ECMWF_sample) 
  res(AOD_40km)
  plot(AOD_40km)


### Extract points from raster

AOD_40km_pts <- rasterToPoints(AOD_40km)
head(AOD_40km_pts)
colnames(AOD_40km_pts) <- c("Lon", "Lat", "AOD_40km")
AOD_40km_pts <- as.data.frame(AOD_40km_pts)
# add a column for the date
AOD_40km_pts$Date <- date
# AOD_40km_pts$Date <- paste0(date,"-", sprintf("%02d", i))


# save csv file for each day----------------------- 
write_csv(AOD_40km_pts, paste0(path = new_dir,"/", date,".csv"))

}

BBB <- lapply(filenames, regrid_40km)

