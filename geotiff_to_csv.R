
library(raster) 
library(stringr)
library(threadr)
library(lubridate)
library(readr)

# setwd("D:/website_MODIS/AOD_web/2015_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/2013_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/2014_AOD_tiff_1km")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/2015_AOD_tiff_1km")
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/2016_AOD_tiff_1km")


filenames <- list.files(pattern = "\\.tif$")
# file <- filenames[1]
file <- filenames

dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_Terra_Aqua"
DAYS <- str_sub(list.files(dir), start = 1, end = -1)


hour_avg <- 0
minutes_avg <- 0

for (i in 1:length(DAYS)) {
# i = 1

### Function to extract points from raster AOD ########################
# extract_points_to_csv <- function (file) {      ## this is the filenames 
  
# extract time from MODIS files~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date <- str_sub(file[i], start = 1, end = -5)

AOD_pts <- rasterToPoints(raster(file[i]))

colnames(AOD_pts) <- c("Lon", "Lat", "AOD_1km")
AOD_pts <- as.data.frame (AOD_pts)
AOD_pts <- subset(AOD_pts, !is.na(AOD_1km) & AOD_1km>0)


dir_DAY <- paste0(dir,"/",DAYS[i])  
HOURS <- str_sub(list.files(dir_DAY), start = 19, end = -25)
MINUTES <- str_sub(list.files(dir_DAY), start = 21, end = -23)
# make average of hours
HOUR_AVG <- mean(as.numeric(HOURS))
MINUTES_AVG <- mean(as.numeric(MINUTES))

AVG_HOUR <- mean(HOUR_AVG, hour_avg) + 4
AVG_MIN <- mean(MINUTES_AVG, minutes_avg)
AVG_time <- cbind(AVG_HOUR, AVG_MIN)


HOUR_AVG <- round(HOUR_AVG, digit = 0) +4  # (from GMT to GST)
HOUR_AVG <- paste0(HOUR_AVG, ":00")
AOD_pts$Date <- paste(date,HOUR_AVG)

AOD_pts <- AOD_pts %>%
  mutate(Date = ymd_hm(Date)) 

# create a csv folder inside the 201i_AOD_tiff_1km directory
write_csv(AOD_pts, file = paste("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/csv/",date,".csv", sep = ""), row.names=FALSE)

write.csv(AVG_time, file = paste("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/csv/AVG_time.csv"))
# print <- i

print <- AVG_time

}

# extract_points_to_csv(filenames[1])
# apply the function to exctract point from raster from all the list of geotif files

# BBB <- lapply(filenames, extract_points_to_csv)


  







