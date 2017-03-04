
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


# script to re-organise all hdf file downloaded from the website of LADADS NASA into folders by day number
# setwd("E:/MASDAR_FK/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_Terra_Aqua")
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2012_Terra_Aqua")


# folder_day <- as.character("002")

wd <- getwd()
Sys.time()
current_date <- str_sub(Sys.time(), start = 1, end = -10)
str(current_date)

# make a list of all .hdf files and create a folder for each day
filenames <- list.files(pattern = "\\.hdf$") 

# make a list of the day (from the .hdf files)
for (i in 1:length(filenames)) {
  year <- str_sub(filenames[i], start = 11, end = -31)
  day <- str_sub(filenames[i], start = 15, end = -28)
  dir.create(day)
}



# list folder names
DAYS <- str_sub(list.dirs(), start = 3, end = -1)
DAYS <- DAYS[-1]

# match folder names with specific filenames
# copy files into respective folder names
for (i in 1:length(DAYS)) {
  # AAA <- list.files(pattern = paste0(year,"001"))
  files_day <- list.files(pattern = paste0("A",year,DAYS[i]))
for (j in 1: length(files_day))
    file.copy(files_day[j], DAYS[i])
}

########################################################################
