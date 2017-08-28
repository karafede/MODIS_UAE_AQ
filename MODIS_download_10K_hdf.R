
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

setwd("/disk3/fkaragulian/MODIS_AOD/")
wd <- getwd()
Sys.time()
current_date <- str_sub(Sys.time(), start = 1, end = -10)
str(current_date)

## 3km
# url_list_data = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/"
## 10km
url_list_data = "ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/" 

list_data = getURL(url_list_data, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
list_data = paste(url_list_data, strsplit(list_data, "\r*\n")[[1]], sep = "") 
# last element (date) of the list in order of time 
list_data <- sort(list_data)
n <- as.numeric(length(list_data))
list_data[n]
folder <- unlist(list_data[n])
# isolate folder name
folder_year <- str_sub(folder, start = 74, end = -5)
folder_day <- str_sub(folder, start = 79, end = -1) 

# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/",folder_day, "/", sep = "")
# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/","349", "/", sep = "")

url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/",folder_day, "/", sep = "")
# url = paste("ftp://karafede:Password08@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/","350", "/", sep = "")


filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 

# select only files in the time range as in the Overpass time

Overpass_times <- read.csv("/disk3/fkaragulian/MODIS_AOD/Overpass_times.csv")

matches <- unique (grep(paste(Overpass_times$x,collapse="|"), 
                        filenames, value=TRUE))


# create a new working directory for each download date
dir.create(folder_year)
# setwd(folder)

# download data
setwd(paste0(wd,"/",folder_year))
dir.create(folder_day)
setwd(paste0(wd,"/",folder_year,"/", folder_day))

filenames_MODIS_10k <- unlist(str_extract_all(matches, ".+(.hdf$)"))
filenames_MODIS_10k <- sort(filenames_MODIS_10k)

# start downloading data in the main directory -----------------------
mapply(download.file, filenames_MODIS_10k,basename(filenames_MODIS_10k)) 

