
library(stringr)

setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed")
wd <- getwd()

# list folder names
DAYS <- str_sub(list.files(), start = 1, end = -1)
# DAYS <- DAYS[-1]
# DAYS <- as.numeric(DAYS)

dir.create("2013_AOD_tiff_1km")


# match folder names with specific filenames and extention .tif
# copy tif files into a new directory
for (i in 1:length(DAYS)) {
  # enter into the directory
  setwd(paste0(wd,"/",DAYS[i]))
  files_day <- list.files(pattern = paste0("_1km_UAE_",DAYS[i],".tif"))
  for (j in 1: length(files_day))
  file.copy(files_day[j], paste0(wd,"/2013_AOD_tiff_1km"))
}

# delete the 10km file
# change day into data on each tiff filename

setwd(paste0(wd,"/2013_AOD_tiff_1km"))
filelist <- list.files(pattern = "\\.tif$")


# for (i in 61:length(DAYS)) {
#   # replace name of a file
#   file.rename(list.files(pattern = paste0(DAYS[i],".tif")),
#               paste0("2015-",format(strptime(as.numeric(DAYS[i])-1, format="%j"), format="%m-%d"), ".tif"))
# }




for (i in 1:(length(DAYS))) {
  # replace name of a file
  file.rename(list.files(pattern = paste0(DAYS[i],".tif")),
             paste0("2013-",format(strptime(DAYS[i], format="%j"), format="%m-%d"), ".tif"))
}

  
# 
# # reset correct working directory
# setwd("C:/MI Drive/MODIS_AOD/2016")
# dir.create("2016_AOD_csv_interp_5km")
# 
# # match folder names with specific filenames and extention .tif
# # copy csv files into a new directory
# for (i in 1:length(DAYS)) {
#   # enter into the directory
#   setwd(paste0(wd,"/",DAYS[i]))
#   files_day <- list.files(pattern = paste0("_UAE_",DAYS[i],"_interp.csv"))
#   for (j in 1: length(files_day))
#     file.copy(files_day[j], paste0(wd,"/2016_AOD_csv_interp_5km"))
# }
# 


