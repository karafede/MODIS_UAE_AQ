
library(stringr)

setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed")
wd <- getwd()

# list folder names
DAYS <- str_sub(list.files(), start = 1, end = -1)
# DAYS <- DAYS[-1]
# DAYS <- as.numeric(DAYS)

dir.create("2016_AOD_csv_interp_5km")


# match folder names with specific filenames and extention _interp_csv
# copy tif files into a new directory
for (i in 1:length(DAYS)) {
  # enter into the directory
  setwd(paste0(wd,"/",DAYS[i]))
  files_day <- list.files(pattern = paste0("_UAE_",DAYS[i],"_interp.csv"))
  for (j in 1: length(files_day))
  file.copy(files_day[j], paste0(wd,"/2016_AOD_csv_interp_5km"))
}


# change day into data on each .csv filename

# setwd(paste0(wd,"/2016_AOD_tiff_5km"))
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_AOD_csv_interp_5km")
filelist <- list.files(pattern = "\\.csv$")


# this is only for the leap year 2016
# run this for the first 60 days and replace 1st march with 29th of February 2016

# for (i in 1:length(DAYS)) {
#   # replace name of a file
#   file.rename(list.files(pattern = paste0(DAYS[i],"_interp.csv")),
#               paste0("2016-",format(strptime(as.numeric(DAYS[i]), format="%j"), format="%m-%d"), ".csv"))
# }


# this is only for the leap year 2016
# run this for days from 61 and replace 1st march with 29th of February 2016

for (i in 61:length(DAYS)) {
  # replace name of a file
  file.rename(list.files(pattern = paste0(DAYS[i],"_interp.csv")),
              paste0("2016-",format(strptime(as.numeric(DAYS[i]), format="%j"), format="%m-%d"), ".csv"))
}



# for (i in 1:(length(DAYS))) {
#   # replace name of a file
#   file.rename(list.files(pattern = paste0(DAYS[i],".tif")),
#              paste0("2016-",format(strptime(DAYS[i], format="%j"), format="%m-%d"), ".tif"))
# }

  
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


