
library(stringr)

setwd("C:/MI Drive/MODIS_AOD/2016")
wd <- getwd()

# list folder names
DAYS <- str_sub(list.files(), start = 1, end = -1)
# DAYS <- DAYS[-1]
# DAYS <- as.numeric(DAYS)


dir.create("2016_AOD_tiff_5km")


# match folder names with specific filenames and extention .tif
# copy tif files into a new directory
for (i in 1:length(DAYS)) {
  # enter into the directory
  setwd(paste0(wd,"/",DAYS[i]))
  files_day <- list.files(pattern = paste0("_UAE_",DAYS[i],".tif"))
  for (j in 1: length(files_day))
  file.copy(files_day[j], paste0(wd,"/2016_AOD_tiff_5km"))
}


# change day into data on each tiff filename

setwd(paste0(wd,"/2016_AOD_tiff_5km"))
filelist <- list.files(pattern = "\\.tif$")


for (i in 1:length(filelist)) {
  # replace name of a file
  file.rename(list.files(pattern = paste0(DAYS[i],".tif")),
              paste0("2016-",format(strptime(DAYS[i], format="%j"), format="%m-%d"), ".tif"))
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


