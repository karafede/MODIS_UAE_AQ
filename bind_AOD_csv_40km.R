
library(readr)

# bind all regridded satellite data from MODIS Terra & Aqua

setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/csv_40km")


filenames <- list.files(pattern = "\\.csv$")

### Function to bind all AOD ########################

LAT = NULL
LON = NULL
aod = NULL


# filenames <- filenames[3]

  for (i in 1:length(filenames)) {
    
    lon <- read_csv(filenames[i])[,1]
    lat <- read_csv(filenames[i])[,2]
    AOD <- read_csv(filenames[i])[,3]
    
    LON = rbind(LON, data.frame(lon))
    LAT = rbind(LAT, data.frame(lat))
    AOD = rbind(aod, data.frame(AOD))
}
  

AOD_data <- cbind(LON, LAT, AOD)

write.csv(AOD_data, "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/AOD_2015_MODIS_all_40km.csv")
