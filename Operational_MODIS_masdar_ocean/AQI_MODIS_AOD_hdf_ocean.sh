#!/bin/bash

mainPath=/home/mariners/MODIS_AOD/  

# source /home/mariners/MODIS_AOD/export.sh 

# setup data and Julian day
year=`date +%Y`  ; jday=`date +%j`

cd $mainPath

Rscript ${mainPath}/AQI_shp_file_ocean.R

#################################################################################################################################
file_geojson=${mainPath}/${year}/${jday}/AQI_${jday}.geojson

###################################################################################################################################

cd ${mainPath}/${year}/${jday}

rsync -avz ${mainPath}/${year}/${jday}/AQI_${jday}.geojson pvernier@atlas-prod.minet.ae:/home/pvernier/scripts_cron/modis_data/
rsync -avz ${mainPath}/${year}/${jday}/AQI_${jday}.geojson fkaragulian@cesam-uat:/home/pvernier/scripts_cron/modis_data/
rsync -avz ${mainPath}/${year}/${jday}/AQI_${jday}.geojson fkaragulian@cesam-web-prod:/data/scripts_cron/modis_data/


echo "fine"

exit
