#!/bin/bash

mainPath=/home/mariners/MODIS_AOD/  

# source /home/mariners/MODIS_AOD/export.sh 

# setup data and Julian day
year=`date +%Y`  ; jday=`date +%j`

cd $mainPath


# Rscript ${mainPath}/Overpass_MODIS_Terra.R
# Rscript ${mainPath}/Overpass_MODIS_Aqua.R
# sleep 1m
# Rscript ${mainPath}/MODIS_download_10K_ocean.R
Rscript ${mainPath}/download_MODIS_https.R

#################################################################################################################################
file_tiff=${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif

if [ -f ${file_tiff} ] ; then 
   echo "file exist : PM2.5 map created"
else
   echo "PM2.5 Creation failed"
   /usr/bin/python ${mainPath}/sendMail.py "PM2.5 Creation failed: Date: ${year}  day: ${jday}" 
fi
###################################################################################################################################

cd ${mainPath}/${year}/${jday}

rsync -avz ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif pvernier@atlas-prod.minet.ae:/home/pvernier/scripts_cron/modis_data/
rsync -avz ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif fkaragulian@cesam-uat:/home/pvernier/scripts_cron/modis_data/
rsync -avz ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif fkaragulian@cesam-web-prod:/data/scripts_cron/modis_data/

echo "fine"

exit
