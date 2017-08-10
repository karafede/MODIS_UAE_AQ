#!/bin/bash


mainPath=/home/fkaragulian/MODIS_AOD/  ; dest=/research/cesam/AirQuality/MODIS_AOD

source /home/fkaragulian/export.sh 

year=`date +%Y`  ; jday=`date +%j`

mkdir -p ${dest}/${year}/${jday}

cd $mainPath


/apps/R/R-3.3.2/bin/Rscript ${mainPath}/Overpass_MODIS_Terra.R
/apps/R/R-3.3.2/bin/Rscript ${mainPath}/Overpass_MODIS_Aqua.R
sleep 1m
 /apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/MODIS_download_10K.R

# /apps/R/R-3.3.2/bin/Rscript ${mainPath}/MODIS_download_10K_NEW_29June2015_HPC.R


cd ${mainPath}/${year}/${jday}


#cp -r ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif   ${dest}/${year}/${jday}   

scp ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif fkaragulian@10.102.14.39:/home/fkaragulian/mtg/_SHARED_FOLDERS/AQ_Website/MODIS_satellite/PM25_MODIS/

rsync -avz ${mainPath}/${year}/${jday}/PM25_MODIS_1km_UAE_${jday}.tif pvernier@atlas-prod.minet.ae:/home/pvernier/scripts_cron/modis_data/


echo "fine"

exit
