#!/bin/bash


mainPath=/home/fkaragulian/MODIS_AOD/
source /home/fkaragulian/export.sh 
cd $mainPath


# /apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/MODIS_download_10K_historical.R
/apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/MODIS_kriege_10K_historical.R


echo "fine"

exit
