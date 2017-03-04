#!/bin/bash

#export PATH=/apps/gcc/gcc-4.9.2/bin/:/apps/R/R-3.3.2/bin:/apps/netcdf/netcdf-4.3.2/bin:/apps/gdal/gdal-2.1.2/bin:/apps/curl/curl-7.37.1/bin:$PATH
#export LD_LIBRARY_PATH=/apps/gcc/gcc-4.9.2/lib64:/apps/R/R-3.3.2/lib64:/apps/netcdf/netcdf-4.3.2/lib:/apps/gdal/gdal-2.1.2/lib:/apps/curl/curl-7.37.1/lib:$LD_LIBRARY_PATH

mainPath=/home/fkaragulian/MODIS_AOD/
source /home/fkaragulian/export.sh 
cd $mainPath

#module load R/3.3.2
#module load netcdf/4.3.2
#module load gdal/2.1.2

/apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/Overpass_MODIS_Terra.R
/apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/Overpass_MODIS_Aqua.R
sleep 1m
/apps/R/R-3.3.2/bin/Rscript /home/fkaragulian/MODIS_AOD/MODIS_download_10K.R
#/apps/R/R-3.3.2/bin/Rscript  /home/fkaragulian/ECMWF_forecasts/CAMS_NRT_download.R

#Rscript /home/fkaragulian/MODIS_AOD/Overpass_MODIS_Terra.R
#Rscript /home/fkaragulian/MODIS_AOD/MODIS_download_10K_hdf.R

echo "fine"

exit
