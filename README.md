# MODIS_UAE_AQ
## An ensemble of scripts to process satellite data from MODIS TERRA and AQUA (also MAIAC data) with the R software. 
## These algorithms generate operational maps from near-real time satellite data through the NASA website. The scripts are designed for Windows and LINUX environments. You can customize these scripts by changing the extension of your geographical area that are you going to choose starting from the shape files of your region.
## Bash codes to setup the operational scripts in LINUX are also reported

## The repositoriy contains code to process historical data of Aerosols Optical Depths from MODIS. Data were previously downloaded from https://ladsweb.modaps.eosdis.nasa.gov/search/

See more etails here below:
1)	<b> MODIS_reorganization.R <b>, order .hdf file by Julian day and create a folder for each day
2)	<b> MODIS_kriege_10K_historical.R <b>  and/or <b> MODIS_download_10K_historical.R <b> extracts points from each .hdf files (from raster to points), it combine tiles and create a raster file only for UAE.
3)	gather_tiff.R and geotiff_to_csv.R reoder all .tif file and extract points into a unique .csv file
4)	In order the extract satellite data (from the above generated txt file) at the locations of monitoring stations, we used a Matlab code: Extract_SAT_Data_auto_FK_DG.m that you can find in the folder Z:\_SHARED_FOLDERS\Air Quality\Phase 1\Pathflow of Phase I_DG\Sat_AOD_Correlation
5)	Correlation between satellite data points and ground measurements is performed through the script corr_AOD_MODIS_vs_PM25.R. The correlation factor can be used to “approximately” rescale all the AOD maps obtained from satellite data


Please hit me up if you need further clarifications.
