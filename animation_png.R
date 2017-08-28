
library(webshot)
library(leaflet)
library(htmlwidgets)
# library(ggmap)
library(RColorBrewer)
library(raster)
library(classInt)
library(stringr)

setwd("D:/website_MODIS/AOD_web/2016_AOD_tiff_1km")

# Load data

list_tiff <- list.files(pattern = "\\.tif$")
list_days <- str_sub(list_tiff, start = 1, end = -5)

df_AOD = data.frame(dates = list_days, image_tiff = list_tiff)
df_AOD$dates <- as.character(df_AOD$dates)
df_AOD$image_tiff <- as.character(df_AOD$image_tiff)


 for (i in 1:nrow(df_AOD)) {
  
#  for (i in 1:5) {
  # AOD_raster_tiff <- raster(as.character(filter(df_AOD, dates == input$dates))[2])
AOD_raster_tiff <- raster(as.character(df_AOD$image_tiff)[i]) 

values(AOD_raster_tiff)[values(AOD_raster_tiff) < 0.05]<- 0.05
values(AOD_raster_tiff)[values(AOD_raster_tiff) > 1.6]<- 1.6


MIN_AOD <- 0.05
MAX_AOD <- 1.6

# pal_AOD <- colorNumeric(c("#9999FF", "#ffd699", "#FFFF00", "#ffbf00", "#ffc700", "#FF0000", "#994c00"),
#                         getValues(AOD_raster_tiff),na.color = "transparent")

pal_AOD <- colorNumeric(c("#9999FF", "#ffd699", "#FFFF00", "#ffbf00", "#ffc700", "#FF0000", "#994c00"),
                        c(MIN_AOD, MAX_AOD), na.color = "transparent")

# define popup for AOD 
 "h1 { font-size: 4px;}"
content <- paste('<h1><strong>', (df_AOD$dates)[i],'', sep = "")



  map <- leaflet() %>% 
    addTiles() %>% 
    addTiles(group = "OSM (default)") %>%
    addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
    addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
    
    addPopups(53.37, 24.8, content,
              options = popupOptions(closeButton = FALSE)) %>%

    addRasterImage(AOD_raster_tiff, 
                   colors = pal_AOD, 
                   opacity = 0.5, group = "AOD") %>%
    # addLegend("bottomright", pal = pal_AOD, values = getValues(AOD_raster_tiff),
    #           title = "<br><strong>AOD - MODIS: </strong>",
    #           labFormat = labelFormat(prefix = ""),
    #           opacity = 0.5) %>%
    addLegend("bottomright", pal = pal_AOD, values = c(MIN_AOD, MAX_AOD),
              title = "<br><strong>AOD - MODIS: </strong>",
              labFormat = labelFormat(prefix = ""),
              opacity = 0.5) %>%
    addLayersControl(
      baseGroups = c("Road map", "Topographical", "Satellite", "Toner Lite"),
      overlayGroups = "AOD",
      options = layersControlOptions(collapsed = TRUE))
  
  
  ## This is the png creation part
  saveWidget(map, 'temp.html', selfcontained = FALSE)
  webshot('temp.html', file = paste0((df_AOD$dates)[i],".png"), vwidth = 900, vheight = 900,
          cliprect = 'viewport')
  
}


# to use with ImageMagik using the commnad line cmd in windows

# magick -delay 100 -loop 0 *.png MODIS_AOD_2016.gif


###############################################################################
###############################################################################
###############################################################################
