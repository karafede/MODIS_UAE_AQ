
library(rsconnect)

rsconnect::setAccountInfo(name='testaq',
                          token='6B39489E91B60D1F34FC7F23EC15599F',
                          secret='+gBT3YxY3rboG+9Le3JXP38HIMlX5YePE990LAAl')

# rsconnect::deployApp("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/R_FK")
rsconnect::deployApp("D:/website_MODIS/AOD_web")
