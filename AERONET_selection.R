
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
library(ggplot2)


setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/AERONET_MODIS")


# 2013-------------------------------------------------------------------------
AERONET_2013_MASDAR <- read.csv("130101_131231_Masdar_Institute.csv")
AERONET_2013_MASDAR <- AERONET_2013_MASDAR %>%
  mutate(Date = mdy(Date.dd.mm.yy.))

# extract hours
AERONET_2013_MASDAR$hour <- str_sub(AERONET_2013_MASDAR$Time.hh.mm.ss., start = 1, end = -7)

# select AOT 500 nm at 12pm
AERONET_2013_MASDAR <- AERONET_2013_MASDAR %>%
  dplyr::select(Date,
                hour,
                AOT_500) %>%
  group_by(Date,
           hour) %>%
  summarise(AOT = mean(AOT_500, na.rm = TRUE)) %>%
  filter(hour == "12")


# 2014-------------------------------------------------------------------------
AERONET_2014_MASDAR <- read.csv("140101_141231_Masdar_Institute.csv")
AERONET_2014_MASDAR <- AERONET_2014_MASDAR %>%
  mutate(Date = mdy(Date.dd.mm.yy.))

# extract hours
AERONET_2014_MASDAR$hour <- str_sub(AERONET_2014_MASDAR$Time.hh.mm.ss., start = 1, end = -7)

# select AOT 500 nm at 12pm
AERONET_2014_MASDAR <- AERONET_2014_MASDAR %>%
  dplyr::select(Date,
                hour,
                AOT_500) %>%
  group_by(Date,
           hour) %>%
  summarise(AOT = mean(AOT_500, na.rm = TRUE)) %>%
  filter(hour == 12)



# 2015-------------------------------------------------------------------------
AERONET_2015_MASDAR <- read.csv("150101_151231_Masdar_Institute.csv")
AERONET_2015_MASDAR <- AERONET_2015_MASDAR %>%
  mutate(Date = mdy(Date.dd.mm.yy.))

# extract hours
AERONET_2015_MASDAR$hour <- str_sub(AERONET_2015_MASDAR$Time.hh.mm.ss., start = 1, end = -7)


# select AOT 500 nm at 12pm
AERONET_2015_MASDAR <- AERONET_2015_MASDAR %>%
  dplyr::select(Date,
                hour,
                AOT_500) %>%
  group_by(Date,
           hour) %>%
  summarise(AOT = mean(AOT_500, na.rm = TRUE)) %>%
  filter(hour == 12)


# 2016-------------------------------------------------------------------------
AERONET_2016_MASDAR <- read.csv("160101_161231_Masdar_Institute.csv")
AERONET_2016_MASDAR <- AERONET_2016_MASDAR %>%
  mutate(Date = mdy(Date.dd.mm.yy.))

# extract hours
AERONET_2016_MASDAR$hour <- str_sub(AERONET_2016_MASDAR$Time.hh.mm.ss., start = 1, end = -7)

# select AOT 500 nm at 12pm
AERONET_2016_MASDAR <- AERONET_2016_MASDAR %>%
  dplyr::select(Date,
                hour,
                AOT_500) %>%
  group_by(Date,
           hour) %>%
  summarise(AOT = mean(AOT_500, na.rm = TRUE)) %>%
  filter(hour == 11)


# bind all years togehter
AERONET_MASDAR <- rbind(AERONET_2013_MASDAR,
                        AERONET_2014_MASDAR,
                        AERONET_2015_MASDAR,
                        AERONET_2016_MASDAR)


####################################################################################

# load MODIS data exctracted at the location of the AERONET station in MASDAR


MODIS_2013_MASDAR <- read.csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2013_MODIS_processed/csv/extracted_MASDAR.csv")
MODIS_2014_MASDAR <- read.csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2014_MODIS_processed/csv/extracted_MASDAR.csv")
MODIS_2015_MASDAR <- read.csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2015_MODIS_processed/csv/extracted_MASDAR.csv")
MODIS_2016_MASDAR <- read.csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/MODIS_LAADS_NASA/2016_MODIS_processed/csv/extracted_MASDAR.csv")

# bind all years together
MODIS_MASDAR <- rbind(MODIS_2013_MASDAR,
                      MODIS_2014_MASDAR,
                      MODIS_2015_MASDAR,
                      MODIS_2016_MASDAR)

# make date for of MODIS data as the one of AERONET
MODIS_MASDAR$Date <- as.Date(parse_date_time(MODIS_MASDAR$Date,"mdy"))



# select only AERONET station site
MODIS_MASDAR <- MODIS_MASDAR %>%
  filter(Site == "MI_AERONET")


##################################################################################
###################################################################################

# merge AERONET data with MODIS data
MODIS_AERONET <- MODIS_MASDAR %>%
  left_join(AERONET_MASDAR, "Date")

# remove all lines with NA
MODIS_AERONET <- na.omit(MODIS_AERONET)

## get the months of observations
MODIS_AERONET$month <- factor(format(MODIS_AERONET$Date, format = "%b"), levels = month.abb)

## Define seasons
MODIS_AERONET$season <- character(length = nrow(MODIS_AERONET))
MODIS_AERONET$season[MODIS_AERONET$month %in% month.abb[c(1:2)]] <- "winter"
MODIS_AERONET$season[MODIS_AERONET$month %in% month.abb[c(12)]] <- "winter"
MODIS_AERONET$season[MODIS_AERONET$month %in% month.abb[c(3:5)]] <- "spring"
MODIS_AERONET$season[MODIS_AERONET$month %in% month.abb[c(6:8)]] <- "summer"
MODIS_AERONET$season[MODIS_AERONET$month %in% month.abb[c(9:11)]] <- "fall"
MODIS_AERONET$season <- factor(MODIS_AERONET$season, levels = c("winter","spring","summer","fall"))

MODIS_AERONET <- MODIS_AERONET %>%
  select(Date,
         Value,
         AOT,
         season)
colnames(MODIS_AERONET) <- c("Date","MODIS_AOD", "AERONET_AOD", "season")


## correlation between MODIS and AERONET~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(MODIS_AERONET$MODIS_AOD, MODIS_AERONET$AERONET_AOD)


# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(MODIS_AERONET, aes(season, MODIS_AOD)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE) +   # no legend
  ylim(0, 2) 
plot


# check your data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot <- ggplot(MODIS_AERONET, aes(season, AERONET_AOD)) +
  theme_bw() +
  geom_boxplot() + 
  guides(fill=FALSE) +   # no legend
  ylim(0, 2) 
plot


#### fit function and label for AOD ----------------------------------------------
#### this funtion FORCE regression to pass through the origin #######################

lm_eqn <- function(df){
  m <- lm(MODIS_AOD ~ -1 + AERONET_AOD, df);
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b = format(coef(m)[1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


# # this function includes the intercept~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lm_eqn <- function(df){
#   m <- lm(MODIS_AOD ~  AERONET_AOD, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
#                    list(a = format(coef(m)[2], digits = 2),
#                         b = format(coef(m)[1], digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
#####################################################################################

# plot with regression line-----


jpeg('Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/AERONET_MODIS/MODIS_vs_AERONET_MASDAR_2013_2016.jpg',    
     quality = 100, bg = "white", res = 200, width = 7, height = 7, units = "in")
par(mar=c(4, 10, 9, 2) + 0.3)
oldpar <- par(las=1)

library(plyr)

# define regression equation for each season
eq <- ddply(MODIS_AERONET, .(season),lm_eqn)


ggplot(MODIS_AERONET, aes(y=MODIS_AOD, x=AERONET_AOD)) +
  theme_bw() +
  geom_jitter(colour=alpha("black",0.15) ) +
  facet_grid(season ~ .) +
  theme( strip.text = element_text(size = 18)) + 
  scale_color_manual(values = c("#ff0000", "#0000ff", "#000000", "#ffb732")) + 
  geom_smooth(method = "lm", formula = y ~ -1 + x) +  # force fit through the origin
  ylab(expression("AOD (MODIS)")) +
  xlab(expression("AOD (AERONET)")) +
  ylim(c(0,1.8)) + 
  xlim(c(0,2)) +
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=13)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
  geom_text(data = eq, aes(x = 1.7, y = 0.25, label = V1),
            parse = TRUE, inherit.aes=FALSE, size = 4, color = "black" ) + 
  ggtitle("AOD-MODIS vs AERONET @ MASDAR") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20, hjust = 0.5)) 
 

par(oldpar)
dev.off()

#############################################################################
