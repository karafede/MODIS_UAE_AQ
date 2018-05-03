

library(stringr)
library(lubridate)
library(plyr)

# read table about time when MODIS Terra overpass Abu Dhabi
# https://oceandata.sci.gsfc.nasa.gov/cgi/overpass_pred
# lat 24.4539
# long 54.3773

setwd("/home/mariners/MODIS_AOD/")
# setwd("D:/MODIS_AOD/")

Sys.time()
current_date <- str_sub(Sys.time(), start = 1, end = -10)
str(current_date)

Overpass <- read.table("satpasses_Aqua.txt", skip = 10, sep="\t")
# current_overpass_time <- str_sub(Overpass[1,1], start = 13, end = -66)

# make table with date and overpass date
Overpass_date <- as.data.frame(str_sub(Overpass[,1], start = 1, end = -72))
colnames(Overpass_date) <- "overpass_date"


# make table with date and overpass time
Overpass_time <- as.data.frame(str_sub(Overpass[,1], start = 13, end = -63))
colnames(Overpass_time) <- "overpass_time"

# generate a sequence of dates
days_1 <- seq(ISOdate(2018,01,01), by = "day", length.out = 1)
days_1 <- as.data.frame(str_sub(days_1, start = 1, end = -10))
colnames(days_1) <- "date"

days_2 <- seq(ISOdate(2018,01,03), by = "day", length.out = 15)
days_2 <- as.data.frame(str_sub(days_2, start = 1, end = -10))
colnames(days_2) <- "date"


days_3 <- seq(ISOdate(2018,01,19), by = "day", length.out = 15)
days_3 <- as.data.frame(str_sub(days_3, start = 1, end = -10))
colnames(days_3) <- "date"


days_4 <- seq(ISOdate(2018,02,04), by = "day", length.out = 15)
days_4 <- as.data.frame(str_sub(days_4, start = 1, end = -10))
colnames(days_4) <- "date"


days_5 <- seq(ISOdate(2018,02,20), by = "day", length.out = 102)
days_5 <- as.data.frame(str_sub(days_5, start = 1, end = -10))
colnames(days_5) <- "date"


days_6 <- seq(ISOdate(2018,06,03), by = "day", length.out = 15)
days_6 <- as.data.frame(str_sub(days_6, start = 1, end = -10))
colnames(days_6) <- "date"


days_7 <- seq(ISOdate(2018,06,19), by = "day", length.out = 42)
days_7 <- as.data.frame(str_sub(days_7, start = 1, end = -10))
colnames(days_7) <- "date"

# days_8 <- seq(ISOdate(2017,05,24), by = "day", length.out = 8)
# days_8 <- as.data.frame(str_sub(days_8, start = 1, end = -10))
# colnames(days_8) <- "date"


days <- rbind(days_1, days_2, days_3, days_4, days_5, days_6, days_7)

Overpass <- cbind(days, Overpass_date, Overpass_time)


########################################################################
######## start overpass matching with today date and time ##############
########################################################################


# match the same date of the sys time
match_date <- unique (grep(paste(current_date,collapse="|"), 
                           Overpass$date, value=TRUE))

# select only the current date
Overpass <- Overpass[grep(match_date, Overpass$date), ]

hour <- str_sub(Overpass$overpass_time[1], start = 1, end = -7)
minutes <- str_sub(Overpass$overpass_time[1], start = 4, end = -4)

# add +/- 5 minutes to the tile download


if (minutes == "00") {
  
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  minutes <- as.numeric(minutes) 
  current_overpass_time_plus <- paste0(".", hour,"0",minutes+5,".061.NRT")
  hour <- as.numeric(hour) 
  current_overpass_time_minus <- paste0(".","0", hour-1,"55",".061.NRT")
  
}


if ((minutes == "05")) {
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  minutes <- as.numeric(minutes) 
   current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,"0",minutes-5,".061.NRT")
  
}


if ( (minutes == "10") ) {
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  minutes <- as.numeric(minutes) 
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,"0",minutes-5,".061.NRT")
  
}


if ((minutes == "15")
    || (minutes == "20") || (minutes == "20") || (minutes == "25") || (minutes == "30")
    || (minutes == "35") || (minutes == "40") || (minutes == "45") || (minutes == "50")
    || (minutes == "55"))  {
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")

}


if (minutes == "01" || (minutes == "02") || (minutes == "03") || (minutes == "04")) {
  minutes <- "05"
  minutes <- as.numeric(minutes)
  current_overpass_time <- paste0(".",hour,"0",minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,"0", minutes-5,".061.NRT")
}


if (minutes == "06" || (minutes == "07") || (minutes == "08") || (minutes == "09")) {
  minutes <- "10"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,"0", minutes-5,".061.NRT")
}



if (minutes == "11" || (minutes == "12") || (minutes == "13") || (minutes == "14")) {
  minutes <- "15"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "16" || (minutes == "17") || (minutes == "18") || (minutes == "19")) {
  minutes <- "20"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "21" || (minutes == "22") || (minutes == "23") || (minutes == "24")) {
  minutes <- "25"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "26" || (minutes == "27") || (minutes == "28") || (minutes == "29")) {
  minutes <- "30"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "31" || (minutes == "32") || (minutes == "33") || (minutes == "34")) {
  minutes <- "35"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "36" || (minutes == "37") || (minutes == "38") || (minutes == "39")) {
  minutes <- "40"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "41" || (minutes == "42") || (minutes == "43") || (minutes == "44")) {
  minutes <- "45"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "46" || (minutes == "47") || (minutes == "48") || (minutes == "49")) {
  minutes <- "50"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_plus <- paste0(".",hour,minutes+5,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
}


if (minutes == "51" || (minutes == "52") || (minutes == "53") || (minutes == "54")) {
  minutes <- "55"
  minutes <- as.numeric(minutes) 
  current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
  current_overpass_time_minus <- paste0(".",hour,minutes-5,".061.NRT")
  hour <- as.numeric(hour)
  current_overpass_time_plus <- paste0(".","0",hour+1,"00",".061.NRT")
}

# replace 56, 57, 58, 59 with 00 and hour ==> + 1 !!!

if (minutes == "56" || (minutes == "57") || (minutes == "58") || (minutes == "59")) {
minutes <- "55"

 hour <- as.character(hour)
# hour <- paste0("0",hour)
 minutes <- as.numeric(minutes) 

 current_overpass_time <- paste0(".",hour,minutes,".061.NRT")
 hour <- as.numeric(hour)
 current_overpass_time_plus <- paste0(".","0",hour+1,"00",".061.NRT")
 current_overpass_time_minus <- paste0(".","0",hour,"50",".061.NRT")
}


Overpass_times <- c(current_overpass_time,current_overpass_time_plus,
                    current_overpass_time_minus)
write.csv(Overpass_times, paste0("/home/mariners/MODIS_AOD/Overpass_times_Aqua_",
                                 current_date,".csv"))

#### --------- ########## ------------ ############## -------------- ######

