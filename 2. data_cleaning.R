# 2. load clean data from EIA_Cleaned_Hourly_Electricity_Demand_Data-v1.1
# June 4, 2021
# read the interconnects data
# merge BAs and aggregate with the clean data
# adjust clean data to fit the regression and combine with weather data from Kaholo

rm(list = ls())

library(tidyverse)
library(data.table)
library(sp)
library(rgdal)
library(lutz)
library(lubridate)

##### load clean data and combine with merged ba
all_files <- list.files(path = "balancing_authorities/", pattern = "*.csv", full.names = TRUE)
all_files_s <- list.files(path = "balancing_authorities/", pattern = "*.csv")

# read file content
all_content <- all_files %>% lapply(read.csv, header = TRUE)
# read file name
all_filenames <- all_files_s %>% basename() %>% as.list()
# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
all_result <- rbindlist(all_lists, fill = T)
names(all_result)[6] <- "File.Path"
all_result$ba <- sub("\\..*", "", all_result$File.Path)

demand.dt <- all_result[, c('date_time', 'cleaned.demand..MW.', 'ba')]
names(demand.dt)[2] <- "value"
demand.dt$date_time <- as.POSIXct(as.character(demand.dt$date_time), format="%Y-%m-%d %H:%M:%S", tz = "UTC")

saveRDS(demand.dt, "demand_raw_newdt.RDS")

#combine BAs to match aggregated BA shapefiles
matching.ids <- fread(file = "eia_sf_match.tsv")
setkey(matching.ids, id)
east.ids <- readRDS("BA_SF/ids/east_all.RDS")
setnames(east.ids, "sf_id", "id")
setkey(east.ids, "id")
west.ids <- readRDS("BA_SF/ids/west_all.RDS")
setnames(west.ids, "sf_id", "id")
setkey(west.ids, "id")
#match BA name to a rowid
east.key <- east.ids[matching.ids, nomatch = 0][, list(rowid, id, ba)]
west.key <- west.ids[matching.ids, nomatch = 0][, list(rowid, id, ba)]

# function to reshape the rowid
reshape.ids <- function(f_path, sf_num){
  dt.ids <- readRDS(f_path)
  dt.ids <- cbind(dt.ids, dt.ids[, tstrsplit(rowid, split = ".", fixed = T, type.convert = T)])
  dt.ids[, rowid := NULL]
  dt.ids <- melt(dt.ids, id.vars = "sf_id", na.rm = T, value.name = "rowid")[, variable := NULL]
  setkey(dt.ids, sf_id, rowid)
  dt.ids[, sf := sf_num]
  return(dt.ids)
}

east.ids <- reshape.ids("BA_SF/ids/east_15.RDS", 15)
setkey(east.ids, sf, sf_id)
west.ids <- reshape.ids("BA_SF/ids/west_15.RDS", 15)
setkey(west.ids, sf, sf_id)

# combine reshaped rowids and ba names
east.ids <- east.ids[east.key, on = "rowid"]
setkey(east.ids, ba)
west.ids <- west.ids[west.key, on = "rowid"]
setkey(west.ids, ba)

#merge and calc total demand by sf_id
sum.by.sf_id <- function(data.dt, id.dt, sf_num){
  dt <- data.dt[id.dt[sf == sf_num], nomatch = 0][, list(demand = sum(value)), by = list(sf_id, date_time)]
  dt[, sf := sf_num]
  setkey(dt, sf, sf_id, date_time)
  return(dt)
}
setkey(demand.dt, ba)

demand.east <- sum.by.sf_id(demand.dt, east.ids, 15)
setkey(demand.east, sf, sf_id, date_time)

demand.west <- sum.by.sf_id(demand.dt, west.ids, 15)
setkey(demand.west, sf, sf_id, date_time)

saveRDS(demand.east, "demande_cleandata.RDS")
saveRDS(demand.west, "demandw_cleandata.RDS")

####################################################
###### change BA time zone to local time zone ######
####################################################
# copy from BA_tz.R

## load 15 BAs shapefile in the east/west and east/west demand data
east15_shp <- readOGR("BA_SF/aggregated/east/east_15.shp")
west15_shp <- readOGR("BA_SF/aggregated/west/west_15.shp")
demand_e <- readRDS("demande_cleandata.RDS")
demand_w <- readRDS("demandw_cleandata.RDS")

# transform the coordinates to long/lat coordinate in order to match time zones
east15_trans <- spTransform(east15_shp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
latlong <- coordinates(east15_trans) # check coordinate [long, lat]
latlong <- data.frame(latlong) %>% 
  mutate(ID = seq(1:nrow(latlong)))

# plot(east15_trans[1,])

## choose the 15 BAs join with the center point (lat, long) of the BA area and corresponding tz
demand_e15 <- demand_e %>%      
  filter(sf == '15') %>% 
  left_join(latlong, by = c("sf_id" = "ID")) %>% # X1=long X2=lat
  mutate(tzname = tz_lookup_coords(X2, X1, method = "accurate")) 

df_e <- data.frame(demand_e15, stringsAsFactors=FALSE)
tz_v <- Vectorize(function(x,y) {format(x, tz=y, usetz=TRUE)})
df_e$newtz <- tz_v(df_e$date_time, df_e$tzname)   # newtz is character
df_e <- data.table(df_e)

saveRDS(df_e, "de15tz_new.RDS")

###======================
## repeat for the west ##
##=======================

# transform the coordinates to long/lat coordinate in order to match time zones
west15_trans <- spTransform(west15_shp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
latlong <- coordinates(west15_trans) # check coordinate [long, lat]
latlong <- data.frame(latlong) %>% 
  mutate(ID = seq(1:nrow(latlong)))

# plot(west15_trans[1,])

## choose the 15 BAs join with the center point (lat, long) of the BA area and corresponding tz
# tzlist <- tz_list() # list of tzname matches with tz and daytime saving indication
demand_w15 <- demand_w %>%      
  filter(sf == '15') %>% 
  left_join(latlong, by = c("sf_id" = "ID")) %>% # X1=long X2=lat
  mutate(tzname = tz_lookup_coords(X2, X1, method = "accurate")) 

# split into different time zones for time changes (calculation purpose)
df_w <- data.frame(demand_w15, stringsAsFactors=FALSE)
tz_v <- Vectorize(function(x,y) {format(x, tz=y, usetz=TRUE)})
df_w$newtz <- tz_v(df_w$date_time, df_w$tzname)   # newtz is character
df_w <- data.table(df_w)

saveRDS(df_w, "dw15tz_new.RDS")


##=============================
#### repeat for ERCOT ########
##=============================

# load raw demend and filter ercot data
D_raw <- readRDS("demand_raw_newdt.RDS")          # demand data from clean data
ercot_d <- D_raw %>% filter(ba == "ERCO") 
ercot_dt <- as.data.table(ercot_d)
ercot_dt$ba <- 1
ercot_dt$ba <- as.numeric(ercot_dt$ba)
summary(ercot_dt)

## load 15 BAs shapefile in the ercot
ercot_shp <- readOGR("BA_SF/aggregated/ercot/ercot.shp")

# transform the coordinates to long/lat coordinate in order to match time zones
ercot_trans <- spTransform(ercot_shp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
latlong <- coordinates(ercot_trans) # check coordinate [long, lat]
latlong <- data.frame(latlong) %>% 
  mutate(ID = seq(1:nrow(latlong)))

plot(ercot_trans[1,])

## join with the center point (lat, long) of the BA area and corresponding tz
demand_ercot <- ercot_dt %>%      
  left_join(latlong, by = c("ba" = "ID")) %>% # X1=long X2=lat
  mutate(tzname = tz_lookup_coords(X2, X1, method = "accurate")) 

# split into different time zones for time changes (calculation purpose)
n <- length(unique(demand_ercot$tzname)) # n = 2
d_ntz <- list()
for (i in 1:n){
  d_ntz[[i]] <- demand_ercot %>% 
    filter(tzname == unique(tzname)[i]) %>% 
    mutate(newtz = with_tz(date_time, tzone = unique(tzname)[i])) # considering daylight savings
}
# merge the list into df
d_ntz_df <- do.call(rbind.data.frame, d_ntz)

saveRDS(d_ntz_df, "detz_ercot_new.RDS")



#############################################################
########## East IC: adjust the demand for regression ########
############################################################
# copy from 3b_Data_cleaning.R

demand15_east <- readRDS("de15tz_new.RDS")

# adjust data for the regressions
de15east <- demand15_east %>% 
  select(sf_id, date_time, demand, tzname, newtz) %>% 
  mutate(weekday = weekdays(as.Date(newtz))) %>% 
  rename(date_time_ntz = newtz) %>% 
  mutate(year = year(date_time_ntz)) %>% 
  mutate(month = month(date_time_ntz)) %>% 
  mutate(day = day(date_time_ntz)) %>% 
  mutate(yday = yday(date_time_ntz)) %>% 
  mutate(hour = as.numeric(substr(date_time_ntz, 12,13))) %>% 
  filter(year > 2015)

de15east <- data.table(de15east)
de15east[is.na(hour), hour := 0]
summary(de15east)

de15east_df <- de15east %>% 
  mutate(weekdaynum = as.integer(car::recode(weekday,
    "'Monday'='1';'Tuesday'='2';'Wednesday'='3';'Thursday'='4';
    'Friday'='5';'Saturday'='6';'Sunday'='7'")))
de15east_dt <- data.table(de15east_df)
de15east_dt$trend <- de15east_dt$year - 2015 # add yearly trend to the data set

saveRDS(de15east_dt, "demand15e_new.RDS") # use for combining with weather data

##==============================
#### run the same for west IC ##
##==============================

demand15_west <- readRDS("dw15tz_new.RDS")

# adjust data for the regressions
de15west <- demand15_west %>% 
  select(sf_id, date_time, demand, tzname, newtz) %>% 
  mutate(weekday = weekdays(as.Date(newtz))) %>% 
  rename(date_time_ntz = newtz) %>% 
  mutate(year = year(date_time_ntz)) %>% 
  mutate(month = month(date_time_ntz)) %>% 
  mutate(day = day(date_time_ntz)) %>% 
  mutate(yday = yday(date_time_ntz)) %>% 
  mutate(hour = as.numeric(substr(date_time_ntz, 12,13))) %>% 
  filter(year > 2015)

de15west <- data.table(de15west)
de15west[is.na(hour), hour := 0]
summary(de15west)

de15west_df <- de15west %>% 
  mutate(weekdaynum = as.integer(car::recode(weekday,
    "'Monday'='1';'Tuesday'='2';'Wednesday'='3';'Thursday'='4';
    'Friday'='5';'Saturday'='6';'Sunday'='7'")))
de15west_dt <- data.table(de15west_df)
de15west_dt$trend <- de15west_dt$year - 2015 # add yearly trend to the data set

saveRDS(de15west_dt, "demand15w_new.RDS") # use for combining with weather data


#############################################
### combine with weather data from Kaholo ###
#############################################

# load climate weather data 
wdh_clim.dt <- readRDS("wdh_clim12.RDS")
# load demand data for east 15
demande15 <- readRDS("demand15e_new.RDS")
wdh_clim_e15 <- wdh_clim.dt[sf == 'e.15'] # select only the east 15 BAs

DWe15 <- demande15 %>% 
  left_join(wdh_clim_e15, by = c("date_time", "sf_id" = "ID")) %>% 
  dplyr::select(sf_id, date_time, demand, date_time_ntz, year, month, day, yday, hour, 
                weekdaynum, trend, temp, hdh_18, cdh_18, climH, climC) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) %>% 
  filter(year != 2019)  # weather data till 2018 for now

DWe15 <- data.table(DWe15)
DWe15[is.na(datentz), datentz := format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d", tz="")), 
                                       format = "%Y-%m-%d")]
summary(DWe15)          # did not remove NA values yet

DWe15_clim <- data.table(DWe15)
DWe15_clim <- DWe15_clim[year > 2015, ]
DWe15_clim[, trend1 := 1:length(date_time_ntz), by = "sf_id"]
DWe15_clim[, yhour := (yday-1)*24 + hour]

DWe15_clim <- DWe15_clim[complete.cases(DWe15_clim[, ])]   # remove 96 NAs from temp
summary(DWe15_clim)
saveRDS(DWe15_clim, "DWe15_clim12_newdt.RDS")

##==============================
#### run the same for west IC ##
##==============================

# load climate weather data 
wdh_clim.dt <- readRDS("wdh_clim12.RDS")
# load demand data for east 15
demandw15 <- readRDS("demand15w_new.RDS")
wdh_clim_w15 <- wdh_clim.dt[sf == 'w.15'] # select only the west 15 BAs

DWw15 <- demandw15 %>% 
  left_join(wdh_clim_w15, by = c("date_time", "sf_id" = "ID")) %>% 
  dplyr::select(sf_id, date_time, demand, date_time_ntz, year, month, day, yday, hour, 
                weekdaynum, trend, temp, hdh_18, cdh_18, climH, climC) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) %>% 
  filter(year != 2019)  # weather data till 2018 for now

DWw15 <- data.table(DWw15)
DWw15[is.na(datentz), datentz := format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d", tz="")), 
                                        format = "%Y-%m-%d")]
summary(DWw15)          # did not remove NA values yet

DWw15_clim <- data.table(DWw15)
DWw15_clim <- DWw15_clim[year > 2015, ]
DWw15_clim[, trend1 := 1:length(date_time_ntz), by = "sf_id"]
DWw15_clim[, yhour := (yday-1)*24 + hour]

DWw15_clim <- DWw15_clim[complete.cases(DWw15_clim[, ])]   # remove 124 NAs from temp
summary(DWw15_clim)
saveRDS(DWw15_clim, "DWw15_clim12_newdt.RDS")

