# climate change case - rerun the model result
# Aug 20, 2020
# temperature rise by 2C everywhere
# need to run in koholo
# added yc (demand under climate change) on Feb 28, 2021
# update Jun 6, 2021 with clean data for demand

rm(list = ls())

library(tidyverse)
library(lubridate)
library(data.table)
library(splines)
library(tseries)
library(effects)
library(foreach)
library(doParallel)
library(gridExtra)

source("functions.R")

# read all temp per cell data, and add 2C to value, adjust cdh and hdh
t1 <- readRDS("temp_by_cell.RDS")
t2 <- readRDS("temp_by_cell1.RDS")
t3 <- readRDS("temp_by_cell2.RDS")

t1[, value := value + 2]
t1[, cdh_18 := pmax(value - 18, 0)]
t1[, hdh_18 := pmax(18 - value, 0)]

t2[, value := value + 2]
t2[, cdh_18 := pmax(value - 18, 0)]
t2[, hdh_18 := pmax(18 - value, 0)]

t3[, value := value + 2]
t3[, cdh_18 := pmax(value - 18, 0)]
t3[, hdh_18 := pmax(18 - value, 0)]

# only climate variables
t1c <- t1[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
mean(t1c$climH) # 6.26
mean(t1c$climC) # 3.57
t2c <- t2[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
t3c <- t3[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
t_c <- rbind(t1c, t2c, t3c)
setkey(t_c, cell)

# 12 year climate avg per cell
tc_avg <- t_c[, list(climH = mean(climH), climC = mean(climC)), by = list(cell)]

# merge the climate avg with temp per cell data t1 by cell
t_merge <- t1[tc_avg, ]
t_merge <- t_merge[, list(value, date_time, cdh_18, hdh_18, climH, climC, 
                              climH_H = climH*hdh_18, climC_C = climC*cdh_18), 
                       by = list(cell)]
saveRDS(t_merge, "temp_cell_clim2C.RDS")

## continue with the same code as in 3_calcBAWDD.R
#merge in temp values to pop
pop.dt <- readRDS("popweight.RDS")
pop.cells <- pop.dt[, unique(cell)]
pop.dt <- split(pop.dt, by = c("sf", "ID"), drop = T)
lapply(pop.dt, function(dt) setkey(dt, cell))
setkey(t_merge, cell)

wdh_clim.dt <- lapply(pop.dt, calc.wdh.by.ID.clim, full.temp.dt = t_merge)
wdh_clim.dt <- rbindlist(wdh_clim.dt) # climH and climC are climate and hdh/cdh interaction term
setkey(wdh_clim.dt, sf, ID, date_time)
saveRDS(wdh_clim.dt, "wdh_clim_2C.RDS")

# below same code as climavg_kaholo.R
##########################################################
## merge climate data with demand data                 ##
## join the two dataset based on date_time and sf_id   ##
##########################################################

wdh_clim.dt <- readRDS("wdh_clim_2C.RDS")
demande15 <- readRDS("demand15e_new.RDS")  # updated with clean data
wdh_clim_e15 <- wdh_clim.dt[sf == 'e.15'] # select only the east 15 BAs

DWe15 <- demande15 %>% 
  left_join(wdh_clim_e15, by = c("date_time", "sf_id" = "ID")) %>% 
  select(sf_id, date_time, demand, date_time_ntz, year, month, day, yday, hour, 
         weekdaynum, trend, temp, hdh_18, cdh_18, climH, climC) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) %>% 
  filter(year != 2019) # weather data till 2018 for now
summary(DWe15) # did not remove NA values yet
saveRDS(DWe15, "DWe15_clim12_2C_new.RDS")

## load demand data for west 15 
demandw15 <- readRDS("demand15w_new.RDS")
wdh_clim_w15 <- wdh_clim.dt[sf == 'w.15'] # select only the west 15 BAs

DWw15 <- demandw15 %>% 
  left_join(wdh_clim_w15, by = c("date_time", "sf_id" = "ID")) %>% 
  select(sf_id, date_time, demand, date_time_ntz, year, month, day, yday, hour, 
         weekdaynum, trend, temp, hdh_18, cdh_18, climH, climC) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) %>% 
  filter(year != 2019) # weather data till 2018 for now
summary(DWw15) # did not remove NA values yet
saveRDS(DWw15, "DWw15_clim12_2C_new.RDS")

# below code is the same as 5. predictNAs - adjusted April 28, 2021: remove NAs, just predict filled data
#############################################################################
#########   run regression for predicted demand with climate change   #######
#############################################################################
# switch east and west readRDS and run the same code
# load demand weather data without removing NAs for prediction (+2C climate change)
DWe15_clim <- readRDS("DWe15_clim12_2C_new.RDS")     # east
DWe15_clim <- readRDS("DWw15_clim12_2C_new.RDS")     # west

DWe15_clim <- DWe15_clim[complete.cases(DWe15_clim[, "temp"]), ]   # remove only the temp NAs
summary(DWe15_clim)
DWe15_clim <- data.table(DWe15_clim)

# load demand weather data without NAs for regression
DWe15_clim_new <- readRDS("DWe15_clim12_newdt.RDS")     # east
DWe15_clim_new <- readRDS("DWw15_clim12_newdt.RDS")     # west
summary(DWe15_clim_new)
DWe15_clim_new <- data.table(DWe15_clim_new)

# load optimal splines
splines_opt <- readRDS("splines_opt_e.RDS")        # east
splines_opt <- readRDS("splines_opt_w.RDS")        # west

# run in parallel 
detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

DW1_all <- list()
results_reg <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim_new[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                    hdh_18 + climH + cdh_18 + climC, 
                  data = DW1_rmna)
  # predict demand (including NAs) for each BA
  DW1       <- DWe15_clim[sf_id == i, ] %>% 
    filter(year > 2015) %>% 
    mutate(trend1 = 1:length(date_time_ntz)) %>% 
    mutate(yhour = (yday-1)*24 + hour)  
  DW1_0     <- DW1 %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0) %>% mutate(climH = 0) %>% mutate(climC = 0)
  logy_hat0 <- predict(fit_final, DW1_0)
  logy_hat1 <- predict(fit_final, DW1)
  # calculate resid to take exponential
  e <- resid(fit_final)
  yc <- exp(logy_hat1+e)   # Mar 4 2021: e has diff length to logy_hat1 - check again
  mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
  y_hat1 <- exp(logy_hat1)*mean_e
  y_hat0 <- exp(logy_hat0)*mean_e
  s <- (y_hat1 - y_hat0)/y_hat1
  DW1$yhat0 <- y_hat0
  DW1$yhat1 <- y_hat1
  DW1$flex <- s
  DW1$yc <- yc
  
  # save DW1_rmna with flex load for each BA
  DW1_all[[i]] <- DW1
}

# stop cluster
stopCluster(cl)

summary(results_reg[[15]])
saveRDS(results_reg, "predicte15_2C_new.RDS")    # added yc column
saveRDS(results_reg, "predictw15_2C_new.RDS")


