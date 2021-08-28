## ERCOT
# Aug 1st, 2020
# read demand data from BA_tz.R
# read weather data from climavg_kaholo.R
# train splines and calculate flex load
# plot flex load and map coef
# update Feb 19, 2021
# add model comparison and coef range

rm(list = ls())

library(tidyverse)
library(data.table)
library(lubridate)
library(splines)
library(tseries)
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)
library(maptools)
library(sp)
library(ggspatial)
library(ggrepel)
library(sf)
library(lattice)
library(sandwich)
library(lmtest)
library(foreach)
library(doParallel)

demand_ercot <- readRDS("detz_ercot_new.RDS")
wdh_clim.dt <- readRDS("wdh_clim12.RDS")
wdh_clim.dt_2c <- readRDS("wdh_clim_2C.RDS")    # climate 2c

# data cleaning
colnames(demand_ercot)[2:3] <- c("demand", "sf_id")
de_ercot <- demand_ercot %>% 
  dplyr::select("sf_id", "date_time", "demand", "tzname", "newtz") %>% 
  mutate(weekday = weekdays(newtz)) %>% 
  rename(date_time_ntz = newtz) %>% 
  mutate(year = year(date_time_ntz)) %>% 
  mutate(month = month(date_time_ntz)) %>% 
  mutate(day = day(date_time_ntz)) %>% 
  mutate(yday = yday(date_time_ntz)) %>% 
  mutate(hour = as.numeric(substr(as.character(date_time_ntz), 12,13))) %>% 
  mutate(weekdaynum = as.integer(car::recode(weekday,
    "'Monday'='1';'Tuesday'='2';'Wednesday'='3';'Thursday'='4';
    'Friday'='5';'Saturday'='6';'Sunday'='7'"))) %>% 
  filter(year > 2015 & year < 2019) %>% 
  mutate(trend1 = 1:length(date_time_ntz)) %>% 
  mutate(yhour = (yday-1)*24 + hour) 

dercot_dt <- data.table(de_ercot)
dercot_dt$trend <- dercot_dt$year - 2015 # add yearly trend to the data set

summary(dercot_dt)

# combine with weather data (no climate case)
wdh_clim_e15 <- wdh_clim.dt[sf == 'm.t'] # select only ercot

DWercot <- dercot_dt %>% 
  left_join(wdh_clim_e15, by = c("date_time")) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) 
summary(DWercot) # did not remove NA values yet
DWercot <- data.table(DWercot)
DWercot <- DWercot[complete.cases(DWercot[, 'temp'])] 

saveRDS(DWercot, "DWercot_clim12_newdt.RDS")
#DW_rmna <- DWercot[complete.cases(DWercot[, 'demand'])]   # remove NAs
#saveRDS(DW_rmna, "Sisi/DWercot_clim12_new.RDS")

# climate 2c case
wdh_clim_2c <- wdh_clim.dt_2c[sf == 'm.t'] # climate 2c
DWercot2c <- dercot_dt %>% 
  left_join(wdh_clim_2c, by = c("date_time")) %>% 
  mutate(datentz = format(as.POSIXct(strptime(date_time_ntz,"%Y-%m-%d %H:%M:%S",tz="")), 
                          format = "%Y-%m-%d")) 
summary(DWercot2c) # did not remove NA values yet

DWercot2c <- data.table(DWercot2c)  
DWercot2c <- DWercot2c[complete.cases(DWercot2c[, 'temp'])]     # remove nonvalid temp value
#DW_rmna2c <- DWercot2c[complete.cases(DWercot2c[, 'demand'])]   # remove NAs

saveRDS(DWercot2c, "DWercot_clim12_2c_newdt.RDS")


###########################################
###### train opt splines ##################
###########################################

DWercot <- readRDS("DWercot_clim12_newdt.RDS")

source("functions.R")

# tune the first parameter: hour of the day
cv <- tune_new1(BAdata = DWercot)
cv_r <- data.table(data.frame(cv))
cv_rmse_r <- cv_r[, list(X2,X5,X8)]
avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])   # 19

# tune the second parameter
cv1 <- tune_new2(hour_df = hour_df, BAdata = DWercot)
cv_r1 <- data.table(data.frame(cv1))
cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])   # 6

#===============================================#
# train again with no climate interaction model #
#===============================================#
# tune the first parameter: hour of the day
cv <- tune_noclim1(BAdata = DWercot)
cv_r <- data.table(data.frame(cv))
cv_rmse_r <- cv_r[, list(X2,X5,X8)]
avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])   # 19

# tune the second parameter
cv1 <- tune_noclim2(hour_df = hour_df, BAdata = DWercot)
cv_r1 <- data.table(data.frame(cv1))
cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])   # 6

#############################################
#### run model comparison and coef range ####
#############################################

DW_rmna <- readRDS("DWercot_clim12_newdt.RDS")
hour_df <- 19
yhour_df <- 6

baseline <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum), 
               data = DW_rmna)

weathermodel <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                   hdh_18 + cdh_18, 
                   data = DW_rmna)

fullmodel <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
     hdh_18 + climH + cdh_18 + climC, 
     data = DW_rmna)


# calc rmse for each model and its reduction
s1 <- summary(baseline)
rmse_b <- sqrt(mean(s1$residuals^2))

s2 <- summary(weathermodel)
rmse_w <- sqrt(mean(s2$residuals^2))

s3 <- summary(fullmodel)
rmse_f <- sqrt(mean(s3$residuals^2))

r1 <- (rmse_b-rmse_w)/rmse_b*100
r2 <- (rmse_b-rmse_f)/rmse_b*100


### coef of cdh/hdh with nwsd
# NW sd function
nwsd.fn <- function(model){
  TT <- 24   # tried diff values - result did not change in this case
  m <- 0.75 * (length(TT))^(1/3)
  nw_sd <- NeweyWest(model,
                     lag = m - 1, prewhite = F,
                     adjust = T)
  return(nw_sd)
}

# coef for weather model
b0_c_w <- coef(weathermodel)["cdh_18"]
b0_h_w <- coef(weathermodel)["hdh_18"]
nwsd_w <- nwsd.fn(weathermodel)
c0_sd_w <- coeftest(weathermodel, vcov = nwsd_w)['cdh_18', 'Std. Error']
h0_sd_w <- coeftest(weathermodel, vcov = nwsd_w)['hdh_18', 'Std. Error']
coef_w <- c(b0_c_w, c0_sd_w, b0_h_w, h0_sd_w)

# coef for full model
beta0_c <- coef(fullmodel)["cdh_18"]
beta0_h <- coef(fullmodel)["hdh_18"]
beta1_c <- coef(fullmodel)["climC"]
beta1_h <- coef(fullmodel)["climH"]
nwsd_f <- nwsd.fn(fullmodel)
c0_sd <- coeftest(fullmodel, vcov = nwsd_f)['cdh_18', 'Std. Error']
h0_sd <- coeftest(fullmodel, vcov = nwsd_f)['hdh_18', 'Std. Error']
c1_sd <- coeftest(fullmodel, vcov = nwsd_f)['climC', 'Std. Error']
h1_sd <- coeftest(fullmodel, vcov = nwsd_f)['climH', 'Std. Error']
c0c1_cov <- nwsd_f['cdh_18', 'climC']
h0h1_cov <- nwsd_f['hdh_18', 'climH']

coef_f <- c( beta0_c, c0_sd, beta0_h, h0_sd, beta1_c, c1_sd, beta1_h, h1_sd, c0c1_cov, h0h1_cov)


### load clim per cell data and matching cells to BA regoin with pop per cell data
clim <- readRDS("clim_cell12.RDS")
pop.dt <- readRDS("popweight.RDS")

pop_ercot <- pop.dt[sf == 'm.t']
climcell_ba <- clim[pop_ercot, on = 'cell']

## calc pop weighted avg for cdh/hdh
clim_avg <- climcell_ba %>% group_by(ID) %>% 
  summarise(climC_avg = sum(climC*w), climH_avg = sum(climH*w))

# calc the range of beta_i
ba1 <- climcell_ba
beta_ci_max <- beta0_c + beta1_c*max(ba1$climC)
beta_ci_min <- beta0_c + beta1_c*min(ba1$climC)
beta_ci_avg <- beta0_c + beta1_c*clim_avg$climC_avg
beta_hi_max <- beta0_h + beta1_h*max(ba1$climH)
beta_hi_min <- beta0_h + beta1_h*min(ba1$climH)
beta_hi_avg <- beta0_h + beta1_h*clim_avg$climH_avg

beta_ci <- c(beta_ci_min, beta_ci_max, beta_ci_avg)
beta_hi <- c(beta_hi_min, beta_hi_max, beta_hi_avg)
beta_ci <- data.frame(beta_ci)
beta_hi <- data.frame(beta_hi)
rownames(beta_ci) <- c("Max", "Min", "Avg")
rownames(beta_hi) <- c("Max", "Min", "Avg")


# calc sd for the range of beta_i
climCi <- ba1$climC
beta_ci_all <- beta0_c + beta1_c*climCi
var_bi_c <- (c0_sd)^2 + (c1_sd)^2*(climCi)^2 + 2*climCi*c0c1_cov
sd_bi_c <- sqrt(var_bi_c)

climHi <- ba1$climH
beta_hi_all <- beta0_h + beta1_h*climHi
var_bi_h <- (h0_sd)^2 + (h1_sd)^2*(climHi)^2 + 2*climHi*h0h1_cov
sd_bi_h <- sqrt(var_bi_h)

# calc 95% CI
ci <- 0.95
a <- qnorm(1-(1-ci)/2)

# weather model CI
l0w_c <- b0_c_w - a*c0_sd_w
h0w_c <- b0_c_w + a*c0_sd_w
CI_w_c <- cbind(l0w_c, h0w_c)

l0w_h <- b0_h_w - a*h0_sd_w
h0w_h <- b0_h_w + a*h0_sd_w
CI_w_h <- cbind(l0w_h, h0w_h)

# full model CI
l1_c <- beta_ci_all - a*sd_bi_c
h1_c <- beta_ci_all + a*sd_bi_c
beta_c_ci <- cbind(min(l1_c), max(h1_c))

l1_h <- beta_hi_all - a*sd_bi_h
h1_h <- beta_hi_all + a*sd_bi_h
beta_h_ci <- cbind(min(l1_h), max(h1_h))

###==================================###
#### calc out of sample R2 (CV) ########
###===================================##
source("functions.R")

# baseline cv
itr=3     # total number of fold cv
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

cv <- matrix(NA, itr, 3)
trainR2 <- matrix(NA, itr, 1)
trainR2_adj <- matrix(NA, itr, 1)
cv_result <- foreach(k=1:itr, .combine=rbind, .packages= c("data.table", "splines")) %dopar% { 
  train <- DW_rmna[year!=(2015+k), ]
  test <- DW_rmna[year==(2015+k), ]
  fit_cv <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum), 
               data = train)
  trainR2_adj[k,] <- var(exp(predict(fit_cv, train)))/var(train$demand)   # adjust R2 to compare linear case
  trainR2[k,] <- summary(fit_cv)$r.squared
  cv[k,] <- RMSE_R2_log_adj(fit_cv, test, train)
  c(trainR2[k,], trainR2_adj[k,], cv[k,])
}

#stop cluster
stopCluster(cl)

baseline_r2 <- colMeans(cv_result)
baseline_r2 <- data.frame(baseline_r2)
rownames(baseline_r2) <- c("R2", "adjR2", "RMSE_test", "R2_test", "adjR2_test")


# weather model cv
itr=3     # total number of fold cv
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

cv <- matrix(NA, itr, 3)
trainR2 <- matrix(NA, itr, 1)
trainR2_adj <- matrix(NA, itr, 1)
cv_result_w <- foreach(k=1:itr, .combine=rbind, .packages= c("data.table", "splines")) %dopar% { 
  train <- DW_rmna[year!=(2015+k), ]
  test <- DW_rmna[year==(2015+k), ]
  fit_cv <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum)
               + cdh_18 + hdh_18, 
               data = train)
  trainR2_adj[k,] <- var(exp(predict(fit_cv, train)))/var(train$demand)   # adjust R2 to compare linear case
  trainR2[k,] <- summary(fit_cv)$r.squared
  cv[k,] <- RMSE_R2_log_adj(fit_cv, test, train)
  c(trainR2[k,], trainR2_adj[k,], cv[k,])
}

#stop cluster
stopCluster(cl)

weather_r2 <- colMeans(cv_result_w)
weather_r2 <- data.frame(weather_r2)
rownames(weather_r2) <- c("R2", "adjR2", "RMSE_test", "R2_test", "adjR2_test")

# full model
itr=3     # total number of fold cv
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

cv <- matrix(NA, itr, 3)
trainR2 <- matrix(NA, itr, 1)
trainR2_adj <- matrix(NA, itr, 1)
cv_result_f <- foreach(k=1:itr, .combine=rbind, .packages= c("data.table", "splines")) %dopar% { 
  train <- DW_rmna[year!=(2015+k), ]
  test <- DW_rmna[year==(2015+k), ]
  fit_cv <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum)
               + cdh_18 + hdh_18 + climC + climH, 
               data = train)
  trainR2_adj[k,] <- var(exp(predict(fit_cv, train)))/var(train$demand)   # adjust R2 to compare linear case
  trainR2[k,] <- summary(fit_cv)$r.squared
  cv[k,] <- RMSE_R2_log_adj(fit_cv, test, train)
  c(trainR2[k,], trainR2_adj[k,], cv[k,])
}

#stop cluster
stopCluster(cl)

full_r2 <- colMeans(cv_result_f)
full_r2 <- data.frame(full_r2)
rownames(full_r2) <- c("R2", "adjR2", "RMSE_test", "R2_test", "adjR2_test")

R2_results <- cbind(baseline_r2, weather_r2, full_r2)
(R2_results$baseline_r2[3]-R2_results$weather_r2[3])/R2_results$baseline_r2[3]

# use the adj R2 from full data i/o avg training set

# baseline_r2 weather_r2    full_r2
# R2           0.8311424 0.94260998 0.94304457
# adjR2        0.8314298 0.94672451 0.94736901
# RMSE_test    0.1061496 0.06317725 0.06319041
# R2_test      0.7692611 0.91849000 0.91845030
# adjR2_test   0.8369853 0.95016696 0.95570187

## result with clean data
# baseline_r2 weather_r2    full_r2
# R2           0.8357990 0.94461576 0.94506245
# adjR2        0.8417005 0.95519749 0.95521155
# RMSE_test    0.1013921 0.05793838 0.05806965
# R2_test      0.7905517 0.93160687 0.93130693
# adjR2_test   0.8195554 0.94552258 0.95305135


############################################
######## predict flex load #################
############################################

fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                  hdh_18 + climH + cdh_18 + climC, 
                data = DW_rmna)

# calculate flex load for each BA
DW_0 <- DW_rmna %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0) %>% mutate(climH = 0) %>% mutate(climC = 0)
logy_hat0 <- predict(fit_final, DW_0)
logy_hat1 <- predict(fit_final, DW_rmna)
e <- resid(fit_final)
mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
y_hat1 <- exp(logy_hat1)*mean_e
y_hat0 <- exp(logy_hat0)*mean_e
s <- (y_hat1 - y_hat0)/y_hat1
DW_rmna$yhat0 <- y_hat0
DW_rmna$yhat1 <- y_hat1
DW_rmna$flex <- s

saveRDS(DW_rmna, "flexload_ercot_new.RDS")

# save the reg coef
beta0_c <- coef(fit_final)["cdh_18"]
beta0_c_std <- coef(summary(fit_final))['cdh_18', 'Std. Error']
beta1_c <- coef(fit_final)["climC"]
beta1_c_std <- coef(summary(fit_final))['climC', 'Std. Error']
beta0_h <- coef(fit_final)["hdh_18"]
beta0_h_std <- coef(summary(fit_final))['hdh_18', 'Std. Error']
beta1_h <- coef(fit_final)["climH"]
beta1_h_std <- coef(summary(fit_final))['climH', 'Std. Error']
r2adj <- summary(fit_final)$adj.r.squared
c(beta0_c, beta0_c_std, beta1_c, beta1_c_std, beta0_h, beta0_h_std, beta1_h, beta1_h_std, r2adj)

# group for summer hours
Hour_summer <- DW_rmna[month == 6| month == 7| month == 8, ]
# group for winter hours
Hour_winter <- DW_rmna[month == 12| month == 1| month == 2, ]

# calculate the hourly avg for each BA and the min and max of all flex load
# winter
houravg_w <- Hour_winter %>% group_by(hour) %>% summarise(avgflex = mean(flex))
min(houravg_w)    # 0.01268
max(houravg_w)    # 0.28256

# summer
houravg_s <- Hour_summer %>% group_by(hour) %>% summarise(avgflex = mean(flex))
min(houravg_s)    # 0.028
max(houravg_s)    # 0.4678


####################################################
##### climate 2c prediction (same as in 7) #########
####################################################

# load demand weather data with NAs for prediction
DW_clim <- readRDS("DWercot_clim12_2c_newdt.RDS")

# load demand weather data without NAs for regression
DW_clim_new <- readRDS("DWercot_clim12_newdt.RDS")
summary(DW_clim_new)

# check the optimal spline from Ercot.R file
splines_opt <- list()
splines_opt[[1]] <- c(19,6)

## run the above code without loop
i <- 1
DW1_rmna <- DW_clim_new                   
# choose the opt splines for each BA
hour_df <- splines_opt[[i]][1]
yhour_df <- splines_opt[[i]][2]
fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                  hdh_18 + climH + cdh_18 + climC, 
                data = DW1_rmna)
# predict demand (including NAs) for each BA
DW1       <- DW_clim %>% 
  filter(year > 2015) %>% 
  mutate(trend1 = 1:length(date_time_ntz)) %>% 
  mutate(yhour = (yday-1)*24 + hour)  
DW1_0     <- DW1 %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0) %>% mutate(climH = 0) %>% mutate(climC = 0)
logy_hat0 <- predict(fit_final, DW1_0)
logy_hat1 <- predict(fit_final, DW1)
# calculate resid to take exponential
e <- resid(fit_final)
yc <- exp(logy_hat1+e)
mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
y_hat1 <- exp(logy_hat1)*mean_e
y_hat0 <- exp(logy_hat0)*mean_e
s <- (y_hat1 - y_hat0)/y_hat1
DW1$yhat0 <- y_hat0
DW1$yhat1 <- y_hat1
DW1$flex <- s
DW1$yc <- yc

# DW1 is the result
summary(DW1)
saveRDS(DW1, "predictErcot_2C_yc_new.RDS")


####################################################
######## climate 2c remodel (same as in 8) #########
####################################################

# load demand weather data with NAs for prediction
DW_clim <- readRDS("DWercot_clim12_2c_newdt.RDS")

# load demand weather data without NAs for regression
DW_clim_new <- readRDS("DWercot_clim12_newdt.RDS")

# check the optimal spline from Ercot.R file
splines_opt <- list()
splines_opt[[1]] <- c(19,6)

## run the above code without loop
i <- 1
DW1_rmna <- DW_clim_new                   
# choose the opt splines for each BA
hour_df <- splines_opt[[i]][1]
yhour_df <- splines_opt[[i]][2]
fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                  hdh_18 + cdh_18, 
                data = DW1_rmna)
# predict demand (including NAs) for each BA
DW1       <- DW_clim %>% 
  filter(year > 2015) %>% 
  mutate(trend1 = 1:length(date_time_ntz)) %>% 
  mutate(yhour = (yday-1)*24 + hour)  
DW1_0     <- DW1 %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0) 
logy_hat0 <- predict(fit_final, DW1_0)
logy_hat1 <- predict(fit_final, DW1)
# calculate resid to take exponential
e <- resid(fit_final)
mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
y_hat1 <- exp(logy_hat1)*mean_e
y_hat0 <- exp(logy_hat0)*mean_e
s <- (y_hat1 - y_hat0)/y_hat1
DW1$yhat0 <- y_hat0
DW1$yhat1 <- y_hat1
DW1$flex <- s

# DW1 is the result
summary(DW1)
saveRDS(DW1, "flexload_Ercot_2C_new.RDS")



