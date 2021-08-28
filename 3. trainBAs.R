## 3. train splines with clean data
# Mar 7, 2021
# run in Kaholo

rm(list = ls())

library(tidyverse)
library(lubridate)
library(data.table)
library(splines)
library(tseries)
library(effects)
library(foreach)
library(doParallel)

source("functions.R")

# read demand weather data 
DWe15_clim <- readRDS("DWe15_clim12_newdt.RDS")
summary(DWe15_clim)
DWe15_clim <- data.table(DWe15_clim)

# loop for training all the eastern 15 BAs: returns optimal splines for each BA
# run in parallel loops
detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

optimal_sp <- matrix(NA, 15, 2)

results_sp <- foreach(i=1:15, .packages = c("data.table", "splines", "dplyr")) %dopar% {
  DW1 <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  
  # tune the first parameter: hour of the day
  cv <- tune_new1(BAdata = DW1)
  cv_r <- data.table(data.frame(cv))
  cv_rmse_r <- cv_r[, list(X2,X5,X8)]
  avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
  hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])
  
  # tune the second parameter
  cv <- tune_new2(hour_df = hour_df, BAdata = DW1)
  cv_r1 <- data.table(data.frame(cv))
  cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
  avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
  yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])
  optimal_sp[i, ] <- c(hour_df, yhour_df)
}

# stop cluster
stopCluster(cl)

saveRDS(results_sp, "splines_opt_e.RDS")   

#================================================#
# run again to train without climate interaction #
#================================================#

detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

optimal_sp <- matrix(NA, 15, 2)

results_sp <- foreach(i=1:15, .packages = c("data.table", "splines", "dplyr")) %dopar% {
  DW1 <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  
  # tune the first parameter: hour of the day
  cv <- tune_noclim1(BAdata = DW1)
  cv_r <- data.table(data.frame(cv))
  cv_rmse_r <- cv_r[, list(X2,X5,X8)]
  avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
  hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])
  
  # tune the second parameter
  cv <- tune_noclim2(hour_df = hour_df, BAdata = DW1)
  cv_r1 <- data.table(data.frame(cv))
  cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
  avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
  yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])
  optimal_sp[i, ] <- c(hour_df, yhour_df)
}

# stop cluster
stopCluster(cl)

saveRDS(results_sp, "splines_noclim_e.RDS")   

###########################################################
########## run the same for west IC #######################
###########################################################


DWw15_clim <- readRDS("DWw15_clim12_newdt.RDS")
summary(DWw15_clim)
DWw15_clim <- data.table(DWw15_clim)

detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

optimal_sp <- matrix(NA, 15, 2)

results_sp <- foreach(i=1:15, .packages = c("data.table", "splines", "dplyr")) %dopar% {
  DW1_rmna <- DWw15_clim[sf_id == i, ] 
  # tune the first parameter: hour of the day
  cv <- tune_new1(BAdata = DW1_rmna)
  cv_r <- data.table(data.frame(cv))
  cv_rmse_r <- cv_r[, list(X2,X5,X8)]
  avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
  hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])
  
  # tune the second parameter
  cv <- tune_new2(hour_df = hour_df, BAdata = DW1_rmna)
  cv_r1 <- data.table(data.frame(cv))
  cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
  avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
  yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])
  optimal_sp[i, ] <- c(hour_df, yhour_df)
}

# stop cluster
stopCluster(cl)

results_sp

saveRDS(results_sp, "splines_opt_w.RDS")   


#================================================#
# run again to train without climate interaction #
#================================================#

detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

optimal_sp <- matrix(NA, 15, 2)

results_sp <- foreach(i=1:15, .packages = c("data.table", "splines", "dplyr")) %dopar% {
  DW1_rmna <- DWw15_clim[sf_id == i, ] 
  # tune the first parameter: hour of the day
  cv <- tune_noclim1(BAdata = DW1_rmna)
  cv_r <- data.table(data.frame(cv))
  cv_rmse_r <- cv_r[, list(X2,X5,X8)]
  avg_rmse_r <- data.frame(knots = cv_r[,1], avgcv = rowMeans(cv_rmse_r))
  hour_df <- as.numeric(avg_rmse_r[which.min(avg_rmse_r$avgcv),][1])
  
  # tune the second parameter
  cv <- tune_noclim2(hour_df = hour_df, BAdata = DW1_rmna)
  cv_r1 <- data.table(data.frame(cv))
  cv_rmse_r1 <- cv_r1[, list(X2,X5,X8)]
  avg_rmse_r1 <- data.frame(knots = cv_r1[,1], avgcv = rowMeans(cv_rmse_r1))
  yhour_df <- as.numeric(avg_rmse_r1[which.min(avg_rmse_r1$avgcv),][1])
  optimal_sp[i, ] <- c(hour_df, yhour_df)
}

# stop cluster
stopCluster(cl)

results_sp

saveRDS(results_sp, "splines_noclim_w.RDS")   


