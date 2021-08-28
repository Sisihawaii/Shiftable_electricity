## change the model without interactions
# sep 29, 2020
# remove the interactions for BAs that have negative flex load
# rerun both the model and prediction with no clim interaction trained splines
# Jun 9, 2021 update with clean data

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

# load east and west data separately
# east
flexload <- readRDS("predicte15_2C_new.RDS")              # neg: ba6 ba13
summary(flexload[[1]])
DWe15_clim <- readRDS("DWe15_clim12_2C_new.RDS")          # east
DWe15_clim_new <- readRDS("DWe15_clim12_newdt.RDS")     # east
splines_opt <- readRDS("splines_noclim_e.RDS")        # east
ba <- c(6, 13)         # east

#=================================================
# or west
flexload <- readRDS("predictw15_2C_new.RDS")              # neg: ba1 ba3 ba8 ba12 
DWe15_clim <- readRDS("DWw15_clim12_2C_new.RDS")          # west
DWe15_clim_new <- readRDS("DWw15_clim12_newdt.RDS")     # west
splines_opt <- readRDS("splines_noclim_w.RDS")        # west
ba <- c(1, 3, 8, 12)   # west

#=================================================

DWe15_clim <- data.table(DWe15_clim)
DWe15_clim_new <- data.table(DWe15_clim_new)

# run in parallel 
# n <- length(ba)
n=15
detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

DW1_all <- list()
results_reg <- foreach(j=1:n, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  #i <- ba[j]
  i <- j
  DW1_rmna <- DWe15_clim_new[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                    hdh_18 + cdh_18, 
                  data = DW1_rmna)
  # predict demand (including NAs) for each BA
  DW1       <- DWe15_clim[sf_id == i, ] %>% 
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
  
  # save DW1_rmna with flex load for each BA
  DW1_all[[i]] <- DW1
}

# stop cluster
stopCluster(cl)

summary(results_reg[[12]])   

# for (i in 1:n){
#   flexload[[ba[i]]] <- results_reg[[i]]
# }

# saveRDS(flexload, "Sisi/flexload_east2.RDS")
# saveRDS(flexload, "Sisi/flexload_west4.RDS")
# all flex are positive now for both the east and west

# tried again with removing all interactions
saveRDS(results_reg, "flexload_e_2c_new.RDS")
saveRDS(results_reg, "flexload_w_2c_new.RDS")


