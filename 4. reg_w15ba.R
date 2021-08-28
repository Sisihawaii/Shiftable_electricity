## 5. regressions for all BAs in the west
# July 8, 2020
# train splines results from Kaholo trainBAs
# load result from splines_opt_w.RDS and splines_noclim_w.RDS
# output flexload and regression coefficients

# update Jun 6, 2021: use cleandata results as input
# update Aug 3, 2021: rerun for no climate splines

rm(list = ls())

library(tidyverse)
library(lubridate)
library(data.table)
library(splines)
library(foreach)
library(doParallel)
library(RColorBrewer)
library(gridExtra)

# load optimal splines
splines_opt <- readRDS("splines_opt_w.RDS")
# load demand climate data for west 15 BAs
DWe15_clim <- readRDS("DWw15_clim12_newdt.RDS")
summary(DWe15_clim)
DWe15_clim <- data.table(DWe15_clim)

# run in parallel
detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

DW1_rmna_all <- list()
results_reg <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                  hdh_18 + climH + cdh_18 + climC, 
                data = DW1_rmna)
  # calculate flex load for each BA
  DW1_0 <- DW1_rmna %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0) %>% mutate(climH = 0) %>% mutate(climC = 0)
  logy_hat0 <- predict(fit_final, DW1_0)
  logy_hat1 <- predict(fit_final, DW1_rmna)
  e <- resid(fit_final)
  mean_e <- mean(exp(e))             # same error term for both yhat and yhat0
  y_hat1 <- exp(logy_hat1)*mean_e    # adjustment for delog
  y_hat0 <- exp(logy_hat0)*mean_e
  s <- (y_hat1 - y_hat0)/y_hat1
  DW1_rmna$yhat0 <- y_hat0
  DW1_rmna$yhat1 <- y_hat1
  DW1_rmna$flex <- s
  
  # save DW1_rmna with flex load for each BA
  DW1_rmna_all[[i]] <- DW1_rmna
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg, "flexload_BAw15_new.RDS")

#=====================================#
# results without climate interaction #
#=====================================#

# load optimal splines
splines_opt1 <- readRDS("splines_noclim_w.RDS")

detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

DW1_rmna_all <- list()
results_reg <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt1[[i]][1]
  yhour_df <- splines_opt1[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                    hdh_18 + cdh_18, 
                  data = DW1_rmna)
  # calculate flex load for each BA
  DW1_0 <- DW1_rmna %>% mutate(hdh_18 = 0) %>% mutate(cdh_18 = 0)
  logy_hat0 <- predict(fit_final, DW1_0)
  logy_hat1 <- predict(fit_final, DW1_rmna)
  e <- resid(fit_final)
  mean_e <- mean(exp(e))             # same error term for both yhat and yhat0
  y_hat1 <- exp(logy_hat1)*mean_e    # adjustment for delog
  y_hat0 <- exp(logy_hat0)*mean_e
  s <- (y_hat1 - y_hat0)/y_hat1
  DW1_rmna$yhat0 <- y_hat0
  DW1_rmna$yhat1 <- y_hat1
  DW1_rmna$flex <- s
  
  # save DW1_rmna with flex load for each BA
  DW1_rmna_all[[i]] <- DW1_rmna
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg, "flexload_BAw15_noclim.RDS")

#=====================================================

# run in parallel again for coefficients
cl <- makeCluster(4)
registerDoParallel(cl)

beta_w <- matrix(NA, 15, 9)
results_reg1 <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                    hdh_18 + climH + cdh_18 + climC, 
                  data = DW1_rmna)
  
  # save the regression coefficients for weather varialbes
  beta0_c <- coef(fit_final)["cdh_18"]
  beta0_c_std <- coef(summary(fit_final))['cdh_18', 'Std. Error']
  beta1_c <- coef(fit_final)["climC"]
  beta1_c_std <- coef(summary(fit_final))['climC', 'Std. Error']
  beta0_h <- coef(fit_final)["hdh_18"]
  beta0_h_std <- coef(summary(fit_final))['hdh_18', 'Std. Error']
  beta1_h <- coef(fit_final)["climH"]
  beta1_h_std <- coef(summary(fit_final))['climH', 'Std. Error']
  r2adj <- summary(fit_final)$adj.r.squared
  beta_w[i, ] <- c(beta0_c, beta0_c_std, beta1_c, beta1_c_std, beta0_h, beta0_h_std, beta1_h, beta1_h_std, r2adj)
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg1, "regcoef_BAw15_new.RDS")

#===============================#
# rerun with no climate splines #
#===============================#

# run in parallel again for coefficients
cl <- makeCluster(4)
registerDoParallel(cl)

beta_w <- matrix(NA, 15, 9)
results_reg1 <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt1[[i]][1]
  yhour_df <- splines_opt1[[i]][2]
  fit_final <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) + 
                    hdh_18 + climH + cdh_18 + climC, 
                  data = DW1_rmna)
  
  # save the regression coefficients for weather varialbes
  beta0_c <- coef(fit_final)["cdh_18"]
  beta0_c_std <- coef(summary(fit_final))['cdh_18', 'Std. Error']
  beta1_c <- coef(fit_final)["climC"]
  beta1_c_std <- coef(summary(fit_final))['climC', 'Std. Error']
  beta0_h <- coef(fit_final)["hdh_18"]
  beta0_h_std <- coef(summary(fit_final))['hdh_18', 'Std. Error']
  beta1_h <- coef(fit_final)["climH"]
  beta1_h_std <- coef(summary(fit_final))['climH', 'Std. Error']
  r2adj <- summary(fit_final)$adj.r.squared
  beta_w[i, ] <- c(beta0_c, beta0_c_std, beta1_c, beta1_c_std, beta0_h, beta0_h_std, beta1_h, beta1_h_std, r2adj)
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg1, "regcoef_BAw15_new1.RDS")



