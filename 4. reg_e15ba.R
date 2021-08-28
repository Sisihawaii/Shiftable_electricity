## 5. regressions for all BAs in the east
# May 27, 2020
# train splines results from Kaholo trainBAs
# load result from splines_noclim_e.RDS and splines_opt_e.RDS 
# output flexload and regression coefficients 

# update Jun 6, 2021: use new clean data results as input
# update Jun 17, 2021: added BA1 east IC CV R2 at the end
# update Aug 3, 2021: rerun everything with noclimate interaction opt splines

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
splines_opt <- readRDS("splines_opt_e.RDS")
# load demand climate data for east 15 BAs
DWe15_clim <- readRDS("DWe15_clim12_newdt.RDS")
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
  # e <- resid(fit_final)
  # mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
  y_hat1 <- exp(logy_hat1)
  y_hat0 <- exp(logy_hat0)
  s <- (y_hat1 - y_hat0)/y_hat1
  DW1_rmna$yhat0 <- y_hat0
  DW1_rmna$yhat1 <- y_hat1
  DW1_rmna$flex <- s

  # save DW1_rmna with flex load for each BA
  DW1_rmna_all[[i]] <- DW1_rmna
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg, "flexload_BAe15_new.RDS")

#=====================================#
# results without climate interaction #
#=====================================#
# load optimal splines
splines_opt1 <- readRDS("splines_noclim_e.RDS")

# run in parallel
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
  # e <- resid(fit_final)
  # mean_e <- mean(exp(e))   # same error term for both yhat and yhat0
  y_hat1 <- exp(logy_hat1)
  y_hat0 <- exp(logy_hat0)
  s <- (y_hat1 - y_hat0)/y_hat1
  DW1_rmna$yhat0 <- y_hat0
  DW1_rmna$yhat1 <- y_hat1
  DW1_rmna$flex <- s
  
  # save DW1_rmna with flex load for each BA
  DW1_rmna_all[[i]] <- DW1_rmna
}

# stop cluster
stopCluster(cl)
saveRDS(results_reg, "flexload_BAe15_noclim.RDS")

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
saveRDS(results_reg1, "regcoef_BAe15_new.RDS")


#==============================================#
# run again for no climate interaction splines #
#==============================================#
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
saveRDS(results_reg1, "regcoef_BAe15_new1.RDS")   # results with no climate opt splines


###################################
#### 3 fold cv for ba1 east IC ####
###################################

splines_opt <- readRDS("splines_noclim_e.RDS")       # use no climate opt splines for all models
DWe15_clim <- readRDS("DWe15_clim12_newdt.RDS")
DWe15_clim <- data.table(DWe15_clim)
DW_rmna <- DWe15_clim[sf_id == 1, ]

hour_df <- splines_opt[[1]][1]
yhour_df <- splines_opt[[1]][2]

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
  trainR2[k,] <- summary(fit_cv)$adj.r.squared
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
  trainR2[k,] <- summary(fit_cv)$adj.r.squared
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
  trainR2[k,] <- summary(fit_cv)$adj.r.squared
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

