# compute BIC for AR model objects of class 'dynlm'
BIC <- function(model) {
  
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  npar <- length(model$coef)
  
  return(
    round(c("p" = npar - 1,
            "BIC" = log(ssr/t) + npar * log(t)/t,
            "R2" = summary(model)$r.squared), 4)
  )
}

######################################
# function to calculate RMSE and R2 ##
######################################
RMSE_R2 <- function(model, testdata, train){
  p <- predict(model, testdata)
  res_p <- testdata$demand - p
  RMSE_p <- sqrt(mean(res_p^2))
  SSR_p <- sum(res_p^2)
  SST_p <- sum((testdata$demand - mean(train$demand))^2)
  R2_p <- 1 - SSR_p/SST_p 
  rmse_r2 <- c(RMSE_p, R2_p)
  return(rmse_r2)
}

RMSE_R2_log <- function(model, testdata, train){
  p <- predict(model, testdata)
  res_p <- log(testdata$demand) - p
  RMSE_p <- sqrt(mean(res_p^2))
  SSR_p <- sum(res_p^2)
  SST_p <- sum((log(testdata$demand) - mean(log(train$demand)))^2)
  R2_p <- 1 - SSR_p/SST_p
  rmse_r2 <- c(RMSE_p, R2_p)
  return(rmse_r2)
}

RMSE_R2_log_adj <- function(model, testdata, train){
  p <- predict(model, testdata)
  res_p <- log(testdata$demand) - p
  RMSE_p <- sqrt(mean(res_p^2))
  SSR_p <- sum(res_p^2)
  SST_p <- sum((log(testdata$demand) - mean(log(train$demand)))^2)
  # R2_p <- 1 - SSR_p/SST_p
  R2_p <- (cor(log(testdata$demand), p))^2
  R2_adj <- var(exp(p))/var(testdata$demand)    # adjust to compare with linear R2
  rmse_r2 <- c(RMSE_p, R2_p, R2_adj)
  return(rmse_r2)
}

RMSE_R2_std <- function(model, testdata, train){
  p <- predict(model, testdata)
  res_p <- testdata$std_logD - p
  RMSE_p <- sqrt(mean(res_p^2))
  SSR_p <- sum(res_p^2)
  SST_p <- sum((testdata$std_logD - mean(train$std_logD))^2)
  R2_p <- 1 - SSR_p/SST_p
  rmse_r2 <- c(RMSE_p, R2_p)
  return(rmse_r2)
}


##################################################
# function for parameter tuning train_splines.R ##
##################################################


#-------------------------------------------------------------------------
# use log(demand) rmse function and 3 fold - one fold per year
#-------------------------------------------------------------------------
tune_new1 <- function(BAdata, yhour_df=5){
  cv <- matrix(NA, 17, 9)
  for (k in 1:3){
    train <- BAdata[year!=(2015+k), ]
    valid <- BAdata[year==(2015+k), ]
    for (i in 4:20){
      fit <- lm(log(demand) ~ trend1 + ns(hour, df = i)*ns(yhour, df = yhour_df)*factor(weekdaynum) +
                   hdh_18 + climH + cdh_18 + climC, 
                 data = train)
      cv[(i-3), ((k-1)*3+1):(k*3)] <- c(i, RMSE_R2_log(fit, valid, train))
    }
  }
  return(cv)
}

tune_new2 <- function(hour_df=5, BAdata){
  cv <- matrix(NA, 17, 9)
  for (k in 1:3){
    train <- BAdata[year!=(2015+k), ]
    valid <- BAdata[year==(2015+k), ]
    for (i in 4:20){
      fit <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = i)*factor(weekdaynum) + 
                  hdh_18 + climH + cdh_18 + climC, 
                data = train)
      cv[(i-3), ((k-1)*3+1):(k*3)] <- c(i, RMSE_R2_log(fit, valid, train))
    }
  }
  return(cv)
}

#==================================#
# tune without climate interaction #
#==================================#
tune_noclim1 <- function(BAdata, yhour_df=5){
  cv <- matrix(NA, 17, 9)
  for (k in 1:3){
    train <- BAdata[year!=(2015+k), ]
    valid <- BAdata[year==(2015+k), ]
    for (i in 4:20){
      fit <- lm(log(demand) ~ trend1 + ns(hour, df = i)*ns(yhour, df = yhour_df)*factor(weekdaynum) +
                  hdh_18 + cdh_18, 
                data = train)
      cv[(i-3), ((k-1)*3+1):(k*3)] <- c(i, RMSE_R2_log(fit, valid, train))
    }
  }
  return(cv)
}

tune_noclim2 <- function(hour_df=5, BAdata){
  cv <- matrix(NA, 17, 9)
  for (k in 1:3){
    train <- BAdata[year!=(2015+k), ]
    valid <- BAdata[year==(2015+k), ]
    for (i in 4:20){
      fit <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = i)*factor(weekdaynum) + 
                  hdh_18 + cdh_18, 
                data = train)
      cv[(i-3), ((k-1)*3+1):(k*3)] <- c(i, RMSE_R2_log(fit, valid, train))
    }
  }
  return(cv)
}



#######################################################
###### function for clim_map.R script #################
#######################################################
# matching raster cells for electricity load prediction and flexible shares

raster.flex <- function(pop_temp, flex_ba, celln_ba){
  # matching raster cells for electricity load prediction and flexible shares for each hour
  # pop_temp: pop raster projected on temp 
  # flex_ba: dataset for selected ba and period - hourly per cell data
  # celln_ba: unique cell values for selected ba
  # returns rasterstack with 24 layers corresponding to 24 hours
  
  # setup the initial rasters
  yhat0_ba <- raster(pop_temp)
  yhat1_ba <- raster(pop_temp)
  s_ba <- raster(pop_temp)
  
  yhat0_ms <- raster(pop_temp)                                     # set an empty raster for stacking in the loop
  yhat1_ms <- raster(pop_temp)
  s_ms <- raster(pop_temp)
  for (hi in 0:23){
    flex_ba_h <- flex_ba[hour == hi, ]                             # for each hour, we will have a raster layer
    for (i in 1:length(celln_ba)){
      yhat0_ba[celln_ba[i]] <- as.matrix(flex_ba_h[i, 3])          # assign each cell with values
      yhat1_ba[celln_ba[i]] <- as.matrix(flex_ba_h[i, 4])  
      s_ba[celln_ba[i]] <- as.matrix(flex_ba_h[i, 5])  
    }
    yhat0_ms <- stack(yhat0_ms, yhat0_ba)                       # stack layers from each loop
    yhat1_ms <- stack(yhat1_ms, yhat1_ba)
    s_ms <- stack(s_ms, s_ba)
  }
  return(list(yhat0_ms, yhat1_ms, s_ms))
}






