# regression coefficient visualization
# feb 8, 2021
# run baseline model for all BAs save the RMSE
# run model with just weather variables for all BAs save the coef sd and RMSE
# calc grid cell cdh/hdh range in full model
# organize the data to fit the visualization function
# apply the schart function for visualization
# save all the variable results in separate folder for east and west IC (run saveRDS only once)

# update April 23: use cleandata rerun for east and west IC
# update Aug 3: use new opt splines from cleandata for all the models rerun results (with climate interaction)

rm(list = ls())

library(tidyverse)
library(lubridate)
library(data.table)
library(splines)
library(foreach)
library(doParallel)
library(sandwich)
library(lmtest)

# NW sd function
nwsd.fn <- function(model){
  TT <- 24   # tried diff values - result did not change in this case
  m <- 0.75 * (length(TT))^(1/3)
  nw_sd <- NeweyWest(model,
                     lag = m - 1, prewhite = F,
                     adjust = T)
  return(nw_sd)
}


# new data set east
# splines_opt <- readRDS("splines_opt_e.RDS")
# DWe15_clim <- readRDS("DWe15_clim12_newdt.RDS")
# regcoef_e15 <- readRDS("regcoef_BAe15_new.RDS")   # used the coef, rerun the model below, calc NW sd

# new data set west
splines_opt <- readRDS("splines_opt_w.RDS")
DWe15_clim <- readRDS("DWw15_clim12_newdt.RDS")
regcoef_e15 <- readRDS("regcoef_BAw15_new.RDS")  

summary(DWe15_clim)
DWe15_clim <- data.table(DWe15_clim)


###########################################################
#### run baseline and with weather variable regression ####
###########################################################

# run for baseline model and save the entired model
cl <- makeCluster(4)
registerDoParallel(cl)

fit_base <- list()
results_reg1 <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_base[[i]] <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum),
                    data = DW1_rmna)
}

# stop cluster
stopCluster(cl)

# calc rmse for each BA baseline
rmse_b <- matrix(NA, 15, 2)
for (i in 1:15){
  rmse_b[i, 1] <- i
  s1 <- summary(results_reg1[[i]])
  rmse_b[i, 2] <- sqrt(mean(s1$residuals^2))
}

colnames(rmse_b) <- c("BA", "rmse_b")

# saveRDS(rmse_b, "Regvis/East15ba/baseline_rmse.RDS")  
saveRDS(rmse_b, "Regvis/West15ba/baseline_rmse.RDS")  

##===============================================##

# remove previous model result (due to out of memory) and run weather model
rm(results_reg1)

cl <- makeCluster(4)
registerDoParallel(cl)

fit_w <- list()
results_reg2 <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_w[[i]] <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) +
                     hdh_18 + cdh_18,
                   data = DW1_rmna)
}

# stop cluster
stopCluster(cl)

# calc rmse and extract reg coef for CDH/HDH and their Robust SD
rmse_w <- matrix(NA, 15, 2)
coef_c <- matrix(NA, 15, 6)
coef_h <- matrix(NA, 15, 6)

for (i in 1:15){
  fit <- results_reg2[[i]]
  # calc rmse
  rmse_w[i, 1] <- i
  s2 <- summary(fit)
  rmse_w[i, 2] <- sqrt(mean(s2$residuals^2))
  # coef and nw sd
  beta0_c <- coef(fit)["cdh_18"]
  beta0_h <- coef(fit)["hdh_18"]
  nwsd <- nwsd.fn(fit)
  c0_sd <- coeftest(fit, vcov = nwsd)['cdh_18', 'Std. Error']
  h0_sd <- coeftest(fit, vcov = nwsd)['hdh_18', 'Std. Error']
  coef_c[i, ] <- c(beta0_c, c0_sd, 1, 1, 0, 0)
  coef_h[i, ] <- c(beta0_h, h0_sd, 1, 1, 0, 0)
}

colnames(coef_c) <- c("coef", "se", "cdh", "hdh", "climC_cdh", "climH_hdh")
colnames(coef_h) <- c("coef", "se", "cdh", "hdh", "climC_cdh", "climH_hdh")
colnames(rmse_w) <- c("BA", "rmse_w")

# saveRDS(rmse_w, "Regvis/East15ba/weather_rmse.RDS")  
# saveRDS(coef_c, "Regvis/East15ba/weathercoef_cdh.RDS")  
# saveRDS(coef_h, "Regvis/East15ba/weathercoef_hdh.RDS")  

saveRDS(rmse_w, "Regvis/West15ba/weather_rmse.RDS")  
saveRDS(coef_c, "Regvis/West15ba/weathercoef_cdh.RDS")  
saveRDS(coef_h, "Regvis/West15ba/weathercoef_hdh.RDS")  

#====================================================#

rm(results_reg2)
# run for weather clim interaction model and save the entired model again 
cl <- makeCluster(4)
registerDoParallel(cl)

fit_f <- list()
results_reg3 <- foreach(i=1:15, .packages = c("data.table", "splines", "tidyverse")) %dopar% {
  DW1_rmna <- DWe15_clim[sf_id == i, ]                      # select only one BA regoin
  # choose the opt splines for each BA
  hour_df <- splines_opt[[i]][1]
  yhour_df <- splines_opt[[i]][2]
  fit_f[[i]] <- lm(log(demand) ~ trend1 + ns(hour, df = hour_df)*ns(yhour, df = yhour_df)*factor(weekdaynum) +
                     hdh_18 + climH + cdh_18 + climC,
                   data = DW1_rmna)
}

# stop cluster
stopCluster(cl)

# save once - takes long to save
# saveRDS(results_reg3, "Regvis/East15ba/modelresults.RDS")
saveRDS(results_reg3, "Regvis/West15ba/modelresults.RDS")

# calc rmse and extract reg coef for CDH/HDH and their Robust SD
# need to remove all variables before and just run the first 50 lines and this part due to my pc capacity
# results_reg3 <- readRDS("Regvis/East15ba/modelresults.RDS")
results_reg3 <- readRDS("Regvis/West15ba/modelresults.RDS")

rmse_f <- matrix(NA, 15, 2)
coef_c1 <- matrix(NA, 15, 5)
coef_h1 <- matrix(NA, 15, 5)

for (i in 1:15){
  fit <- results_reg3[[i]]
  # calc rmse
  rmse_f[i, 1] <- i
  s3 <- summary(fit)
  rmse_f[i, 2] <- sqrt(mean(s3$residuals^2))
  # coef and nw sd
  beta0_c <- coef(fit)["cdh_18"]
  beta0_h <- coef(fit)["hdh_18"]
  beta1_c <- coef(fit)["climC"]
  beta1_h <- coef(fit)["climH"]
  nwsd <- nwsd.fn(fit)
  c0_sd <- coeftest(fit, vcov = nwsd)['cdh_18', 'Std. Error']
  h0_sd <- coeftest(fit, vcov = nwsd)['hdh_18', 'Std. Error']
  c1_sd <- coeftest(fit, vcov = nwsd)['climC', 'Std. Error']
  h1_sd <- coeftest(fit, vcov = nwsd)['climH', 'Std. Error']
  c0c1_cov <- nwsd['cdh_18', 'climC']
  h0h1_cov <- nwsd['hdh_18', 'climH']
  coef_c1[i, ] <- c(beta0_c, c0_sd, beta1_c, c1_sd, c0c1_cov)
  coef_h1[i, ] <- c(beta0_h, h0_sd, beta1_h, h1_sd, h0h1_cov)
}

colnames(coef_c1) <- c("coef0", "se0", "coef1", "se1", "cov")
colnames(coef_h1) <- c("coef0", "se0", "coef1", "se1", "cov")
colnames(rmse_f) <- c("BA", "rmse_f")

coef_c1 <- data.frame(coef_c1)
coef_h1 <- data.frame(coef_h1)
rmse_f <- data.frame(rmse_f)

# saveRDS(rmse_f, "Regvis/East15ba/model_rmse.RDS")  
# saveRDS(coef_c1, "Regvis/East15ba/modelcoef_cdh.RDS")  
# saveRDS(coef_h1, "Regvis/East15ba/modelcoef_hdh.RDS")  

saveRDS(rmse_f, "Regvis/West15ba/model_rmse.RDS")  
saveRDS(coef_c1, "Regvis/West15ba/modelcoef_cdh.RDS")  
saveRDS(coef_h1, "Regvis/West15ba/modelcoef_hdh.RDS")  

##############################################
### calc cdh/hdh range with all grid cells ###
##############################################

# load clim per cell data and matching cells to BA regoin with pop per cell data
clim <- readRDS("clim_cell12.RDS")
pop.dt <- readRDS("popweight.RDS")

# pop_e15 <- pop.dt[sf == 'e.15']
pop_e15 <- pop.dt[sf == 'w.15']
climcell_ba <- clim[pop_e15, on = 'cell']

## calc pop weighted avg for cdh/hdh
clim_avg <- climcell_ba %>% group_by(ID) %>% 
  summarise(climC_avg = sum(climC*w), climH_avg = sum(climH*w))

# calc coef for cdh/hdh for climate interaction model
beta_i_list <- list()
for (i in 1:15){
  ba1 <- climcell_ba[ID == i]
  b0_c <- regcoef_e15[[i]]['cdh_18']
  b1_c <- regcoef_e15[[i]]['climC']
  b0_h <- regcoef_e15[[i]]['hdh_18']
  b1_h <- regcoef_e15[[i]]['climH']
  beta_ci_max <- b0_c + b1_c*max(ba1$climC)
  beta_ci_min <- b0_c + b1_c*min(ba1$climC)
  beta_ci_avg <- b0_c + b1_c*clim_avg$climC_avg[i]
  beta_hi_max <- b0_h + b1_h*max(ba1$climH)
  beta_hi_min <- b0_h + b1_h*min(ba1$climH)
  beta_hi_avg <- b0_h + b1_h*clim_avg$climH_avg[i]
  # combine range and average 
  beta_i <- matrix(0, 2, 3)
  beta_i[1,] <- cbind(beta_ci_min, beta_ci_max, beta_ci_avg)
  beta_i[2,] <- cbind(beta_hi_min, beta_hi_max, beta_hi_avg)
  beta_i_list[[i]] <- cbind(beta_i, i)   # first row cdh, second row hdh
}

beta_e15 <- do.call(rbind, beta_i_list)
beta_e15 <- data.frame(beta_e15)
# rearrange min/max
beta_e15$min <- apply(beta_e15[,1:2], 1, min)
beta_e15$max <- apply(beta_e15[,1:2], 1, max)
beta_e15$avg <- beta_e15$V3
beta_e15 <- beta_e15[,4:7]

# saveRDS(beta_e15, "Regvis/East15ba/coef_climdh_range.RDS")
saveRDS(beta_e15, "Regvis/West15ba/coef_climdh_range.RDS")

# calc se for cdh/hdh in clim interaction model for each grid cell i
# var(beta_i) = var(beta0) + climCi^2*var(beta1) + 2*climCi*cov(beta0,beta1)

# load model coefs and se
# coef_c1 <- readRDS("Regvis/East15ba/modelcoef_cdh.RDS")
# coef_h1 <- readRDS("Regvis/East15ba/modelcoef_hdh.RDS")

coef_c1 <- readRDS("Regvis/West15ba/modelcoef_cdh.RDS")
coef_h1 <- readRDS("Regvis/West15ba/modelcoef_hdh.RDS")

# calc beta_i for all BA15 - cdh
sd_range_c <- matrix(NA, 15, 2)
for (i in 1:15){
  b0_c <- coef_c1$coef0[i]
  b1_c <- coef_c1$coef1[i]
  se0_c <- coef_c1$se0[i]
  se1_c <- coef_c1$se1[i]
  ba1 <- climcell_ba[ID == i]
  climCi <- ba1$climC
  cov_c <- coef_c1$cov[i]
  
  beta_ci <- b0_c + b1_c*climCi
  var_bi <- (se0_c)^2 + (se1_c)^2*(climCi)^2 + 2*climCi*cov_c
  sd_bi <- sqrt(var_bi)
  
  # combine coef and se and calc ci
  df_c <- cbind.data.frame(beta_ci, sd_bi)
  ci <- 0.95
  a <- qnorm(1-(1-ci)/2)
  df_c$l1 <- df_c$beta_ci - a*df_c$sd_bi
  df_c$h1 <- df_c$beta_ci + a*df_c$sd_bi
  sd_range_c[i, ] <- cbind(min(df_c$l1), max(df_c$h1))
}

# do the same for hdh
sd_range_h <- matrix(NA, 15, 2)
for (i in 1:15){
  b0_c <- coef_h1$coef0[i]
  b1_c <- coef_h1$coef1[i]
  se0_c <- coef_h1$se0[i]
  se1_c <- coef_h1$se1[i]
  ba1 <- climcell_ba[ID == i]
  climCi <- ba1$climH
  cov_c <- coef_h1$cov[i]
  
  beta_ci <- b0_c + b1_c*climCi
  var_bi <- (se0_c)^2 + (se1_c)^2*(climCi)^2 + 2*climCi*cov_c
  sd_bi <- sqrt(var_bi)
  
  # combine coef and se and calc ci
  df_c <- cbind.data.frame(beta_ci, sd_bi)
  ci <- 0.95
  a <- qnorm(1-(1-ci)/2)
  df_c$l1 <- df_c$beta_ci - a*df_c$sd_bi
  df_c$h1 <- df_c$beta_ci + a*df_c$sd_bi
  sd_range_h[i, ] <- cbind(min(df_c$l1), max(df_c$h1))
}

# adjust cdh and hdh sd range 
sd_range_c <- data.frame(sd_range_c)
sd_range_h <- data.frame(sd_range_h)
colnames(sd_range_c) <- c("l2", "h2")
colnames(sd_range_h) <- c("l2", "h2")

###########################################
######### applying schart function ########
###########################################

# load downloaded function
source("spec_chart_function.R")   

# coef_c <- readRDS("Regvis/East15ba/weathercoef_cdh.RDS")
# coef_h <- readRDS("Regvis/East15ba/weathercoef_hdh.RDS")

coef_c <- readRDS("Regvis/West15ba/weathercoef_cdh.RDS")
coef_h <- readRDS("Regvis/West15ba/weathercoef_hdh.RDS")

df_c <- data.frame(coef_c)
df_c$cdh <- TRUE
df_c$hdh <- TRUE
df_c$climC_cdh <- FALSE
df_c$climH_hdh <- FALSE

df_h <- data.frame(coef_h)
df_h$cdh <- TRUE
df_h$hdh <- TRUE
df_h$climC_cdh <- FALSE
df_h$climH_hdh <- FALSE

# combine cdh and hdh
ba_e15 <- list()
for (i in 1:15){
  ba <- rbind(df_c[i, ], df_h[i,])   # first row cdh, second row hdh
  ba$ba <- i
  ba_e15[[i]] <- ba
}

data_ba <- do.call(rbind, ba_e15)
rownames(data_ba) <- 1:30

# calc ci using se - create new variables l1, h1
data1 <- data_ba[, 1:6]
data1$se <- NULL
ci <- .95
a <- qnorm(1-(1-ci)/2)
est <- data_ba$coef
se <- data_ba$se
data1$l1 <- est - a*se
data1$h1 <- est + a*se

#=================================#

### adjust data format with clim interaction cdh/hdh range as se
# beta_e15 <- readRDS("Regvis/East15ba/coef_climdh_range.RDS")
beta_e15 <- readRDS("Regvis/West15ba/coef_climdh_range.RDS")

data2 <- beta_e15[,2:4]
data2[,4:7] <- data_ba[, 3:6]
data2$climC_cdh <- TRUE
data2$climH_hdh <- TRUE
data2$min <- NULL
data2$max <- NULL
data2$l1 <- beta_e15$min
data2$h1 <- beta_e15$max
colnames(data2)[1] <- "coef"

# combine weather coef and interaction coef
data_list <- list()
for (i in 1:15){
  data_cb <- data1[(2*i-1):(2*i),]
  data_cb1 <- data2[(2*i-1):(2*i),]
  data_list[[i]] <- rbind(data_cb, data_cb1)
}

data_all <- do.call(rbind, data_list)
rownames(data_all) <- 1:60

# index.ci <- match(c("l1","h1"), names(data_all))

# combine sd range as l2 h2
data_all$l2 <- data_all$l1
data_all$h2 <- data_all$l2

for (i in 1:15){
  data_all$l2[4*i-1] <- sd_range_c$l2[i]
  data_all$h2[4*i-1] <- sd_range_c$h2[i]
  data_all$l2[4*i] <- sd_range_h$l2[i]
  data_all$h2[4*i] <- sd_range_h$h2[i]
}


index.ci <- match(c("l1","h1","l2","h2"), names(data_all))

### calc daily avg demand per region
DD <- DWe15_clim %>% group_by(sf_id, year, month, day) %>% summarise(dailydemand = sum(demand))
DD_avg <- DD %>% group_by(sf_id) %>% summarise(daily_avg = mean(dailydemand))

### add RMSE changes
# rmse_b <- readRDS("Regvis/East15ba/baseline_rmse.RDS")
# rmse_w <- readRDS("Regvis/East15ba/weather_rmse.RDS")
# rmse_f <- readRDS("Regvis/East15ba/model_rmse.RDS")

rmse_b <- readRDS("Regvis/West15ba/baseline_rmse.RDS")
rmse_w <- readRDS("Regvis/West15ba/weather_rmse.RDS")
rmse_f <- readRDS("Regvis/West15ba/model_rmse.RDS")

rmse_b_df <- data.frame(rmse_b)
rmse_w_df <- data.frame(rmse_w)
rmse_df <- merge(rmse_b_df, rmse_w_df)
rmse_df1 <- merge(rmse_df, rmse_f)
rmse_df1$diff <- (rmse_df1$rmse_b - rmse_df1$rmse_w)/rmse_df1$rmse_b*100
rmse_df1$diff1 <- (rmse_df1$rmse_b - rmse_df1$rmse_f)/rmse_df1$rmse_b*100

# add zero for the same model cdh rows
rmse_reduction <- matrix(0, 75, 1)
for (i in 1:15){
  rmse_reduction[2+5*(i-1)] <- rmse_df1$diff[i]
  rmse_reduction[3+5*(i-1)] <- rmse_df1$diff1[i]
}

#=======================================#
# apply schart for e15 BAs
# schart(df_c, leftmargin = 5, heights = c(4,1), lwd.symbol = 2, pch.est = 21, lwd.ref = 1, 
#        adj = c(0.5,0.5), offset = c(3, 2))

# ylim <- c(-0.15,0.2)   # east
ylim <- c(-0.3,0.3)
schart(data_all, highlight = c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60),
       leftmargin = 5, heights = c(4,1), n = 4, ylim = ylim, index.ci = index.ci, axes = F)
abline(v = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70))
# y <- 0.2   # east
y <- 0.3   # west
text(2,y,"BA1")
text(7.5,y,"BA2")
text(12.5,y,"BA3")
text(17.5,y,"BA4")
text(22.5,y,"BA5")
text(27.5,y,"BA6")
text(32.5,y,"BA7")
text(37.5,y,"BA8")
text(42.5,y,"BA9")
text(47.5,y,"BA10")
text(52.5,y,"BA11")
text(57.5,y,"BA12")
text(62.5,y,"BA13")
text(67.5,y,"BA14")
text(73,y,"BA15")
# add avg daily MWh in each region
# y1 <- 0.19 # east
y1 <- 0.27 # west
text(2,y1, paste0(round(DD_avg[1,2]/1000), "GWh"), cex = 0.7)
text(7.5,y1, paste0(round(DD_avg[2,2]/1000), "GWh"), cex = 0.7)
text(12.5,y1, paste0(round(DD_avg[3,2]/1000), "GWh"), cex = 0.7)
text(17.5,y1, paste0(round(DD_avg[4,2]/1000), "GWh"), cex = 0.7)
text(22.5,y1, paste0(round(DD_avg[5,2]/1000), "GWh"), cex = 0.7)
text(27.5,y1, paste0(round(DD_avg[6,2]/1000), "GWh"), cex = 0.7)
text(32.5,y1, paste0(round(DD_avg[7,2]/1000), "GWh"), cex = 0.7)
text(37.5,y1, paste0(round(DD_avg[8,2]/1000), "GWh"), cex = 0.7)
text(42.5,y1, paste0(round(DD_avg[9,2]/1000), "GWh"), cex = 0.7)
text(47.5,y1, paste0(round(DD_avg[10,2]/1000), "GWh"), cex = 0.7)
text(52.5,y1, paste0(round(DD_avg[11,2]/1000), "GWh"), cex = 0.7)
text(57.5,y1, paste0(round(DD_avg[12,2]/1000), "GWh"), cex = 0.7)
text(62.5,y1, paste0(round(DD_avg[13,2]/1000), "GWh"), cex = 0.7)
text(67.5,y1, paste0(round(DD_avg[14,2]/1000), "GWh"), cex = 0.7)
text(73,y1, paste0(round(DD_avg[15,2]/1000), "GWh"), cex = 0.7)

# rmse
# h0 <- -0.15   # east
# h1 <- -0.09   # east
# ytest <- -.15-.005   # east

h0 <- -0.3   # west
h1 <- -0.2   # west
ytest <- -0.3-0.01   # west

abline(h=h0)
abline(h=h1)
a <- (h1-h0)/100   # scale to 100
lapply(1:length(rmse_reduction), function(i) {
  rect(xleft=i-.4, ybottom=min(ylim), xright=i+.4, ytop=min(ylim)+rmse_reduction[i]*a, border=NA, col="royalblue")
})
text(x=mean(1:length(rmse_reduction)), y=ytest, "RMSE Reduction", col="royalblue", font=2)

# add legend axis ticks
# legend(65, -0.03, lwd=1:2, col=c("grey", "red"), c("cdh", "hdh"),
# box.lty = 0, cex = 0.6)   # east

legend(65, -0.1, lwd=1:2, col=c("grey", "red"), c("cdh", "hdh"),
       box.lty = 0, cex = 0.6)   # west

# east
# ticks1 <- seq(-0.08, 0.2, 0.04)
# ticks2 <- seq(-0.15, -0.08, 0.015)
# west
ticks1 <- round(seq(-0.15, 0.3, 0.05), digits = 2)
ticks2 <- seq(-0.3, -0.2, 0.025)

axis(side=2, at=ticks1, labels=ticks1, las = 2)
axis(side=2, at=ticks2, labels=c("0%", "25%", "50%", "75%", "100%"), 
     tck=-.01, las = 2, col.axis="royalblue", cex.axis = 0.7)




