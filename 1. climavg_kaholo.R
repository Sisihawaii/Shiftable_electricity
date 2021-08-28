## 1. take average of the temperature for each cell per hour avg for 12 years
# use results from 0.3 for temperature per cell
# run in kaholo for the climateavg
# output 'wdh_clim12.RDS' 

rm(list = ls())

library(data.table)
library(foreach)
library(parallel)
library(ggplot2)
library(rgdal)
library(raster)
library(rgeos)
library(maptools)
library(sp)
library(tidyr)
library(dplyr)

detectCores()
# read all temp per cell data
t1 <- readRDS("temp_by_cell.RDS")
t2 <- readRDS("temp_by_cell1.RDS")
t3 <- readRDS("temp_by_cell2.RDS")
t4 <- readRDS("temp_by_cell3.RDS")
t5 <- readRDS("temp_by_cell4.RDS")
summary(t1)

# only climate variables
t1c <- t1[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
mean(t1c$climH) # 7.433
mean(t1c$climC) # 2.736
# saveRDS(t1c, "clim_cell.RDS")
t2c <- t2[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
t3c <- t3[, list(climH = mean(hdh_18), climC = mean(cdh_18)), by = list(cell)]
t_c <- rbind(t1c, t2c, t3c)
setkey(t_c, cell)
# 12 year climate avg per cell
tc_avg <- t_c[, list(climH = mean(climH), climC = mean(climC)), by = list(cell)]
saveRDS(tc_avg, "clim_cell12.RDS")

# take avg for each file
t1_avg <- t1[, list(value, date_time, cdh_18, hdh_18, 
                    climH = mean(hdh_18), climC = mean(cdh_18),
                    climH_H = mean(hdh_18)*hdh_18, climC_C = mean(cdh_18)*cdh_18), 
              by = list(cell)] # working 
# saveRDS(t1_avg, "clim_dh_cell.RDS")

# merge the climate avg with temp per cell data t1 by cell
t_merge <- t1[tc_avg, ]
t_merge_avg <- t_merge[, list(value, date_time, cdh_18, hdh_18, climH, climC, 
                              climH_H = climH*hdh_18, climC_C = climC*cdh_18), 
                       by = list(cell)]

# stats for US overall without weight
climH_avg <- mean(t1_avg$climH)  # equal weight mean/no pop weight 7.43
climC_avg <- mean(t1_avg$climC)  # 2.74
climH_sd <- sd(t1_avg$climH)     # 3.77
climC_sd <- sd(t1_avg$climC)     # 1.77

#################################################################

# modifiy the function for climate variable aggregation
calc.wdh.by.ID.clim <- function(ba.pop.dt, full.temp.dt){
  #pop.dt should be the population weights for all of the cells in BA i
  #full.temp.dt should be the full temperature history for all cells in US
  ba.temp <- full.temp.dt[ba.pop.dt, nomatch = 0] 
  setkey(ba.temp, date_time)
  
  wdh <- ba.temp[, list(temp = sum(value * w), hdh_18 = sum(hdh_18 * w), cdh_18 = sum(cdh_18 * w),
                        climH = sum(climH_H * w), climC = sum(climC_C * w)), 
                 by = list(sf, ID, date_time)]
  setkey(wdh, date_time)
  return(wdh)
}

##################################################################

#merge in temp values to pop
pop.dt <- readRDS("popweight.RDS")
pop.cells <- pop.dt[, unique(cell)]
pop.dt <- split(pop.dt, by = c("sf", "ID"), drop = T)
lapply(pop.dt, function(dt) setkey(dt, cell))
setkey(t_merge_avg, cell)

wdh_clim.dt <- lapply(pop.dt, calc.wdh.by.ID.clim, full.temp.dt = t_merge_avg)
wdh_clim.dt <- rbindlist(wdh_clim.dt) # climH and climC are climate and hdh/cdh interaction term
setkey(wdh_clim.dt, sf, ID, date_time)
saveRDS(wdh_clim.dt, "wdh_clim12.RDS")





