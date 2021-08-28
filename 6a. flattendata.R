
#script to flatten data
#last edit: 22 June 2021

library(tidyverse)
library(parallel)
library(data.table)


mc = 4

setwd("~/Desktop/datasets")
# load("06apr_preFL_new.RData")
    #loads from Sisi raw files:
    # easti <- readRDS("~/Downloads/flexload_BAe15_noclim.RDS")
    # westi <- readRDS("~/Downloads/flexload_BAw15_noclim.RDS")
    # ercotOG <- readRDS("~/Downloads/flexload_ercot_new.RDS")
# load("06apr_preFL_CLIM.RData")
    # westi <- readRDS("~/Downloads/flexload_BAw15_new.RDS")
    # easti <- readRDS("~/Downloads/flexload_BAe15_new.RDS")
    # ercotOG <- readRDS("~/Downloads/flexload_ercot_new.RDS")

#-----------------------------
#LOAD FUNCTIONS

flatten <- function(h,f) {
  fill.order = order(h) #save original order
  flex.total = sum(f)
  x = as.numeric(length(f)) #find number of hours in each day
  cum.mean = (cumsum(sort(h)) + flex.total)/(1:x)
  fill.no  = which.min(cum.mean) #find the smallest cumulative mean
  smooth.load = c(rep(cum.mean[fill.no], fill.no), sort(h)[(fill.no+1):x]) #fill with with smalles cumulative mean (in the order of least to greatest hard demand) until that observation
  smooth.load = smooth.load[order(fill.order)] #return to original order
  return(smooth.load)
}

#-----------------------------
#filling & flattening/smoothing
ercotOG$sf_id <- 1
ercoti <- list(ercotOG)
c.times <- data.frame(date_time = seq(from=as.POSIXct("2016-01-01 0:00", tz="UTC"),
                                to=as.POSIXct("2018-12-31 23:00", tz="UTC"),
                                by="hour"))
mega.og <- list(easti, westi, ercoti)
data5 <- lapply(mega.og, function(g){ #function variable HAS to be i
  #within each interconnect
  ax <- mclapply(g, function(i){ #i is a BA
    #within each region
    df1 <- i
    
    #hardd, alpha0.5, alpha0.25
    df1$flexd = df1$demand*(df1$yhat1 - df1$yhat0)/ df1$yhat1 #in real demand values #alpha flex by multiplying alpha
    df1$hardd = df1$demand*df1$yhat0/df1$yhat1 #hard load is inflexible remaining load
    df1$hardd0.5 <- df1$hardd + 0.5*df1$flexd
    df1$hardd0.25 <- df1$hardd + (1-0.25)*df1$flexd
    
    df1 <- df1[!is.na(df1$hardd),] #removing rows with no hardd (needs demand & yhat1 & yhat0 values)
    df1 <- right_join(df1, c.times, by = "date_time") #leaves NA spaces in missing hours
    # setkey(df1, date_time) #assure correct order before flattening
    df1$datentz[is.na(df1$datentz)] <- 0 #fill all NAs with 0
    
    #flattening
    for(d in unique(df1$datentz)){ #for each date
      
      if(sum(is.na(df1$hardd[df1$datentz == d]))==0){ #see if all hours of the day are present
      
      #flattening 1
      h = na.omit(df1$hardd[df1$datentz == d])
      f = na.omit(df1$demand[df1$datentz == d] - h)
      df1$smooth.load[df1$datentz == d & !is.na(df1$hardd)] <- flatten(h,f)
      
      #flattening 0.5
      h = na.omit(df1$hardd0.5[df1$datentz == d])
      f = na.omit(df1$demand[df1$datentz == d] - h)
      df1$smooth.load0.5[df1$datentz == d & !is.na(df1$hardd)] <- flatten(h,f)
      
      #flattening 0.25
      h = na.omit(df1$hardd0.25[df1$datentz == d])
      f = na.omit(df1$demand[df1$datentz == d] - h)
      df1$smooth.load0.25[df1$datentz == d & !is.na(df1$hardd)] <- flatten(h,f)
      
      }
    }
    
    df1 <- df1[!is.na(df1$hardd),] #taking out NAs
    
    df1$date_time_ntz <- as.character(df1$date_time_ntz)
    df1$date_time_ntz[str_length(df1$date_time_ntz)==14] <- 
      paste(substr(df1$date_time_ntz[str_length(df1$date_time_ntz)==14], 1, 10), 
            "00:00:00", substr(df1$date_time_ntz[str_length(df1$date_time_ntz)==14], 12, 14))
    
    df1$lresid <- log(df1$demand) - log(df1$yhat1)
    
    df1 <- rename(df1, "UTC" = "date_time")
    
    return(df1)
  }, mc.cores = mc)
  return(ax)
})

x <- colnames(data5[[1]][[1]])
data5[[3]][[1]] <- select(data5[[3]][[1]], c(x))
eastINTFL <- data5[[1]]
westINTFL <- data5[[2]]
ercotINTFL <- data5[[3]]

#-----------------------------
#finding residuals
residuals <- mclapply(data5, function(g){
  #within each interconnect
  t.resid <- NULL
  for(i in 1:length(g)){
    df1 <- select(g[[i]], UTC, sf_id, lresid, date_time_ntz) #ln residual
    t.resid <- rbind(t.resid, df1)
  }
  return(t.resid)
}, mc.cores = 3)

#-----------------------------
#SAVING
save(eastINTFL, westINTFL, ercotINTFL, file = "06apo_INTFLs_new.RData")
saveRDS(residuals, file = "INTFL_residuals.RDS")


# save(eastINTFL, westINTFL, ercotINTFL, file = "06apo_INTFLs_CLIM.RData")

rm(x, data5, mega.og)




