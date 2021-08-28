
#smoothing data for 2c

library(tidyverse)
library(parallel)

load("08pr_2cs_new.RData") #flexload_w_2C_new.RDS, flexload_e_2C_new.RDS, flexload_Ercot_2C_new.RDS

mc = 4

ptm <- proc.time()

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
#INSERT UTC (temp fix)

# names <- c("east", "west")
# for(a in 1:2){
#   data1 <- get(paste0(names[a], "2c"))
#   data7 <- get(paste0(names[a], "INTFL"))
#   for(i in 1:length(data1)){
#     df1 <- data1[[i]]
#     df2 <- select(data7[[i]], UTC, d = date_time_ntz)
#     df2$merge <- substr(df2$d, 1, 18)
#     df1$merge <- substr(df1$date_time_ntz, 1,18)
#     df1 <- left_join(df1, df2, by = "merge")
#     data1[[i]] <- select(df1, -merge, -d)
#   }
#   assign(paste0(names[a], "2c"), data1)
# }

# ercot2c <- rename(ercot2c, "UTC" = "date_time")
#-----------------------------
#adding in residuals (NO LONGER NEEDED, use yc instead)

names <- colnames(east2c[[1]])
#EAST
data1 <- east2c
data6 <- residuals[[1]]

east2c <- mclapply(data1, function(i){
  df1 <- as.data.frame(i) #data1 is all BA regions with demand
  idn = first(df1$sf_id)
  df6 <- filter(data6, sf_id == idn)
  df6 <- select(df6, UTC, lresid)
  df1 <- left_join(df1, df6, by = c("date_time" = "UTC")) #adding residuals row

  return(df1)
}, mc.cores = mc)

#WEST
data1 <- west2c
data6 <- residuals[[2]]

west2c <- mclapply(data1, function(i){
  df1 <- as.data.frame(i) #data1 is all BA regions with demand
  idn = first(df1$sf_id)
  df6 <- filter(data6, sf_id == idn)
  df6 <- select(df6, UTC, lresid)
  df1 <- left_join(df1, df6, by = c("date_time" = "UTC")) #adding residuals row
  return(df1)
}, mc.cores = mc)

#ERCOT
data1 <- list()
data1[[1]] <- ercot2c #ALREADY HAS UTC
data6 <- residuals[[3]]
ercot2c <- list()

#ercot needs:
df1 <- select(data1[[1]], names)
df6 <- select(data6, UTC, lresid)
df1 <- left_join(df1, df6, by = c("date_time" = "UTC")) #adding residuals row
ercot2c[[1]] <- df1


#-----------------------------
#flattening and smoothing
c.times <- data.frame(date_time = seq(from=as.POSIXct("2016-01-01 0:00", tz="UTC"),
                                      to=as.POSIXct("2018-12-31 23:00", tz="UTC"),
                                      by="hour"))
mega.2c <- list(east2c, west2c, ercot2c)
data5 <- lapply(mega.2c, function(a){ #function variable HAS to be i
  #within each interconnect
  ax <- mclapply(a, function(i){ #i is a BA
    #within each region
    df1 <- i
    df1$datentz <- as.character(substr(df1$date_time_ntz, 1, 10))

    # #adding residual on
    # df1$demand2 <- df1$yhat1+df1$resid
    # df1$demand2 <- df1$yc #yc is with interaction terms
    df1$demand2 <- exp(log(df1$yhat1)+df1$lresid)
    
    #hardd, alpha0.5, alpha0.25
    df1$flexd = df1$demand2*(df1$yhat1 - df1$yhat0)/ df1$yhat1 #in real demand values #alpha flex by multiplying alpha
    df1$hardd = df1$demand2*df1$yhat0/df1$yhat1 #hard load is inflexible remaining load
    df1$hardd0.5 <- df1$hardd + 0.5*df1$flexd
    df1$hardd0.25 <- df1$hardd + (1-0.25)*df1$flexd

    df1 <- df1[!is.na(df1$hardd),] #removing rows with no hardd (needs yhat1, yhat0, demand2 values)
    df1 <- right_join(df1, c.times, by = "date_time") #leaves NA spaces in missing hours
    # setkey(df1, date_time) #assure correct order before flattening
    df1$datentz[is.na(df1$datentz)] <- 0 #fill all NAs with 0
    
    for(d in unique(df1$datentz)){ #for each date
      
      if(sum(is.na(df1$hardd[df1$datentz == d]))==0){ #see if all hours of the day are present
    
    #flattening
      h = na.omit(df1$hardd[df1$datentz == d])
      f = na.omit(df1$demand2[df1$datentz == d] - h)
      df1$smooth.load[df1$datentz == d & !is.na(df1$hardd)] = flatten(h,f)
    #flattening 0.5
      h = na.omit(df1$hardd0.5[df1$datentz == d])
      f = na.omit(df1$demand2[df1$datentz == d] - h)
      df1$smooth.load0.5[df1$datentz == d & !is.na(df1$hardd)] <- flatten(h,f)

    # #flattening 0.25
    #   h = na.omit(df1$hardd0.25[df1$datentz == d])
    #   f = na.omit(df1$demand2[df1$datentz == d] - h)
    #     df1$smooth.load0.25[df1$datentz == d & !is.na(df1$hardd)] <- flatten(h,f)
      }
    }
      
    df1 <- df1[!is.na(df1$hardd),] #removing un-smoothable days
    
    df1$date_time_ntz <- as.character(df1$date_time_ntz)  #making values consistent length
    df1$date_time_ntz[str_length(df1$date_time_ntz)==14] <- 
      paste(substr(df1$date_time_ntz[str_length(df1$date_time_ntz)==14], 1, 10), 
            "00:00:00", substr(df1$date_time_ntz[str_length(df1$date_time_ntz)==14], 12, 14))
    
    df1 <- rename(df1, "UTC" = "date_time")
    return(df1)
  }, mc.cores = mc)
  return(ax)
})

x <- colnames(data5[[1]][[1]])
data5[[3]][[1]] <- select(data5[[3]][[1]], c(x))
eastINTFL2c <- data5[[1]]
westINTFL2c <- data5[[2]]
ercotINTFL2c <- data5[[3]]


#-----------------------------
#SAVING
save(eastINTFL2c, westINTFL2c, ercotINTFL2c, file = "08po_INTFLS2c_new.RData")
# save(east2c, west2c, ercot2c, file = "08pr_2cs_new.RData")

rm(df1, data5, data1, names, c.times, data6,
   east2c, west2c, ercot2c, residuals,
   mega.2c)

