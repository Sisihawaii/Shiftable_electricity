

#Script for building table of peak % reduction, base % increase, cv or sd % reduction, % of flat days
  #peak and base calculated on index (divide by mean demand). Important for interconnect & national
  #for BA values: season & flattening calculated based on ntz;each BA equally weighted to get final
  #for interc & nation values: season & flattening calculated based on PST
    #1. each BA w/in interconnect or national is average weighted for each day; use as baseline comparison (x0)
    #2. peakred/baseinc/cvred % calculated per day: (x0 - x1)/x0 where x1 is aggdemand or aggsmooth.load, etc
    #3. 365 days*3 yrs equally weighted to get final value
  #for overall version: using single peak/base/cv spread of all years
#dataset from EY; result of 6a.flattendata1.R

#manually enter: type ("Daily"/"Overall) & mc (number of cores on computer for parallel processing)

library(tidyverse)
library(parallel)
library(xtable)
library(data.table)
library(lubridate)
mc = 1 #number of cores

#-----------------------------
#load functions
flatten <- function(h,f){
  fill.order = order(h) #save original order
  flex.total = sum(f)
  x = as.numeric(length(f)) #find number of hours in each day
  cum.mean = (cumsum(sort(h)) + flex.total)/(1:x)
  fill.no  = which.min(cum.mean) #find the sma llest cumulative mean
  smooth.load = c(rep(cum.mean[fill.no], fill.no), sort(h)[(fill.no+1):x]) #fill with with smalles cumulative mean (in the order of least to greatest hard demand) until that observation
  smooth.load = smooth.load[order(fill.order)] #return to original order
  return(smooth.load)
}

getSeason <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,319,0620,0921,1220,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(cuts)
}

cv <- function(x) sd(x) / mean(x)


#-----------------------------
#LOAD & PREP
dailytable2 <- list()
overalltable2 <- list()
setwd("C:/Users/Eleanor/Desktop/SWITCH")
load("09pr_INTFLn2cs_new.RData")
for(r in 1:3){
  if(r == 2) { eastINTFL <- eastINTFL2c; westINTFL <- westINTFL2c; ercotINTFL <- ercotINTFL2c
    rm(eastINTFL2c, westINTFL2c, ercotINTFL2c)
  }
  if(r == 3){ rm(eastINTFL,westINTFL, ercotINTFL)
  load("06apo_INTFLs_CLIM.RData") }

mega <- list(eastINTFL, westINTFL, ercotINTFL)
mega <- lapply(mega, function(a){
  ax <- mclapply(a, function(i){
    
    # Find complete days with non-missing synchronized demand across all regions
    # convert from UTC to PST. Hours closer to actual day & begin on 1/1 all regions
    i$PST <- format(i$UTC, usetz = TRUE,  tz="Etc/GMT+8") #standard PST, no PDT
    i$datePST <- as.Date(substr(i$PST, 1,10))
    i$season <- getSeason(i$datePST) #PST season
    i <- as.data.table(i)
    return(i)
  }, mc.cores = mc)
  return(ax)
})



dailytable <- NULL
overalltable <- NULL
for(type in c("Daily", "Overall")){
  # type <- "Overall" #"Daily" or "Overall" type in
  ptm <- proc.time()

for(filter.s in c("ALL", "Winter", "Summer")){

if(filter.s == "Winter" | filter.s == "Summer"){
  mega.s <- lapply(mega, function(a){
    ax <- mclapply(a, function(i){
      i <- filter(i, season==filter.s)
      return(i)
    }, mc.cores = mc)
    return(ax)
  })
}else{mega.s <- mega}

#-----------------------------
#REGIONS flattening counted by regional day; final value not weighted by day or region
Xx <- lapply(mega, function(a){ #a is an interconnect
  #within each interconnect
  ax <- mclapply(a, function(i){ #i is a BA
    #within each region
    peakstats <- NULL
    i$season <- getSeason(i$datentz) #ntz season
    if(filter.s == "Winter" | filter.s == "Summer"){
      i <- filter(i, season==filter.s)}
    df1 <- i
    
    #compare peakstats against daily demand values (pre-smoothing) of peak, base, cv
    peakstats <- df1 %>% 
      {if(type == "Daily")group_by(., datentz) else group_by(., NULL)}%>%
      summarise(peakred0 = 0, peakred.25 = (max(demand) - max(smooth.load0.25))/max(demand), 
                peakred.5 = (max(demand) - max(smooth.load0.5))/max(demand), peakred = (max(demand) - max(smooth.load))/max(demand),
                base0 = 0, baseinc.25 = -1*(min(demand) - min(smooth.load0.25))/min(demand), 
                baseinc.5 = -1*(min(demand) - min(smooth.load0.5))/min(demand), baseinc = -1*(min(demand) - min(smooth.load))/min(demand),
                cvred0 = 0, cvred.25 = (cv(demand) - cv(smooth.load0.25))/cv(demand),
                cvred.5 = (cv(demand) - cv(smooth.load0.5))/cv(demand), cvred = (cv(demand) - cv(smooth.load))/cv(demand), 
                flat0 = ifelse(sd(demand)==0, 1, 0), flat.25 = ifelse(sd(smooth.load0.25)==0, 1, 0), 
                flat.5 = ifelse(sd(smooth.load0.5)==0, 1, 0), flat = ifelse(sd(smooth.load)==0, 1, 0),
                sf_id = first(sf_id), season = first(season),
                sdred0 = 0, sdred.25 = (sd(demand) - sd(smooth.load0.25))/sd(demand),
                sdred.5 = (sd(demand) - sd(smooth.load0.5))/sd(demand), sdred = (sd(demand) - sd(smooth.load))/sd(demand))
    
    return(peakstats) #return peakstats df for each BA region
  }, mc.cores = mc)
  bx <- rbindlist(ax)
  return(bx) #return one peakstats df for each interconnect
})

#rbind and label all BA peakstats:
Xx[[1]]$interc <- "East"
Xx[[2]]$interc <- "West"
Xx[[3]]$interc <- "Ercot"
peakstats <- rbindlist(Xx)
rm(Xx)
if(type == "Overall") peakstats$flat0 = peakstats$flat.25 = peakstats$flat.5 = peakstats$flat = rep(NA, nrow(peakstats))

#save for Sisi's map
if(filter.s == "ALL" & type == "Daily" & r == 1){
  savesd <- select(peakstats, datentz, peakred, baseinc, season, sdred, cvred, flat, sf_id, interc)
  peakstatss <- list()
  peakstatss[[1]] <- filter(savesd, interc == "East")
  peakstatss[[2]] <- filter(savesd, interc == "West")
  peakstatss[[3]] <- filter(savesd, interc == "Ercot")
  rm(savesd)
}

#constructing table 
f.table <- data.frame(NULL)
if(type == "Daily") peakstats <- as.data.frame(t(apply(peakstats[,2:17], 2, mean))) *100 #*100 = percentage form. find mean
if(type == "Overall") peakstats <- as.data.frame(t(apply(peakstats[,1:12], 2, mean))) *100
peakstats$level <- "BA"
f.table <- rbind(f.table, peakstats)


#-----------------------------
#INTERCONNECT flattening counted by PST. right now IS weighted according to day; but not region
peakstats <- NULL
for(n in 1:3){ #New Calculations, ##EDIT FOR: Ercot needs to be recalculted according to PST
 a <- mega.s[[n]]
  #find synchronous days
  date_times <- a[[1]]$PST 
  if(n==1|n==2){
  for(i in 2:15) date_times <- c(date_times, a[[i]]$PST)}
  dates <- substr(date_times, 1,10)
  dateCount <- tapply(dates, dates, length)
  dates <- as.Date(names(dateCount[dateCount==length(a)*24]))  # 360 = 24*15 (hours/day * # regions)
  
  a <- mclapply(a, function(i) i <- i[datePST %in% dates], mc.cores = mc) #only synchronized dates
  
  calcs <- mclapply(a, function(i){
    df1 <- i[datePST %in% dates] 
    calc <- i %>% #to get base values of each BA, to be mean-weighted (only looking at original demand)
      {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>% #peak and base values are indexed to individual mean demand
      summarise(daily.cv = cv(demand), daily.mean = mean(demand),
                daily.peak = max(demand)/mean(demand), daily.base = min(demand/mean(demand)),
                sf_id = first(sf_id))
    return(calc)
  }, mc.cores = mc)


  calcs <- rbindlist(calcs)

  OGvalues <- calcs %>%  #indexed peak and base are weighted (now 1 value for each interconnect per day); can compare to indexed agg value now
    {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>%  #15 or 1 values to weight
    summarise(weightcv = weighted.mean(daily.cv, daily.mean), #OR overall mean
              weightpeak = weighted.mean(daily.peak, daily.mean), weightbase = weighted.mean(daily.base, daily.mean))
  
  agg <- a[[1]][,c("datePST","PST","demand","flexd","hardd","hardd0.5","hardd0.25")]#,"smooth.load","smooth.load0.5","smooth.load0.25")]
  if(n==1|n==2){ #agg demand
    for(i in 2:15){
      agg[,3:7] = agg[,3:7] + a[[i]][,c("demand","flexd","hardd","hardd0.5","hardd0.25")]#,"smooth.load","smooth.load0.5","smooth.load0.25")]
    }} #3:10
  
  agg$aggSmooth = agg$aggSmooth.50 = agg$aggSmooth.25 = rep(NA, nrow(agg)) #initialize
  
  for(d in unique(dates)){
    agg$aggSmooth[agg$datePST==d] = flatten(agg$hardd[agg$datePST==d], agg$flexd[agg$datePST==d])
    agg$aggSmooth.50[agg$datePST==d] = flatten(agg$hardd0.5[agg$datePST==d], agg$flexd[agg$datePST==d]*0.5)
    agg$aggSmooth.25[agg$datePST==d] = flatten(agg$hardd0.25[agg$datePST==d], agg$flexd[agg$datePST==d]*0.25)
  }
  agg$year <- as.numeric(substr(agg$datePST, 1, 4))
  
  aggvalues <- agg %>% #getting indexed peak and base & cv for agg for each day
    {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>%  #base, peak, sd indexed to avg demand of the day (for the total interconnect)
    summarise(peak0 = max(demand)/mean(demand), peak.25 = max(aggSmooth.25)/mean(aggSmooth.25),
              peak.5 = max(aggSmooth.50)/mean(aggSmooth.50), peak = max(aggSmooth)/mean(aggSmooth),
              base0 = min(demand)/mean(demand), base.25 = min(aggSmooth.25)/mean(aggSmooth.25),
              base.5 = min(aggSmooth.50)/mean(aggSmooth.50), base = min(aggSmooth)/mean(aggSmooth),
              cv0 = cv(demand), cv.25 = cv(aggSmooth.25), cv.5 = cv(aggSmooth.50), cv = cv(aggSmooth),
              flat0 = ifelse(sd(demand)==0,1,0), flat.25 = ifelse(sd(aggSmooth.25)==0,1,0),
              flat.5 = ifelse(sd(aggSmooth.50)==0,1,0), flat = ifelse(sd(aggSmooth)==0,1,0),
              sumdem = sum(demand))
  if(type == "Overall") aggvalues$flat0 = aggvalues$flat.25 = aggvalues$flat.5 = aggvalues$flat = rep(NA, nrow(aggvalues))
  
  #percentage reduction/increase form:
  if(type == "Daily") aggcalc <- left_join(aggvalues, OGvalues, by = "datePST")
  if(type == "Overall") aggcalc <- data.frame(aggvalues, OGvalues) #picking single base/peak
    
  aggcalc[,c("peak0", "peak.25", "peak.5", "peak")] <-  (aggcalc$weightpeak - aggcalc[,c("peak0", "peak.25", "peak.5", "peak")])/aggcalc$weightpeak
  aggcalc[,c("base0", "base.25", "base.5", "base")] <-  -1*(aggcalc$weightbase - aggcalc[,c("base0", "base.25", "base.5", "base")])/aggcalc$weightbase
  aggcalc[,c("cv0", "cv.25", "cv.5", "cv")] <-  (aggcalc$weightcv - aggcalc[,c("cv0", "cv.25", "cv.5", "cv")])/aggcalc$weightcv
  if(type == "Daily") aggcalc <- select(aggcalc, -datePST)
  colnames(aggcalc)[1:16] <- c("peakred0","peakred.25","peakred.5","peakred","base0","baseinc.25","baseinc.5","baseinc",
                               "cvred0", "cvred.25", "cvred.5","cvred","flat0","flat.25","flat.5","flat") #FLAT VALUES DON'T MEAN ANYTHING FOR OVERALL
  if(type == "Overall") aggcalc$flat0 = aggcalc$flat.25 = aggcalc$flat.5 = aggcalc$flat = rep(NA, nrow(aggcalc))
  
  #final weighted mean reduction/increase %; weighted across the days; no change for overall
  ps <- data.frame(t(apply(aggcalc[1:12], 2, mean)), #peak/base/cv
                             t(apply(aggcalc[13:16], 2, function(x) sum(x)/sum(!is.na(x)) )) ) #flat
  peakstats <- rbind(peakstats, ps)
  rm(ps, aggcalc, aggvalues, OGvalues, agg, a, date_times, dateCount, dates, calcs)
}

df3 <- as.data.frame(t(apply(peakstats, 2, mean))) *100 #percentage value; unweighted mean among interconnects
df3$level <- "Interc"
if(type == "Overall") df3 <- select(df3, -flat0, -flat.25, -flat.5, -flat)
f.table <- rbind(f.table,df3)

#-----------------------------
#NATIONAL
date_times <- NULL
for(n in 1:3){
  for(i in 1:length(mega.s[[n]])) date_times <- c(date_times, mega.s[[n]][[i]]$PST)
}
dates <- substr(date_times, 1,10)
dateCount <- tapply(dates, dates, length)
dates <- as.Date(names(dateCount[dateCount==744])) # 744 = 24*31 (hours/day * # regions)

calcs <- NULL
for(n in 1:3){ #synchronized dates only
  for(i in 1:length(mega.s[[n]])) mega.s[[n]][[i]] <- mega.s[[n]][[i]][datePST %in% dates]

  calcs2 <- mclapply(mega.s[[n]], function(i){
    calc <- i %>% #to get base values of each BA, to be mean-weighted
      {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>% #peak and base values are indexed to individual mean demand
      summarise(daily.cv = cv(demand), daily.mean = mean(demand),
                daily.peak = max(demand)/mean(demand), daily.base = min(demand/mean(demand)),
                sf_id = first(sf_id))
    return(calc)
  }, mc.cores = mc)
  calcs[[n]] <- rbindlist(calcs2)
}

calcs <- rbindlist(calcs) #base and peak index & cv for each BA calculated at interconnect

OGvalues <- calcs %>%  #peak and base are indexed to individual BA demand
  {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>% #31 values to weight
  summarise(weightcv = weighted.mean(daily.cv, daily.mean), 
            weightpeak = weighted.mean(daily.peak, daily.mean), weightbase = weighted.mean(daily.base, daily.mean))


agg <- mega.s[[3]][[1]][,c("PST", "datePST", "demand","flexd","hardd","hardd0.5","hardd0.25","smooth.load","smooth.load0.5","smooth.load0.25")]
for(n in 1:2){ #East and West join Ercot
  for(i in 1:15){
    agg[,3:10] = agg[,3:10] + mega.s[[n]][[i]][,c("demand","flexd","hardd","hardd0.5","hardd0.25","smooth.load","smooth.load0.5","smooth.load0.25")]
  }}

agg$aggSmooth = agg$aggSmooth.50 = agg$aggSmooth.25 = rep(NA, nrow(agg)) #initialize

for(d in unique(dates)){
  agg$aggSmooth[agg$datePST==d] = flatten(agg$hardd[agg$datePST==d], agg$flexd[agg$datePST==d])
  agg$aggSmooth.50[agg$datePST==d] = flatten(agg$hardd0.5[agg$datePST==d], agg$flexd[agg$datePST==d]*0.5)
  agg$aggSmooth.25[agg$datePST==d] = flatten(agg$hardd0.25[agg$datePST==d], agg$flexd[agg$datePST==d]*0.25)
}
agg$year <- as.numeric(substr(agg$datePST, 1, 4))

aggvalues <- agg %>% #getting indexed peak and base & cv for agg
  {if(type == "Daily")group_by(., datePST) else group_by(., NULL)}%>%
  summarise(peak0 = max(demand)/mean(demand), peak.25 = max(aggSmooth.25)/mean(aggSmooth.25),
            peak.5 = max(aggSmooth.50)/mean(aggSmooth.50), peak = max(aggSmooth)/mean(aggSmooth),
            base0 = min(demand)/mean(demand), base.25 = min(aggSmooth.25)/mean(aggSmooth.25),
            base.5 = min(aggSmooth.50)/mean(aggSmooth.50), base = min(aggSmooth)/mean(aggSmooth),
            cv0 = cv(demand), cv.25 = cv(aggSmooth.25), cv.5 = cv(aggSmooth.50), cv = cv(aggSmooth),
            flat0 = ifelse(sd(demand)==0,1,0), flat.25 = ifelse(sd(aggSmooth.25)==0,1,0),
            flat.5 = ifelse(sd(aggSmooth.50)==0,1,0), flat = ifelse(sd(aggSmooth)==0,1,0),
            sumdem = sum(demand))

#percentage reduction/increase form:
if(type == "Daily") aggcalc <- left_join(aggvalues, OGvalues, by = "datePST")
if(type == "Overall") aggcalc <- data.frame(aggvalues, OGvalues) #picking single base/peak
aggcalc[,c("peak0", "peak.25", "peak.5", "peak")] <-  (aggcalc$weightpeak - aggcalc[,c("peak0", "peak.25", "peak.5", "peak")])/aggcalc$weightpeak
aggcalc[,c("base0", "base.25", "base.5", "base")] <-  -1*(aggcalc$weightbase - aggcalc[,c("base0", "base.25", "base.5", "base")])/aggcalc$weightbase
aggcalc[,c("cv0", "cv.25", "cv.5", "cv")] <-  (aggcalc$weightcv - aggcalc[,c("cv0", "cv.25", "cv.5", "cv")])/aggcalc$weightcv
if(type == "Daily") aggcalc <- select(aggcalc, -datePST)
colnames(aggcalc)[1:16] <- c("peakred0","peakred.25","peakred.5","peakred","base0","baseinc.25","baseinc.5","baseinc",
                       "cvred0", "cvred.25", "cvred.5","cvred","flat0","flat.25","flat.5","flat")
if(type == "Overall") aggcalc$flat0 = aggcalc$flat.25 = aggcalc$flat.5 = aggcalc$flat = rep(NA, nrow(aggcalc))

#final weighted mean reduction/increase %; not weighted across days; #same for overall
peakstats <- data.frame(t(apply(aggcalc[1:12], 2, mean)), #peak/base/cv
                 t(apply(aggcalc[13:16], 2, function(x) sum(x)/sum(!is.na(x)) )) ) #flat
rm(aggcalc, aggvalues, OGvalues, agg, date_times, dateCount, dates, calcs, calcs2)

df3 <- peakstats *100 #percentage value; unweighted mean across days
df3$level <- "National"
if(type == "Overall") df3 <- select(df3, -flat0, -flat.25, -flat.5, -flat)
f.table <- rbind(f.table,df3)


#FINAL
f.table$season = filter.s
if(type == "Daily") dailytable <- rbind(dailytable, f.table)
if(type == "Overall") overalltable <- rbind(overalltable, f.table)
rm(f.table, df3, d, filter.s, i, n, peakstats)
}

#END
if(type == "Daily") dailytable <- data.frame(round(dailytable[,1:16], digits = 1), level = dailytable[,17], season = dailytable$season)
if(type == "Overall") overalltable <- data.frame(round(overalltable[,1:12], digits = 1), level = overalltable[,13], season = overalltable$season)

proc.time() - ptm
}
dailytable2[[r]] <- dailytable
overalltable2[[r]] <- overalltable
}

condition <- c("W/O INTERACTION", "+2C", "W/INTERACTION")
for(r in 1:3){
  print(condition[r])
  print("DAILY")
  print(xtable(dailytable2[[r]]))
  print("OVERALL")
  print(xtable(overalltable2[[r]]))
}
# 
# saveRDS(peakstatss, "peakstatss.RDS") #save for Sisi's sd map

