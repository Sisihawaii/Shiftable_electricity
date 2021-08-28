

#Script for making bar graph peak/base as a percent of (daily or Overall) demand mean (1 = 100%)
#Also makes bar graph that includes +2c data (compare side by side, red is +2c)
#Also makes bar graph that includes climate interaction data (purple is present w/interaction)
#whiskers are set at 1% & 99% percentile values; change "avgs" function for different percentile
#dataset from EY; result of 6a.flattendata1.R & 8a.flattendata2c.R

#manually enter: mc (number of cores on computer for parallel processing)
#end results: 2 ggplots & final dataset tot.agg2 
#& latex table for daily/overall base/peak 1%/avg/99% Present/+2C

library(ggplot2)
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)
library(data.table)
library(lubridate)
library(xtable)


load("06apo_INTFLs_CLIM.RData")
mega.clim <- list(eastINTFL, westINTFL, ercotINTFL)
rm(list = setdiff(ls(), "mega.clim"))
load("09pr_INTFLn2cs_new.RData") #made of 08po_INTFLS2c_new.RData & 06apo_INTFLs_new.RData
# alternatively: load("08po_INTFLS2c_new.RData"); load("06apo_INTFLs_new.RData")

mc = 1 #number of cores for mclapply


#-----------------------------
#LOAD FUNCTIONS

pstats2 <- function(df2){ #switch between dividing by daily mean and overall mean
  colnames(df2) [3:4] = c("dem", "og") #fourth column will always be original demand
  daily <- df2 %>%
    group_by(date)%>%
    summarise(mean = mean(og), 
              base = min(dem), peak = max(dem))
  daily$baseperc <- daily$base/daily$mean
  daily$peakperc <- daily$peak/daily$mean
  return(select(daily, date, baseperc, peakperc))
}

pstats2ov <- function(df2){ #dividing by overall mean #df2 is data of df for each BA 
  colnames(df2) [3:4] = c("dem", "og") #fourth column will always be original demand
  daily <- df2 %>%
    group_by(date)%>%
    summarise(base = min(dem), peak = max(dem))
  meanov <- mean(df2$og)
  daily$basepercyr <- daily$base/meanov
  daily$peakpercyr <- daily$peak/meanov
  return(select(daily, date, basepercyr, peakpercyr))
}

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

avgs <- function(df1){ #average values for the graph + whiskers
  avg <- mean(df1$baseperc, na.rm = TRUE)
  top <- quantile(df1$baseperc, .99, na.rm = TRUE) #99% percentile
  bottom <- quantile(df1$baseperc, .01,na.rm = TRUE) #1% percentile
  # top <- max(df1$baseperc)
  # bottom <- min(df1$baseperc)
  perctype <- "Base"
  a <- data.frame(avg, top, bottom, perctype)
  avg <- mean(df1$peakperc, na.rm = TRUE)
  top <- quantile(df1$peakperc, .99, na.rm = TRUE) #99% percentile
  bottom <- quantile(df1$peakperc, .01, na.rm = TRUE) #1% percentile
  # top <- max(df1$peakperc)
  # bottom <- min(df1$peakperc)
  perctype <- "Peak"
  b <- data.frame(avg, top, bottom, perctype)
  return(as.data.frame(rbind(a,b)))
}

#-----------------------------
#MAKING DATASET LISTS

mega.clim <- lapply(mega.clim, function(a){
  ax <- mclapply(a, function(i){
    i$PST <- format(i$UTC, usetz = TRUE,  tz="Etc/GMT+8") #standard PST, no PDT
    i$datePST <- as.Date(substr(i$PST, 1,10))
    i$raw <- i$demand
    i <- as.data.table(i)
    return(i)
  }, mc.cores = mc)
  return(ax)
})

mega <- list(eastINTFL, westINTFL, ercotINTFL)
mega <- lapply(mega, function(a){
  ax <- mclapply(a, function(i){
    i$PST <- format(i$UTC, usetz = TRUE,  tz="Etc/GMT+8") #standard PST, no PDT
    i$datePST <- as.Date(substr(i$PST, 1,10))
    i$raw <- i$demand
    i <- as.data.table(i)
    return(i)
  }, mc.cores = mc)
  return(ax)
})
mega.pr <- mega

#MAKING mega.co dataset (need mega)
#checking if they are the same length, deleting unmatched times between present and 2c
#joining dataframes to use mclapply

east.co <- list()
for(i in 1:length(eastINTFL2c)){
  df1 <- eastINTFL2c[[i]]
  df2 <- mega[[1]][[i]]
  if(nrow(df1) != nrow(df2)){
    df1 <- df1[df1$UTC %in% df2$UTC,]
    df2 <- df2[df2$UTC %in% df1$UTC,]
  }
  df1 <- as.data.table(select(df1, UTC, date_time_ntz, year, datentz, demand = demand2, flex, 
                flexd, hardd, hardd0.5, hardd0.25, smooth.load, smooth.load0.5, #smooth.load0.25,
                sf_id))
  df2 <- select(df2, UTC, raw, PST, datePST)
  east.co[[i]] <- as.data.table(left_join(df1, df2, by = "UTC"))
}

west.co <- list()
for(i in 1:length(westINTFL2c)){
  df1 <- westINTFL2c[[i]]
  df2 <- mega[[2]][[i]]
  if(nrow(df1) != nrow(df2)){
    df1 <- df1[df1$UTC %in% df2$UTC,]
    df2 <- df2[df2$UTC %in% df1$UTC,]
  }
  df1 <- as.data.table(select(df1, UTC, date_time_ntz, year, datentz, demand = demand2, flex, 
                flexd, hardd, hardd0.5, hardd0.25, smooth.load, smooth.load0.5, #smooth.load0.25,
                sf_id))
  df2 <- select(df2, UTC, raw, PST, datePST)
  west.co[[i]] <- as.data.table(left_join(df1, df2, by = "UTC"))
}

ercot.co <- list()
i = 1
df1 <- ercotINTFL2c[[i]]
df2 <- mega[[3]][[i]]
if(nrow(df1) != nrow(df2)){
  df1 <- df1[df1$UTC %in% df2$UTC,]
  df2 <- df2[df2$UTC %in% df1$UTC,]
}
df1 <- as.data.table(select(df1, UTC, date_time_ntz, year, datentz, demand = demand2, flex,
              flexd, hardd, hardd0.5, hardd0.25, smooth.load, smooth.load0.5, #smooth.load0.25,
              sf_id))
df2 <- select(df2, UTC, raw, PST, datePST)
ercot.co[[i]] <- as.data.table(left_join(df1, df2, by = "UTC"))

mega.co <- list(east.co, west.co, ercot.co)
rm(east.co, west.co, ercot.co)

#-----------------------------
#MASSIVE FOR LOOPS: Daily & Overall; Present & +2c & Climate interacitons
  #difference for Daily & Overall: use pstats2ov instead (overall is second loop)
  #difference for Present & +2c: use mega.co instead of mega
tot.agg2_climm <- list()
for(type in c("Daily", "Overall")){ #Daily & Overall base and peak
  if(type=="Overall")pstats2 <- pstats2ov #replace orig function with overall version
    
finalvalues <- list()
for( r in 1:3){ #Present & +2c
  if(r==1)mega <- mega.pr
  if(r==2)mega <- mega.co #use mega.co dataset instead; mega.co already calculated w/mega
  if(r==3)mega <- mega.clim #use mega.clim dataset; mega w/climate interactions

#-----------------------------
#REGIONS

  regionsagg <- NULL
  Xx <- lapply(mega, function(a){
    #within each interconnect
    ax <- mclapply(a, function(i){ #i is a BA
      #within each region
      regionsagg <- NULL
      df2 <- select(i, date = datentz, year, smooth.load, raw) #a=1
      df3 <- pstats2(df2)
      df2 <- select(i, date = datentz, year, smooth.load0.5, raw) #a=0.05
      df3 <- left_join(df3, pstats2(df2), by="date")
      df2 <- select(i, date = datentz, year, demand, raw) #a=0
      df3 <- left_join(df3, pstats2(df2), by="date")
      colnames(df3) <- c("date", "baseperc.sm", "peakperc.sm", "baseperc.sm0.5", "peakperc.sm0.5", "baseperc.d",  "peakperc.d")
      df3$date <- as.Date(df3$date)
      regionsagg <- rbind(regionsagg, df3)
      return(regionsagg)
    }, mc.cores = mc)
    bx <- rbindlist(ax)
    return(bx)
  })
  
  regionsagg <- rbindlist(Xx)
  rm(Xx)

#-----------------------------
#INTERCONNECT
  interconnectagg <- NULL
  for(n in 1:3){ #New Calculations, ##EDIT FOR: Ercot needs to be recalculted according to PST
    a <- mega[[n]]
    #find synchronous days
    date_times <- a[[1]]$PST 
    if(n==1|n==2){
      for(i in 2:15) date_times <- c(date_times, a[[i]]$PST)}
    dates <- substr(date_times, 1,10)
    dateCount <- tapply(dates, dates, length)
    dates <- as.Date(names(dateCount[dateCount==length(a)*24]))  # 360 = 24*15 (hours/day * # regions)
    
    ax <- mclapply(a, function(i){
      df1 <- i[datePST %in% dates] #only synchronized dates
      setkey(df1, PST)
      return(df1)
    }, mc.cores = mc)
    
    agg <- ax[[1]][,c("datePST","PST","demand","flexd","hardd","hardd0.5","raw")]
    if(n==1|n==2){
      for(i in 2:15){
        agg[,3:7] = agg[,3:7] + ax[[i]][,c("demand","flexd","hardd","hardd0.5","raw")]
      }}
    
    agg$aggSmooth = agg$aggSmooth.50 = rep(NA, nrow(agg)) #initialize
    
    for(d in unique(dates)){
      agg$aggSmooth[agg$datePST==d] = flatten(agg$hardd[agg$datePST==d], agg$flexd[agg$datePST==d])
      agg$aggSmooth.50[agg$datePST==d] = flatten(agg$hardd0.5[agg$datePST==d], agg$flexd[agg$datePST==d]*0.5)
    }
    
    agg$year <- as.numeric(substr(agg$datePST, 1, 4))
    dem.sel <- select(agg, date = datePST, year, aggSmooth, raw) #if days were max smoothed, how smooth? look at base and peak vs its avg
    stats <- pstats2(dem.sel)
    dem.sel <- select(agg, date = datePST, year, aggSmooth.50, raw)
    stats <- left_join(stats, pstats2(dem.sel), by="date")
    dem.sel <- select(agg, date = datePST, year, demand, raw)
    stats <- left_join(stats, pstats2(dem.sel), by="date")
    colnames(stats) <- c("date", "baseperc.sm", "peakperc.sm", "baseperc.sm0.5", "peakperc.sm0.5", "baseperc.d",  "peakperc.d")
    interconnectagg <- rbind(interconnectagg, stats) 
    rm(a, agg, stats, dem.sel, ax)
  }

#-----------------------------
#NATIONAL
  date_times <- NULL
  for(n in 1:3){
    for(i in 1:length(mega[[n]])) date_times <- c(date_times, mega[[n]][[i]]$PST)
  }
  dates <- substr(date_times, 1,10)
  dateCount <- tapply(dates, dates, length)
  dates <- as.Date(names(dateCount[dateCount==744])) # 744 = 24*31 (hours/day * # regions)
  
  for(n in 1:3){ #synchronized dates only
    for(i in 1:length(mega[[n]])) mega[[n]][[i]] <- mega[[n]][[i]][datePST %in% dates]; setkey(mega[[n]][[i]], PST)
  }
  
  agg <- mega[[3]][[1]][,c("PST", "datePST", "demand","flexd","hardd","hardd0.5", "raw")]
  for(n in 1:2){ #East and West join Ercot
    for(i in 1:15){
      agg[,3:7] = agg[,3:7] + mega[[n]][[i]][,c("demand","flexd","hardd","hardd0.5", "raw")]
    }}
  
  agg$aggSmooth = agg$aggSmooth.50 = rep(NA, nrow(agg)) #initialize
  
  for(d in unique(dates)){
    agg$aggSmooth[agg$datePST==d] = flatten(agg$hardd[agg$datePST==d], agg$flexd[agg$datePST==d])
    agg$aggSmooth.50[agg$datePST==d] = flatten(agg$hardd0.5[agg$datePST==d], agg$flexd[agg$datePST==d]*0.5)
  }
  
  #remember, this is within each interconnect
  agg$year <- as.numeric(substr(agg$datePST, 1, 4))
  dem.sel <- select(agg, date = datePST, year, aggSmooth, raw) #if days were max smoothed, how smooth? look at base and peak vs its avg
  stats <- pstats2(dem.sel)
  dem.sel <- select(agg, date = datePST, year, aggSmooth.50, raw)
  stats <- left_join(stats, pstats2(dem.sel), by="date") #not actually ntz (just for pstats2 function)
  dem.sel <- select(agg, date = datePST, year, demand, raw)
  stats <- left_join(stats, pstats2(dem.sel), by="date")
  colnames(stats) <- c("date", "baseperc.sm", "peakperc.sm", "baseperc.sm0.5", "peakperc.sm0.5", "baseperc.d",  "peakperc.d")
  nationalagg <- stats
  rm(agg, stats, dem.sel, date_times, dates, dateCount, d, n)

#-----------------------------
#COMBINE
  regionsagg$aggregate <- "Region"
  interconnectagg$aggregate <- "Interconnect"
  nationalagg$aggregate <- "Nation"
  
  tot.agg <- rbind(regionsagg,interconnectagg)
  tot.agg <- rbind(tot.agg, nationalagg)
  
  tot.agg2 <- NULL #one dataframe for entire graph that can be filtered
  for(e in c("Region", "Interconnect", "Nation")){
    df3 <- filter(tot.agg, aggregate == e)
    df1 <- select(df3, date, baseperc.d, peakperc.d, aggregate)
    colnames(df1)[2:3] <- c("baseperc", "peakperc")
    df2 <- avgs(df1)
    df2$aggregate <- first(df1$aggregate)
    df2$demandtype <- as.factor("Raw, alpha = 0")
    tot.agg2 <- rbind(tot.agg2, df2)
    df1 <- select(df3, date, baseperc.sm0.5, peakperc.sm0.5, aggregate)
    colnames(df1)[2:3] <- c("baseperc", "peakperc")
    df2 <- avgs(df1)
    df2$aggregate <- first(df1$aggregate)
    df2$demandtype <- as.factor("Smooth, alpha = 0.5")
    tot.agg2 <- rbind(tot.agg2, df2)
    df1 <- select(df3, date, baseperc.sm, peakperc.sm, aggregate)
    colnames(df1)[2:3] <- c("baseperc", "peakperc")
    df2 <- avgs(df1)
    df2$aggregate <- first(df1$aggregate)
    df2$demandtype <- as.factor("Smooth, alpha = 1")
    tot.agg2 <- rbind(tot.agg2, df2)
  }
  
  
  
  tot.agg2$aggregate = factor(tot.agg2$aggregate, levels = c("Region", "Interconnect", "Nation"))
  tot.agg2$demandtype = factor(tot.agg2$demandtype, levels = c("Raw, alpha = 0", "Smooth, alpha = 0.5", "Smooth, alpha = 1")) #naming for visual
  tot.agg2$perctype = as.character(tot.agg2$perctype)
  
  finalvalues[[r]] <- tot.agg2
  
  rm(tot.agg, tot.agg2, df1, df2, df3, e, i, regionsagg, interconnectagg, nationalagg)
}

# # -----------------------------
#VISUAL
tot.agg2_orig <- finalvalues[[1]]
tot.agg2_2c <- finalvalues[[2]]
tot.agg2_clim <- finalvalues[[3]]
rm(finalvalues, r)
tot.agg2_orig$aggregate = factor(tot.agg2_orig$aggregate, levels = c("Region", "Interconnect", "Nation"))
tot.agg2_orig$demandtype = factor(tot.agg2_orig$demandtype, levels = c("Raw, alpha = 0", "Smooth, alpha = 0.5", "Smooth, alpha = 1")) #naming for visual
tot.agg2_orig$perctype = as.character(tot.agg2_orig$perctype)
tot.agg2_2c$aggregate = factor(tot.agg2_2c$aggregate, levels = c("Region", "Interconnect", "Nation"))
tot.agg2_2c$demandtype = factor(tot.agg2_2c$demandtype, levels = c("Raw, alpha = 0", "Smooth, alpha = 0.5", "Smooth, alpha = 1")) #naming for visual
tot.agg2_2c$perctype = as.character(tot.agg2_2c$perctype)
tot.agg2_clim$aggregate = factor(tot.agg2_clim$aggregate, levels = c("Region", "Interconnect", "Nation"))
tot.agg2_clim$demandtype = factor(tot.agg2_clim$demandtype, levels = c("Raw, alpha = 0", "Smooth, alpha = 0.5", "Smooth, alpha = 1")) #naming for visual
tot.agg2_clim$perctype = as.character(tot.agg2_clim$perctype)


if(type == "Daily") bcolor = "#00BFC4"
if(type == "Overall") bcolor = "#008D91"

#VISUAL ONLY ORIG
g <- ggplot(tot.agg2_orig)+
  geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, fill = bcolor)+
  geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.6, color = "darkorchid4")+
  facet_grid(aggregate ~ demandtype)+
  theme_bw()+
  # coord_cartesian(ylim=c(0.6,1.5))+
  theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank())
if(type == "Daily"){
  g <- g + labs(x = "Daily Base or Peak Demand",
                y = "(Flattened Peak or Base Demand) / (Regional Mean of Daily Demand)",
                title = "Average Values of Base and Peak as Percentage of Mean Daily Demand")}
if(type == "Overall"){
  g <- g + labs(x = "Overall Base or Peak Demand",
                y = "(Flattened Overall Peak or Base Demand) / (Regional Mean of Overall Demand)",
                title = "Average Values of Overall Base and Peak as Percentage of Mean Overall Demand")}
  
# ggsave(plot = g, file = paste0("4.trans_flex_", type, ".pdf"), width = 7.41, height = 6.21)
rm(g)

#VISUAL WITH 2C (combo)
tot.agg2_orig$deg <- factor("Present")
tot.agg2_2c$deg <- factor("+2C")
tot.agg2_clim$deg <- factor("Present")
tot.agg2 <- rbind(tot.agg2_orig, tot.agg2_2c)

g <- ggplot(tot.agg2, aes(fill = deg))+
  geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, position = "dodge")+
  geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.5,
                color = "darkorchid4", position = position_dodge(width = 0.9))+
  facet_grid(aggregate ~ demandtype)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        legend.position = "bottom")#+
  # coord_cartesian(ylim=c(0.25, 1.75))
if(type == "Daily"){
  g <- g + labs(x = "Daily Base or Peak Demand",
                y = "(Flattened Peak or Base Demand) / (Regional Mean of Daily Demand)",
                title = "Average Values of Base and Peak as Percentage of Mean Daily Demand",
                fill = "Climate Scenario")+
    scale_fill_manual(values = c("Present" = "#00BFC4", "+2C" = "#F8766D"))}
if(type == "Overall"){
  g <- g + labs(x = "Daily Base or Peak Demand",
                y = "(Flattened Overall Peak or Base Demand) / (Regional Mean of Overall Demand)",
                title = "Average Values of Overall Base and Peak as Percentage of Mean Overall Demand",
                fill = "Climate Scenario")+
    scale_fill_manual(values = c("Present" = "#008D91", "+2C" = "#F53224"))}

# ggsave(plot = g, file = paste0("4b.trans_flex_", type, "_2c.pdf"), width = 7.41, height = 6.21)
rm(g)

tot.agg2$type <- type
if(type == "Daily") dailytot.agg2 <- tot.agg2
if(type == "Overall") overalltot.agg2 <- tot.agg2

#VISUAL WITH CLIM INTERACTIONS
tot.agg2_climm[[type]] <- tot.agg2_clim %>% mutate(type = type)
tot.agg2_orig$interact <- factor("Without Interaction")
tot.agg2_clim$interact <- factor("With Interaction")
tot.agg2 <- rbind(tot.agg2_orig, tot.agg2_clim)

# g <- ggplot(tot.agg2, aes(fill = interact))+
#   geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, position = "dodge")+
#   geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.5,
#                 color = "darkorchid4", position = position_dodge(width = 0.9))+
#   facet_grid(aggregate ~ demandtype)+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
#         axis.title = element_text(size = 11),
#         axis.title.x = element_blank(),
#         legend.position = "bottom")#+
# # coord_cartesian(ylim=c(0.25, 1.75))
# if(type == "Daily"){
#   g <- g + labs(x = "Daily Base or Peak Demand",
#                 y = "(Flattened Peak or Base Demand) / (Regional Mean of Daily Demand)",
#                 title = "Average Values of Base and Peak as Percentage of Mean Daily Demand",
#                 fill = "Climate Interaction")+
#     scale_fill_manual(values = c("Without Interaction" = "#00BFC4", "With Interaction" = "#B27BD7"))}
# if(type == "Overall"){
#   g <- g + labs(x = "Daily Base or Peak Demand",
#                 y = "(Flattened Overall Peak or Base Demand) / (Regional Mean of Overall Demand)",
#                 title = "Average Values of Overall Base and Peak as Percentage of Mean Overall Demand",
#                 fill = "Climate Interaction")+
#     scale_fill_manual(values = c("Without Interaction" = "#008D91", "With Interaction" = "#652E8B"))}

# ggsave(plot = g, file = paste0("4e.trans_flex_", type, "_interact.pdf"), width = 7.41, height = 6.21)

}

rm(mega.co, mega, mega.pr, tot.agg2, tot.agg2_2c, tot.agg2_orig, tot.agg2_clim, type, bcolor)
# # -----------------------------
#VISUAL W/BOTH OVERALL & DAILY

alltot.agg2 <- rbind(dailytot.agg2, overalltot.agg2)
present <-  filter(alltot.agg2, deg == "Present")
present$type <- factor(present$type, levels = c("Daily", "Overall"))
climchange <- filter(alltot.agg2, deg == "+2C")
climchange$type <- factor(climchange$type, levels = c("Daily", "Overall"))


for(r in c("Present", "+2C")){
  if(r=="Present"){
    g <- ggplot(present, aes(fill = type))+
      scale_fill_manual(values = c("Daily" = "#00BFC4", "Overall" = "#008D91"))
  }
  if(r=="+2C"){
    g <- ggplot(climchange, aes(fill = type))+
      scale_fill_manual(values = c("Daily" = "#F8766D", "Overall" = "#F53224"))
  }
  g <- g + geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, position = "dodge")+
    geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.5,
                  color = "darkorchid4", position = position_dodge(width = 0.9))+
    facet_grid(aggregate ~ demandtype)+
    theme_bw()+
    labs(x = "Base or Peak Demand",
         y = "(Flattened Peak or Base Demand) / (Regional Mean Demand)",
         title = paste(r,"Overall or Daily Average Values of Base and Peak as Percentage of Mean Demand"),
         fill = "Time Period")+
    theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
          axis.title = element_text(size = 11),
          axis.title.x = element_blank(),
          legend.position = "bottom")
  
  # ggsave(plot = g, file = paste0("4c.trans_flex_", r, "dayov.pdf"), width = 7.41, height = 6.21)
  
  rm(g)
}

# # -----------------------------
#VISUAL WITH ALL (Present & +2c, Daily & Overall)

alltot.agg2$degXtype <- alltot.agg2$type
alltot.agg2$degXtype[alltot.agg2$degXtype == "Overall"] <- "Overall"
alltot.agg2$degXtype[alltot.agg2$deg != "Present"] <- paste(alltot.agg2$degXtype[alltot.agg2$deg != "Present"], "+2c")
alltot.agg2$degXtype[alltot.agg2$degXtype == "Overall"] <- "Overall (3 Year)"
alltot.agg2$degXtype <- as.factor(alltot.agg2$degXtype)

g <- ggplot(alltot.agg2, aes(fill = degXtype))+
    scale_fill_manual(values = c("#00BFC4", "#F8766D",
                                 "#008D91", "#F53224"))+
  geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, position = "dodge")+
  geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.5,
                color = "darkorchid4", position = position_dodge(width = 0.9))+
  facet_grid(aggregate ~ demandtype)+
  theme_bw()+
  labs(x = "Base or Peak Demand",
       y = "(Flattened Peak or Base Demand) / (Regional Mean Demand)",
       title = paste("Average Values of Base and Peak as Percentage of Mean Demand"),
       fill = "Time Period")+
  theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank())

# ggsave(plot = g, file = paste0("4d.trans_flex_ALL.pdf"), width = 7.41, height = 6.21)


#VISUAL WITH INTERACTIONS (4 bars)
tot.agg2_climm <- rbindlist(tot.agg2_climm)
tot.agg2_climm$interact <- "With Interaction"
alltot.agg2$interact <- "Without Interaction"
alltot.agg2 <- rbind(select(alltot.agg2, -degXtype), tot.agg2_climm)

alltot.agg2f <- filter(alltot.agg2, deg != "+2C")
alltot.agg2f$interactXtype <- paste(alltot.agg2f$type, alltot.agg2f$interact)

g <- ggplot(alltot.agg2f, aes(fill = interactXtype))+
  scale_fill_manual(values = c("#00BFC4", "#B27BD7",
                               "#008D91", "#652E8B"))+
  geom_bar(aes(x = perctype, y = avg), stat = "identity", alpha = 0.7, position = "dodge")+
  geom_errorbar(aes(x = perctype, ymax = top, ymin = bottom), width = 0.5,
                color = "darkorchid4", position = position_dodge(width = 0.9))+
  facet_grid(aggregate ~ demandtype)+
  theme_bw()+
  labs(x = "Base or Peak Demand",
       y = "(Flattened Peak or Base Demand) / (Regional Mean Demand)",
       title = paste("Average Values of Base and Peak as Percentage of Mean Demand"),
       fill = "Interaction")+
  theme(plot.title = element_text(hjust = 0.6, size = 11, face = "bold"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank())

# ggsave(plot = g, file = paste0("4e.trans_flex_ALL_interact.pdf"), width = 7.41, height = 6.21)


rm(present, climchange, overalltot.agg2, dailytot.agg2,r,g, tot.agg2_climm, alltot.agg2f)


#MAKING PEAK & BASE TABLE
bgtable <- select(filter(alltot.agg2, interact == "Without Interaction"), -interact)
bgtable <- as.data.table(bgtable)
bgtable$demandtype <- as.factor(bgtable$demandtype)

bgtable2 <-dcast(melt(bgtable, id.vars = c("aggregate","perctype", "deg", "type", "demandtype")), aggregate+demandtype+deg~variable+perctype+type)
bgtable2$aggregate <- factor(bgtable2$aggregate, levels = c("Region", "Interconnect","Nation"))
bgtable2 <- bgtable2[order(bgtable2$aggregate),]
bgtable2[,4:15] <- round((bgtable2[,4:15]*100), digits = 1)
bgtable2 <- select(bgtable2, aggregate, demandtype, deg, 
                   "Base1D" = "bottom_Base_Daily", "Base1O" = "bottom_Base_Overall",
                   "BaseAD" = "avg_Base_Daily", "BaseAO" = "avg_Base_Overall",
                   "Base99D" = "top_Base_Daily", "Base99O" = "top_Base_Overall",
                   "Peak1D" = "bottom_Peak_Daily", "Peak1O" = "bottom_Peak_Overall",
                   "PeakAD" = "avg_Peak_Daily", "PeakAO" = "avg_Peak_Overall",
                   "Peak99D" = "top_Peak_Daily", "Peak99O" = "top_Peak_Overall")
xtable(bgtable2)




