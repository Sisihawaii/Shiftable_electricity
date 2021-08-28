
library(tidyverse)
library(parallel)
library(ggplot2)

load("06apo_INTFLs_new_rev.RData")

ptm <- proc.time()
mc = 16  #number of cores for mclapply
#16 core 40 GB ram ran everything in 63 min

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
#add seasons
mega <- list(eastINTFL, westINTFL, ercotINTFL)
mega.s <- lapply(mega, function(a){
  ax <- mclapply(a, function(i){
    i$season <- NA
    i$season[i$month == "12" | i$month == "1" | i$month == "2"] <- "Winter"
    i$season[i$month == "3" | i$month == "4" | i$month == "5"] <- "Spring"
    i$season[i$month == "6" | i$month == "7" | i$month == "8"] <- "Summer"
    i$season[i$month == "9" | i$month == "10" | i$month == "11"] <- "Fall"
    # i <- i[i$season == filter.s,]
    return(i)
  }, mc.cores = mc)
  return(ax)
})

#-----------------------------
#DAILY CALCULATIONS + ANNUAL/SEASON PT 1 CALCULATIONS
Xx <- lapply(mega.s, function(c){
  #within each interconnect
  ax <- mclapply(c, function(i){ #i is a BA
    #within each region
    df2 <- select(i, UTC, date_time_ntz, sf_id, datentz, month,
                  demand, flexd, hardd, smooth.load, season)
    df3 <- NULL
    df4 <- NULL
    ptm <- proc.time()
    for(a in seq(0, 1, by = 0.05)){
      for(d in unique(df2$datentz)){ #set df to each day to calculate flattening
        h <- na.omit(df2$hardd[df2$datentz == d] + (1-a)*df2$flexd[df2$datentz == d]) #hard demand after a
        f <- na.omit(df2$demand[df2$datentz == d]) - h
        dem <- na.omit(df2$demand[df2$datentz == d])
        if(length(h) > 22 & length(f) > 22){ #remove if less than 23 observations
          sm <- flatten(h,f)
          b <- data.frame(datentz = first(df2$datentz[df2$datentz == d]), season = first(df2$season[df2$datentz == d]),
                          sdred = (sd(dem) - sd(sm)) / sd(dem), alpha = a, sf_id = first(df2$sf_id[df2$datentz == d]))
          c <- data.frame(season = first(df2$season[df2$datentz == d]), year = substr(first(df2$datentz[df2$datentz == d]), 1, 4),
                          dem, sm, alpha = a) #TO CALCULATE ANNUAL (later)
          df3 <- rbind(df3, b)
          df4 <- rbind(df4, c)
        }
      }
    }
    both <- list(df3, df4)
    proc.time() - ptm
    return(both)
  }, mc.cores = mc)
  return(ax)
})

sdred.t <- NULL
sdred.t.an <- NULL

for(i in 1:length(Xx[[1]])){ 
  df1 <- Xx[[1]][[i]][[1]] #daily
  df1$interc <- "East"
  df1$sf_id <- i
  df1 <- select(df1, sdred, alpha, sf_id, interc)
  sdred.t <- rbind(sdred.t, df1)
  df2 <- Xx[[1]][[i]][[2]] #annual pt1
  df2$interc <- "East"
  df2$sf_id <- i
  sdred.t.an <- rbind(sdred.t.an, df2)
}
for(i in 1:length(Xx[[2]])){ 
  df1 <- Xx[[2]][[i]][[1]]
  df1$interc <- "West"
  df1$sf_id <- i
  df1 <- select(df1, sdred, alpha, sf_id, interc)
  sdred.t <- rbind(sdred.t, df1)
  df2 <- Xx[[2]][[i]][[2]] #annual pt1
  df2$interc <- "West"
  df2$sf_id <- i
  sdred.t.an <- rbind(sdred.t.an, df2)
}
i = 1
df1 <- Xx[[3]][[i]][[1]]
df1$interc <- "Ercot"
df1$sf_id <- i
df1 <- select(df1, sdred, alpha, sf_id, interc)
sdred.t <- rbind(sdred.t, df1)
df2 <- Xx[[3]][[i]][[2]] #annual pt1
df2$interc <- "Ercot"
df2$sf_id <- i
sdred.t.an <- rbind(sdred.t.an, df2)

sdred.t$sdred[is.infinite(sdred.t$sdred)] <- 0 #infinite is 0/0; these results are from alpha = 0

#mean for each BA
sdred.BA <- sdred.t %>%
  group_by(alpha, sf_id, interc)%>%
  summarise(sdmean = mean(sdred, na.rm = TRUE))
#general mean
sdred.agg <- sdred.t %>%
  group_by(alpha)%>%
  arrange(alpha, .by_group = TRUE) %>%
  summarise(sdmean = mean(sdred, na.rm = TRUE))

#-----------------------------
#DAILY GRAPH
sdred.BA$interc[sdred.BA$interc == "Ercot"] <- "ERCOT (1 region)" #use all caps bc acronym
sdred.BA$interc[sdred.BA$interc == "East"] <- "East (15 regions)"
sdred.BA$interc[sdred.BA$interc == "West"] <- "West (15 regions)"
ggplot()+
  geom_line(sdred.BA, mapping = aes(x = alpha, y = sdmean, group = interaction(sf_id, interc), 
                                    color = interc), alpha = 0.7)+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), y = "SD Reduction Mean", title = expression("Daily SD Reduction Mean according to "*alpha*" value"), color = "Area of region")+
  theme_minimal()+
  geom_line(sdred.agg, mapping = aes(x = alpha, y = sdmean, color = "Average overall"), size = 1)+
  scale_color_manual(values = c(
    "Average overall" = "black",
    "East (15 regions)" = "#619CFF",
    "West (15 regions)" = "#F8766D",
    "ERCOT (1 region)" = "#00BA38"))+
  theme(legend.position = c(.85,.15))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
ggsave("92.daily_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")
#-----------------------------
#ANNUAL/SEASON PT2 CALCULATIONS

sdred.BA.an <- NULL
sdred.BA.an <- sdred.t.an %>% 
  group_by(alpha, year, sf_id, interc) %>% 
  summarise(sdred.an = (sd(dem) - sd(sm))/sd(dem)) #annual sd
sdred.BA.an <- sdred.BA.an %>% 
  group_by(alpha, sf_id, interc) %>% #leaving out year
  summarise(sdmean = mean(sdred.an)) #mean of yearly SD red spread according to BA; use for graph
sdred.BAagg.an <- sdred.BA.an %>% 
  group_by(alpha) %>% 
  summarise(sdmean = mean(sdmean)) #mean of yearly SD red across (BA weighted equally); black line on graph

sdred.BAseas.an <- NULL
sdred.BAseas.an <- sdred.t.an %>%
  group_by(alpha, season, year, sf_id, interc) %>% 
  summarise(sdred.an = (sd(dem) - sd(sm))/sd(dem)) #seasonal sd (every season, every year) according to BA
sdred.BAseas.an <- sdred.BAseas.an %>% 
  group_by(alpha, season) %>% 
  summarise(sdmean = mean(sdred.an)) #mean of seasonal spread (BA weighted equally); use for graph
sdred.BAaggseas.an <- sdred.BAseas.an %>% 
  group_by(alpha) %>% 
  summarise(sdmean = mean(sdmean)) #mean of all seasons; black line on graph

#-----------------------------
#ANNUAL/SEASON GRAPH

sdred.BA.an$interc[sdred.BA.an$interc == "Ercot"] <- "ERCOT (1 region)" #use all caps bc acronym
sdred.BA.an$interc[sdred.BA.an$interc == "East"] <- "East (15 regions)"
sdred.BA.an$interc[sdred.BA.an$interc == "West"] <- "West (15 regions)"
ggplot()+
  geom_line(sdred.BA.an, mapping = aes(x = alpha, y = sdmean, group = interaction(sf_id, interc), 
                                    color = interc), alpha = 0.7)+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), y = "SD Reduction Mean", title = expression("Annual SD Reduction Mean according to "*alpha*" value"), color = "Area of region")+
  theme_minimal()+
  geom_line(sdred.BAagg.an, mapping = aes(x = alpha, y = sdmean, color = "Average overall"), size = 1)+
  scale_color_manual(values = c(
    "Average overall" = "black",
    "East (15 regions)" = "#619CFF",
    "West (15 regions)" = "#F8766D",
    "ERCOT (1 region)" = "#00BA38"))+
  theme(legend.position = c(.85,.15))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
ggsave("92.annual_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")

ggplot()+ 
  geom_line(sdred.BAseas.an, mapping = aes(x = alpha, y = sdmean, group = season, color = season), alpha = 0.5) + 
  geom_line(sdred.BAaggseas.an, mapping = aes(x = alpha, y = sdmean, color = "Average"), size = 1)+
  theme_minimal()+
  scale_color_manual(values = c(
    "Winter" = "#C77CFF",
    "Spring" = "#7CAE00",
    "Summer" = "#00BFC4",
    "Fall" = "#F8766D",
    "Average" = "black"))+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), y = "SD Reduction Mean", title = expression("Average Seasonal SD Reduction Mean according to "*alpha*" value across all regions"), color = "Season")+
  theme(legend.position = c(0.87, 0.18))
ggsave("92.season_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")
#-----------------------------
proc.time() - ptm

rm(df1,df2,mega,mega.s, sdred.t, sdred.t.an,Xx)




