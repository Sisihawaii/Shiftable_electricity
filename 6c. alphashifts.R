
library(tidyverse)
library(parallel)
library(ggplot2)
library(data.table)
library(lubridate)
library(grid)
library(gridExtra)

setwd("~/SWITCHuploads")
load("06apo_INTFLs_new.RData")

ptm <- proc.time()
mc = 18  #number of cores for mclapply

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
getSeason <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,319,0620,0921,1220,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(cuts)
}
#-----------------------------
#add seasons & ensure all hours day
c.times <- data.frame(UTC = seq(from=as.POSIXct("2016-01-01 0:00", tz="UTC"),
                                to=as.POSIXct("2018-12-31 23:00", tz="UTC"),
                                by="hour"))

mega <- list(eastINTFL, westINTFL, ercotINTFL)
mega.s <- lapply(mega, function(a){
  ax <- mclapply(a, function(i){
    i <- as.data.table(right_join(i, c.times, by = "UTC")) #leaves NA spaces in missing hours
    keep <- unique(i$datentz[!is.na(i$datentz)]) #should have none from INTFL dataset (safety)
    i <- i[datentz %in% keep]
    setkey(i, UTC) #assure correct order before flattening
    i$season <- getSeason(as.Date(i$datentz))
    return(i)
  }, mc.cores = mc)
  return(ax)
})

#-----------------------------
#DAILY CALCULATIONS + OVERALL/SEASON PT 1 CALCULATIONS
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
        sm <- flatten(h,f)
        b <- data.frame(datentz = first(df2$datentz[df2$datentz == d]), season = first(df2$season[df2$datentz == d]),
                        sdred = (sd(dem) - sd(sm)) / sd(dem), alpha = a, sf_id = first(df2$sf_id[df2$datentz == d]))
        c <- data.frame(season = first(df2$season[df2$datentz == d]), year = substr(first(df2$datentz[df2$datentz == d]), 1, 4),
                        dem, sm, alpha = a) #TO CALCULATE OVERALL (later)
        df3 <- rbind(df3, b)
        df4 <- rbind(df4, c)
      }
    }
    both <- list(df3, df4)
    proc.time() - ptm
    return(both)
  }, mc.cores = mc)
  return(ax)
})
ptm - proc.time()

sdred.t <- NULL
sdred.t.ov <- NULL

for(i in 1:length(Xx[[1]])){ 
  df2 <- Xx[[1]][[i]][[2]] #overall pt1
  df2$interc <- "East"
  df2$sf_id <- i
  sdred.t.ov <- rbind(sdred.t.ov, df2) 
  df1 <- Xx[[1]][[i]][[1]] #daily
  df1$interc <- "East"
  df1$sf_id <- i
  df1 <- select(df1, sdred, alpha, sf_id, interc)
  df1$avgdem <- mean(df2$dem)
  sdred.t <- rbind(sdred.t, df1) 
}
for(i in 1:length(Xx[[2]])){ 
  df2 <- Xx[[2]][[i]][[2]] #overall pt1
  df2$interc <- "West"
  df2$sf_id <- i
  sdred.t.ov <- rbind(sdred.t.ov, df2)
  df1 <- Xx[[2]][[i]][[1]] #daily
  df1$interc <- "West"
  df1$sf_id <- i
  df1 <- select(df1, sdred, alpha, sf_id, interc)
  df1$avgdem <- mean(df2$dem)
  sdred.t <- rbind(sdred.t, df1)
  
}
i = 1
df2 <- Xx[[3]][[i]][[2]] #overall pt1
df2$interc <- "Ercot"
df2$sf_id <- i
sdred.t.ov <- rbind(sdred.t.ov, df2)
df1 <- Xx[[3]][[i]][[1]] #daily
df1$interc <- "Ercot"
df1$sf_id <- i
df1 <- select(df1, sdred, alpha, sf_id, interc)
df1$avgdem <- mean(df2$dem)
sdred.t <- rbind(sdred.t, df1)


# sdred.t$sdred[is.infinite(sdred.t$sdred)] <- 0 #infinite is 0/0; these results are from alpha = 0
rm(df1, df2, c.times, mega, mega.s)

#mean for each BA
sdred.BA <- sdred.t %>%
  group_by(alpha, sf_id, interc)%>%
  summarise(sdmean = mean(sdred, na.rm = TRUE), .groups = "drop_last")
#general weighted mean
sdred.agg <- sdred.t %>%
  group_by(alpha)%>%
  arrange(alpha, .by_group = TRUE) %>%
  summarise(sdmean = weighted.mean(sdred, w = avgdem, na.rm = TRUE), .groups = "drop_last")

#-----------------------------
#DAILY GRAPH
sdred.BA$interc[sdred.BA$interc == "Ercot"] <- "ERCOT (1 region)" #use all caps bc acronym
sdred.BA$interc[sdred.BA$interc == "East"] <- "East (15 regions)"
sdred.BA$interc[sdred.BA$interc == "West"] <- "West (15 regions)"
dayplot <- ggplot()+
  geom_line(sdred.BA, mapping = aes(x = alpha, y = sdmean, group = interaction(sf_id, interc), 
                                    color = interc), alpha = 0.7)+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), 
       y = "Proportional Reduction in SD of Load", title = expression("Daily SD Reduction Mean according to "*alpha*" value"), color = "Area of region")+
  theme_minimal()+
  geom_line(sdred.agg, mapping = aes(x = alpha, y = sdmean, color = "Region-Weighted Average"), size = 1)+
  scale_color_manual(values = c(
    "Region-Weighted Average" = "black",
    "East (15 regions)" = "#619CFF",
    "West (15 regions)" = "#F8766D",
    "ERCOT (1 region)" = "#00BA38"))+
  theme(legend.position = "none",
        # c(.85,.15))+
        title = element_blank(),
        axis.text = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  ylim(-0.1, 1)
# ggsave(plot = dayplot, file = "92.daily_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")
#-----------------------------
#OVERALL/SEASON PT2 CALCULATIONS

sdred.BA.ov <- NULL
sdred.BA.ov <- sdred.t.ov %>% 
  group_by(alpha, sf_id, interc) %>% 
  summarise(sdmean = (sd(dem) - sd(sm))/sd(dem), avgdem = mean(dem), .groups = "drop_last") #overall sd
sdred.BAagg.ov <- sdred.BA.ov %>% 
  group_by(alpha) %>% 
  summarise(sdmean = weighted.mean(sdmean, w = avgdem), .groups = "drop_last") #mean of yearly SD red across (BA weighted); black line on graph

sdred.BAseas.ov <- NULL
sdred.BAseas.ov <- sdred.t.ov %>%
  group_by(alpha, season, sf_id, interc) %>% 
  summarise(sdmean = (sd(dem) - sd(sm))/sd(dem), avgdem = mean(dem), .groups = "drop_last") #seasonal sd (every season, every year) according to BA
sdred.BAseas.ov <- sdred.BAseas.ov %>% 
  group_by(alpha, season) %>% 
  summarise(sdmean = weighted.mean(sdmean, w = avgdem), .groups = "drop_last") #mean of seasonal spread (BA weighted avg); use for graph
sdred.BAaggseas.ov <- sdred.BAseas.ov %>% 
  group_by(alpha) %>% 
  summarise(sdmean = mean(sdmean), .groups = "drop_last") #unweighted mean of all seasons; black line on graph

#-----------------------------
#OVERALL/SEASON GRAPH

sdred.BA.ov$interc[sdred.BA.ov$interc == "Ercot"] <- "ERCOT (1 region)" #use all caps bc acronym
sdred.BA.ov$interc[sdred.BA.ov$interc == "East"] <- "East (15 regions)"
sdred.BA.ov$interc[sdred.BA.ov$interc == "West"] <- "West (15 regions)"
ovplot <- ggplot()+
  geom_line(sdred.BA.ov, mapping = aes(x = alpha, y = sdmean, group = interaction(sf_id, interc), 
                                       color = interc), alpha = 0.7)+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), y = "Overall SD Reduction Mean", title = expression("Overall SD Reduction Mean according to "*alpha*" value"), color = "Area of region")+
  theme_minimal()+
  geom_line(sdred.BAagg.ov, mapping = aes(x = alpha, y = sdmean, color = "Region-Weighted\nAverage"), size = 1)+
  scale_color_manual(values = c(
    "Region-Weighted\nAverage" = "black",
    "East (15 regions)" = "#619CFF",
    "West (15 regions)" = "#F8766D",
    "ERCOT (1 region)" = "#00BA38"))+
  # theme(legend.position = c(.85,.15),
  theme(legend.position = c(.8,.85),
        title = element_blank(),
        axis.text = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  ylim(-0.1, 1)
# ggsave(plot = ovplot, file = "92.overall_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")

ggplot()+ 
  geom_line(sdred.BAseas.ov, mapping = aes(x = alpha, y = sdmean, group = season, color = season), alpha = 0.5) + 
  geom_line(sdred.BAaggseas.ov, mapping = aes(x = alpha, y = sdmean, color = "Average\nAcross Seasons"), size = 1)+
  theme_minimal()+
  scale_color_manual(values = c(
    "Average\nAcross Seasons" = "black",
    "Winter" = "#C77CFF",
    "Spring" = "#7CAE00",
    "Summer" = "#00BFC4",
    "Fall" = "#F8766D"))+
  labs(x = expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), y = "Seasonal SD Reduction Mean", title = expression("Seasonal SD Reduction Mean according to "*alpha*" value across all regions"), color = "Season")+
  theme(legend.position = c(0.87, 0.18))
ggsave("92.season_sdred.pdf", width = 7.14, height = 6.46, dpi = 300, units = "in")
#-----------------------------
proc.time() - ptm

rm(sdred.t, sdred.t.ov,Xx)

combo <- cowplot::plot_grid(dayplot, ovplot, labels = "AUTO", ncol = 2, align = "h", 
                            label_size = 18, label_x = 0.1, label_y = 0.95)
x.grob <- textGrob(expression("Share of Temperature-Sensitive Load that is Flexible ("*alpha*")"), 
                   gp=gpar(fontsize = 16))
y.grob <- textGrob("Proportional Reduction in SD of Load", 
                   gp=gpar(fontsize = 16), rot=90)
title.grob <- textGrob(expression("Daily or Overall SD Reduction Mean according to "*alpha*" value"), 
                       gp=gpar(fontface="bold", fontsize = 20))
final <- grid.arrange(arrangeGrob(combo, bottom = x.grob), left = y.grob, top = title.grob)
ggsave(plot = final, file = "92.combo_sdred.pdf", width = 10.28, height = 5.81, dpi = 300, units = "in")



