
#Script for making East/West/ERCOT interconnect graphs of demand VS temp (2 types)
#1. season_temp.pdf: line graph of seasonal avg demand for each BA region VS histogram count of temp occurrences
#2. cdhhdh_temp.pdf: line graph of CDH/HDH avg demand for each BA region VS histogram count of CDH/HDH magnitude occurrence (overlayed)
#for loop: saves 3 graphs for each type
#dataset from EY & SZ; result of 6a.flattendata1.R

library(sf)
library(ggplot2) 
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(egg)
library(parallel)

mc = 1 #cores for parallel processing

#load data
eastSH <- st_read("east_15/east_15.shp") #from Sisi
westSH <- st_read("west_15/west_15.shp") #from Sisi
load("06apo_INTFLs_new.RData")

#-----------------------------
#SETUP
qe <- list()
for(i in 1:15){
  qe[[i]] <- ggplot()+
    geom_sf(data = eastSH, col = "gray87", fill = "gray87", alpha = 0.8)+
    geom_sf(data = eastSH[i,], col = "black", fill = "black", alpha = 0.8)+
    theme_void()
}
qw <- list()
for(i in 1:15){
  qw[[i]] <- ggplot()+
    geom_sf(data = westSH, col = "gray87", fill = "gray87", alpha = 0.8)+
    geom_sf(data = westSH[i,], col = "black", fill = "black", alpha = 0.8)+
    theme_void()
}


mega <- list(eastINTFL, westINTFL, ercotINTFL)
sorted <- mclapply(mega, function(a){
  order <- NULL
  for(i in 1:length(a)){
    df1 <- a[[i]]
    df <- data.frame(sf_id = i, avgtemp = mean(df1$temp, na.rm = TRUE), n = sum(!is.na(df1$temp)))
    order <- rbind(order, df)
  }
  n.order <- order(order$avgtemp) #lowest to highest temp
  return(n.order)
}, mc.cores = mc)
ex <- sorted[[1]]
wx <- sorted[[2]]
labelss <- LETTERS[1:15]

mega.w <- lapply(mega, function(a){
  ax <- mclapply(a, function(i){
    x = mean(i$demand, na.rm = TRUE)
    i$weightdem <- i$demand / x
    return(i)
  }, mc.cores = mc)
  return(ax)
})

rm(eastSH, westSH, mega, sorted) #save memory

names <- c("East", "West")

#-----------------------------
#***AVG SEASON DEMAND VS TEMP + COUNT***
#East & West
for(n in 1:2){ #1 is East; 2 is West
  aa <- mega.w[[n]]
  p <- list()
  
  for(i in 1:length(aa)){
    if(n==1) j = ex[i]
    if(n==2) j = wx[i]
    
    df1 = as.data.frame(aa[[j]])
    df1 <- df1[!is.na(df1$weightdem),] #remove NAs
    
    df1$roundtemp<- round(df1$temp)
    winter <- filter(df1, df1$month == "12" | df1$month == "1" | df1$month == "2")
    spring <- filter(df1, df1$month == "3" | df1$month == "4" | df1$month == "5")
    summer <- filter(df1, df1$month == "6" | df1$month == "7" | df1$month == "8")
    fall <- filter(df1, df1$month == "9" | df1$month == "10" | df1$month == "11")
    
    wintdemandtemp <- winter %>%
      group_by(roundtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
    sprdemandtemp <- spring %>%
      group_by(roundtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
    sumdemandtemp <- summer %>%
      group_by(roundtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
    falldemandtemp <- fall %>%
      group_by(roundtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
    
    a <- ggplot()+
      geom_line(wintdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Winter"))+
      geom_line(sprdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Spring"))+
      geom_line(sumdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Summer"))+
      geom_line(falldemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Fall"))+
      geom_vline(xintercept=c(18), linetype="dotted", color = "blue")+
      labs(x = "Temperature", y = "Demand Index")+
      theme_minimal()+
      theme(axis.title.x = element_blank(), #all no x axis or marker labels
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = unit(c(0.1,0.1,0,0.1), "cm"))+
      scale_color_manual(values = c(
        "Winter" = "#C77CFF",
        "Spring" = "#7CAE00",
        "Summer" = "#00BFC4",
        "Fall" = "#F8766D"))+
      theme(plot.title = element_text(size = 10, face = "bold")) +
      xlim(-27.5, 46.5)+ #xlim(-27.5, 46.5)+ -27.5, 46.5; -8.5,42.5
      ylim(0.5, 2.5)+
      # labs(title = "Ercot 2016 - 2018 Demand Index VS Temperature According to Season") #HERE
      annotate(geom = "text", x = -25, y = 2, label = labelss[i], size = 5)
    
    if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
       i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){ #right 2 columns no y axis or marker labels
      a <- a + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank())
    }
    if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){ #left column no y axis, but maintain marker labels
      a <- a + theme(axis.title.y = element_blank())
    }
    
    b <- ggplot(df1, aes(x = roundtemp))+
      geom_histogram(binwidth = 1, color = "blue", fill = "blue", alpha = 0.5)+
      geom_vline(xintercept=c(18), linetype="dotted", color = "blue")+
      labs(y = "Freq.", x = "Temperature")+
      theme_minimal()+
      theme(axis.title.y = element_text(colour = "blue", size = 8),
            axis.title.x = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(color = "blue", size = 8),
            plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) +
      xlim(-27.5, 46.5)+ # -8.5, 42.5
      ylim(0, 3000) #ylim(0, 3000)
    
    if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
       i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){ #right 2 columns no y axis or marker labels
      b <- b + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank())
    }
    if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){ #left column no y axis, but maintain marker labels
      b <- b + theme(axis.title.y = element_blank())
    }
    
    if(i=="1"|i=="2"|i=="3"|i=="4"|
       i=="6"|i=="7"|i=="8"|i=="9"|
       i=="11"|i=="12"|i=="13"|i=="14"){ #above bottom row no x axis or marker labels
      b <- b + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank())
    }
    if(i=="5"|i=="10"|i=="15"){ #bottom row no x axis, but maintain marker labels
      b <- b + theme(axis.title.x = element_blank())
    }
    
    
    p[[i]] <- plot_grid(a, b, align = "v",  ncol = 1, rel_heights = c(0.75, 0.5))
    
    if(n==1){ #insert east maps
      if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qe[[j]], x = 0.10, y = 0.15, width = 0.4, height = 0.39)
      }
      if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
         i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qe[[j]], x = 0.025, y = 0.15, width = 0.4, height = 0.39)
      }
    }
    
    if(n==2){ #insert west maps
      if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qw[[j]], x = 0.10, y = 0.15, width = 0.4, height = 0.39)
      }
      if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
         i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qw[[j]], x = 0.025, y = 0.15, width = 0.4, height = 0.39)
      }
      
    }
  }
  
  clump <- plot_grid( p[[1]],  p[[6]],  p[[11]], 
                      p[[2]],  p[[7]],  p[[12]],
                      p[[3]],  p[[8]],  p[[13]],
                      p[[4]],  p[[9]],  p[[14]],
                      p[[5]],  p[[10]], p[[15]],
                      ncol = 3, rel_widths = c(1.1,1,1), 
                      rel_heights = c(1,1,1,1,1.07))
  
  a <- a + theme(axis.text.x = element_text(size = 8),
                 axis.text.y = element_text(size = 8))
  ggdraw(cowplot::get_y_axis(a, position = "right"))
  
  #title & axis labeling (change title)
  title.grob <- textGrob(paste(names[n], "Interconnect 2016 - 2018 Demand Index VS Temperature According to Season"), 
                         gp=gpar(fontface="bold", fontsize=10))
  x.grob <- textGrob("Temperature", gp=gpar(fontsize=8))
  space <- textGrob(".", gp = gpar(col = "white")) #adding a space for padding for the y axis and vertical arrow
  clump <- grid.arrange(arrangeGrob(clump, top = title.grob, left = space, bottom = x.grob, right = space))
  clump1 <- ggdraw(clump) + draw_label("Demand Index", x = unit(0, "npc"), y = unit(0.54, "npc"), vjust = 1.5, angle = 90, color = "black", size = 8) +
    draw_label("Freq.", x = unit(0, "npc"), y = unit(0.45, "npc") , vjust = 1.5,  angle = 90, color = "blue", size = 8)
  
  #adding legend to bottom
  a <- a + theme(legend.position = "bottom", legend.title = element_text(size = 8)) + 
    labs(color = "Season")+
    theme(legend.box.margin = margin(0, 0, 10, 25))
  legendd <- cowplot::get_legend(a)
  clump2 <- plot_grid(clump1, legendd, ncol = 1, rel_heights = c(3, 0.2))
  plot_grid(p[[1]], legendd, ncol = 1, rel_heights = c(3, 0.1))
  
  #adding arrows & arrow labels to bottom and side
  aa1 <- linesGrob(x = unit(-0.3, "npc"), y = unit(c(0.6,0.4), "npc"), arrow = arrow(length = unit(3,"mm")),
                   gp=gpar(lwd = 1.5))
  aa2 <- linesGrob(x = unit(c(0.4,0.6), "npc"), y = unit(1.3, "npc"), arrow = arrow(length = unit(3,"mm")),
                   gp=gpar(lwd = 1.5))
  clump2 <- grid.arrange(arrangeGrob(clump2, bottom = aa2, right = aa1))
  arrow1 <- textGrob(label = "Hotter average temperature" , rot = -90, hjust = 0.6, vjust = 0.85, gp = gpar(fontsize = 8))
  arrow2 <- textGrob(label = "Hotter average temperature" , hjust = 0.55, vjust = 0.4, gp = gpar(fontsize = 8))
  
  #final graph
  final <- grid.arrange(arrangeGrob(clump2, bottom = arrow2, right = arrow1))
  
  ggsave(paste0("9.", names[n], "_season_temp.pdf"), plot = final, width = 8.43, height = 8.96, dpi = 300, units = "in", device = "pdf")
  
  rm(clump,clump1,clump2,a,legendd,final, title.grob,x.grob,space,aa1,aa2,arrow1,arrow2, p, b) #clearing up memory
}
#ERCOT
{
  df1 = as.data.frame(mega.w[[3]])
  df1 <- df1[!is.na(df1$weightdem),] #remove NAs
  
  df1$roundtemp<- round(df1$temp)
  winter <- filter(df1, df1$month == "12" | df1$month == "1" | df1$month == "2")
  spring <- filter(df1, df1$month == "3" | df1$month == "4" | df1$month == "5")
  summer <- filter(df1, df1$month == "6" | df1$month == "7" | df1$month == "8")
  fall <- filter(df1, df1$month == "9" | df1$month == "10" | df1$month == "11")
  
  wintdemandtemp <- winter %>%
    group_by(roundtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
  sprdemandtemp <- spring %>%
    group_by(roundtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
  sumdemandtemp <- summer %>%
    group_by(roundtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
  falldemandtemp <- fall %>%
    group_by(roundtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) #HERE
  
  a <- ggplot()+
    geom_line(wintdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Winter"))+
    geom_line(sprdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Spring"))+
    geom_line(sumdemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Summer"))+
    geom_line(falldemandtemp, mapping = aes(x = roundtemp, y = mean, color = "Fall"))+
    geom_vline(xintercept=c(18), linetype="dotted", color = "blue")+
    labs(x = "Temperature", y = "Demand Index", color = "Season")+
    theme_minimal()+
    theme(axis.title.x = element_blank(), #all no x axis or marker labels
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          plot.margin = unit(c(0.1,0.1,0,0.1), "cm"))+
    scale_color_manual(values = c(
      "Winter" = "#C77CFF",
      "Spring" = "#7CAE00",
      "Summer" = "#00BFC4",
      "Fall" = "#F8766D"))+
    theme(plot.title = element_text(size = 10, face = "bold")) +
    xlim(-8.5, 42.5)+ #xlim(-27.5, 46.5)+ -27.5, 46.5; -8.5,42.5
    # ylim(0.5, 2.5)+
    labs(title = "Ercot 2016 - 2018 Demand Index VS Temperature According to Season") #HERE
  # annotate(geom = "text", x = -25, y = 2, label = labelss[i], size = 5)
  
  b <- ggplot(df1, aes(x = roundtemp))+
    geom_histogram(binwidth = 1, color = "blue", fill = "blue", alpha = 0.5)+
    geom_vline(xintercept=c(18), linetype="dotted", color = "blue")+
    theme_minimal()+
    theme(axis.title.y = element_text(colour = "blue", size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(color = "blue", size = 8),
          plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) +
    xlim(-8.5, 42.5)+ #-8.5, 42.5 #-27.5, 46.5
    # ylim(0, 3000)
    labs(y = "Freq.", x = "Temperature")
  
  final <- plot_grid(a, b, align = "v",  ncol = 1, rel_heights = c(0.75, 0.5))
  a <- a + theme(legend.position = "bottom", legend.title = element_text(size = 8))+
    theme(legend.box.margin = margin(0, 0, 10, 25))
  legendd <- cowplot::get_legend(a)
  final <- plot_grid(final, legendd, ncol = 1, rel_heights = c(3, 0.2))
  
  ggsave("9.ERCOT_season_temp.pdf", plot = final, width = 5.88, height = 5.07, dpi = 300, units = "in", device = "pdf")
  rm(final,a,b,df1, winter,spring,summer,fall,wintdemandtemp,sprdemandtemp,sumdemandtemp,falldemandtemp)
}

#-----------------------------
gc()
#***AVG CDH/HDH VS TEMP + COUNT***

#East & West
for(n in 1:2){ #1 is East; 2 is West
  aa <- mega.w[[n]]
  p <- list()
  
  for(i in 1:length(aa)) {
    if(n==1) j = ex[i] #i is letter number; j is BA number
    if(n==2) j = wx[i]
    df1 = as.data.frame(aa[[j]])
    df1 <- df1[!is.na(df1$weightdem), ] #remove demand NAs
    df1 <- df1[!is.na(df1$temp), ] #remove temp NAs
    
    df1$ch <- ifelse(df1$temp > 18, "CDH", "HDH")
    df1$dhtemp <- round(ifelse(df1$ch == "CDH", df1$temp - 18, 18 - df1$temp))
    tempCDH <- filter(df1, ch == "CDH") %>%
      group_by(dhtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) ##########HERE
    tempHDH <- filter(df1, ch == "HDH") %>%
      group_by(dhtemp)%>%
      summarise(mean = mean(weightdem, na.rm = TRUE)) ##########HERE
    
    a2 <- ggplot()+
      geom_line(tempCDH, mapping = aes(x = dhtemp, y = mean, group = 1, color = "CDH" ))+
      geom_line(tempHDH, mapping = aes(x = dhtemp, y = mean, group = 1, color = "HDH"))+
      scale_color_manual(values = c(
        'CDH' = '#F8766D',
        'HDH' = '#00BFC4')) +
      labs(y = "Demand Index", color = "Type")+
      theme_minimal()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(face = "bold"),
            plot.margin = unit(c(0.1,0.1,0,0.1), "cm"))+
      # labs(title = "Ercot 2016 - 2018 Demand Index VS Cooling or Heating Degree Hours")+
      # labs(title = paste("BA Region", j, sep = " "))+
      theme(plot.title = element_text(size = 10)) +
      xlim(-0.5,  45.5)+ #-0.5, 45.5 #-0.5, 28.5
      ylim(0.5, 2.5)+
      annotate(geom = "text", x = 0.5, y = 2.2, label = labelss[i], size = 5)
    
    if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
       i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
      a2 <- a2 + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank())
    }
    if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
      a2 <- a2 + theme(axis.title.y = element_blank())
    }
    
    b2 <- ggplot()+
      geom_histogram(filter(df1, ch=="CDH"), mapping = aes(x = dhtemp), binwidth = 1, color = "firebrick4", alpha = 0.5, fill = "#F8766D")+
      geom_histogram(filter(df1, ch=="HDH"), mapping = aes(x = dhtemp), binwidth = 1, color = "firebrick4", alpha = 0.5, fill = "#00BFC4")+
      labs(y = "Freq.", x = "CDH or HDH: Temperature Measured Distance from 18 C")+
      theme_minimal()+
      theme(legend.position = "none",
            axis.title.y = element_text(color = "firebrick4", size = 8),
            axis.text.y = element_text(color = "firebrick4", size = 8),
            axis.title.x = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            plot.margin = unit(c(0,0.1,0.1,0.1), "cm"))+
      xlim(-0.5, 45.5)+ #45.5 #28.5
      ylim(0, 3050) #3050 #2050
    
    if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
       i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
      b2 <- b2 + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank())
    }
    if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
      b2 <- b2 + theme(axis.title.y = element_blank())
    }
    
    if(i=="1"|i=="2"|i=="3"|i=="4"|
       i=="6"|i=="7"|i=="8"|i=="9"|
       i=="11"|i=="12"|i=="13"|i=="14"){
      b2 <- b2 + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank())
    }
    if(i=="5"|i=="10"|i=="15"){
      b2 <- b2 + theme(axis.title.x = element_blank())
    }
    
    p[[i]] <- plot_grid(a2, b2, align = "v", ncol = 1, rel_heights = c(0.75, 0.5))
    
    
    if(n==1){ #EAST
      if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qe[[j]], x = 0.613, y = 0.15, width = 0.4, height = 0.39) #x = 0.10
      }
      if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
         i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qe[[j]], x = 0.6, y = 0.15, width = 0.4, height = 0.39) #x = 0.03
      }
    }
    if(n==2){ #WEST
      if(i=="1"|i=="2"|i=="3"|i=="4"|i=="5"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qw[[j]], x = 0.613, y = 0.15, width = 0.4, height = 0.39) #x = 0.10
      }
      if(i=="6"|i=="7"|i=="8"|i=="9"|i=="10"|
         i=="11"|i=="12"|i=="13"|i=="14"|i=="15"){
        p[[i]] <- ggdraw()+
          draw_plot(p[[i]])+
          draw_plot(qw[[j]], x = 0.6, y = 0.15, width = 0.4, height = 0.39) #x = 0.03
      }
    }
  }
  
  
  clump <- plot_grid( p[[1]], p[[6]], p[[11]], 
                      p[[2]], p[[7]], p[[12]],
                      p[[3]], p[[8]], p[[13]],
                      p[[4]], p[[9]], p[[14]],
                      p[[5]], p[[10]], p[[15]],
                      ncol = 3, rel_widths = c(1.1,1,1), 
                      rel_heights = c(1,1,1,1,1.07))
  
  
  x.grob <- textGrob("CDH or HDH: Temperature Measured Distance from 18 C", 
                     gp=gpar(fontsize=8))
  title.grob <- textGrob(paste(names[n], "2016 - 2018 Average Weighted Demand VS Cooling or Heating Degree Hours"), 
                         gp=gpar(fontface="bold", fontsize=10))
  space <- textGrob(".", gp = gpar(col = "white"))
  clump <- grid.arrange(arrangeGrob(clump, bottom = x.grob, top = title.grob, right = space, left = space))
  #yaxis
  clump1 <- ggdraw(clump) + draw_label("Demand Index", x = unit(0, "npc"), y = unit(0.54, "npc"), vjust = 1.5, angle = 90, color = "black", size = 8) +
    draw_label("Freq.", x = unit(0, "npc"), y = unit(0.45, "npc" ) , vjust = 1.5, angle = 90, color = "firebrick4", size = 8)
  
  #adding legend to bottom
  a2 <- a2 + theme(legend.position = "bottom", legend.title = element_text(size = 10))#+
  #theme(legend.box.margin = margin(0, 0, 10, 25))
  legendd <- cowplot::get_legend(a2)
  clump2 <- plot_grid(clump1, legendd, ncol = 1, rel_heights = c(3, 0.2))
  plot_grid(p[[1]], legendd, ncol = 1, rel_heights = c(3, 0.1))
  
  #adding arrows & arrow labels to bottom and side
  aa1 <- linesGrob(x = unit(-0.3, "npc"), y = unit(c(0.6,0.4), "npc"), arrow = arrow(length = unit(3,"mm")),
                   gp=gpar(lwd = 1.5))
  aa2 <- linesGrob(x = unit(c(0.4,0.6), "npc"), y = unit(1.3, "npc"), arrow = arrow(length = unit(3,"mm")),
                   gp=gpar(lwd = 1.5))
  clump2 <- grid.arrange(arrangeGrob(clump2, bottom = aa2, right = aa1))
  arrow1 <- textGrob(label = "Hotter average temperature" , rot = -90, hjust = 0.6, vjust = 0.85, gp = gpar(fontsize = 8))
  arrow2 <- textGrob(label = "Hotter average temperature" , hjust = 0.55, vjust = 0.4, gp = gpar(fontsize = 8))
  #final graph
  final <- grid.arrange(arrangeGrob(clump2, bottom = arrow2, right = arrow1))
  
  ggsave(paste0("9.", names[n], "_cdhhdh_temp.pdf"), plot = final, width = 8.43, height = 8.96, dpi = 300, units = "in", device = "pdf")
  
  rm(clump,clump1,clump2, p, a2, b2,legendd,final, title.grob,x.grob,space,aa1,aa2, arrow1, arrow2) #clearing up memory
}
#ERCOT
{
  df1 = as.data.frame(mega.w[[3]])
  df1 <- df1[!is.na(df1$weightdem), ] #remove demand NAs
  df1 <- df1[!is.na(df1$temp), ] #remove temp NAs
  
  df1$ch <- ifelse(df1$temp > 18, "CDH", "HDH")
  df1$dhtemp <- round(ifelse(df1$ch == "CDH", df1$temp - 18, 18 - df1$temp))
  tempCDH <- filter(df1, ch == "CDH") %>%
    group_by(dhtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) ##########HERE
  tempHDH <- filter(df1, ch == "HDH") %>%
    group_by(dhtemp)%>%
    summarise(mean = mean(weightdem, na.rm = TRUE)) ##########HERE
  
  a2 <- ggplot()+
    geom_line(tempCDH, mapping = aes(x = dhtemp, y = mean, group = 1, color = "CDH" ))+
    geom_line(tempHDH, mapping = aes(x = dhtemp, y = mean, group = 1, color = "HDH"))+
    scale_color_manual(values = c(
      'CDH' = '#F8766D',
      'HDH' = '#00BFC4')) +
    labs(y = "Demand Index", color = "Type")+
    theme_minimal()+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold"),
          plot.margin = unit(c(0.1,0.1,0,0.1), "cm"))+
    labs(title = "Ercot 2016 - 2018 Demand Index VS Cooling or Heating Degree Hours")+
    # labs(title = paste("BA Region", j, sep = " "))+
    theme(plot.title = element_text(size = 10)) +
    xlim(-0.5, 28.5) #-0.5, 45.5 #-0.5, 28.5
  # ylim(0.5, 2.5)+
  # annotate(geom = "text", x = 0.5, y = 2.2, label = labelss[i], size = 5)
  
  b2 <- ggplot()+
    geom_histogram(filter(df1, ch=="CDH"), mapping = aes(x = dhtemp), binwidth = 1, color = "firebrick4", alpha = 0.5, fill = "#F8766D")+
    geom_histogram(filter(df1, ch=="HDH"), mapping = aes(x = dhtemp), binwidth = 1, color = "firebrick4", alpha = 0.5, fill = "#00BFC4")+
    labs(y = "Freq.", x = "CDH or HDH: Temperature Measured Distance from 18 C")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_text(color = "firebrick4", size = 8),
          axis.text.y = element_text(color = "firebrick4", size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          plot.margin = unit(c(0,0.1,0.1,0.1), "cm"))+
    xlim(-0.5, 28.5) #45.5 #28.5
  # ylim(0, 3050) #3050 #2050
  
  final <- plot_grid(a2, b2, align = "v", ncol = 1, rel_heights = c(0.75, 0.5))
  a2 <- a2 + theme(legend.position = "bottom", legend.title = element_text(size = 10))#+
  #theme(legend.box.margin = margin(0, 0, 10, 25))
  legendd <- cowplot::get_legend(a2)
  final <- plot_grid(final, legendd, ncol = 1, rel_heights = c(3, 0.2))
  ggsave("91.ERCOT_cdhhdh_temp.pdf", plot = final, width = 5.88, height = 5.07, dpi = 300, units = "in", device = "pdf")
  rm(final,a2,b2, df1, tempCDH, tempHDH)
}
rm(aa, i, j, n, names, labelss, qe, qw, mega.w)
#-----------------------------


