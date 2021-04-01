## logmeanNe2

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmeanNe2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logmeanNe2",".RData"))

oob.logmeanNe2.Ava <- data.frame(x = list.logmeanNe2$Ava,
                                 y = list.reg.logmeanNe2$Ava$model.rf$predictions)

oob.logmeanNe2.Hum <- data.frame(x = list.logmeanNe2$Hum,
                                 y = list.reg.logmeanNe2$Hum$model.rf$predictions)

oob.logmeanNe2.Dav <- data.frame(x = list.logmeanNe2$Dav,
                                 y = list.reg.logmeanNe2$Dav$model.rf$predictions)

oob.logmeanNe2.Sta <- data.frame(x = list.logmeanNe2$Sta,
                                 y = list.reg.logmeanNe2$Sta$model.rf$predictions)

oob.logmeanNe2.Ste <- data.frame(x = list.logmeanNe2$Ste,
                                 y = list.reg.logmeanNe2$Ste$model.rf$predictions)

oob.logmeanNe2.Riv <- data.frame(x = list.logmeanNe2$Riv,
                                 y = list.reg.logmeanNe2$Riv$model.rf$predictions)

oob.logmeanNe2.Pla <- data.frame(x = list.logmeanNe2$Pla,
                                 y = list.reg.logmeanNe2$Pla$model.rf$predictions)
## StaLON
my_breaks_1 <- c(54.598150, 7.389056, 1)
ne.Ava   <- ggplot(oob.logmeanNe2.Sta, aes(x,y))
ne2d.Ava <- ne.Ava + 
            stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
            scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                                 ,
                                 breaks = round(my_breaks_1),labels = round(my_breaks_1)
            ) +
            xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
            ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
            ggtitle("Avalon") +
            xlim(-0.5,3.5) +
            ylim(-0.5,3.5) +
            geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
            theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                               legend.position = "right", axis.text = element_text(size = 12),
                               axis.title=element_text(size=12))
ne2d.Ava

## HUMBOLDT POPULATION
#my_breaks_1 <- c(54.598150, 7.389056, 1)
ne.Hum   <- ggplot(oob.logmeanNe2.Hum, aes(x,y))
ne2d.Hum <- ne.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Hum

## DAVIS POPULATION
ne.Dav   <- ggplot(oob.logmeanNe2.Dav, aes(x,y))
ne2d.Dav <- ne.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Dav

## STANISLAUS POPULATION
ne.Sta   <- ggplot(oob.logmeanNe2.Sta, aes(x,y))
ne2d.Sta <- ne.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Sta

## STEBBINS POPULATION
ne.Ste   <- ggplot(oob.logmeanNe2.Ste, aes(x,y))
ne2d.Ste <- ne.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Ste

## RIVERSIDE POPULATION
ne.Riv   <- ggplot(oob.logmeanNe2.Riv, aes(x,y))
ne2d.Riv <- ne.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Riv

## PLACERITA POPULATION
my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
ne.Pla   <- ggplot(oob.logmeanNe2.Pla, aes(x,y))
ne2d.Pla <- ne.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.Pla

cowplot::plot_grid(ne2d.Ava, ne2d.Hum,
                   ne2d.Dav, ne2d.Sta,
                   ne2d.Ste, ne2d.Riv,
                   ne2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
