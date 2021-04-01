## loggammamean

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_loggammamean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_loggammamean",".RData"))

oob.loggammamean.Ava <- data.frame(x = list.loggammamean$Ava,
                                   y = list.reg.loggammamean$Ava$model.rf$predictions)

oob.loggammamean.Hum <- data.frame(x = list.loggammamean$Hum,
                                   y = list.reg.loggammamean$Hum$model.rf$predictions)

oob.loggammamean.Dav <- data.frame(x = list.loggammamean$Dav,
                                   y = list.reg.loggammamean$Dav$model.rf$predictions)

oob.loggammamean.Sta <- data.frame(x = list.loggammamean$Sta,
                                   y = list.reg.loggammamean$Sta$model.rf$predictions)

oob.loggammamean.Ste <- data.frame(x = list.loggammamean$Ste,
                                   y = list.reg.loggammamean$Ste$model.rf$predictions)

oob.loggammamean.Riv <- data.frame(x = list.loggammamean$Riv,
                                   y = list.reg.loggammamean$Riv$model.rf$predictions)

oob.loggammamean.Pla <- data.frame(x = list.loggammamean$Pla,
                                   y = list.reg.loggammamean$Pla$model.rf$predictions)
## AVALON
my_breaks_1 <- c(7.389056, 2.718282, 1)
gammamean.Ava   <- ggplot(oob.loggammamean.Ava, aes(x,y))
gammamean2d.Ava <- gammamean.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Ava

## HUMBOLDT
my_breaks_2 <- c(7.389056, 2.718282, 1)
gammamean.Hum   <- ggplot(oob.loggammamean.Hum, aes(x,y))
gammamean2d.Hum <- gammamean.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Hum

## DAVIS
my_breaks_3 <- c(7.389056, 2.718282, 1)
gammamean.Dav   <- ggplot(oob.loggammamean.Dav, aes(x,y))
gammamean2d.Dav <- gammamean.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Dav

## STANISLAUS
my_breaks_4 <- c(7.389056, 2.718282, 1)
gammamean.Sta   <- ggplot(oob.loggammamean.Sta, aes(x,y))
gammamean2d.Sta <- gammamean.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Sta

## STEBBINS
my_breaks_5 <- c(7.389056, 2.718282, 1)
gammamean.Ste   <- ggplot(oob.loggammamean.Ste, aes(x,y))
gammamean2d.Ste <- gammamean.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Ste

## RIVERSIDE
my_breaks_6 <- c(7.389056, 2.718282, 1)
gammamean.Riv   <- ggplot(oob.loggammamean.Riv, aes(x,y))
gammamean2d.Riv <- gammamean.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Riv

## PLACERITA
my_breaks_7 <- c(7.389056, 2.718282, 1)
gammamean.Pla   <- ggplot(oob.loggammamean.Pla, aes(x,y))
gammamean2d.Pla <- gammamean.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7),
  ) +
  xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
  ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-3,0) +
  ylim(-3,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammamean2d.Pla

cowplot::plot_grid(gammamean2d.Ava, gammamean2d.Hum,
                   gammamean2d.Dav, gammamean2d.Sta,
                   gammamean2d.Ste, gammamean2d.Riv,
                   gammamean2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
