## logmeanNe2ncs

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logmeanNe2ncs",".RData"))

oob.logmeanNe2ncs.Ava <- data.frame(x = list.logmeanNe2ncs$Ava,
                                    y = list.reg.logmeanNe2ncs$Ava$model.rf$predictions)

oob.logmeanNe2ncs.Hum <- data.frame(x = list.logmeanNe2ncs$Hum,
                                    y = list.reg.logmeanNe2ncs$Hum$model.rf$predictions)

oob.logmeanNe2ncs.Dav <- data.frame(x = list.logmeanNe2ncs$Dav,
                                    y = list.reg.logmeanNe2ncs$Dav$model.rf$predictions)

oob.logmeanNe2ncs.Sta <- data.frame(x = list.logmeanNe2ncs$Sta,
                                    y = list.reg.logmeanNe2ncs$Sta$model.rf$predictions)

oob.logmeanNe2ncs.Ste <- data.frame(x = list.logmeanNe2ncs$Ste,
                                    y = list.reg.logmeanNe2ncs$Ste$model.rf$predictions)

oob.logmeanNe2ncs.Riv <- data.frame(x = list.logmeanNe2ncs$Riv,
                                    y = list.reg.logmeanNe2ncs$Riv$model.rf$predictions)

oob.logmeanNe2ncs.Pla <- data.frame(x = list.logmeanNe2ncs$Pla,
                                    y = list.reg.logmeanNe2ncs$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(403.428793, 54.598150, 7.389056, 1)
nencs.Ava   <- ggplot(oob.logmeanNe2ncs.Ava, aes(x,y))
nencs2d.Ava <- nencs.Ava + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Ava

## HUMBOLDT POPULATION
my_breaks_2 <- c(2980.957987 ,403.428793, 54.598150, 7.389056, 1)
nencs.Hum   <- ggplot(oob.logmeanNe2ncs.Hum, aes(x,y))
nencs2d.Hum <- nencs.Hum + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Hum

## DAVIS POPULATION
nencs.Dav   <- ggplot(oob.logmeanNe2ncs.Dav, aes(x,y))
nencs2d.Dav <- nencs.Dav + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Dav

## STANISLAUS POPULATION
nencs.Sta   <- ggplot(oob.logmeanNe2ncs.Sta, aes(x,y))
nencs2d.Sta <- nencs.Sta + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Sta

## STEBBINS POPULATION
nencs.Ste   <- ggplot(oob.logmeanNe2ncs.Ste, aes(x,y))
nencs2d.Ste <- nencs.Ste + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Ste

## RIVERSIDE POPULATION
nencs.Riv   <- ggplot(oob.logmeanNe2ncs.Riv, aes(x,y))
nencs2d.Riv <- nencs.Riv + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Riv

## PLACERITA POPULATION
nencs.Pla   <- ggplot(oob.logmeanNe2ncs.Pla, aes(x,y))
nencs2d.Pla <- nencs.Pla + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.Pla

cowplot::plot_grid(nencs2d.Ava, nencs2d.Hum,
                   nencs2d.Dav, nencs2d.Sta,
                   nencs2d.Ste, nencs2d.Riv,
                   nencs2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
