## ThetaPS

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logthetaPS",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logthetaPS",".RData"))

oob.logthetaPS.Ava <- data.frame(x = list.logthetaPS$Ava,
                                 y = list.reg.logthetaPS$Ava$model.rf$predictions)

oob.logthetaPS.Hum <- data.frame(x = list.logthetaPS$Hum,
                                 y = list.reg.logthetaPS$Hum$model.rf$predictions)

oob.logthetaPS.Dav <- data.frame(x = list.logthetaPS$Dav,
                                 y = list.reg.logthetaPS$Dav$model.rf$predictions)

oob.logthetaPS.Sta <- data.frame(x = list.logthetaPS$Sta,
                                 y = list.reg.logthetaPS$Sta$model.rf$predictions)

oob.logthetaPS.Ste <- data.frame(x = list.logthetaPS$Ste,
                                 y = list.reg.logthetaPS$Ste$model.rf$predictions)

oob.logthetaPS.Riv <- data.frame(x = list.logthetaPS$Riv,
                                 y = list.reg.logthetaPS$Riv$model.rf$predictions)

oob.logthetaPS.Pla <- data.frame(x = list.logthetaPS$Pla,
                                 y = list.reg.logthetaPS$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Ava   <- ggplot(oob.logthetaPS.Ava, aes(x,y))
thetaPS2d.Ava <- thetaPS.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Ava

## HUMBOLDT
my_breaks_2 <- c(54.598150, 7.389056, 1)
thetaPS.Hum   <- ggplot(oob.logthetaPS.Hum, aes(x,y))
thetaPS2d.Hum <- thetaPS.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Hum

## DAVIS
my_breaks_3 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Dav   <- ggplot(oob.logthetaPS.Dav, aes(x,y))
thetaPS2d.Dav <- thetaPS.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Dav

## STANISLAUS
my_breaks_4 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Sta   <- ggplot(oob.logthetaPS.Sta, aes(x,y))
thetaPS2d.Sta <- thetaPS.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Sta

## STEBBINS
my_breaks_5 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Ste   <- ggplot(oob.logthetaPS.Ste, aes(x,y))
thetaPS2d.Ste <- thetaPS.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Ste

## RIVERSIDE
my_breaks_6 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Riv   <- ggplot(oob.logthetaPS.Riv, aes(x,y))
thetaPS2d.Riv <- thetaPS.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Riv

## PLACERITA
my_breaks_7 <- c(20.085537, 7.389056, 2.718282 , 1)
thetaPS.Pla   <- ggplot(oob.logthetaPS.Pla, aes(x,y))
thetaPS2d.Pla <- thetaPS.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.Pla

cowplot::plot_grid(thetaPS2d.Ava, thetaPS2d.Hum,
                   thetaPS2d.Dav, thetaPS2d.Sta,
                   thetaPS2d.Ste, thetaPS2d.Riv,
                   thetaPS2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
