## logmu

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmu",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logmu",".RData"))

oob.logmu.Ava <- data.frame(x = list.logmu$Ava,
                            y = list.reg.logmu$Ava$model.rf$predictions)

oob.logmu.Hum <- data.frame(x = list.logmu$Hum,
                            y = list.reg.logmu$Hum$model.rf$predictions)

oob.logmu.Dav <- data.frame(x = list.logmu$Dav,
                            y = list.reg.logmu$Dav$model.rf$predictions)

oob.logmu.Sta <- data.frame(x = list.logmu$Sta,
                            y = list.reg.logmu$Sta$model.rf$predictions)

oob.logmu.Ste <- data.frame(x = list.logmu$Ste,
                            y = list.reg.logmu$Ste$model.rf$predictions)

oob.logmu.Riv <- data.frame(x = list.logmu$Riv,
                            y = list.reg.logmu$Riv$model.rf$predictions)

oob.logmu.Pla <- data.frame(x = list.logmu$Pla,
                            y = list.reg.logmu$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(54.598150, 7.389056, 1)
mu.Ava   <- ggplot(oob.logmu.Ava , aes(x,y))
mu2d.Ava <- mu.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Ava

## HUMBOLDT POPULATION
mu.Hum   <- ggplot(oob.logmu.Hum , aes(x,y))
mu2d.Hum <- mu.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Hum

## DAVIS POPULATION
mu.Dav   <- ggplot(oob.logmu.Dav , aes(x,y))
mu2d.Dav <- mu.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Dav

## STANISLAUS POPULATION
mu.Sta   <- ggplot(oob.logmu.Sta , aes(x,y))
mu2d.Sta <- mu.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Sta

## STEBBINS POPULATION
mu.Ste   <- ggplot(oob.logmu.Ste , aes(x,y))
mu2d.Ste <- mu.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Ste

## RIVERSIDE POPULATION
mu.Riv   <- ggplot(oob.logmu.Riv , aes(x,y))
mu2d.Riv <- mu.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Riv

## PLACERITA POPULATION
mu.Pla   <- ggplot(oob.logmu.Pla , aes(x,y))
mu2d.Pla <- mu.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(mu))))) +
  ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-10,-6.5) +
  ylim(-10,-6.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
mu2d.Pla

cowplot::plot_grid(mu2d.Ava, mu2d.Hum,
                   mu2d.Dav, mu2d.Sta,
                   mu2d.Ste, mu2d.Riv,
                   mu2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
