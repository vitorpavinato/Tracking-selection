## logrr

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logrr",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logrr",".RData"))

oob.logrr.Ava <- data.frame(x = list.logrr$Ava,
                            y = list.reg.logrr$Ava$model.rf$predictions)

oob.logrr.Hum <- data.frame(x = list.logrr$Hum,
                            y = list.reg.logrr$Hum$model.rf$predictions)

oob.logrr.Dav <- data.frame(x = list.logrr$Dav,
                            y = list.reg.logrr$Dav$model.rf$predictions)

oob.logrr.Sta <- data.frame(x = list.logrr$Sta,
                            y = list.reg.logrr$Sta$model.rf$predictions)

oob.logrr.Ste <- data.frame(x = list.logrr$Ste,
                            y = list.reg.logrr$Ste$model.rf$predictions)

oob.logrr.Riv <- data.frame(x = list.logrr$Riv,
                            y = list.reg.logrr$Riv$model.rf$predictions)

oob.logrr.Pla <- data.frame(x = list.logrr$Pla,
                            y = list.reg.logrr$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(7.389056, 2.718282, 1)
rr.Ava   <- ggplot(oob.logrr.Ava, aes(x,y))
rr2d.Ava <- rr.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Ava

## HUMBOLDT POPULATION
rr.Hum   <- ggplot(oob.logrr.Hum, aes(x,y))
rr2d.Hum <- rr.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Hum

## DAVIS POPULATION
rr.Dav   <- ggplot(oob.logrr.Dav, aes(x,y))
rr2d.Dav <- rr.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Dav

## STANISLAUS POPULATION
rr.Sta   <- ggplot(oob.logrr.Sta, aes(x,y))
rr2d.Sta <- rr.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Sta

## STEBBINS POPULATION
rr.Ste   <- ggplot(oob.logrr.Ste, aes(x,y))
rr2d.Ste <- rr.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Ste

## RIVERSIDE POPULATION
rr.Riv   <- ggplot(oob.logrr.Riv, aes(x,y))
rr2d.Riv <- rr.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Riv

## PLACERITA POPULATION
rr.Pla   <- ggplot(oob.logrr.Pla, aes(x,y))
rr2d.Pla <- rr.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
  ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-8,-4) +
  ylim(-8,-4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
rr2d.Pla

cowplot::plot_grid(rr2d.Ava, rr2d.Hum,
                   rr2d.Dav, rr2d.Sta,
                   rr2d.Ste, rr2d.Riv,
                   rr2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)