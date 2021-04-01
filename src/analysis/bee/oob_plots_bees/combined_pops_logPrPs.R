## logPrPs

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logPrPs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logPrPs",".RData"))

oob.logPrPs.Ava <- data.frame(x = list.logPrPs$Ava,
                              y = list.reg.logPrPs$Ava$model.rf$predictions)

oob.logPrPs.Hum <- data.frame(x = list.logPrPs$Hum,
                              y = list.reg.logPrPs$Hum$model.rf$predictions)

oob.logPrPs.Dav <- data.frame(x = list.logPrPs$Dav,
                              y = list.reg.logPrPs$Dav$model.rf$predictions)

oob.logPrPs.Sta <- data.frame(x = list.logPrPs$Sta,
                              y = list.reg.logPrPs$Sta$model.rf$predictions)

oob.logPrPs.Ste <- data.frame(x = list.logPrPs$Ste,
                              y = list.reg.logPrPs$Ste$model.rf$predictions)

oob.logPrPs.Riv <- data.frame(x = list.logPrPs$Riv,
                              y = list.reg.logPrPs$Riv$model.rf$predictions)

oob.logPrPs.Pla <- data.frame(x = list.logPrPs$Pla,
                              y = list.reg.logPrPs$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(7.389056, 2.718282, 1)
prps.Ava   <- ggplot(oob.logPrPs.Ava, aes(x,y))
prps2d.Ava <- prps.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Ava

## HUMBOLDT
my_breaks_2 <- c(7.389056, 2.718282, 1)
prps.Hum   <- ggplot(oob.logPrPs.Hum, aes(x,y))
prps2d.Hum <- prps.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Hum

## DAVIS
my_breaks_3 <- c(7.389056, 2.718282, 1)
prps.Dav   <- ggplot(oob.logPrPs.Dav, aes(x,y))
prps2d.Dav <- prps.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Dav

## STANISLAUS
my_breaks_4 <- c(7.389056, 2.718282, 1)
prps.Sta   <- ggplot(oob.logPrPs.Sta, aes(x,y))
prps2d.Sta <- prps.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Sta

## STEBBINS
my_breaks_5 <- c(7.389056, 2.718282, 1)
prps.Ste   <- ggplot(oob.logPrPs.Ste, aes(x,y))
prps2d.Ste <- prps.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Ste

## RIVERSIDE
my_breaks_7 <- c(7.389056, 2.718282, 1)
prps.Riv   <- ggplot(oob.logPrPs.Riv, aes(x,y))
prps2d.Riv <- prps.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Riv

## PLACERITA
my_breaks_8 <- c(7.389056, 2.718282, 1)
prps.Pla   <- ggplot(oob.logPrPs.Pla, aes(x,y))
prps2d.Pla <- prps.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_8),labels = round(my_breaks_8)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
  ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-6,0) +
  ylim(-6,0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
prps2d.Pla

cowplot::plot_grid(prps2d.Ava, prps2d.Hum,
                   prps2d.Dav, prps2d.Sta,
                   prps2d.Ste, prps2d.Riv,
                   prps2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
