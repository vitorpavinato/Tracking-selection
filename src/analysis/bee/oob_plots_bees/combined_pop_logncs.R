## logncs

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logncs",".RData"))

oob.logncs.Ava <- data.frame(x = list.logncs$Ava,
                             y = list.reg.logncs$Ava$model.rf$predictions)

oob.logncs.Hum <- data.frame(x = list.logncs$Hum,
                             y = list.reg.logncs$Hum$model.rf$predictions)

oob.logncs.Dav <- data.frame(x = list.logncs$Dav,
                             y = list.reg.logncs$Dav$model.rf$predictions)

oob.logncs.Sta <- data.frame(x = list.logncs$Sta,
                             y = list.reg.logncs$Sta$model.rf$predictions)

oob.logncs.Ste <- data.frame(x = list.logncs$Ste,
                             y = list.reg.logncs$Ste$model.rf$predictions)

oob.logncs.Riv <- data.frame(x = list.logncs$Riv,
                             y = list.reg.logncs$Riv$model.rf$predictions)

oob.logncs.Pla <- data.frame(x = list.logncs$Pla,
                             y = list.reg.logncs$Pla$model.rf$predictions)

## AVALON
my_breaks_1 <- c(54.598150, 7.389056, 1)
ncs.Ava   <- ggplot(oob.logncs.Ava, aes(x,y))
ncs2d.Ava <- ncs.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Ava

## HUMBOLDT POPULATION
ncs.Hum   <- ggplot(oob.logncs.Hum, aes(x,y))
ncs2d.Hum <- ncs.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Hum

## DAVIS POPULATION
ncs.Dav   <- ggplot(oob.logncs.Dav, aes(x,y))
ncs2d.Dav <- ncs.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Dav

## STANISLAUS POPULATION
ncs.Sta   <- ggplot(oob.logncs.Sta, aes(x,y))
ncs2d.Sta <- ncs.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Sta

## STEBBINS POPULATION
ncs.Ste   <- ggplot(oob.logncs.Ste, aes(x,y))
ncs2d.Ste <- ncs.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Ste

## RIVERSIDE POPULATION
ncs.Riv   <- ggplot(oob.logncs.Riv, aes(x,y))
ncs2d.Riv <- ncs.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Riv

## PLACERITA POPULATION
ncs.Pla   <- ggplot(oob.logncs.Pla, aes(x,y))
ncs2d.Pla <- ncs.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N))))) +
  ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.Pla

cowplot::plot_grid(ncs2d.Ava, ncs2d.Hum,
                   ncs2d.Dav, ncs2d.Sta,
                   ncs2d.Ste, ncs2d.Riv,
                   ncs2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
