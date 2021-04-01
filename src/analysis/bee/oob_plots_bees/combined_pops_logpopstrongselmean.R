## logpopstrongselmean

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logpopstrongselmean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logpopstrongselmean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_zero_popstrongselmean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_posterior_zero_popstrongselmean",".RData"))

oob.logpopstrongselmean.Ava <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Ava)), list.logpopstrongselmean$Ava),
                                          y = c(list.posterior.zero.popstrongselmean$Ava$expectation, list.reg.logpopstrongselmean$Ava$model.rf$predictions))

oob.logpopstrongselmean.Hum <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Hum)), list.logpopstrongselmean$Hum),
                                          y = c(list.posterior.zero.popstrongselmean$Hum$expectation, list.reg.logpopstrongselmean$Hum$model.rf$predictions))

oob.logpopstrongselmean.Dav <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Dav)), list.logpopstrongselmean$Dav),
                                          y = c(list.posterior.zero.popstrongselmean$Dav$expectation, list.reg.logpopstrongselmean$Dav$model.rf$predictions))

oob.logpopstrongselmean.Sta <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Sta)), list.logpopstrongselmean$Sta),
                                          y = c(list.posterior.zero.popstrongselmean$Sta$expectation, list.reg.logpopstrongselmean$Sta$model.rf$predictions))

oob.logpopstrongselmean.Ste <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Ste)), list.logpopstrongselmean$Ste),
                                          y = c(list.posterior.zero.popstrongselmean$Ste$expectation, list.reg.logpopstrongselmean$Ste$model.rf$predictions))

oob.logpopstrongselmean.Riv <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Riv)), list.logpopstrongselmean$Riv),
                                          y = c(list.posterior.zero.popstrongselmean$Riv$expectation, list.reg.logpopstrongselmean$Riv$model.rf$predictions))

oob.logpopstrongselmean.Pla <- data.frame(x = c(rep(-4, length(list.zero.popstrongselmean$Pla)), list.logpopstrongselmean$Pla),
                                          y = c(list.posterior.zero.popstrongselmean$Pla$expectation, list.reg.logpopstrongselmean$Pla$model.rf$predictions))

## AVALON
my_breaks_1 <- c(54.598150, 7.389056, 1)
gammas.Ava   <- ggplot(oob.logpopstrongselmean.Ava, aes(x,y))
gammas2d.Ava <- gammas.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Ava

## HUMBOLDT
my_breaks_2 <- c(54.598150, 7.389056, 1)
gammas.Hum   <- ggplot(oob.logpopstrongselmean.Hum, aes(x,y))
gammas2d.Hum <- gammas.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Hum

## DAVIS
my_breaks_3 <- c(54.598150, 7.389056, 1)
gammas.Dav   <- ggplot(oob.logpopstrongselmean.Dav, aes(x,y))
gammas2d.Dav <- gammas.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Dav

## STANISLAUS
my_breaks_4 <- c(54.598150, 7.389056, 1)
gammas.Sta   <- ggplot(oob.logpopstrongselmean.Sta, aes(x,y))
gammas2d.Sta <- gammas.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Sta

## STEBBINS
my_breaks_5 <- c(54.598150, 7.389056, 1)
gammas.Ste   <- ggplot(oob.logpopstrongselmean.Ste, aes(x,y))
gammas2d.Ste <- gammas.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Ste

## RIVERSIDE
my_breaks_6 <- c(54.598150, 7.389056, 1)
gammas.Riv   <- ggplot(oob.logpopstrongselmean.Riv, aes(x,y))
gammas2d.Riv <- gammas.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Riv

## PLACERITA
my_breaks_7 <- c(403.428793, 54.598150, 7.389056, 1)
gammas.Pla   <- ggplot(oob.logpopstrongselmean.Pla, aes(x,y))
gammas2d.Pla <- gammas.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
  ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-4.5,1) +
  ylim(-3,1) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gammas2d.Pla

cowplot::plot_grid(gammas2d.Ava, gammas2d.Hum,
                   gammas2d.Dav, gammas2d.Sta,
                   gammas2d.Ste, gammas2d.Riv,
                   gammas2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)