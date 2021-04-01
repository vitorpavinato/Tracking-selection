## AVERAGE GENETIC LOAD

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_averageGenLoad",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_averageGenLoad",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_zero_averageGenLoad",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_posterior_zero_averageGenLoad",".RData"))

oob.lgtgl.Ava <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Ava)), list.averageGenLoad$Ava),
                            y = c(list.posterior.zero.averageGenLoad$Ava$expectation, list.reg.averageGenLoad$Ava$model.rf$predictions))

oob.lgtgl.Hum <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Hum)), list.averageGenLoad$Hum),
                            y = c(list.posterior.zero.averageGenLoad$Hum$expectation, list.reg.averageGenLoad$Hum$model.rf$predictions))

oob.lgtgl.Dav <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Dav)), list.averageGenLoad$Dav),
                            y = c(list.posterior.zero.averageGenLoad$Dav$expectation, list.reg.averageGenLoad$Dav$model.rf$predictions))

oob.lgtgl.Sta <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Sta)), list.averageGenLoad$Sta),
                            y = c(list.posterior.zero.averageGenLoad$Sta$expectation, list.reg.averageGenLoad$Sta$model.rf$predictions))

oob.lgtgl.Ste <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Ste)), list.averageGenLoad$Ste),
                            y = c(list.posterior.zero.averageGenLoad$Ste$expectation, list.reg.averageGenLoad$Ste$model.rf$predictions))

oob.lgtgl.Riv <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Riv)), list.averageGenLoad$Riv),
                            y = c(list.posterior.zero.averageGenLoad$Riv$expectation, list.reg.averageGenLoad$Riv$model.rf$predictions))

oob.lgtgl.Pla <- data.frame(x = c(rep(-37.5, length(list.zero.averageGenLoad$Pla)), list.averageGenLoad$Pla),
                            y = c(list.posterior.zero.averageGenLoad$Pla$expectation, list.reg.averageGenLoad$Pla$model.rf$predictions))

## AVALON
my_breaks_1 <- c(54.598150, 7.389056, 1)
gl.Ava   <- ggplot(oob.lgtgl.Ava, aes(x,y))
gl2d.Ava <- gl.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Ava

## HUMBOLDT
my_breaks_2 <- c(54.598150, 7.389056, 1)
gl.Hum   <- ggplot(oob.lgtgl.Hum, aes(x,y))
gl2d.Hum <- gl.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Hum

## Davis
my_breaks_3 <- c(54.598150, 7.389056, 1)
gl.Dav   <- ggplot(oob.lgtgl.Dav, aes(x,y))
gl2d.Dav <- gl.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Dav

## Stanislaus
my_breaks_4 <- c(54.598150, 7.389056, 1)
gl.Sta   <- ggplot(oob.lgtgl.Sta, aes(x,y))
gl2d.Sta <- gl.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Sta

## STEBBINS
my_breaks_5 <- c(54.598150, 7.389056, 1)
gl.Ste   <- ggplot(oob.lgtgl.Ste, aes(x,y))
gl2d.Ste <- gl.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Ste

## RIVERSIDE
my_breaks_6 <- c(54.598150, 7.389056, 1)
gl.Riv   <- ggplot(oob.lgtgl.Riv, aes(x,y))
gl2d.Riv <- gl.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Riv

## Placerita
my_breaks_7 <- c(54.598150, 7.389056, 1)
gl.Pla   <- ggplot(oob.lgtgl.Pla, aes(x,y))
gl2d.Pla <- gl.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-38,7) +
  ylim(-20,7) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.Pla

cowplot::plot_grid(gl2d.Ava, gl2d.Hum,
                   gl2d.Dav,gl2d.Sta,
                   gl2d.Ste,gl2d.Riv,
                   gl2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)

