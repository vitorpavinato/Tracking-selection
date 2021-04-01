## Logit Proportion of Strongly Selected Mutations

#### Load librarires
library(ggplot2)
library(cowplot)

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logitpopstrongmsel",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_reg_logitpopstrongmsel",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_zero_popstrongmsel",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_pods/list_posterior_zero_popstrongmsel",".RData"))

oob.lgtps.Ava <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Ava)), list.logitpopstrongmsel$Ava),
                            y = c(list.posterior.zero.popstrongmsel$Ava$expectation, list.reg.logitpopstrongmsel$Ava$model.rf$predictions))

oob.lgtps.Hum <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Hum)), list.logitpopstrongmsel$Hum),
                            y = c(list.posterior.zero.popstrongmsel$Hum$expectation, list.reg.logitpopstrongmsel$Hum$model.rf$predictions))

oob.lgtps.Dav <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Dav)), list.logitpopstrongmsel$Dav),
                            y = c(list.posterior.zero.popstrongmsel$Dav$expectation, list.reg.logitpopstrongmsel$Dav$model.rf$predictions))

oob.lgtps.Sta <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Sta)), list.logitpopstrongmsel$Sta),
                            y = c(list.posterior.zero.popstrongmsel$Sta$expectation, list.reg.logitpopstrongmsel$Sta$model.rf$predictions))

oob.lgtps.Ste <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Ste)), list.logitpopstrongmsel$Ste),
                            y = c(list.posterior.zero.popstrongmsel$Ste$expectation, list.reg.logitpopstrongmsel$Ste$model.rf$predictions))

oob.lgtps.Riv <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Riv)), list.logitpopstrongmsel$Riv),
                            y = c(list.posterior.zero.popstrongmsel$Riv$expectation, list.reg.logitpopstrongmsel$Riv$model.rf$predictions))

oob.lgtps.Pla <- data.frame(x = c(rep(-17.5, length(list.zero.popstrongmsel$Pla)), list.logitpopstrongmsel$Pla),
                            y = c(list.posterior.zero.popstrongmsel$Pla$expectation, list.reg.logitpopstrongmsel$Pla$model.rf$predictions))

## AVALON
my_breaks_1 <- c(54.598150, 7.389056, 1)
ps.Ava   <- ggplot(oob.lgtps.Ava, aes(x,y))
ps2d.Ava <- ps.Ava + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Avalon") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Ava

## HUMBOLDT
my_breaks_2 <- c(54.598150, 7.389056, 1)
ps.Hum   <- ggplot(oob.lgtps.Hum, aes(x,y))
ps2d.Hum <- ps.Hum + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Humboldt") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Hum

## DAVIS
my_breaks_3 <- c(54.598150, 7.389056, 1)
ps.Dav   <- ggplot(oob.lgtps.Dav, aes(x,y))
ps2d.Dav <- ps.Dav + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Davis") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Dav

## STANISLAUS
my_breaks_4 <- c(54.598150, 7.389056, 1)
ps.Sta   <- ggplot(oob.lgtps.Sta, aes(x,y))
ps2d.Sta <- ps.Sta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Stanislaus") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Sta

## STEBBINS
my_breaks_5 <- c(54.598150, 7.389056, 1)
ps.Ste   <- ggplot(oob.lgtps.Ste, aes(x,y))
ps2d.Ste <- ps.Ste + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Stebbins") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Ste

## RIVERSIDE
my_breaks_6 <- c(54.598150, 7.389056, 1)
ps.Riv   <- ggplot(oob.lgtps.Riv, aes(x,y))
ps2d.Riv <- ps.Riv + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Riverside") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Riv

## PLACERITA
my_breaks_7 <- c(54.598150, 7.389056, 1)
ps.Pla   <- ggplot(oob.lgtps.Pla, aes(x,y))
ps2d.Pla <- ps.Pla + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  ggtitle("Placerita") +
  xlim(-18,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.Pla

cowplot::plot_grid(ps2d.Ava, ps2d.Hum,
                   ps2d.Dav, ps2d.Sta,
                   ps2d.Ste, ps2d.Riv,
                   ps2d.Pla,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)
