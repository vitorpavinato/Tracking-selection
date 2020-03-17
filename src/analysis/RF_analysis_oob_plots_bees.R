##### manuscript plots

#### Load librarires
library(ggplot2)
library(cowplot)


## PLACERITA
##----------------------

### Load plot dataframes
## Selection
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.lgtgl.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.lgtps.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logthetaPS.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logNs.pla",".RData"))

## Demography
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logmeanNe2.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logncs.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logmeanNe2ncs.pla",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logtheta2.pla",".RData"))

### ALL OOB Plots (demography and selection)
###-----------------------------------------

## logit genetic load
my_breaks_1 <- c(20.085537, 7.389056, 2.718282, 1)
gl   <- ggplot(oob.lgtgl.pla, aes(x,y))
gl2d <- gl + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
                       ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  xlim(-22,6) +
  ylim(-22,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d

## logit Proportion of strongly selected mutations
my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
ps   <- ggplot(oob.lgtps.pla, aes(x,y))
ps2d <- ps + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
                       ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  xlim(-15,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d

## log10 Theta selected mutations
my_breaks_3 <- c(20.085537, 7.389056, 2.718282, 1)
thetaPS   <- ggplot(oob.logthetaPS.pla, aes(x,y))
thetaPS2d <- thetaPS + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d

## log10 Ns
my_breaks_4 <- c(20.085537, 7.389056, 2.718282, 1)
ns   <- ggplot(oob.logNs.pla, aes(x,y))
ns2d <- ns + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]*italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(OOB prediction)")))+
  xlim(-4,4) +
  ylim(-4,4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns2d

## log mean Ne2
my_breaks_5 <- c(54.598150, 7.389056, 1)
ne   <- ggplot(oob.logmeanNe2.pla, aes(x,y))
ne2d <- ne + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  xlim(-0.5,4.5) +
  ylim(-0.5,4.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d

## log NCS2
my_breaks_6 <- c(54.598150, 7.389056, 1)
ncs   <- ggplot(oob.logncs.pla, aes(x,y))
ncs2d <- ncs + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(OOB prediction)")))+
  xlim(-0.5,4.5) +
  ylim(-0.5,4.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d

## log Ne/NCS
my_breaks_7 <- c(2980.9557987, 403.42879, 54.598150, 7.389056, 1)
nencs   <- ggplot(oob.logmeanNe2ncs.pla, aes(x,y))
nencs2d <- nencs + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(OOB prediction)")))+
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d

## log theta2
my_breaks_8 <- c(54.598150, 7.389056, 1)
theta   <- ggplot(oob.logtheta2.pla, aes(x,y))
theta2d <- theta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_8),labels = round(my_breaks_8)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(OOB prediction)")))+
  xlim(-1.5,6) +
  ylim(-1.5,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta2d

cowplot::plot_grid(gl2d,
                   ps2d,
                   thetaPS2d,
                   ns2d,
                   ne2d,
                   ncs2d,
                   nencs2d,
                   theta2d,
                   nrow = 4,
                   labels = NULL)

## AVALON
##----------------------

### Load plot dataframes
## Selection
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.lgtgl.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.lgtps.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logthetaPS.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logNs.ava",".RData"))

## Demography
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logmeanNe2.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logncs.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logmeanNe2ncs.ava",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/oob.logtheta2.ava",".RData"))

### ALL OOB Plots (demography and selection)
###-----------------------------------------

## logit genetic load
my_breaks_1 <- c(20.085537, 7.389056, 2.718282, 1)
gl   <- ggplot(oob.lgtgl.ava, aes(x,y))
gl2d <- gl + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_1),labels = round(my_breaks_1)
  ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
  xlim(-22,6) +
  ylim(-22,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d

## logit Proportion of strongly selected mutations
my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
ps   <- ggplot(oob.lgtps.ava, aes(x,y))
ps2d <- ps + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
  ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
  xlim(-15,2) +
  ylim(-15,2) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d

## log10 Theta selected mutations
my_breaks_3 <- c(20.085537, 7.389056, 2.718282, 1)
thetaPS   <- ggplot(oob.logthetaPS.ava, aes(x,y))
thetaPS2d <- thetaPS + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")))+
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d

## log10 Ns
my_breaks_4 <- c(20.085537, 7.389056, 2.718282, 1)
ns   <- ggplot(oob.logNs.ava, aes(x,y))
ns2d <- ns + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]*italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(OOB prediction)")))+
  xlim(-4,4) +
  ylim(-4,4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns2d

## log mean Ne2
my_breaks_5 <- c(54.598150, 7.389056, 1)
ne   <- ggplot(oob.logmeanNe2.ava, aes(x,y))
ne2d <- ne + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
  xlim(-0.5,4.5) +
  ylim(-0.5,4.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d

## log NCS2
my_breaks_6 <- c(54.598150, 7.389056, 1)
ncs   <- ggplot(oob.logncs.ava, aes(x,y))
ncs2d <- ncs + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(OOB prediction)")))+
  xlim(-0.5,4.5) +
  ylim(-0.5,4.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d

## log Ne/NCS
my_breaks_7 <- c(403.42879, 54.598150, 7.389056, 1)
nencs   <- ggplot(oob.logmeanNe2ncs.ava, aes(x,y))
nencs2d <- nencs + 
  stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(OOB prediction)")))+
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d

## log theta2
my_breaks_8 <- c(54.598150, 7.389056, 1)
theta   <- ggplot(oob.logtheta2.ava, aes(x,y))
theta2d <- theta + 
  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_8),labels = round(my_breaks_8)
  ) +
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(OOB prediction)")))+
  xlim(-1.5,6) +
  ylim(-1.5,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta2d

cowplot::plot_grid(gl2d,
                   ps2d,
                   thetaPS2d,
                   ns2d,
                   ne2d,
                   ncs2d,
                   nencs2d,
                   theta2d,
                   nrow = 4,
                   labels = NULL)
