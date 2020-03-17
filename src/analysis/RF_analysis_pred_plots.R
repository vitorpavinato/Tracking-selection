##### manuscript plots

rm(list=ls())
ls()

#### Load librarires
library(ggplot2)
library(cowplot)

### Load plot dataframes -- RANDOM PODS
## Selection
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logNs",".RData"))

## Demography
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logtheta2",".RData"))

# Estimated Ne - WF-ABC and Temporal FST
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.randompods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.randompods",".RData"))

### Load plot dataframes -- FIXED PODS (ORIGINAL + ADDITIONAL)
## Selection
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logNs.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logNs.fixedpods_add",".RData"))

## Demography
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logtheta2.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logtheta2.fixedpods_add",".RData"))

# Estimated Ne - WF-ABC and Temporal FST
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.fixedpods_add",".RData"))

load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.fixedpods",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.fixedpods_add",".RData"))

### ALL Prediction Plots (demography and selection) - RANDOM PODS
###--------------------------------------------------------------

#par(mfrow=c(4, 2))
#par(mar=c(5,5,4,1)+.1)

## logit genetic load
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.lgtgl, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(L)))),
#       ylab = expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")),
#       xlim = c(-17,6),
#       ylim = c(-17,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_1 <- c(54.598150, 7.389056, 1)
gl.random   <- ggplot(pred.lgtgl, aes(x,y))
gl2d.random <- gl.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",  name="counts",
                       breaks = round(my_breaks_1), labels = round(my_breaks_1)
                       ) +
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")))+
  xlim(-17.5,6) +
  ylim(-17.5,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl2d.random

## logit Proportion of strongly selected mutations
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.lgtps, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(P)))),
#       ylab = expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")),
#       xlim = c(-15,2),
#       ylim = c(-15,2),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
ps.random   <- ggplot(pred.lgtps, aes(x,y))
ps2d.random <- ps.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
                       ) +
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")))+
  xlim(-15.5,0.5) +
  ylim(-15.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ps2d.random

## log10 Theta selected mutations
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logthetaPS, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(θ) * italic(P)[S]))),
#       ylab = expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")),
#       xlim = c(-8,6),
#       ylim = c(-8,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_3 <- c(7.389056, 2.718282,1)
thetaPS.random   <- ggplot(pred.logthetaPS, aes(x,y))
thetaPS2d.random <- thetaPS.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")))+
  xlim(-8,6) +
  ylim(-8,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS2d.random

## log10 Ns
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logNs, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(N)[e]*italic(s)))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(ABC-RF)")),
#       xlim = c(-4,4),
#       ylim = c(-4,4),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_4 <- c(7.389056, 2.718282,1)
ns.random   <- ggplot(pred.logNs, aes(x,y))
ns2d.random <- ns.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]*italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(ABC-RF)")))+
  xlim(-4,4) +
  ylim(-4,4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns2d.random

## log mean Ne2
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logmeanNe2, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log(italic(hat(N))[e]), " ","(ABC-RF)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_5 <- c(54.598150, 7.389056, 1)
ne.random   <- ggplot(pred.logmeanNe2, aes(x,y))
ne2d.random <- ne.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne2d.random

## log NCS2
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logncs, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[cs]))),
#       ylab = expression(paste(log(italic(hat(N))[cs]), " ","(ABC-RF)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_6 <- c(54.598150, 7.389056, 1)
ncs.random   <- ggplot(pred.logncs, aes(x,y))
ncs2d.random <- ncs.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(ABC-RF)")))+
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs2d.random

## log Ne/NCS
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logmeanNe2ncs, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")),
#       xlim = c(-4,1),
#       ylim = c(-4,1),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_7 <- c(403.42879, 54.598150, 7.389056,1)
nencs.random   <- ggplot(pred.logmeanNe2ncs, aes(x,y))
nencs2d.random <- nencs.random + 
  stat_bin2d(bins=50,na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")))+
  xlim(-3.5,0.5) +
  ylim(-3.5,0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs2d.random

## log theta2
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(pred.logtheta2, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(θ)))),
#       ylab = expression(paste(log(italic(hat(θ))), " ","(ABC-RF)")),
#       xlim = c(-1.5,6),
#       ylim = c(-1.5,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_8 <- c(20.085537, 7.389056, 2.718282, 1)
theta.random   <- ggplot(pred.logtheta2, aes(x,y))
theta2d.random <- theta.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                       ,
                       breaks = round(my_breaks_8),labels = round(my_breaks_8)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(ABC-RF)")))+
  xlim(-1.5,6) +
  ylim(-1.5,6) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta2d.random

cowplot::plot_grid(gl2d.random,
                   ps2d.random,
                   thetaPS2d.random,
                   ns2d.random,
                   ne2d.random,
                   ncs2d.random,
                   nencs2d.random,
                   theta2d.random,
                   nrow = 4,
                   labels = NULL)

### ALL Prediction Plots (demography and selection) - FIXED PODS
###--------------------------------------------------------------

## SELECTION
##--------------------------------

## logit genetic load
## ORIGINAL s=0.1
gl.fix <- ggplot(pred.lgtgl.fixedpods, aes(x=x , y=y, color=c))
gl.fix <- gl.fix + 
  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")))+
  xlim(-17, 4) +
  ylim(-17, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl.fix

## ADDITIONAL s=0.01
gl.fix_add1 <- ggplot(pred.lgtgl.fixedpods_add[c(1:100, 201:300), ], aes(x=x , y=y, color=c))
gl.fix_add1 <- gl.fix_add1 + 
  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")))+
  xlim(-17, 4) +
  ylim(-17, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl.fix_add1

## ADDITIONAL s=0.25
gl.fix_add2 <- ggplot(pred.lgtgl.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
gl.fix_add2 <- gl.fix_add2 + 
  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(L))))) +
  ylab(expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")))+
  xlim(-17, 4) +
  ylim(-17, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
gl.fix_add2

## logit Proportion of strongly selected mutations
## ORIGINAL s=0.1
ps.fix <- ggplot(pred.lgtps.fixedpods, aes(x=x , y=y, color=c))
ps.fix <- ps.fix + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")))+
  xlim(-15, 0) +
  ylim(-15, 0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))

ps.fix 

## ADDITIONAL s=0.01
ps.fix_add1 <- ggplot(pred.lgtps.fixedpods_add[c(1:100, 201:300), ], aes(x=x , y=y, color=c))
ps.fix_add1 <- ps.fix_add1 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")))+
  xlim(-15, 0) +
  ylim(-15, 0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))

ps.fix_add1 

## ADDITIONAL s=0.25
ps.fix_add2 <- ggplot(pred.lgtps.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
ps.fix_add2 <- ps.fix_add2 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", logit(italic(P))))) +
  ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")))+
  xlim(-15, 0) +
  ylim(-15, 0) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))

ps.fix_add2 

## log10 Theta selected mutations
## ORIGINAL s=0.1
thetaPS.fix <- ggplot(pred.logthetaPS.fixedpods, aes(x=x , y=y, color=c))
thetaPS.fix <- thetaPS.fix + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")))+
  xlim(-4, 4) +
  ylim(-4, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS.fix

## ADDITIONAL s=0.01
thetaPS.fix_add1 <- ggplot(pred.logthetaPS.fixedpods_add[c(1:100, 201:300), ], aes(x=x , y=y, color=c))
thetaPS.fix_add1 <- thetaPS.fix_add1 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")))+
  xlim(-4, 4) +
  ylim(-4, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS.fix_add1

## ADDITIONAL s=0.25
thetaPS.fix_add2 <- ggplot(pred.logthetaPS.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
thetaPS.fix_add2 <- thetaPS.fix_add2 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("limited", "unlimited"),
                      labels = ""##c("mutation limited", "mutation unlimited"),
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")))+
  xlim(-4, 4) +
  ylim(-4, 4) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
thetaPS.fix_add2

## log10 Ns
## ORIGINAL s=0.1
ns.fix <- ggplot(pred.logNs.fixedpods, aes(x=x , y=y, color=c))
ns.fix <- ns.fix + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A", "B"),
                      labels = ""#c("mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e] * italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e] * italic(s)), " ","(ABC-RF)")))+
  xlim(0, 3) +
  ylim(0, 3) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns.fix

## ADDITIONAL s=0.01
ns.fix_add1 <- ggplot(pred.logNs.fixedpods_add[c(1:100, 201:300), ], aes(x=x , y=y, color=c))
ns.fix_add1 <- ns.fix_add1 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A", "B"),
                      labels = ""#c("mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e] * italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e] * italic(s)), " ","(ABC-RF)")))+
  xlim(0, 3) +
  ylim(0, 3) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns.fix_add1

## ADDITIONAL s=0.25
ns.fix_add2 <- ggplot(pred.logNs.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
ns.fix_add2 <- ns.fix_add2 + geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A", "B"),
                      labels = ""#c("mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e] * italic(s))))) +
  ylab(expression(paste(log[10](italic(hat(N))[e] * italic(s)), " ","(ABC-RF)")))+
  xlim(0, 3) +
  ylim(0, 3) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ns.fix_add2

## FIXED PODS SELECTION
##-------------------------------------------------
cowplot::plot_grid(gl.fix_add1, gl.fix, gl.fix_add2,
                   ps.fix_add1, ps.fix, ps.fix_add2,
                   thetaPS.fix_add1, thetaPS.fix, thetaPS.fix_add2,
                   ns.fix_add1, ns.fix, ns.fix_add2,
                   nrow = 4,
                   labels = NULL)


## DEMOGRAPHY
##--------------------------------

## log mean Ne2
## ORIGINAL s=0.1
ne.fix <- ggplot(pred.logmeanNe2.fixedpods, aes(x=x , y=y, color=c))
ne.fix <- ne.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = "",#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne.fix

## ADDITIONAL s=0.01
ne.fix_add1 <- ggplot(pred.logmeanNe2.fixedpods_add[c(1:100,101:200, 301:400), ], aes(x=x , y=y, color=c))
ne.fix_add1 <- ne.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = "",#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne.fix_add1

## ADDITIONAL s=0.25
ne.fix_add2 <- ggplot(pred.logmeanNe2.fixedpods_add[c(1:100,201:300, 401:500), ], aes(x=x , y=y, color=c))
ne.fix_add2 <- ne.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = "",#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ne.fix_add2

## log NCS2
## ORIGINAL s=0.1
ncs.fix <- ggplot(pred.logncs.fixedpods, aes(x=x , y=y, color=c))
ncs.fix <- ncs.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs.fix

## ADDITIONAL s=0.01
ncs.fix_add1 <- ggplot(pred.logncs.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
ncs.fix_add1 <- ncs.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs.fix_add1

## ADDITIONAL s=0.25
ncs.fix_add2 <- ggplot(pred.logncs.fixedpods_add[c(201:300, 401:500), ], aes(x=x , y=y, color=c))
ncs.fix_add2 <- ncs.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
ncs.fix_add2

## log Ne/NCS
## ORIGINAL s=0.1
nencs.fix <- ggplot(pred.logmeanNe2ncs.fixedpods, aes(x=x , y=y, color=c)) 
nencs.fix <- nencs.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")))+
  xlim(-0.9,0.9) +
  ylim(-0.9,0.9) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs.fix

## ADDITIONAL s=0.01
nencs.fix_add1 <- ggplot(pred.logmeanNe2ncs.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c)) 
nencs.fix_add1 <- nencs.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")))+
  xlim(-0.9,0.9) +
  ylim(-0.9,0.9) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs.fix_add1

## ADDITIONAL s=0.25
nencs.fix_add2 <- ggplot(pred.logmeanNe2ncs.fixedpods_add[c(201:300, 401:500), ], aes(x=x , y=y, color=c)) 
nencs.fix_add2 <- nencs.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")))+
  xlim(-0.9,0.9) +
  ylim(-0.9,0.9) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
nencs.fix_add2

## log theta2
## ORIGINAL s=0.1
theta.fix <- ggplot(pred.logtheta2.fixedpods, aes(x=x , y=y, color=c))
theta.fix <- theta.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta.fix

## ADDITIONAL s=0.01
theta.fix_add1 <- ggplot(pred.logtheta2.fixedpods_add[c(101:200, 301:400), ], aes(x=x , y=y, color=c))
theta.fix_add1 <- theta.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta.fix_add1

## ADDITIONAL s=0.25
theta.fix_add2 <- ggplot(pred.logtheta2.fixedpods_add[c(201:300, 401:500), ], aes(x=x , y=y, color=c))
theta.fix_add2 <- theta.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("B","C"),
                      labels = ""#c("mutation limited", "mutation unlimited")
  ) + 
  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
  ylab(expression(paste(log[10](italic(hat(θ))), " ","(ABC-RF)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
theta.fix_add2

## FIXED PODS DEMOGRAPHY
cowplot::plot_grid(ne.fix_add1, ne.fix, ne.fix_add2,
                   ncs.fix_add1, ncs.fix, ncs.fix_add2,
                   nencs.fix_add1, nencs.fix, nencs.fix_add2,
                   theta.fix_add1, theta.fix, theta.fix_add2,
                   nrow = 4,
                   labels = NULL)


### Prediction Demography and selection (joint inference)
###-----------------------------------------

## log Ne2 WF-ABC
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(estimated.wfabc.randompods, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log(italic(hat(N))[e]), " ","(WF-ABC)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_9 <- c(7.389056, 2.718282,1)
wfne.random   <- ggplot(estimated.wfabc.randompods, aes(x,y))
wfne2d.random <- wfne.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_9),labels = round(my_breaks_9)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
wfne2d.random

## log Ne2 FST
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(estimated.fstNe2.randompods, nbins=50, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_10 <- c(20.085537,7.389056, 2.718282,1)
fstne.random   <- ggplot(estimated.fstNe2.randompods, aes(x,y))
fstne2d.random <- fstne.random + 
  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                       ,
                       breaks = round(my_breaks_10),labels = round(my_breaks_10)
                       ) +
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
  xlim(-0.5,3.5) +
  ylim(-0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
fstne2d.random

## log Ne2 WF-ABC
## ORIGINAL s = 0.1
wfne.fix <- ggplot(estimated.wfabc.fixedpods, aes(x=x , y=y, color=c))
wfne.fix <- wfne.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
wfne.fix

## ADDITIONAL s = 0.01
wfne.fix_add1 <- ggplot(estimated.wfabc.fixedpods_add[c(1:100, 101:200, 301:400),], aes(x=x , y=y, color=c))
wfne.fix_add1 <- wfne.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
wfne.fix_add1

## ADDITIONAL s = 0.25
wfne.fix_add2 <- ggplot(estimated.wfabc.fixedpods_add[c(1:100, 201:300, 401:500),], aes(x=x , y=y, color=c))
wfne.fix_add2 <- wfne.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
wfne.fix_add2

## log Ne2 FST
## ORIGINAL s = 0.1
fstne.fix <- ggplot(estimated.fstNe2.fixedpods, aes(x=x , y=y, color=c))
fstne.fix <- fstne.fix +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
fstne.fix

## ADDITIONAL s = 0.01
fstne.fix_add1 <- ggplot(estimated.fstNe2.fixedpods_add[c(1:100, 101:200, 301:400),], aes(x=x , y=y, color=c))
fstne.fix_add1 <- fstne.fix_add1 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
fstne.fix_add1

## ADDITIONAL s = 0.25
fstne.fix_add2 <- ggplot(estimated.fstNe2.fixedpods_add[c(1:100, 201:300, 401:500),], aes(x=x , y=y, color=c))
fstne.fix_add2 <- fstne.fix_add2 +  geom_point(size=3, shape=1) + 
  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                      name   = "",#"Scenario",
                      breaks = "",#c("A","B","C"),
                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
                      ) + 
  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
  xlim(0.5,3.5) +
  ylim(0.5,3.5) +
  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "top", axis.text = element_text(size = 12),
                     axis.title=element_text(size=12))
fstne.fix_add2

cowplot::plot_grid(ne2d.random, wfne2d.random,fstne2d.random,
                   ne.fix_add1, wfne.fix_add1, fstne.fix_add1,
                   ne.fix, wfne.fix, fstne.fix,
                   ne.fix_add2, wfne.fix_add2, fstne.fix_add2,
                   nrow = 4,
                   labels = NULL)
