##### manuscript plots

rm(list=ls())
ls()

#### Load librarires
library(ggplot2)
library(cowplot)

### Load plot dataframes
## Selection
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.lgtgl",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.lgtps",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logthetaPS",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logNs",".RData"))

## Demography
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logmeanNe2",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logncs",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logmeanNe2ncs",".RData"))
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logtheta2",".RData"))

### ALL OOB Plots (demography and selection)
###-----------------------------------------
#par(mfrow=c(4, 2))
#par(mar=c(5,5,4,1)+.1)

## logit genetic load
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.lgtgl, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(L)))),
#       ylab = expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")),
#       xlim = c(-17,6),
#       ylim = c(-17,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_1 <- c(403.428793, 54.598150, 7.389056, 1)
gl   <- ggplot(oob.lgtgl, aes(x,y))
gl2d <- gl + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_1),labels = round(my_breaks_1)
                             ) +
        xlab(expression(paste("true", " ", logit(italic(L))))) +
        ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)")))+
        xlim(-17,6) +
        ylim(-17,6) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gl2d

## logit Proportion of strongly selected mutations
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.lgtps, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(P)))),
#       ylab = expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")),
#       xlim = c(-15,2),
#       ylim = c(-15,2),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_2 <- c(54.598150, 7.389056, 1)
ps   <- ggplot(oob.lgtps, aes(x,y))
ps2d <- ps + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,breaks = round(my_breaks_2),labels = round(my_breaks_2)
                             ) +
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
        xlim(-15,2) +
        ylim(-15,2) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 16),
                           axis.title=element_text(size=16))
ps2d

## log10 Theta selected mutations
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logthetaPS, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(θ) * italic(P)[S]))),
#       ylab = expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(OOB prediction)")),
#       xlim = c(-8,6),
#       ylim = c(-8,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_3 <- c(54.598150, 7.389056, 1)
thetaPS   <- ggplot(oob.logthetaPS, aes(x,y))
thetaPS2d <- thetaPS + 
             stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
             scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                                  breaks = round(my_breaks_3),labels = round(my_breaks_3)) +
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
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logNs, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]*italic(s)))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(OOB prediction)")),
#       xlim = c(-4,4),
#       ylim = c(-4,4),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_4 <- c(54.598150, 7.389056, 1)
ns   <- ggplot(oob.logNs, aes(x,y))
ns2d <- ns + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                             breaks = round(my_breaks_4),labels = round(my_breaks_4)) +
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
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logmeanNe2, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_5 <- c(403.428793, 54.598150, 7.389056, 1)
ne   <- ggplot(oob.logmeanNe2, aes(x,y))
ne2d <- ne + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,breaks = round(my_breaks_5),labels = round(my_breaks_5)
                             ) +
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
        xlim(-0.5,3.5) +
        ylim(-0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 16),
                           axis.title=element_text(size=16))
ne2d

## log NCS2
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logncs, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[cs]))),
#       ylab = expression(paste(log[10](italic(hat(N))[cs]), " ","(OOB prediction)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_6 <- c(403.428793, 54.598150, 7.389056, 1)
ncs   <- ggplot(oob.logncs, aes(x,y))
ncs2d <- ncs + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
                             breaks = round(my_breaks_6),labels = round(my_breaks_6)) +
        xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
        ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(OOB prediction)")))+
        xlim(-0.5,3.5) +
        ylim(-0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
ncs2d

## log Ne/NCS
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logmeanNe2ncs, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(OOB prediction)")),
#       xlim = c(-3.5,0.5),
#       ylim = c(-3.5,0.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_7 <- c(8103.08393, 403.42879, 20.08554, 1)
nencs   <- ggplot(oob.logmeanNe2ncs, aes(x,y))
nencs2d <- nencs + 
        stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts",
                             breaks = round(my_breaks_7),labels = round(my_breaks_7)) +
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
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logtheta2, nbins=100,
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(θ)))),
#       ylab = expression(paste(log[10](italic(hat(θ))), " ","(OOB prediction)")),
#       xlim = c(-1.5,6),
#       ylim = c(-1.5,6),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)
#
# ggplot2
#par(mar=c(5,5,4,1)+.1)
my_breaks_8 <- c(403.428793, 54.598150, 7.389056, 1)
theta   <- ggplot(oob.logtheta2, aes(x,y))
theta2d <- theta + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts",
                             breaks = round(my_breaks_8),labels = round(my_breaks_8)) +
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


### OOA Demography and selection (joint inference)
###-----------------------------------------

### Ne estimates based on temporal FST
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/estimated.fstNe2",".RData"))

#par(mfrow=c(1, 3))
#par(mar=c(5,5,4,1)+.1)

## logit Proportion of strongly selected mutations
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.lgtps, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", logit(italic(P)))),
#       ylab = expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")),
#       xlim = c(-15,2),
#       ylim = c(-15,2),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)

## log mean Ne2
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(oob.logmeanNe2, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)

## NE temporal FST
#hist2d plot
#par(mar=c(5,5,4,1)+.1)
#hist2d(estimated.fstNe2, nbins=100, 
#       col=viridis::viridis(32), 
#       FUN=function(x) log(length(x)),
#       xlab = expression(paste("true", " ", log[10](italic(N)[e]))),
#       ylab = expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")),
#       xlim = c(-0.5,3.5),
#       ylim = c(-0.5,3.5),
#       cex.axis = 1.2,
#       cex.lab  = 1.2)
#abline(a=0, b=1, col = "black", lty = 3)

my_breaks_9 <- c(403.428793, 54.598150, 7.389056, 1)
fstne   <- ggplot(estimated.fstNe2, aes(x,y))
fstne2d <- fstne + 
           stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
           scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                                ,breaks =round(my_breaks_9),labels = round(my_breaks_9)
                                ) +
           xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
           ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
           xlim(-0.5,3.5) +
           ylim(-0.5,3.5) +
           geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
           theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              legend.position = "right", axis.text = element_text(size = 16),
                              axis.title=element_text(size=16))
fstne2d

cowplot::plot_grid(ps2d,
                   ne2d,
                   fstne2d,
                   nrow = 1,
                   labels = NULL)

