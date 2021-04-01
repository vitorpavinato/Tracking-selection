########################################
##          Manuscript plots          ##
##          Proof of Concept          ##
########################################

rm(list=ls())
ls()

#### Load librarires
library(ggplot2)
library(ggpubr)
library(cowplot)

## OOB plots
##---------------------------------------

### Load plot dataframes

## Selection

# Log PrPs - Proportion of beneficial mutation - parameter
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logPrPs_c",".RData"))

# PrPs - Proportion of beneficial mutation - parameter
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logitPrPs_c",".RData"))

# Log gamma mean - parameter
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.loggammamean_c",".RData"))

# Logit Proportion of strongly selected mutation - latent variable + boxplot with neutral simulation predictions
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.lgtps_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps.neutralsims",".RData"))

# Log gamma mean - strongly selected mutation - latent variable + boxplot with neutral simulation predictions
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logpopstrongselmean_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logpopStrongSelMean.zeroValues",".RData"))

# Logit Average Genetic Load - latent variable + boxplot with neutral simulation predictions
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.lgtgl_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl.neutralsims",".RData"))

# Log Theta Selected Mutations
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logthetaPS_c",".RData"))

## Demography

# Log Per site mutation rate mu - parameter
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logmu_c",".RData"))

# Log recombination rate rr - paramter
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logrr_c",".RData"))

# Log Harmonic mean of the effective population sizes of period 2 - latent variable
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logmeanNe2_c",".RData"))

# Log Population size of period 2 - parameter 
load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logncs_c",".RData"))

# Log Ratio Mean effective size / population size of period 2 - MeanNe2/Ncs - latent variable
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/oob.logmeanNe2ncs_c",".RData"))

### Figures 

## Selection - scatterplots

# Log PrPs - Proportion of beneficial mutation - parameter
#par(mar=c(5,5,4,1)+.1)
my_breaks_1 <- c(20.085537, 7.389056, 2.718282, 1)
prps   <- ggplot(oob.logPrPs, aes(x,y))
prps2d <- prps + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_1),labels = round(my_breaks_1),
        ) +
        xlab(expression(paste("true", " ", log[10](italic(P)[R] * italic(P)[B])))) +
        ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
        xlim(-6,0) +
        ylim(-6,0) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
prps2d

# Logit PrPs - Proportion of beneficial mutation - parameter
#par(mar=c(5,5,4,1)+.1)
#my_breaks_1 <- c(54.598150, 7.389056, 1)
#prps   <- ggplot(oob.logitPrPs, aes(x,y))
#prps2d <- prps + 
#        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
#        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                             ,
#                             breaks = round(my_breaks_1),labels = round(my_breaks_1),
#        ) +
#        xlab(expression(paste("true", " ", logit(italic(P)[R] * italic(P)[B])))) +
#        ylab(expression(paste(logit(hat(italic(P)[R] * italic(P)[B])), " ","(OOB prediction)")))+
#        xlim(-15,3) +
#        ylim(-15,3) +
#        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                           panel.grid.major = element_blank(),
#                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                           legend.position = "right", axis.text = element_text(size = 12),
#                           axis.title=element_text(size=12))
#prps2d

# Log gamma mean - parameter
#par(mar=c(5,5,4,1)+.1)
my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
gammamean   <- ggplot(oob.loggammamean, aes(x,y))
gammamean2d <- gammamean + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_2),labels = round(my_breaks_2),
        ) +
        xlab(expression(paste("true", " ", log[10](italic(gamma))))) +
        ylab(expression(paste(log[10](italic(hat(gamma))), " ","(OOB prediction)")))+
        xlim(-3,0) +
        ylim(-3,0) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gammamean2d

# Logit Proportion of strongly selected mutation - latent variable
#par(mar=c(5,5,4,1)+.1)
my_breaks_3 <- c(54.598150, 7.389056, 1)
ps   <- ggplot(oob.lgtps, aes(x,y))
ps2d <- ps + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,breaks = round(my_breaks_3),labels = round(my_breaks_3)
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

# Log gamma mean - strongly selected mutation - latent variable
#par(mar=c(5,5,4,1)+.1)
my_breaks_4 <- c(54.598150, 7.389056, 1)
strongSelMean   <- ggplot(oob.logpopstrongselmean, aes(x,y))
strongSelMean2d <- strongSelMean + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_4),labels = round(my_breaks_4),
        ) +
        xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
        ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
        xlim(-3,1) +
        ylim(-3,1) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
strongSelMean2d


# Logit Average Genetic Load - latent variable
#par(mar=c(5,5,4,1)+.1)
my_breaks_5 <- c(403.428793, 54.598150, 7.389056, 1)
gl   <- ggplot(oob.lgtgl, aes(x,y))
gl2d <- gl + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_5),labels = round(my_breaks_5),
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

# Log Theta Selected Mutations
#par(mar=c(5,5,4,1)+.1)
my_breaks_6 <- c(54.598150, 7.389056, 1)
thetaPS   <- ggplot(oob.logthetaPS, aes(x,y))
thetaPS2d <- thetaPS + 
             stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
             scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                                  ,
                                  breaks = round(my_breaks_6),labels = round(my_breaks_6)
                                  ) +
             xlab(expression(paste("true", " ", log[10](italic(θ)[b])))) +
             ylab(expression(paste(log[10](italic(hat(θ))[b]), " ","(OOB prediction)")))+
             xlim(-8,6) +
             ylim(-8,6) +
             geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
             theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 16),
                           axis.title=element_text(size=16))

thetaPS2d

## Boxplots for the estimated values for the neutral simulations

# Logit Proportion of strongly selected mutations - neutral simulations
neutral.p.plot <- ggplot(pred.lgtps.neutralsims, aes(x=col, y=est_lgtp)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        ylab(expression(paste(logit(italic(hat(P))), " ","(OOB prediction)")))+
        ylim(-15,2) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.p.plot

# Log gamma mean of strongly selected mutations - neutral simulations
zerovalues.strongSelMean.plot <- ggplot(pred.logpopStrongSelMean.zeroValues, aes(x=col, y=est_logsmean)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", log[10](italic(bar(s)))))) +
        ylab(expression(paste(log[10](hat(italic(bar(s)))), " ","(OOB prediction)")))+
        ylim(-3,1) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

zerovalues.strongSelMean.plot


# Logit Average Genetic load - neutral simulations
neutral.l.plot <- ggplot(pred.lgtgl.neutralsims, aes(x=col, y=est_lgtl)) + 
                  geom_boxplot() + 
                  #xlab("") +
                  #ylab("") +
                  xlab(expression(paste("true", " ", logit(italic(L))))) +
                  ylab(expression(paste(logit(italic(hat(L))), " ","(OOB prediction)"))) +
                  ylim(-17,6) +
                  #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
                  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     legend.position = "right", axis.text = element_text(size = 12),
                                     axis.title=element_text(size=12))
neutral.l.plot



box_plots <- cowplot::plot_grid(neutral.p.plot, neutral.l.plot, zerovalues.strongSelMean.plot,
                                nrow = 1, ncol = 3,
                                labels = NULL)

cowplot::plot_grid(prps2d, ps2d,
                   gammamean2d, strongSelMean2d,
                   gl2d, thetaPS2d,
                   box_plots,
                   nrow = 4,
                   ncol = 2,
                   labels = NULL)

## Demography - scatterplots

# Log Per site mutation rate mu - parameter
#par(mar=c(5,5,4,1)+.1)
my_breaks_7 <- c(54.598150, 7.389056, 1)
mu   <- ggplot(oob.logmu, aes(x,y))
mu2d <- mu + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                             ,
                             breaks = round(my_breaks_7),labels = round(my_breaks_7)
        ) +
        xlab(expression(paste("true", " ", log[10](italic(mu))))) +
        ylab(expression(paste(log[10](italic(hat(mu))), " ","(OOB prediction)")))+
        xlim(-10,-6) +
        ylim(-10,-6) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
mu2d

# Log recombination rate rr - paramter
my_breaks_8 <- c(20.085537, 7.389056, 2.718282, 1)
rr   <- ggplot(oob.logrr, aes(x,y))
rr2d <- rr + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                             ,
                             breaks = round(my_breaks_8),labels = round(my_breaks_8)
        ) +
        xlab(expression(paste("true", " ", log[10](italic(c)[0])))) +
        ylab(expression(paste(log[10](hat(italic(c)[0])), " ","(OOB prediction)")))+
        xlim(-9.4,-6.4) +
        ylim(-9.4,-6.4) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
rr2d

# Log Harmonic mean of the effective population sizes of period 2 - latent variable
#par(mar=c(5,5,4,1)+.1)
my_breaks_9 <- c(2980.957987,403.428793, 54.598150, 7.389056, 1)
ne   <- ggplot(oob.logmeanNe2, aes(x,y))
ne2d <- ne + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,breaks = round(my_breaks_9),labels = round(my_breaks_9)
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

# Log Population size of period 2 - parameter 
#par(mar=c(5,5,4,1)+.1)
my_breaks_10 <- c(403.428793, 54.598150, 7.389056, 1)
ncs   <- ggplot(oob.logncs, aes(x,y))
ncs2d <- ncs + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                             ,
                             breaks = round(my_breaks_10),labels = round(my_breaks_10)
                             ) +
        xlab(expression(paste("true", " ", log[10](italic(N))))) +
        ylab(expression(paste(log[10](italic(hat(N))), " ","(OOB prediction)")))+
        xlim(-0.5,3.5) +
        ylim(-0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 16),
                           axis.title=element_text(size=16))
ncs2d

# Log Ratio Mean effective size / population size of period 2 - MeanNe2/Ncs - latent variable
#par(mar=c(5,5,4,1)+.1)
my_breaks_11 <- c(8103.08393, 403.42879, 20.08554, 1)
nencs   <- ggplot(oob.logmeanNe2ncs, aes(x,y))
nencs2d <- nencs + 
        stat_bin2d(bins=100,na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                             ,
                             breaks = round(my_breaks_11),labels = round(my_breaks_11)
                             ) +
        xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)), " ","(OOB prediction)")))+
        xlim(-3.5,0.5) +
        ylim(-3.5,0.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
nencs2d

## OOB Demography and selection (joint inference)
##-----------------------------------------

## Ne estimates based on temporal FST
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob_tables/estimated.fstNe2_c",".RData"))

my_breaks_fstnerf <- c(403.428793, 54.598150, 7.389056, 1)
fstne   <- ggplot(estimated.fstNe2, aes(x,y))
fstne2d <- fstne + 
           stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
           scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                                ,breaks =round(my_breaks_fstnerf),labels = round(my_breaks_fstnerf)
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

# Joint Inference of Demography and Selection
join_demo_sel_plot <- cowplot::plot_grid(thetaPS2d,ncs2d, 
                                         ne2d, fstne2d,
                                         nrow = 2,
                                         ncol = 2,
                                         labels=""
                                         #,
                                         #align = "hv"
                                         )

ggsave(join_demo_sel_plot, filename = "results/pipeline_v5/random_forests/join_demo_sel.pdf", device = cairo_pdf, 
       width = 13, height = 9, units = "in", dpi = "retina")

## Variable importance plots
##-----------------------------------------

par(mfrow=c(2, 4))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
plot(x = reg_averageGenLoad_2, n.var = 20, main=expression(logit(italic(L))), cex.axis = 0.5)

# logit Proportion of strongly selected mutations
plot(x = reg_logitpopstrongmsel, n.var = 20, main=expression(logit(italic(P))), cex.axis = 0.5)

# log10 Theta selected mutations
plot(x = reg_logthetaPS, n.var = 20, main=expression(log[10](italic(θ)[b])), cex.axis = 0.5)

# log mean Ne2
plot(x = reg_logmeanNe2, n.var = 20, main=expression(log[10](italic(N)[e])), cex.axis = 0.5)

# log N
plot(x = reg_logncs, n.var = 20, main=expression(log[10](italic(N))), cex.axis = 0.5)

# log Ne2/N
plot(x = reg_logmeanNe2ncs, n.var = 20, main=expression(log[10](italic(N)[e]/italic(N))), cex.axis = 0.5)

# log mu
plot(x = reg_logmu, n.var = 20, main=expression(log[10](italic(mu))), cex.axis = 0.5)

## Histograms Priors
##-----------------------------------------

par(mfrow=c(4, 2))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
hist(logitaverageGenload, freq = TRUE, xlab = expression(logit(italic(L))), main = "", col = "#bdbdbd")

# logit Proportion of strongly selected mutations
hist(logitpopstrongmsel, freq = TRUE, xlab = expression(logit(italic(P))), main = "", col = "#bdbdbd")

# log10 Theta selected mutations
hist(logthetaPS, freq = TRUE, xlab = expression(log[10](italic(θ)[b])), main = "", col = "#bdbdbd")

# log mean Ne2
hist(logmeanNe2, freq = TRUE, xlab = expression(log[10](italic(N)[e])) , main = "", col = "#bdbdbd")

# log N
hist(logncs, freq = TRUE, xlab = expression(log[10](italic(N))), main = "", col = "#bdbdbd")

# log Ne2/NCS
hist(logmeanNe2ncs, freq = TRUE, xlab = expression(log[10](italic(N)[e]/italic(N))), main = "", col = "#bdbdbd")

# log mu
hist(logmu, freq = TRUE, xlab = expression(log[10](italic(mu))), main = "", col = "#bdbdbd")

## Predition PODS plots
##---------------------------------------

### Load plot dataframes -- True vs predicted for RANDOM PODS
## Selection
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logNs",".RData"))
#
## Demography
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logtheta2",".RData"))
#
# Estimated Ne - WF-ABC and Temporal FST
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.randompods",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.randompods",".RData"))

### Load plot dataframes -- True vs predicted for FIXED PODS - UPDATE [17-03-20] with combined reftable
## Selection
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logPrPs.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logPrPs.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps.fixedpods_c_2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtps.fixedpods_add_c_2",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.loggammamean.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.loggammamean.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logpopStrongSelMean.fixedpods_c_2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logpopStrongSelMean.fixedpods_add_c_2",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS.fixedpods_c_2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logthetaPS.fixedpods_add_c_2",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl.fixedpods_c_2",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.lgtgl.fixedpods_add_c_2",".RData"))

## Demography
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmu.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmu.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logrr.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logrr.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logncs.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2.fixedpods_add_c",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs.fixedpods_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/pred.logmeanNe2ncs.fixedpods_add_c",".RData"))

# Estimated Ne - WF-ABC and Temporal FST
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.fixedpods",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.fixedpods_add",".RData"))

#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.fixedpods",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.fstNe2.fixedpods_add",".RData"))

# Estimated Ne - WF-ABC and Temporal FST
# subset of random PODs - original analysis
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/posterior_pods/estimated.wfabc.randompods",".RData"))

# All simulations - used to train the RF
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/oob.logmeanNe2_c",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v5/random_forests/estimated.fstNe2_c",".RData"))

## ALL Prediction Plots (demography and selection) - RANDOM PODS
#par(mfrow=c(4, 2))
#par(mar=c(5,5,4,1)+.1)

# logit genetic load
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
#my_breaks_1 <- c(54.598150, 7.389056, 1)
#gl.random   <- ggplot(pred.lgtgl, aes(x,y))
#gl2d.random <- gl.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",  name="counts",
#                       breaks = round(my_breaks_1), labels = round(my_breaks_1)
#                       ) +
#  xlab(expression(paste("true", " ", logit(italic(L))))) +
#  ylab(expression(paste(logit(italic(hat(L))), " ","(ABC-RF)")))+
#  xlim(-17.5,6) +
#  ylim(-17.5,6) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#gl2d.random
#
# logit Proportion of strongly selected mutations
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
#my_breaks_2 <- c(20.085537, 7.389056, 2.718282, 1)
#ps.random   <- ggplot(pred.lgtps, aes(x,y))
#ps2d.random <- ps.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts",
#                       breaks = round(my_breaks_2),labels = round(my_breaks_2)
#                       ) +
#  xlab(expression(paste("true", " ", logit(italic(P))))) +
#  ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)")))+
#  xlim(-15.5,0.5) +
#  ylim(-15.5,0.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ps2d.random
#
# log10 Theta selected mutations
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
#my_breaks_3 <- c(7.389056, 2.718282,1)
#thetaPS.random   <- ggplot(pred.logthetaPS, aes(x,y))
#thetaPS2d.random <- thetaPS.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,
#                       breaks = round(my_breaks_3),labels = round(my_breaks_3)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
#  ylab(expression(paste(log[10](italic(hat(θ))*italic(P)[S]), " ","(ABC-RF)")))+
#  xlim(-8,6) +
#  ylim(-8,6) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#thetaPS2d.random
#
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
#my_breaks_4 <- c(7.389056, 2.718282,1)
#ns.random   <- ggplot(pred.logNs, aes(x,y))
#ns2d.random <- ns.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,
#                       breaks = round(my_breaks_4),labels = round(my_breaks_4)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[e]*italic(s))))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]*italic(s)), " ","(ABC-RF)")))+
#  xlim(-4,4) +
#  ylim(-4,4) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ns2d.random
#
# log mean Ne2
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
#my_breaks_5 <- c(54.598150, 7.389056, 1)
#ne.random   <- ggplot(pred.logmeanNe2, aes(x,y))
#ne2d.random <- ne.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,
#                       breaks = round(my_breaks_5),labels = round(my_breaks_5)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ne2d.random
#
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
#my_breaks_6 <- c(54.598150, 7.389056, 1)
#ncs.random   <- ggplot(pred.logncs, aes(x,y))
#ncs2d.random <- ncs.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,
#                       breaks = round(my_breaks_6),labels = round(my_breaks_6)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[cs])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[cs]), " ","(ABC-RF)")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ncs2d.random
#
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
#my_breaks_7 <- c(403.42879, 54.598150, 7.389056,1)
#nencs.random   <- ggplot(pred.logmeanNe2ncs, aes(x,y))
#nencs2d.random <- nencs.random + 
#  stat_bin2d(bins=50,na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
#                       ,
#                       breaks = round(my_breaks_7),labels = round(my_breaks_7)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N)[cs])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]/italic(N)[cs]), " ","(ABC-RF)")))+
#  xlim(-3.5,0.5) +
#  ylim(-3.5,0.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#nencs2d.random
#
# log theta2
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
#my_breaks_8 <- c(20.085537, 7.389056, 2.718282, 1)
#theta.random   <- ggplot(pred.logtheta2, aes(x,y))
#theta2d.random <- theta.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
#                       ,
#                       breaks = round(my_breaks_8),labels = round(my_breaks_8)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(θ))))) +
#  ylab(expression(paste(log[10](italic(hat(θ))), " ","(ABC-RF)")))+
#  xlim(-1.5,6) +
#  ylim(-1.5,6) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#theta2d.random
#
#cowplot::plot_grid(gl2d.random,
#                   ps2d.random,
#                   thetaPS2d.random,
#                   ns2d.random,
#                   ne2d.random,
#                   ncs2d.random,
#                   nencs2d.random,
#                   theta2d.random,
#                   nrow = 4,
#                   labels = NULL)

### FIXED PODS - SELECTION

## log10 PrPb
# ORIGINAL s=0.1
PrPs.fix <- ggplot(pred.logPrPs.fixedpods, aes(x=c, y=y, fill=c))
PrPs.fix <- PrPs.fix +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(ABC-RF)")))+
        ylim(-5.0, -1.0) +
        geom_abline(intercept = c(-4.60206, -1.30103) , slope = 0,color = c("#bdbdbd","black"),  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
PrPs.fix

# ADDITIONAL s=0.01
PrPs.fix_add1 <- ggplot(pred.logPrPs.fixedpods_add[c(1:100, 101:200, 301:400), ], aes(x=c, y=y, fill=c))
PrPs.fix_add1 <- PrPs.fix_add1 +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(P)[R] * italic(P)[B])), " ","(ABC-RF)")))+
        ylim(-5.0, -1.0) +
        geom_abline(intercept = c(-4.60206, -1.30103) , slope = 0,color = c("#bdbdbd","black"),  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
PrPs.fix_add1

## logit Proportion of strongly selected mutations
# ORIGINAL s=0.1
ps.boxplot.values <- which(pred.lgtps.fixedpods$obs_lgt == -12)

ps.fix <- ggplot(pred.lgtps.fixedpods[-ps.boxplot.values,], aes(x=obs_lgt , y=est_lgt, color=c))
ps.fix <- ps.fix + geom_point(size=3, shape=1, alpha=0.65, position = "jitter") + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("neutral","mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        ylab(expression(paste(logit(hat(italic(P))), " ","(ABC-RF)"))) +
        #ylab("") +
        xlim(-12.5, -3) +
        ylim(-12.5, -3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

ps.fix 

# ADDITIONAL s=0.01
pred.lgtps.fixedpods_add1 <- pred.lgtps.fixedpods_add[c(1:200, 301:400), ]
ps.boxplot.values_add1 <- which(pred.lgtps.fixedpods_add1$obs_lgt == -12)

ps.fix_add1 <- ggplot(pred.lgtps.fixedpods_add1[-ps.boxplot.values_add1, ], aes(x=obs_lgt , y=est_lgt, color=c))
ps.fix_add1 <- ps.fix_add1 + geom_point(size=3, shape=1,alpha=0.65,position = "jitter") + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("neutral","mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        ylab(expression(paste(logit(hat(italic(P))), " ","(ABC-RF)"))) +
        #ylab("") +
        xlim(-12.5, -3) +
        ylim(-12.5, -3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

ps.fix_add1 

# BOX-PLOTS NEUTRAL
neutral.ps.plot <- ggplot(pred.lgtps.fixedpods[ps.boxplot.values,], aes(x=c, y=est_lgt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        #ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-12.5, -3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.ps.plot

neutral.ps.add.plot <- ggplot(pred.lgtps.fixedpods_add1[ps.boxplot.values_add1, ], aes(x=c, y=est_lgt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", logit(italic(P))))) +
        #ylab(expression(paste(logit(italic(hat(P))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-12.5, -3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

neutral.ps.add.plot

## log gamma mean
# ORIGINAL s=0.1
gammamean.fix <- ggplot(pred.loggammamean.fixedpods, aes(x=c, y=y, fill=c))
gammamean.fix <- gammamean.fix +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](italic(hat(gamma))), " ","(ABC-RF)")))+
        ylim(-3.0, -0.5) +
        geom_abline(intercept = -1 , slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gammamean.fix

# ADDITIONAL s=0.01
gammamean.fix_add1 <- ggplot(pred.loggammamean.fixedpods_add[c(1:100, 101:200, 301:400), ], aes(x=c, y=y, fill=c))
gammamean.fix_add1 <- gammamean.fix_add1 +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](italic(hat(gamma))), " ","(ABC-RF)")))+
        ylim(-3.0, -0.5) +
        geom_abline(intercept = -2 , slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gammamean.fix_add1

## Log bar gamma
# ORIGINAL s=0.1
selmean.boxplot.values <- which(pred.logpopStrongSelMean.fixedpods$obs_log == -4)

selmean.fix <- ggplot(pred.logpopStrongSelMean.fixedpods[-selmean.boxplot.values, ], aes(x=obs_log , y=est_log, color=c))
selmean.fix <- selmean.fix + geom_point(size=3, shape=1, alpha=0.65, position = "jitter") + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("neutral","mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", italic(bar(s))))) +
        ylab(expression(paste(logit(hat(italic(bar(s)))), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-5, 2) +
        ylim(-5, 2) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

selmean.fix

# ADDITIONAL s=0.01
pred.logpopStrongSelMean.fixedpods_add1 <- pred.logpopStrongSelMean.fixedpods_add[c(1:200, 301:400), ]
selmean.boxplot.values_add1 <- which(pred.logpopStrongSelMean.fixedpods_add1$obs_log == -4)

selmean.fix_add1 <- ggplot(pred.logpopStrongSelMean.fixedpods_add1[-selmean.boxplot.values_add1, ], aes(x=obs_log , y=est_log, color=c))
selmean.fix_add1 <- selmean.fix_add1 + geom_point(size=3, shape=1, alpha=0.65, position = "jitter") + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("neutral","mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", italic(bar(s))))) +
        ylab(expression(paste(logit(hat(italic(bar(s)))), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-5, 2) +
        ylim(-5, 2) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

selmean.fix_add1 

# BOX-PLOTS NEUTRAL
neutral.selmean.plot <- ggplot(pred.logpopStrongSelMean.fixedpods[selmean.boxplot.values, ], aes(x=c, y=est_log)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", bar(italic(s))))) +
        #ylab(expression(paste(log[10](hat(bar(italic(s)))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-5, 2) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.selmean.plot

neutral.selmean.add.plot <- ggplot(pred.logpopStrongSelMean.fixedpods_add1[selmean.boxplot.values_add1, ], aes(x=c, y=est_log)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", bar(italic(s))))) +
        #ylab(expression(paste(log[10](hat(bar(italic(s)))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-5, 2) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.selmean.add.plot

## log10 Theta selected mutations
# ORIGINAL s=0.1
thetaPS.fix <- ggplot(pred.logthetaPS.fixedpods[-c(1:100),], aes(x=obs_logt , y=est_logt, color=c))
thetaPS.fix <- thetaPS.fix + geom_point(size=3, shape=1,alpha=0.65) + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
        ylab(expression(paste(log[10](hat(italic(θ)*italic(P)[S])), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-3, 3) +
        ylim(-3, 3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
thetaPS.fix

# ADDITIONAL s=0.01
thetaPS.fix_add1 <- ggplot(pred.logthetaPS.fixedpods_add[c(101:200, 301:400), ], aes(x=obs_logt , y=est_logt, color=c))
thetaPS.fix_add1 <- thetaPS.fix_add1 + geom_point(size=3, shape=1,alpha=0.65) + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A","B", "C"),
                labels = ""##c("neutral","mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
        ylab(expression(paste(log[10](hat(italic(θ)*italic(P)[S])), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-3, 3) +
        ylim(-3, 3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
thetaPS.fix_add1

# BOX-PLOTS NEUTRAL
neutral.thetaPS.plot <- ggplot(pred.logthetaPS.fixedpods[c(1:100),], aes(x=c, y=est_logt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
        #ylab(expression(paste(log[10](hat(italic(θ)*italic(P)[S])), " ","(ABC-RF)")))+
        ylab("") +
        ylim(-3, 3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.thetaPS.plot

neutral.thetaPS.add.plot <- ggplot(pred.logthetaPS.fixedpods_add[c(1:100), ], aes(x=c, y=est_logt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", log[10](italic(θ) * italic(P)[S])))) +
        #ylab(expression(paste(log[10](hat(italic(θ)*italic(P)[S])), " ","(ABC-RF)")))+
        ylab("") +
        ylim(-3, 3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

neutral.thetaPS.add.plot

## logit genetic load
# ORIGINAL s=0.1
gl.boxplot.values <- which(pred.lgtgl.fixedpods$obs_lgt == -17)

gl.fix <- ggplot(pred.lgtgl.fixedpods[-gl.boxplot.values,], aes(x=obs_lgt , y=est_lgt, color=c))
gl.fix <- gl.fix + 
        geom_point(size=3, shape=1, alpha=0.65) + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A", "B", "C"),
                labels = ""##c("neutral", "mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", logit(italic(L))))) +
        ylab(expression(paste(logit(hat(italic(L))), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-17.5, 3) +
        ylim(-6, 3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gl.fix

# ADDITIONAL s=0.01
pred.lgtgl.fixedpods_add1 <- pred.lgtgl.fixedpods_add[c(1:200, 301:400), ]
gl.boxplot.values_add1 <- which(pred.lgtgl.fixedpods_add1$obs_lgt == -17)

gl.fix_add1 <- ggplot(pred.lgtgl.fixedpods_add1[-gl.boxplot.values_add1, ], aes(x=obs_lgt , y=est_lgt, color=c))
gl.fix_add1 <- gl.fix_add1 + 
        geom_point(size=3, shape=1, alpha=0.65) + 
        scale_colour_manual(values = c(#"black",
                "#2ca25f","#fd8d3c"),
                name   = "",#"Scenario",
                breaks = "",#c("A", "B", "C"),
                labels = ""##c("neutral", "mutation limited", "mutation unlimited"),
        ) + 
        xlab(expression(paste("true", " ", logit(italic(L))))) +
        ylab(expression(paste(logit(hat(italic(L))), " ","(ABC-RF)")))+
        #ylab("") +
        xlim(-17.5, 3) +
        ylim(-6, 3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
gl.fix_add1

# BOX-PLOTS NEUTRAL
neutral.load.plot <- ggplot(pred.lgtgl.fixedpods[gl.boxplot.values,], aes(x=c, y=est_lgt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", logit(italic(L))))) +
        #ylab(expression(paste(logit(hat(italic(L))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-6, 3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
neutral.load.plot

neutral.load.add.plot <- ggplot(pred.lgtgl.fixedpods_add1[gl.boxplot.values_add1, ], aes(x=c, y=est_lgt)) + 
        geom_boxplot() + 
        #xlab("") +
        #ylab("") +
        xlab(expression(paste("true", " ", logit(italic(L))))) +
        #ylab(expression(paste(logit(hat(italic(L))), " ","(ABC-RF)"))) +
        ylab("") +
        ylim(-6, 3) +
        #geom_abline(intercept = 0, slope = 0, color = "darkgrey",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "right", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))

neutral.load.add.plot

## Combined plots
# ORIGINAL s=0.1
box_plots <- cowplot::plot_grid(neutral.ps.plot, neutral.selmean.plot,
                                neutral.thetaPS.plot, neutral.load.plot,
                                nrow = 1, ncol = 4,
                                labels = NULL)

cowplot::plot_grid(PrPs.fix, ps.fix,
                   gammamean.fix, selmean.fix,
                   thetaPS.fix, gl.fix,
                   box_plots,
                   ncol = 2,
                   nrow = 4,
                   labels = NULL)

# ADDITIONAL s=0.01
box_plots_add1 <- cowplot::plot_grid(neutral.ps.add.plot, neutral.selmean.add.plot,
                                     neutral.thetaPS.add.plot, neutral.load.add.plot,
                                     nrow = 1, ncol = 4,
                                     labels = NULL)

cowplot::plot_grid(PrPs.fix_add1, ps.fix_add1,
                   gammamean.fix_add1, selmean.fix_add1,
                   thetaPS.fix_add1, gl.fix_add1,
                   box_plots_add1,
                   ncol = 2,
                   nrow = 4,
                   labels = NULL)

### FIXED PODS - DEMOGRAPHY

## log mu
# ORIGINAL s=0.1
mu.fix <- ggplot(pred.logmu.fixedpods, aes(x=c, y=y, fill=c))
mu.fix <- mu.fix +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(mu))), " ","(ABC-RF)")))+
        ylim(-8.2,-7.70) +
        geom_abline(intercept = -8, slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
mu.fix

# ADDITIONAL s=0.01
mu.fix_add1 <- ggplot(pred.logmu.fixedpods_add[c(1:200, 301:400), ], aes(x=c, y=y, fill=c))
mu.fix_add1 <- mu.fix_add1 +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(mu))), " ","(ABC-RF)")))+
        ylim(-8.2,-7.80) +
        geom_abline(intercept = -8, slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
mu.fix_add1

## log rr
# ORIGINAL s=0.1
rr.fix <- ggplot(pred.logrr.fixedpods_add, aes(x=c, y=y, fill=c))
rr.fix <- rr.fix +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(c))[0]), " ","(ABC-RF)")))+
        ylim(-8.8,-6.40) +
        geom_abline(intercept = -7.30103 , slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
rr.fix

# ADDITIONAL s=0.01
rr.fix_add1 <- ggplot(pred.logrr.fixedpods_add[c(1:200, 301:400), ], aes(x=c, y=y, fill=c))
rr.fix_add1 <- rr.fix_add1 +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(c))[0]), " ","(ABC-RF)")))+
        ylim(-8.8,-6.40) +
        geom_abline(intercept = -7.30103 , slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
rr.fix_add1

## log NCS2
# ORIGINAL s=0.1
ncs.fix <- ggplot(pred.logncs.fixedpods, aes(x=c, y=y, fill=c))
ncs.fix <- ncs.fix +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "", #"Scenario",
                          breaks = "", #c("A","B","C"),
                          labels = "", #c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(N))), " ","(ABC-RF)")))+
        ylim(2.16,3) +
        geom_abline(intercept = 2.69897, slope = 0,color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
ncs.fix

# ADDITIONAL s=0.01
ncs.fix_add1 <- ggplot(pred.logncs.fixedpods_add[c(1:200, 301:400), ], aes(x=c , y=y, fill=c))
ncs.fix_add1 <- ncs.fix_add1 +  geom_boxplot(alpha=1) + 
        scale_fill_manual(values = c("black","#2ca25f","#fd8d3c"),
                          name   = "",#"Scenario",
                          breaks = "",#c("A","B","C"),
                          labels = ""#c("neutral","mutation limited", "mutation unlimited")
        ) + 
        scale_x_discrete(labels = c("neutral","limited", "unlimited")) +
        xlab("") +
        ylab(expression(paste(log[10](hat(italic(N))), " ","(ABC-RF)"))) +
        ylim(2.45,3) +
        geom_abline(intercept = 2.69897, slope = 0, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "top", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
ncs.fix_add1

## log mean Ne2
# ORIGINAL s=0.1
ne.fix <- ggplot(pred.logmeanNe2.fixedpods, aes(x=x , y=y))
ne.fix <- ne.fix +  geom_point(aes(color=c, shape=c, size=c)) + 
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](hat(italic(N))[e]), " ","(ABC-RF)")))+
        xlim(2,2.8) +
        ylim(1.7,3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
ne.fix

# ADDITIONAL s=0.01
ne.fix_add1 <- ggplot(pred.logmeanNe2.fixedpods_add[c(1:100,101:200, 301:400), ], aes(x=x , y=y))
ne.fix_add1 <- ne.fix_add1 +  geom_point(aes(color=c, shape=c, size=c)) + 
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("B","C"),
                            #labels = c("mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](hat(italic(N))[e]), " ","(ABC-RF)")))+
        xlim(2.5,2.8) +
        ylim(2.25,3) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
ne.fix_add1

## log Ne/NCS
# ORIGINAL s=0.1
nencs.fix <- ggplot(pred.logmeanNe2ncs.fixedpods, aes(x=x , y=y)) 
nencs.fix <- nencs.fix + geom_point(aes(color=c, shape=c, size=c)) +
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
        ylab(expression(paste(log[10](hat(italic(N)[e]/italic(N))), " ","(ABC-RF)")))+
        xlim(-0.6,0.1) +
        ylim(-0.8,0.1) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
nencs.fix

# ADDITIONAL s=0.01
nencs.fix_add1 <- ggplot(pred.logmeanNe2ncs.fixedpods_add[c(1:200, 301:400), ], aes(x=x , y=y)) 
nencs.fix_add1 <- nencs.fix_add1 +  geom_point(aes(color=c, shape=c, size=c)) + 
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e]/italic(N))))) +
        ylab(expression(paste(log[10](hat(italic(N)[e]/italic(N))), " ","(ABC-RF)")))+
        xlim(-0.2,0.1) +
        ylim(-0.2,0.1) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
nencs.fix_add1

## COMBINED Plots
# ORIGINAL s=0.1 - main text
cowplot::plot_grid(rr.fix, mu.fix,
                   ncs.fix,ne.fix,
                   nencs.fix,
                   nrow = 3,
                   ncol = 2,
                   labels = NULL)

# ADDITIONAL s=0.01 -supplementary material
cowplot::plot_grid(rr.fix_add1,mu.fix_add1,
                   ncs.fix_add1, ne.fix_add1,
                   nencs.fix_add1,
                   nrow = 3,
                   ncol = 2,
                   labels = NULL)


### Prediction Demography and selection (joint inference)
###-----------------------------------------

## log10 Ne2 - with data used to train the RF
#my_breaks_ne2rf <- c(2980.957987,403.428793, 54.598150, 7.389056, 1)
#ne   <- ggplot(oob.logmeanNe2, aes(x,y))
#ne2d <- ne + 
#  stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,breaks = round(my_breaks_ne2rf),labels = round(my_breaks_ne2rf)
#                        ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(OOB prediction)")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ne2d

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
#my_breaks_wfabc <- c(7.389056, 2.718282,1)
#wfne.random   <- ggplot(estimated.wfabc.randompods, aes(x,y))
#wfne2d.random <- wfne.random + 
#        stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#        scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                             ,
#                             breaks = round(my_breaks_wfabc),labels = round(my_breaks_wfabc)
#        ) +
#        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
#        xlim(-0.5,3.5) +
#        ylim(-0.5,3.5) +
#        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                           panel.grid.major = element_blank(),
#                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                           legend.position = "right", axis.text = element_text(size = 12),
#                           axis.title=element_text(size=12))
#wfne2d.random

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
#my_breaks_10 <- c(20.085537,7.389056, 2.718282,1)
#fstne.random   <- ggplot(estimated.fstNe2.randompods, aes(x,y))
#fstne2d.random <- fstne.random + 
#  stat_bin2d(bins=50, na.rm = T, show.legend=T) + 
#  scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
#                       ,
#                       breaks = round(my_breaks_10),labels = round(my_breaks_10)
#                       ) +
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "right", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#fstne2d.random

## log Ne2 FST - with data used to train the RF
my_breaks_fstrf <- c(403.428793, 54.598150, 7.389056, 1)
fstne   <- ggplot(estimated.fstNe2, aes(x,y))
fstne2d <- fstne + 
        stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
        scale_fill_gradientn(colours=viridis::viridis(32), trans="log",name="counts"
                             ,breaks =round(my_breaks_fstrf),labels = round(my_breaks_fstrf)
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
fstne2d

## log mean Ne2
## ORIGINAL s=0.1
#pred.logmeanNe2.fixedpods$x[1:100] <- pred.logmeanNe2.fixedpods$x[1:100] + 0.068
#ne.fix_2 <- ggplot(pred.logmeanNe2.fixedpods, aes(x=x , y=y, color=c))
#ne.fix_2 <- ne.fix_2 +  geom_point(size=3, shape=1, alpha=0.65, position = "jitter") + 
#  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
#                      name   = "",#"Scenario",
#                      breaks = "",#c("A","B","C"),
#                      labels = "",#c("neutral","mutation limited", "mutation unlimited")
#                      ) + 
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "top", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ne.fix_2
#
## ADDITIONAL s=0.01
#pred.logmeanNe2.fixedpods_add$x[1:100] <- pred.logmeanNe2.fixedpods_add$x[1:100] + 0.072
#pred.logmeanNe2.fixedpods_add$x[101:200] <- pred.logmeanNe2.fixedpods_add$x[101:200] + 0.058
#
#ne.fix_add1_2 <- ggplot(pred.logmeanNe2.fixedpods_add[c(1:100,101:200, 301:400), ], aes(x=x , y=y, color=c))
#ne.fix_add1_2 <- ne.fix_add1_2 +  geom_point(size=3, shape=1, alpha=0.65, position = "jitter") + 
#  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
#                      name   = "",#"Scenario",
#                      breaks = "",#c("B","C"),
#                      labels = "",#c("mutation limited", "mutation unlimited")
#                      ) + 
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(ABC-RF)")))+
#  xlim(-0.5,3.5) +
#  ylim(-0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "top", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#ne.fix_add1_2

## log Ne2 WF-ABC
## ORIGINAL s = 0.1
#estimated.wfabc.fixedpods$x[1:100] <- estimated.wfabc.fixedpods$x[1:100] + 0.068
wfne.fix <- ggplot(estimated.wfabc.fixedpods, aes(x=x , y=y))
wfne.fix <- wfne.fix +  geom_point(aes(color=c, shape=c, size=c)) +  
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
        xlim(-0.5,3.5) +
        ylim(-0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
wfne.fix

## ADDITIONAL s = 0.01
#estimated.wfabc.fixedpods_add$x[1:100] <- estimated.wfabc.fixedpods_add$x[1:100] + 0.068
wfne.fix_add1 <- ggplot(estimated.wfabc.fixedpods_add[c(1:100, 101:200, 301:400),], aes(x=x , y=y))
wfne.fix_add1 <- wfne.fix_add1 +  geom_point(aes(color=c, shape=c, size=c)) + 
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
        xlim(-0.5,3.5) +
        ylim(-0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
wfne.fix_add1

## ADDITIONAL s = 0.25
#wfne.fix_add2 <- ggplot(estimated.wfabc.fixedpods_add[c(1:100, 201:300, 401:500),], aes(x=x , y=y, color=c))
#wfne.fix_add2 <- wfne.fix_add2 +  geom_point(size=3, shape=1) + 
#  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
#                      name   = "",#"Scenario",
#                      breaks = "",#c("A","B","C"),
#                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
#                      ) + 
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(WF-ABC)")))+
#  xlim(0.5,3.5) +
#  ylim(0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "top", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#wfne.fix_add2

## log Ne2 FST
## ORIGINAL s = 0.1
#estimated.fstNe2.fixedpods$x[1:100] <- estimated.fstNe2.fixedpods$x[1:100] + 0.068
fstne.fix <- ggplot(estimated.fstNe2.fixedpods, aes(x=x , y=y))
fstne.fix <- fstne.fix +  geom_point(aes(color=c, shape=c, size=c)) +  
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) +  
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
        xlim(0.5,3.5) +
        ylim(0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
fstne.fix

## ADDITIONAL s = 0.01
#estimated.fstNe2.fixedpods_add$x[1:100] <- estimated.fstNe2.fixedpods_add$x[1:100] + 0.068
fstne.fix_add1 <- ggplot(estimated.fstNe2.fixedpods_add[c(1:100, 101:200, 301:400),], aes(x=x , y=y))
fstne.fix_add1 <- fstne.fix_add1 +  geom_point(aes(color=c, shape=c, size=c)) +
        scale_shape_manual(values=c(1, 4, 2)) +
        scale_size_manual(values=c(6, 3, 3)) +
        scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
                            #name   = "Scenario",
                            #breaks = c("A","B","C"),
                            #labels = c("neutral","mutation limited", "mutation unlimited")
        ) + 
        xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
        ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
        xlim(0.5,3.5) +
        ylim(0.5,3.5) +
        geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
        theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "", axis.text = element_text(size = 12),
                           axis.title=element_text(size=12))
fstne.fix_add1

## ADDITIONAL s = 0.25
#fstne.fix_add2 <- ggplot(estimated.fstNe2.fixedpods_add[c(1:100, 201:300, 401:500),], aes(x=x , y=y, color=c))
#fstne.fix_add2 <- fstne.fix_add2 +  geom_point(size=3, shape=1) + 
#  scale_colour_manual(values = c("black","#2ca25f","#fd8d3c"),
#                      name   = "",#"Scenario",
#                      breaks = "",#c("A","B","C"),
#                      labels = ""#c("neutral","mutation limited", "mutation unlimited")
#                      ) + 
#  xlab(expression(paste("true", " ", log[10](italic(N)[e])))) +
#  ylab(expression(paste(log[10](italic(hat(N))[e]), " ","(from", " ", italic(hat(F))[ST], ")")))+
#  xlim(0.5,3.5) +
#  ylim(0.5,3.5) +
#  geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
#  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                     panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                     legend.position = "top", axis.text = element_text(size = 12),
#                     axis.title=element_text(size=12))
#fstne.fix_add2

## COMBINED Plots
cowplot::plot_grid(wfne2d.random,fstne2d,
                   wfne.fix_add1, fstne.fix_add1,
                   wfne.fix, fstne.fix,
                   nrow = 3,
                   labels = NULL)
