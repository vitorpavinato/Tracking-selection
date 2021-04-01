########################################
##          Manuscript plots          ##
##             Application            ##
########################################

library(gtools)

## DENSITY PLOTS
##---------------------------------------

## JOINT INFERENCE OF DEMOGRAPY AND SELECTION
##--------------------------------------------

# LOAD LIST OF VECTOR OF PRIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logthetaPS",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmeanNe2",".RData"))

# LOAD LIST OF VECTORS OF POSTERIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logthetaPS",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logmeanNe2",".RData"))

## SELECTION ONLY
##-------------------

# LOAD LIST OF VECTOR OF PRIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logPrPs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logitpopstrongmsel",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_loggammamean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logpopstrongselmean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_averageGenLoad",".RData"))

# LOAD LIST OF VECTORS OF POSTERIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logPrPs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_popstrongmsel",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_loggammamean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_popstrongselmean",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_averageGenLoad",".RData"))

population.vector <- c("Avalon","Humboldt","Davis","Stanislaus","Stebbins","Riverside","Platerita")
color.vector <- c("black", "red", "green", "blue", "cyan", "#FFDB58", "orange")

# FUNCTION TO PRODUCE THE PLOTS - LOG and LOGIT SCALE
plot.mult.densities <- function(list.prior=NULL, list.posterior=NULL, par.name=NULL, 
                                col.vect=NULL, pop.vect=NULL, y_lim=c(0,0.5), cex_axis = 1.2, cex_lab  = 1.2,
                                plot.legend=TRUE, legend.side="topright")
{
  xlim.min = min(sapply(list.prior, min))
  xlim.max = min(sapply(list.prior, max))
  
  plot(density(list.prior[[1]]), lty=3, col=col.vect[1], xlim=c(xlim.min, xlim.max), ylim=y_lim, 
       main = "", ylab = "Density", xlab = par.name, cex.axis = cex_axis, cex.lab  = cex_lab)
  lines(density(list.prior[[1]], weights = list.posterior[[1]]$weights), col=col.vect[1])
  
  for (p in 2:length(list.prior))
  {
    lines(density(list.prior[[p]]),lty=3, col=col.vect[p])
    lines(density(list.prior[[p]], weights = list.posterior[[p]]$weights), col=col.vect[p])
  }
  if (plot.legend) {
    legend(legend.side, col = c(rep(col.vect[1],2),col.vect[1:7]), lty = c(3,1,rep(1,7)), cex = 1.4,
           legend = c("prior", "posterior", pop.vect), box.lwd = 0, box.col = NULL, bg = NULL)
  }
}

plot.mult.densities.orig <- function(list.prior=NULL, list.posterior=NULL, par.name=NULL, 
                                     col.vect=NULL, pop.vect=NULL, y_lim=c(0,0.5), cex_axis = 1.2, cex_lab  = 1.2,
                                     plot.legend=TRUE, legend.side="topright", fromScale="logit")
{
  if (fromScale=="logit")
  {
    xlim.min = min(inv.logit(sapply(list.prior, min)))
    xlim.max = min(inv.logit(sapply(list.prior, max)))
    
    plot(density(inv.logit(list.prior[[1]])), lty=3, col=col.vect[1], xlim=c(xlim.min, xlim.max), ylim=y_lim, 
         main = "", ylab = "Density", xlab = par.name, cex.axis = cex_axis, cex.lab  = cex_lab)
    lines(density(inv.logit(list.prior[[1]]), weights = list.posterior[[1]]$weights), col=col.vect[1])
    
    for (p in 2:length(list.prior))
    {
      lines(density(inv.logit(list.prior[[p]])),lty=3, col=col.vect[p])
      lines(density(inv.logit(list.prior[[p]]), weights = list.posterior[[p]]$weights), col=col.vect[p])
    }
  } else {
    xlim.min = min(10^(sapply(list.prior, min)))
    xlim.max = min(10^(sapply(list.prior, max)))
    
    plot(density(10^(list.prior[[1]])), lty=3, col=col.vect[1], xlim=c(xlim.min, xlim.max), ylim=y_lim, 
         main = "", ylab = "Density", xlab = par.name, cex.axis = cex_axis, cex.lab  = cex_lab)
    lines(density(10^(list.prior[[1]]), weights = list.posterior[[1]]$weights), col=col.vect[1])
    
    for (p in 2:length(list.prior))
    {
      lines(density(10^(list.prior[[p]])),lty=3, col=col.vect[p])
      lines(density(10^(list.prior[[p]]), weights = list.posterior[[p]]$weights), col=col.vect[p])
    }
  }
  
  if (plot.legend) {
    legend(legend.side, col = c(rep(col.vect[1],2),col.vect[1:7]), lty = c(3,1,rep(1,7)), cex = 1.4,
           legend = c("prior", "posterior", pop.vect), box.lwd = 0,box.col = NULL, bg = NULL)
  }
}

## MANUSCRIPT PLOT - BEES JOINT INFERENCE
pdf(file = "joint.pdf", height = 5.5, width = 15.85)
par(mar=c(5,5,4,1)+.1, mfrow=c(1,3))

# THETA PS
plot.mult.densities(list.prior = list.logthetaPS,
                    list.posterior = list.posterior.logthetaPS,
                    par.name = expression(log[10](italic(theta)[b])),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,0.4),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = TRUE, legend.side = "topleft")

# N
plot.mult.densities(list.prior = list.logncs,
                    list.posterior = list.posterior.logncs,
                    par.name = expression(log[10](italic(N))),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,1.5),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = FALSE)

# NE
plot.mult.densities(list.prior = list.logmeanNe2,
                    list.posterior = list.posterior.logmeanNe2,
                    par.name = expression(log[10](italic(N)[e])),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,1.5),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = FALSE)
dev.off()

## SELECTION ONLY
##-----------------
pdf(file = "selection.pdf", height = 11.00, width = 8.5)
par(mar=c(5,5,4,1)+.1, mfrow=c(3,2))
# PrPs - parameter
plot.mult.densities(list.prior = list.logPrPs,
                    list.posterior = list.posterior.logPrPs,
                    par.name = expression(log[10](italic(P)[R] * italic(P)[S])),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,0.5),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = TRUE,
                    legend.side = "topleft")

# NUMBER OF STRONGLY SELECTED MUTATIONS
plot.mult.densities(list.prior = list.logitpopstrongmsel,
                    list.posterior = list.posterior.popstrongmsel,
                    par.name = expression(logit(italic(P))),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,0.4),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = FALSE)

# GAMMA MEAN - Parameter
plot.mult.densities(list.prior = list.loggammamean,
                    list.posterior = list.posterior.loggammamean,
                    par.name = expression(log[10](italic(gamma))),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,0.6),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = FALSE)

# GAMMA STRONG SELECTION
#plot.mult.densities(list.prior = list.logpopstrongselmean,
#                    list.posterior = list.posterior.popstrongselmean,
#                    par.name = expression(log[10](italic(bar(s)))),
#                    col.vect = color.vector,
#                    pop.vect = population.vector,
#                    y_lim = c(0,1.2),
#                    cex_axis = 1.6, cex_lab = 1.8,
#                    plot.legend = FALSE)

plot.mult.densities.orig(list.prior = list.logpopstrongselmean,
                         list.posterior = list.posterior.popstrongselmean,
                         par.name = expression(italic(bar(s))),
                         col.vect = color.vector,
                         pop.vect = population.vector,
                         y_lim = c(0,4),
                         cex_axis = 1.6, cex_lab = 1.8,
                         plot.legend = FALSE,
                         fromScale = "log")


# AVERAGE GENETIC LOAD
#plot.mult.densities(list.prior = list.averageGenLoad,
#                    list.posterior = list.posterior.averageGenLoad,
#                    par.name = expression(log[10](italic(L))),
#                    col.vect = color.vector,
#                    pop.vect = population.vector,
#                    y_lim = c(0,0.5),
#                    cex_axis = 1.6, cex_lab = 1.8,
#                    plot.legend = FALSE)

plot.mult.densities.orig(list.prior = list.averageGenLoad,
                         list.posterior = list.posterior.averageGenLoad,
                         par.name = expression(italic(L)),
                         col.vect = color.vector,
                         pop.vect = population.vector,
                         y_lim = c(0,7),
                         cex_axis = 1.6, cex_lab = 1.8,
                         plot.legend = FALSE,
                         fromScale = "logit")
dev.off()

## DEMOGRAPHY ONLY
##--------------------

# LOAD LIST OF VECTOR OF PRIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logmu",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/random_forests/list_vector_logrr",".RData"))

## LOAD LIST OF VECTORS OF POSTERIORS
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logmeanNe2ncs",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logmu",".RData"))
#load(file=paste0("~/My_repositories/Tracking-selection/results/pipeline_v6_bees/posterior_obs/list_posterior_logrr",".RData"))

pdf(file = "demography.pdf", height = 11.00, width = 8.5)
par(mar=c(5,5,4,1)+.1, mfrow=c(2,2))
# mu
plot.mult.densities(list.prior = list.logmu,
                    list.posterior = list.posterior.logmu,
                    par.name = expression(log[10](italic(mu))),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,2.5),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = TRUE,
                    legend.side = "topleft")

# c0
plot.mult.densities(list.prior = list.logrr,
                    list.posterior = list.posterior.logrr,
                    par.name = expression(log[10](italic(c)[0])),
                    col.vect = color.vector,
                    pop.vect = population.vector,
                    y_lim = c(0,1.5),
                    cex_axis = 1.6, cex_lab = 1.8,
                    plot.legend = FALSE)

# Ne/N
#plot.mult.densities(list.prior = list.logmeanNe2ncs,
#                    list.posterior = list.posterior.logmeanNe2ncs,
#                    par.name = expression(log[10](italic(N)[e]/italic(N))),
#                    col.vect = color.vector,
#                    pop.vect = population.vector,
#                    y_lim = c(0,12.0),
#                    cex_axis = 1.6, cex_lab = 1.8,
#                    plot.legend = FALSE)

plot.mult.densities.orig(list.prior = list.logmeanNe2ncs,
                         list.posterior = list.posterior.logmeanNe2ncs,
                         par.name = expression(italic(N)[e]/italic(N)),
                         col.vect = color.vector,
                         pop.vect = population.vector,
                         y_lim = c(0,12),
                         cex_axis = 1.6, cex_lab = 1.8,
                         plot.legend = FALSE,
                         fromScale = "log")
dev.off()
