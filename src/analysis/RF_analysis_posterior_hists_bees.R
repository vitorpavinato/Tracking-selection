##### manuscript plots

### Histograms with Prior and Posterior (demography and selection)

## PLACERITA
##-------------------------------
par(mfrow=c(4, 2))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
hist(x      = logitaverageGenload_placerita,
     breaks = seq(-37,7,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(logit(italic(L))),
     ylab   = "probability density",
     ylim   = c(0,0.25))
wtd.hist(x      = logitaverageGenload_placerita,
         breaks = seq(-37,7,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_averageGenLoad_placerita_2$weights)
abline(v=c(posterior_averageGenLoad_placerita_2$med,
           posterior_averageGenLoad_placerita_2$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# logit Proportion of strongly selected mutations
hist(x      = logitpopstrongmsel_placerita,
     breaks = seq(-14,2,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(logit(italic(P))),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logitpopstrongmsel_placerita,
         breaks = seq(-14,2,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logitpopstrongmsel_placerita$weights)
abline(v=c(posterior_logitpopstrongmsel_placerita$med,
           posterior_logitpopstrongmsel_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log10 Theta selected mutations
hist(x      = logthetaPS_pla,
     breaks = seq(-7,6,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(θ) * italic(P)[S])),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logthetaPS_pla,
         breaks = seq(-7,6,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logthetaPS_placerita$weights)
abline(v=c(posterior_logthetaPS_placerita$med,
           posterior_logthetaPS_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log10 Ns
hist(x      = logNs_pla,
     breaks = seq(-3,4,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[e] * italic(s))),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logNs_pla,
         breaks = seq(-3,4,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logNs_placerita$weights)
abline(v=c(posterior_logNs_placerita$med,
           posterior_logNs_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log mean Ne2
hist(x      = logmeanNe2_pla,
     breaks = seq(-1,5,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[e])),
     ylab   = "probability density",
     ylim   = c(0,2.0))
wtd.hist(x      = logmeanNe2_pla,
         breaks = seq(-1,5,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logmeanNe2_placerita$weights)
abline(v=c(posterior_logmeanNe2_placerita$med,
           posterior_logmeanNe2_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log NCS
hist(x      = logncs_pla,
     breaks = seq(0,4.0,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[cs])),
     ylab   = "probability density",
     ylim   = c(0,2.0))
wtd.hist(x      = logncs_pla,
         breaks = seq(0,4,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logncs_placerita$weights)
abline(v=c(posterior_logncs_placerita$med,
           posterior_logncs_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log Ne2/NCS
hist(x      = logmeanNe2ncs_pla,
     breaks = seq(-4,1,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(italic(N)[e]/italic(N)[cs]),
     ylab   = "probability density",
     ylim   = c(0,4.0))
wtd.hist(x      = logmeanNe2ncs_pla,
         breaks = seq(-4,1,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logmeanNe2ncs_placerita$weights)
abline(v=c(posterior_logmeanNe2ncs_placerita$med,
           posterior_logmeanNe2ncs_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log theta2
hist(x      = logtheta2_pla,
     breaks = seq(-1,6,0.2),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(θ))),
     ylab   = "probability density",
     ylim   = c(0,2))
wtd.hist(x      = logtheta2_pla,
         breaks = seq(-1,6,0.2),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logtheta2_placerita$weights)
abline(v=c(posterior_logtheta2_placerita$med,
           posterior_logtheta2_placerita$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

## AVALON
##-------------------------------
par(mfrow=c(4, 2))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
hist(x      = logitaverageGenload_avalon,
     breaks = seq(-22,6,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(logit(italic(L))),
     ylab   = "probability density",
     ylim   = c(0,0.2))
wtd.hist(x      = logitaverageGenload_avalon,
         breaks = seq(-22,6,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_averageGenLoad_avalon_2$weights)
abline(v=c(posterior_averageGenLoad_avalon_2$med,
           posterior_averageGenLoad_avalon_2$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# logit Proportion of strongly selected mutations
hist(x      = logitpopstrongmsel_avalon,
     breaks = seq(-17,1,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(logit(italic(P))),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logitpopstrongmsel_avalon,
         breaks = seq(-16,1,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logitpopstrongmsel_avalon$weights)
abline(v=c(posterior_logitpopstrongmsel_avalon$med,
           posterior_logitpopstrongmsel_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log10 Theta selected mutations
hist(x      = logthetaPS_ava,
     breaks = seq(-7,6,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(θ) * italic(P)[S])),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logthetaPS_ava,
         breaks = seq(-7,6,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logthetaPS_avalon$weights)
abline(v=c(posterior_logthetaPS_avalon$med,
           posterior_logthetaPS_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log10 Ns
hist(x      = logNs_ava,
     breaks = seq(-3,4,0.5),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[e] * italic(s))),
     ylab   = "probability density",
     ylim   = c(0,0.5))
wtd.hist(x      = logNs_ava,
         breaks = seq(-3,4,0.5),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logNs_avalon$weights)
abline(v=c(posterior_logNs_avalon$med,
           posterior_logNs_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log mean Ne2
hist(x      = logmeanNe2_ava,
     breaks = seq(-2,5,0.1),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[e])),
     ylab   = "probability density",
     ylim   = c(0,2.0))
wtd.hist(x      = logmeanNe2_ava,
         breaks = seq(-2,5,0.1),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logmeanNe2_avalon$weights)
abline(v=c(posterior_logmeanNe2_avalon$med,
           posterior_logmeanNe2_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log NCS
hist(x      = logncs_ava,
     breaks = seq(0,4.0,0.1),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(N)[cs])),
     ylab   = "probability density",
     ylim   = c(0,2.0))
wtd.hist(x      = logncs_ava,
         breaks = seq(0,4,0.1),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logncs_avalon$weights)
abline(v=c(posterior_logncs_avalon$med,
           posterior_logncs_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log Ne2/NCS
hist(x      = logmeanNe2ncs_ava,
     breaks = seq(-6,1,0.1),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(italic(N)[e]/italic(N)[cs]),
     ylab   = "probability density",
     ylim   = c(0,4.0))
wtd.hist(x      = logmeanNe2ncs_ava,
         breaks = seq(-6,1,0.1),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logmeanNe2ncs_avalon$weights)
abline(v=c(posterior_logmeanNe2ncs_avalon$med,
           posterior_logmeanNe2ncs_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))

# log theta2
hist(x      = logtheta2_ava,
     breaks = seq(-1,6,0.1),
     col    = "#969696",
     freq   = FALSE,
     main   = "",
     xlab   = expression(log[10](italic(θ))),
     ylab   = "probability density",
     ylim   = c(0,1.5))
wtd.hist(x      = logtheta2_ava,
         breaks = seq(-1,6,0.1),
         col    = adjustcolor( "#41b6c4", alpha.f = 0.5), freq=FALSE, add=TRUE,
         weight = posterior_logtheta2_avalon$weights)
abline(v=c(posterior_logtheta2_avalon$med,
           posterior_logtheta2_avalon$quantiles[c(1,3)]),
       col="red",
       lty=c(1,3,3))