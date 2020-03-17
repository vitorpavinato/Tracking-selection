# manuscript plots

### PLACERITA
###--------------------

par(mfrow=c(2, 4))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
plot(x = reg_averageGenLoad_placerita_2, n.var = 20, main=expression(logit(italic(L))), cex.axis = 0.5)

# logit Proportion of strongly selected mutations
plot(x = reg_logitpopstrongmsel_placerita, n.var = 20, main=expression(logit(italic(P))), cex.axis = 0.5)

# log10 Theta selected mutations
plot(x = reg_logthetaPS_placerita, n.var = 20, main=expression(log[10](italic(θ) * italic(P)[S])), cex.axis = 0.5)

# log10 Ns
plot(x = reg_logNs_placerita, n.var = 20, main=expression(log[10](italic(N)[e] * italic(s))), cex.axis = 0.5)

# log mean Ne2
plot(x = reg_logmeanNe2_placerita, n.var = 20, main=expression(log[10](italic(N)[e])), cex.axis = 0.5)

# log NCS
plot(x = reg_logncs_placerita, n.var = 20, main=expression(log[10](italic(N)[cs])), cex.axis = 0.5)

# log Ne2/NCS
plot(x = reg_logmeanNe2ncs_placerita, n.var = 20, main=expression(log[10](italic(N)[e]/italic(N)[cs])), cex.axis = 0.5)

# log theta2
plot(x = reg_logtheta2_placerita, n.var = 20, main=expression(log[10](italic(θ))), cex.axis = 0.5)

### AVALON
###--------------------

par(mfrow=c(2, 4))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
plot(x = reg_averageGenLoad_avalon_2, n.var = 20, main=expression(logit(italic(L))), cex.axis = 0.5)

# logit Proportion of strongly selected mutations
plot(x = reg_logitpopstrongmsel_avalon, n.var = 20, main=expression(logit(italic(P))), cex.axis = 0.5)

# log10 Theta selected mutations
plot(x = reg_logthetaPS_avalon, n.var = 20, main=expression(log[10](italic(θ) * italic(P)[S])), cex.axis = 0.5)

# log10 Ns
plot(x = reg_logNs_avalon, n.var = 20, main=expression(log[10](italic(N)[e] * italic(s))), cex.axis = 0.5)

# log mean Ne2
plot(x = reg_logmeanNe2_avalon, n.var = 20, main=expression(log[10](italic(N)[e])), cex.axis = 0.5)

# log NCS
plot(x = reg_logncs_avalon, n.var = 20, main=expression(log[10](italic(N)[cs])), cex.axis = 0.5)

# log Ne2/NCS
plot(x = reg_logmeanNe2ncs_avalon, n.var = 20, main=expression(log[10](italic(N)[e]/italic(N)[cs])), cex.axis = 0.5)

# log theta2
plot(x = reg_logtheta2_avalon, n.var = 20, main=expression(log[10](italic(θ))), cex.axis = 0.5)