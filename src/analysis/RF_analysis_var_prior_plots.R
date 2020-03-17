# manuscript plots

# 1 - var plots
# 2 - prior histograms

### VAR PLOTS
###--------------------

par(mfrow=c(2, 4))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
plot(x = reg_averageGenLoad_2, n.var = 20, main=expression(logit(italic(L))), cex.axis = 0.5)

# logit Proportion of strongly selected mutations
plot(x = reg_logitpopstrongmsel, n.var = 20, main=expression(logit(italic(P))), cex.axis = 0.5)

# log10 Theta selected mutations
plot(x = reg_logthetaPS, n.var = 20, main=expression(log[10](italic(θ) * italic(P)[S])), cex.axis = 0.5)

# log10 Ns
plot(x = reg_logNs, n.var = 20, main=expression(log[10](italic(N)[e] * italic(s))), cex.axis = 0.5)

# log mean Ne2
plot(x = reg_logmeanNe2, n.var = 20, main=expression(log[10](italic(N)[e])), cex.axis = 0.5)

# log NCS
plot(x = reg_logncs, n.var = 20, main=expression(log[10](italic(N)[cs])), cex.axis = 0.5)

# log Ne2/NCS
plot(x = reg_logmeanNe2ncs, n.var = 20, main=expression(log[10](italic(N)[e]/italic(N)[cs])), cex.axis = 0.5)

# log theta2
plot(x = reg_logtheta2, n.var = 20, main=expression(log[10](italic(θ))), cex.axis = 0.5)

### PRIOR HISTOGRAMS
###------------------------

par(mfrow=c(4, 2))
par(mar=c(5,5,4,1)+.1)

# logit average genetic load
hist(logitaverageGenload, freq = TRUE, xlab = expression(logit(italic(L))), main = "", col = "#bdbdbd")

# logit Proportion of strongly selected mutations
hist(logitpopstrongmsel, freq = TRUE, xlab = expression(logit(italic(P))), main = "", col = "#bdbdbd")

# log10 Theta selected mutations
hist(logthetaPS, freq = TRUE, xlab = expression(log[10](italic(θ) * italic(P)[S])), main = "", col = "#bdbdbd")

# log10 Ns
hist(logNs, freq = TRUE, xlab = expression(log[10](italic(N)[e] * italic(s))), main = "", col = "#bdbdbd")

# log mean Ne2
hist(logmeanNe2, freq = TRUE, xlab = expression(log[10](italic(N)[e])) , main = "", col = "#bdbdbd")

# log NCS
hist(logncs, freq = TRUE, xlab = expression(log[10](italic(N)[cs])), main = "", col = "#bdbdbd")

# log Ne2/NCS
hist(logmeanNe2ncs, freq = TRUE, xlab = expression(log[10](italic(N)[e]/italic(N)[cs])), main = "", col = "#bdbdbd")

# log theta2
hist(logtheta2, freq = TRUE, xlab = expression(log[10](italic(θ))), main = "", col = "#bdbdbd")
