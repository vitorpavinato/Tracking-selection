setwd("/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_abc")

install.packages("EasyABC")
library(EasyABC)

runSLiM <- function(x){
  seed <- x[1]
  mu <- x[2]
  #cat('Running SLiM with seed ', seed, ', mu = ', mu, '\n')
  output <- system2("/usr/local/bin/slim", c("-d", paste0("mu=", mu), "-s", seed, "src/abc.slim"), stdout=TRUE)
  as.numeric(output[length(output)])
}

runSLiM(c(100030, 1e-7))

# Set up and run ABC
prior <- list(c("unif", 1e-9, 1e-6))
observed <- 281

abcSLiM <- ABC_sequential(method = "Lenormand", use_seed = TRUE, model = runSLiM, 
                          prior = prior, summary_stat_target = observed, nb_simul = 500)

sum(abcSLiM$param * abcSLiM$weights)

log_parm <- log(abcSLiM$param, 10)
breaks <- seq(from=min(log_parm), to=max(log_parm), length.out = 8)

pdf(file = "results/mu_posterior.pdf")
hist(log_parm, xlim = c(-9,-6), breaks = breaks, col = "grey",
     main="Posterior distribution of mu", xlab = "Estimate of mu", xaxt="n")
axis(side=1, at=-6:-9, labels = c("1e-6", "1e-7", "1e-8", "1e-9"))
dev.off()
