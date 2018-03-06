#################################################
## Tracking Selection in Time-Series Data with ## 
##          Forward-time Simulations           ##
#################################################

## Vitor Pavinato & Miguel Navascues
## vitor.pavinato@supagro.fr
## INRA

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("abc","abcrf","plyr","dplyr"));
#library(abc);
library(abcrf);
library(plyr); # for the convertEggLib function & join function
library(dplyr); # for the convertEggLib function & join function
library(ggplot2);

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;
#sessionInfo()$otherPkgs$plyr$Version;

recordSessionInfo <- sessionInfo();
working_dir <- "/home/pavinato/My_repositories/Tracking-selection/SLiM_in_R";

if (is.na(match(working_dir,getwd()))){
  setwd(working_dir)
} else {
  print("Working directory have already been set")
};

save(recordSessionInfo,file="results/run1-cluster/sessionInfo.RData");
save.image(file = "results/run1-cluster/workspaceFile.RData");
load(file = "results/run1-cluster/workspaceFile.RData");

#################################################
##       Runing Forward-time Simulations       ##
#################################################

### Creating SLiM infiles

## Set the number of simulations

#nsim <- 10000; 

## Set values for the prior - only for Ne

#sim_ne <- as.integer(runif(nsim, 100, 1000));
#sim_seed <- as.integer(runif(nsim, 100000000, 900000000));

#priors <- as.data.frame(cbind(sim=1:nsim, seed=sim_seed,ne=sim_ne));
#hist(priors$ne);

#write.table(priors, file = "results/run1-cluster/ref-tables/sim_priors.txt", quote = FALSE, row.names = FALSE);

priors <- read.table(file = "results/run1-cluster/ref-tables/sim_priors.txt", header = TRUE);
dim(priors);

priors <- priors[1:3412, ];
pdf(file="results/run1-cluster/figures/Ne_prior.pdf");
hist(priors$ne, col = "#f0f0f0", main = "",xlab = expression("Simulated N"[e]));
dev.off();

save.image(file = "results/run1-cluster/workspaceFile.RData");

## Create SLiM infiles using values taken from the prior

#load(file = "results/run1-cluster/workspaceFile.RData");
#source("src/create_slim_infile.r");

#ifelse(file.exists(dir("results/run1-cluster/slim-infiles/", pattern = ".slim$", full.names = TRUE)), 
#       file.remove(dir("results/run1-cluster/slim-infiles/", pattern = ".slim$", full.names = TRUE)), NULL);
      
#createslim <- model1(      ne = sim_ne,
#                     filename = "infile_slim", 
#                       folder = "results/run1-cluster/slim-infiles/",
#                       output = "output_slim",
#                    folderout = "results/run1-cluster/slim-outputs/");

#priors <- cbind(priors, createslim); 
#priors <- priors[,-4];

#save.image(file = "results/run1-cluster/workspaceFile.RData");

### Running SLiM simulations 

#load(file = "results/run1-cluster/workspaceFile.RData");
#source("src/slim.r");

#ifelse(file.exists(dir("results/run1-cluster/slim-outputs/", pattern = "output_slim_", full.names = TRUE)), 
#       file.remove(dir("results/run1-cluster/slim-outputs/", pattern = "output_slim_", full.names = TRUE)), NULL);

#ifelse(file.exists(dir("results/run1-cluster/slim-outputs/", pattern = "outfile_slim_", full.names = TRUE)), 
#       file.remove(dir("results/run1-cluster/slim-outputs/", pattern = "outfile_slim_", full.names = TRUE)), NULL);

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

#slim(     seed = sim_seed,
#       outfile = "outfile_slim",
#     folderout = "results/run1-cluster/slim-outputs/",
#        infile = "infile_slim",
#      folderin = "results/run1-cluster/slim-infiles/");

# SLiM2 runs v2:
# It runs SLiM2 and saves only the outputs;

#save.image(file = "results/run1-cluster/workspaceFile.RData");

#################################################
##     Calculating the Summary Statistics      ##
#################################################

### Converting SLiM outputs to EggLib inputs

load(file = "results/run1-cluster/workspaceFile.RData");
source("src/create_ref_table.r");

#ifelse(file.exists(dir("results/run1-cluster/egglib-inputs/", pattern = "input_egglib_", full.names = TRUE)), 
#       file.remove(dir("results/run1-cluster/egglib-inputs/", pattern = "input_egglib_", full.names = TRUE)), NULL);

dataconvertion <- convertEggLib(    simlist = priors$sim,
                                select_freq = 0.10,
                                     output = "input_egglib" , 
                                  folderout = "results/run1-cluster/egglib-inputs/", 
                                      input = "output_slim", 
                                   folderin = "results/run1-cluster/slim-outputs/",
                                     remove = FALSE);

save.image(file = "results/run1-cluster/workspaceFile.RData");

### Running EggLib

#ifelse(file.exists(dir('results/egglig-outputs/run1-cluster/', pattern = "output_egglib_", full.names = TRUE)), 
#       file.remove(dir('results/egglig-outputs/run1-cluster/', pattern = "output_egglib_", full.names = TRUE)), NULL)

load(file = "results/run1-cluster/workspaceFile.RData");

runEggLib(   simlist = priors$sim,
              output = "output_egglib", 
           folderout = "results/run1-cluster/egglib-outputs/", 
               input = "input_egglib",
            folderin = "results/run1-cluster/egglib-inputs/",
          pythonpath = "/usr/bin/python",
          sourcepath = "src/summstats.py",
                 lss = c("He", "Dj", "WCst"),
                 wss = c("S","thetaW","D", "Da", "ZZ"),
                 gss = c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"),
               wspan = 100, select = "all");

save.image(file = "results/run1-cluster/workspaceFile.RData");

#################################################
##                    ABC-RF                   ##
#################################################

load(file = "results/run1-cluster/workspaceFile.RData");

### Creating the SIMULATIONS Reference Tables

## Global Reference Table - GRT

globalrf <- createGRT( simpriors = priors, 
                          output = "output_egglib", 
                       outfolder = "results/run1-cluster/egglib-outputs/",
                        gss_list = c("GSS.He", "GSS.Dj","GSS.WCst", "GSS.S", "GSS.thetaW","GSS.D", "GSS.Da","GSS.ZZ"));

write.table(globalrf, file = "results/run1-cluster/ref-tables/global_rawrf.txt", quote = FALSE, row.names = FALSE);

globalrf <- read.table(file = "results/run1-cluster/ref-tables/global_rawrf.txt", header = TRUE);

grf <- data.frame(ne=globalrf$ne, He=globalrf$GSS.He, Dj=globalrf$GSS.Dj, 
                  WCst=globalrf$GSS.WCst, S=globalrf$GSS.S, thetaW=globalrf$GSS.thetaW, 
                  D=globalrf$GSS.D, Da=globalrf$GSS.Da);


## Locus-specific Reference Table - LRF

locusrf <- createLRT( simpriors = priors,
                       datalist = dataconvertion,
                          input = "input_egglib",
                       infolder = "results/run1-cluster/egglib-inputs/",
                         output = "output_egglib", 
                      outfolder = "results/run1-cluster/egglib-outputs/",
                         remove = FALSE);

write.table(locusrf, file = "results/run1-cluster/ref-tables/locus_rawrf.txt", quote = FALSE, row.names = FALSE);

locusrf <- read.table(file = "results/run1-cluster/ref-tables/locus_rawrf.txt", header = TRUE);

selected <- as.factor(ifelse(locusrf$ID == "chr1:1", "A", "N"));

lrf <- data.frame(selected=selected, lsHe=locusrf$LSS.He, lsDj=locusrf$LSS.Dj, lsWCst=locusrf$LSS.WCst, 
                  wsS=locusrf$WSS.S, wsthetaW=locusrf$WSS.thetaW, wsD=locusrf$WSS.D, wsDa=locusrf$WSS.Da, 
                  gsHe=locusrf$GSS.He, gsDj=locusrf$GSS.Dj, gsWCst=locusrf$GSS.WCst, gsS=locusrf$GSS.S, 
                  gsthetaW=locusrf$GSS.thetaW, gsD=locusrf$GSS.D, gsDa=locusrf$GSS.Da);

(dim(lrf[which(lrf$selected == "A"),])[1]/3214)*100;

save.image(file = "results/run1-cluster/workspaceFile.RData");

### REGRESSION - Parameter estimation

load(file = "results/run1-cluster/workspaceFile.RData");

regabc <- regAbcrf(formula = ne~., 
                   data = grf,
                   ntree = 2000,
                   paral = TRUE);

pdf(file = "results/run1-cluster/figures/regression_abc_vi.pdf")
plot(regabc)  
dev.off();

regabc$model.rf 

pdf(file = "results/run1-cluster/figures/regression_abc_mse.pdf")
err.regAbcrf(  object = regabc,
               training = grf,
               paral = TRUE)
dev.off();

### CLASSIFICATION

claabc <- abcrf(  formula = selected~.,
                  data    = lrf,
                  lda     = TRUE,
                  ntree   = 2000,
                  paral   = TRUE);

pdf(file = "results/run1-cluster/figures/classification_abc_vi.pdf")
plot(claabc, training=lrf)
dev.off();

claabc$prior.err 

pdf(file = "results/run1-cluster/figures/classification_abc_mse.pdf")
err.abcrf(claabc,
          training=lrf,
          paral=TRUE)
dev.off();

save.image(file = "results/run1-cluster/workspaceFile.RData");

### Creating the OBSERVED Reference Table

load(file = "results/run1-cluster/workspaceFile.RData");

#dataconvertion1 <- convertEggLib(   simlist = 1,
#                                    select_freq = 0.10,
#                                    output = "dataset_egglib" , 
#                                    folderout = "data/egglib-inputs/", 
#                                    input = "dataset", 
#                                    folderin = "data/slim-outputs/",
#                                    remove = FALSE);

dataconvertion2 <- convertEggLib(   simlist = 2,
                                    select_freq = 0.10,
                                    output = "dataset_egglib" , 
                                    folderout = "data/egglib-inputs/", 
                                    input = "dataset", 
                                    folderin = "data/slim-outputs/",
                                    remove = FALSE);

runEggLib(   simlist = 2,
             output = "dataset_output_egglib", 
             folderout = "data/egglib-outputs/", 
             input = "dataset_input_egglib",
             folderin = "data/egglib-inputs/",
             pythonpath = '/usr/bin/python',
             sourcepath = 'src/summstats.py',
             lss = c("He", "Dj", "WCst"),
             wss = c("S","thetaW","D", "Da", "ZZ"),
             gss = c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"),
             wspan = 100, select = "all");


#data <- read.table(file = "data/egglib-outputs/dataset_output_egglib_1.txt", header = TRUE);

#obsgrf <- data.frame(He=data[1,]$GSS.He, Dj=data[1,]$GSS.Dj, WCst=data[1,]$GSS.WCst, S=data[1,]$GSS.S, 
#                     thetaW=data[1,]$GSS.thetaW, D=data[1,]$GSS.D, Da=data[1,]$GSS.Da);

data2 <- read.table(file = "data/egglib-outputs/dataset_output_egglib_2.txt", header = TRUE);

obsgrf2 <- data.frame(He=data2[1,]$GSS.He, Dj=data2[1,]$GSS.Dj, WCst=data2[1,]$GSS.WCst, S=data2[1,]$GSS.S, 
                     thetaW=data2[1,]$GSS.thetaW, D=data2[1,]$GSS.D, Da=data2[1,]$GSS.Da);


## Locus-specific Reference Table - LRF

#obslrf <- data.frame(lsHe=data$LSS.He, lsDj=data$LSS.Dj, lsWCst=data$LSS.WCst, 
#                  wsS=data$WSS.S, wsthetaW=data$WSS.thetaW, wsD=data$WSS.D, wsDa=data$WSS.Da, 
#                  gsHe=data$GSS.He, gsDj=data$GSS.Dj, gsWCst=data$GSS.WCst, gsS=data$GSS.S, 
#                  gsthetaW=data$GSS.thetaW, gsD=data$GSS.D, gsDa=data$GSS.Da);

obslrf2 <- data.frame(lsHe=data2$LSS.He, lsDj=data2$LSS.Dj, lsWCst=data2$LSS.WCst, 
                     wsS=data2$WSS.S, wsthetaW=data2$WSS.thetaW, wsD=data2$WSS.D, wsDa=data2$WSS.Da, 
                     gsHe=data2$GSS.He, gsDj=data2$GSS.Dj, gsWCst=data2$GSS.WCst, gsS=data2$GSS.S, 
                     gsthetaW=data2$GSS.thetaW, gsD=data2$GSS.D, gsDa=data2$GSS.Da);


save.image(file = "results/run1-cluster/workspaceFile.RData");

## Prediction

## Based on Regression ABC - Parameter estimation

#pred.par <- predict(object = regabc,
#                       obs = obsgrf,
#                  training = grf,
#                     paral = TRUE);

#pred.par$expectation
#pred.par$med
#pred.par$quantiles

#pdf(file = "results/run1-cluster/figures/regression_abc_posterior.pdf")
#densityPlot(object    = regabc,
#            obs       = obsgrf,
#            training  = grf,
#            main      = expression("N"[e]),
#            paral     = TRUE)

#lines(density(priors$ne), col="grey")
#abline(v=500, col = "red")
#dev.off();

pred.par2 <- predict(object = regabc,
                    obs = obsgrf2,
                    training = grf,
                    paral = TRUE);

pred.par2$expectation
pred.par2$med
pred.par2$quantiles

pdf(file = "results/run1-cluster/figures/regression_abc_posterior_data2.pdf")
densityPlot(object    = regabc,
            obs       = obsgrf2,
            training  = grf,
            main      = expression("N"[e]),
            paral     = TRUE)

lines(density(priors$ne), col="grey")
abline(v=1500, col = "red")
dev.off();

# Based on Classification ABC - Identification advantage mutation

#selected_locus <- predict(       object  = claabc,
#                                    obs  = obslrf,
#                                training = lrf,
#                                   ntree = 1000,
#                                   paral = TRUE,
#                          paral.predict  = TRUE);

#(selected_locus)

#selected_locus$allocation[which(selected_locus$allocation == "A") ];

#pdf(file = "results/run1-cluster/figures/classification_abc_posterior.pdf")
#plot(selected_locus$vote, main = "Loci Classification",
#     col = ifelse(selected_locus$allocation == "A", "red", "blue"),
#     pch = ifelse(selected_locus$allocation == "A", 19, 1))
#abline(v=0.95, col = "red", lty = 2)
#dev.off();

#post_probA <- ifelse(selected_locus$allocation == "A", 
#                     selected_locus$post.prob, 
#                     (1-selected_locus$post.prob));

#chr_pos <- as.character(data$ID);
#chr_splitted <- strsplit(chr_pos, "\\:"); chr_splitted <- do.call(rbind, chr_splitted);


#genescan <- data.frame(   snp = 1:length(chr_splitted[,2]),
#                          chr = ifelse(chr_splitted[,1]=="chr1", 1, 2), 
#                         pos = as.numeric(chr_splitted[,2]),
#                       classy = selected_locus$allocation,
#                         prob = post_probA, 
#                          fst = obslrf$lsWCst,
#                           he = obslrf$lsHe);

#genescan[which(genescan$classy == "A"), ]

#pdf(file = "results/run1-cluster/figures/classification_abc_manhattan.pdf")

#p_abc <- ggplot(genescan, aes(y = prob, x = pos)) +
#     geom_point(aes(colour = factor(chr))) +
#     scale_color_manual(values=c("black", "grey")) + 
#     facet_wrap("chr") +
#     xlab("Position") + 
#     ylab("Posterior probability")
#p_abc
#dev.off();

#pdf(file = "results/run1-cluster/figures/classification_abc_manhattan_fst.pdf")

#p_fst <- ggplot(genescan, aes(y = fst, x = pos)) +
#  geom_point(aes(colour = factor(chr))) +
#  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
#  facet_wrap("chr") +
#  xlab("Position") + 
#  ylab(expression("F"[ST])) +
#  geom_hline(yintercept = 0.211504, color = "#fb8072")
#p_fst
#dev.off();

#pdf(file = "results/run1-cluster/figures/classification_abc_fstout.pdf")

#p_out <- ggplot(genescan, aes(y = fst, x = he)) +
#  geom_point(aes(colour = factor(chr))) +
#  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
#  facet_wrap("chr") +
#  xlab(expression("H"[E])) + 
#  ylab(expression("F"[ST]))
#p_out
#dev.off();

# Plot snp along the chromosome

#chrstruct <- structure(list(chromosome = 1:2, size = c(50000L, 50000L)), 
#                 .Names = c("chr", "size"), class = "data.frame", row.names = c(NA, -2L))

#classycol <- ifelse(genescan$classy == "N", "#74a9cf", "#e34a33")

#pdf(file = "results/run1-cluster/figures/classification_abc_chr.pdf")
#bp <- barplot(chrstruct$size, border=NA, col="grey80", names.arg = c("chr 1", "chr 2"))
#with(genescan,segments(bp[chr,]-0.5, pos, bp[chr,]+0.5, pos, col=classycol, lwd=2, lend=1))
#dev.off();

selected_locus2 <- predict(      object  = claabc,
                                 obs  = obslrf2,
                                 training = lrf,
                                 ntree = 1000,
                                 paral = TRUE,
                                 paral.predict  = TRUE);

(selected_locus2)

selected_locus2$allocation[which(selected_locus2$allocation == "A") ];

#pdf(file = "results/run1-cluster/figures/classification_abc_posterior_data2.pdf")
#plot(selected_locus2$vote, main = "Loci Classification",
#     col = ifelse(selected_locus2$allocation == "A", "red", "blue"),
#     pch = ifelse(selected_locus2$allocation == "A", 19, 1))
#abline(v=0.95, col = "red", lty = 2)
#dev.off();

post_probA2 <- ifelse(selected_locus2$allocation == "A", 
                     selected_locus2$post.prob, 
                     (1-selected_locus2$post.prob));

chr_pos2 <- as.character(data2$ID);
chr_splitted2 <- strsplit(chr_pos2, "\\:"); chr_splitted2 <- do.call(rbind, chr_splitted2);


genescan2 <- data.frame(  snp = 1:length(chr_splitted2[,2]),
                          chr = ifelse(chr_splitted2[,1]=="chr1", 1, 2), 
                          pos = as.numeric(chr_splitted2[,2]),
                          classy = selected_locus2$allocation,
                          prob = post_probA2, 
                          fst = obslrf2$lsWCst,
                          he = obslrf2$lsHe);

genescan2[which(genescan2$classy == "A"), ]

pdf(file = "results/run1-cluster/figures/classification_abc_manhattan_data2.pdf")

p_abc <- ggplot(genescan2, aes(y = prob, x = pos)) +
  geom_point(aes(colour = factor(chr))) +
  scale_color_manual(values=c("black", "grey")) + 
  facet_wrap("chr") +
  xlab("Position") + 
  ylab("Posterior probability")
p_abc
dev.off();

pdf(file = "results/run1-cluster/figures/classification_abc_manhattan_fst_data2.pdf")

p_fst <- ggplot(genescan2, aes(y = fst, x = pos)) +
  geom_point(aes(colour = factor(chr))) +
  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
  facet_wrap("chr") +
  xlab("Position") + 
  ylab(expression("F"[ST])) +
  geom_hline(yintercept = 0.211504, color = "#fb8072")
p_fst
dev.off();

pdf(file = "results/run1-cluster/figures/classification_abc_fstout_data2.pdf")

p_out <- ggplot(genescan2, aes(y = fst, x = he)) +
  geom_point(aes(colour = factor(chr))) +
  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
  facet_wrap("chr") +
  xlab(expression("H"[E])) + 
  ylab(expression("F"[ST]))
p_out
dev.off();

# Plot snp along the chromosome

chrstruct2 <- structure(list(chromosome = 1:2, size = c(50000L, 50000L)), 
                       .Names = c("chr", "size"), class = "data.frame", row.names = c(NA, -2L))

classycol2 <- ifelse(genescan2$classy == "N", "#74a9cf", "#e34a33")

pdf(file = "results/run1-cluster/figures/classification_abc_chr_data2.pdf")
bp2 <- barplot(chrstruct2$size, border=NA, col="grey80", names.arg = c("chr 1", "chr 2"))
with(genescan2,segments(bp2[chr,]-0.5, pos, bp2[chr,]+0.5, pos, col=classycol2, lwd=2, lend=1))
dev.off();


# Relationship between genetic distance and recombination rate (or probability)

# Distance in centimorgans

d <- function(rrate) {50*log(1/(1-2*rrate))}

d(1e-08)

# 1cM = 1000000 bp in humans

save.image(file = "results/run1-cluster/workspaceFile.RData");
