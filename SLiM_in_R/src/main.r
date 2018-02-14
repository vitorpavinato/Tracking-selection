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
working_dir <- "/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_in_R";

if (is.na(match(working_dir,getwd()))){
  setwd(working_dir)
} else {
  print("Working directory have already been set")
};

#save(recordSessionInfo,file="results/sessionInfo.RData");
#save.image(file = "results/workspaceFile.RData");

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

## Create SLiM infiles using values taken from the prior

#source("src/create_slim_infile.r");

#ifelse(file.exists(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), 
#       file.remove(dir("results/slim-infiles/", pattern = ".slim$", full.names = TRUE)), NULL);
      
#createslim <- model1(      ne = sim_ne,
#                     filename = "infile_slim", 
#                       folder = "results/slim-infiles/",
#                       output = "output_slim",
#                    folderout = "results/slim-outputs");

#priors <- cbind(priors, createslim); 
#priors <- priors[,-4];

#write.table(priors, file = "results/ref-tables/sim_priors.txt", quote = FALSE, row.names = FALSE);

priors <- read.table(file = "results/ref-tables/run1/sim_priors.txt", header = TRUE);
#dim(priors);

priors <- priors[1:1000, ];
#pdf(file="results/figures/run1/Ne_prior_1000.pdf");
#hist(priors$ne, col = "#f0f0f0", main = "",xlab = expression("Simulated N"[e]));
#dev.off();

### Running SLiM simulations 

#source("src/slim.r");

#ifelse(file.exists(dir('results/slim-outputs/', pattern = "output_slim_", full.names = TRUE)), 
#       file.remove(dir('results/slim-outputs/', pattern = "output_slim_", full.names = TRUE)), NULL);

#ifelse(file.exists(dir('results/slim-outputs/', pattern = "outfile_slim_", full.names = TRUE)), 
#       file.remove(dir('results/slim-outputs/', pattern = "outfile_slim_", full.names = TRUE)), NULL);

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

#slim(     seed = sim_seed,
#       outfile = 'outfile_slim',
#     folderout = 'results/slim-outputs/',
#        infile = 'infile_slim',
#      folderin = 'results/slim-infiles/');

# SLiM2 runs v2:
# It runs SLiM2 and saves only the outputs;

# slimclean(    seed = sim_seed, 
#             infile = 'infile_slim',      
#           folderin = 'infiles/');

#save.image(file = "results/workspaceFile.RData");

#################################################
##     Calculating the Summary Statistics      ##
#################################################

#load(file = "results/workspaceFile.RData");

### Converting SLiM outputs to EggLib inputs

#source("src/create_ref_table.r");

#ifelse(file.exists(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), 
#       file.remove(dir('results/egglib-inputs/', pattern = "input_egglib_", full.names = TRUE)), NULL);

#dataconvertion <- convertEggLib(    simlist = priors$sim,
#                                select_freq = 0.10,
#                                     output = "input_egglib" , 
#                                  folderout = "results/egglib-inputs/run1/", 
#                                      input = "output_slim", 
#                                   folderin = "results/slim-outputs/run1/",
#                                     remove = FALSE);

### Running EggLib

#ifelse(file.exists(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), 
#       file.remove(dir('results/egglig-outputs/', pattern = "output_egglib_", full.names = TRUE)), NULL)

#runEggLib(   simlist = priors$sim,
#              output = "output_egglib", 
#           folderout = "results/egglig-outputs/run1/", 
#               input = "input_egglib",
#            folderin = "results/egglib-inputs/run1/",
#          pythonpath = '/Users/vitorpavinato/anaconda/bin/python',
#          sourcepath = 'src/summstats.py',
#                 lss = c("He", "Dj", "WCst"),
#                 wss = c("S","thetaW","D", "Da", "ZZ"),
#                 gss = c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"),
#               wspan = 100, select = "all");

#save.image(file = "results/workspaceFile.RData");

#################################################
##                    ABC-RF                   ##
#################################################

#load(file = "results/workspaceFile.RData");

### Creating the SIMULATIONS Reference Tables

## Global Reference Table - GRT

#globalrf <- createGRT( simpriors = priors, 
#                          output = "output_egglib", 
#                       outfolder = "results/egglig-outputs/run1/",
#                        gss_list = c("GSS.He", "GSS.Dj","GSS.WCst", "GSS.S", "GSS.thetaW","GSS.D", "GSS.Da","GSS.ZZ"));

#write.table(globalrf, file = "results/ref-tables/run1/global_rawrf.txt", quote = FALSE, row.names = FALSE);

globalrf <- read.table(file = "results/ref-tables/run1/global_rawrf.txt", header = TRUE);

grf <- data.frame(ne=globalrf$ne, He=globalrf$GSS.He, Dj=globalrf$GSS.Dj, 
                  WCst=globalrf$GSS.WCst, S=globalrf$GSS.S, thetaW=globalrf$GSS.thetaW, 
                  D=globalrf$GSS.D, Da=globalrf$GSS.Da);

## Locus-specific Reference Table - LRF

#locusrf <- createLRT( simpriors = priors,
#                       datalist = dataconvertion,
#                          input = "input_egglib",
#                       infolder = "results/egglib-inputs/run1/",
#                         output = "output_egglib", 
#                      outfolder = "results/egglig-outputs/run1/",
#                         remove = FALSE);

#write.table(locusrf, file = "results/ref-tables/run1/locus_rawrf.txt", quote = FALSE, row.names = FALSE);

locusrf <- read.table(file = "results/ref-tables/run1/locus_rawrf.txt", header = TRUE);

selected <- as.factor(ifelse(locusrf$ID == "chr1:1", "A", "N"));

lrf <- data.frame(selected=selected, lsHe=locusrf$LSS.He, lsDj=locusrf$LSS.Dj, lsWCst=locusrf$LSS.WCst, 
                  wsS=locusrf$WSS.S, wsthetaW=locusrf$WSS.thetaW, wsD=locusrf$WSS.D, wsDa=locusrf$WSS.Da, 
                  gsHe=locusrf$GSS.He, gsDj=locusrf$GSS.Dj, gsWCst=locusrf$GSS.WCst, gsS=locusrf$GSS.S, 
                  gsthetaW=locusrf$GSS.thetaW, gsD=locusrf$GSS.D, gsDa=locusrf$GSS.Da);

(dim(lrf[which(lrf$selected == "A"),])[1]/1000)*100;

#save.image(file = "results/workspaceFile.RData");

### REGRESSION - Parameter estimation

regabc <- regAbcrf(formula = ne~., 
                   data = grf,
                   ntree = 500,
                   paral = TRUE);

pdf(file = "results/figures/run1/regression_abc_vi_orig.pdf")
plot(regabc)  
dev.off();

regabc$model.rf 

pdf(file = "results/figures/run1/regression_abc_mse_orig.pdf")
err.regAbcrf(  object = regabc,
               training = grf,
               paral = TRUE)
dev.off();

### CLASSIFICATION

claabc <- abcrf(  formula = selected~.,
                  data    = lrf,
                  lda     = FALSE,
                  ntree   = 500,
                  paral   = TRUE);

pdf(file = "results/figures/run1/classification_abc_vi.pdf")
plot(claabc, training=lrf)
dev.off();

claabc$prior.err 

pdf(file = "results/figures/run1/classification_abc_mse.pdf")
err.abcrf(claabc,
          training=lrf,
          paral=TRUE)
dev.off();

#save.image(file = "results/workspaceFile.RData");

### Creating the OBSERVED Reference Table

#dataconvertion1 <- convertEggLib(   simlist = 1,
#                                    select_freq = 0.10,
#                                    output = "dataset_egglib" , 
#                                    folderout = "data/egglib-inputs/", 
#                                    input = "dataset", 
#                                    folderin = "data/slim-outputs/",
#                                    remove = FALSE);


#runEggLib(   simlist = 1,
#             output = "dataset_output_egglib", 
#             folderout = "data/egglib-outputs/", 
#             input = "dataset_egglib",
#             folderin = "data/egglib-inputs/",
#             pythonpath = '/Users/vitorpavinato/anaconda/bin/python',
#             sourcepath = 'src/summstats.py',
#             lss = c("He", "Dj", "WCst"),
#             wss = c("S","thetaW","D", "Da", "ZZ"),
#             gss = c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"),
#             wspan = 100, select = "all");

data <- read.table(file = "data/egglib-outputs/dataset_output_egglib_1.txt", header = TRUE);

obsgrf <- data.frame(He=data[1,]$GSS.He, Dj=data[1,]$GSS.Dj, WCst=data[1,]$GSS.WCst, S=data[1,]$GSS.S, 
                     thetaW=data[1,]$GSS.thetaW, D=data[1,]$GSS.D, Da=data[1,]$GSS.Da);

## Locus-specific Reference Table - LRF

obslrf <- data.frame(lsHe=data$LSS.He, lsDj=data$LSS.Dj, lsWCst=data$LSS.WCst, 
                  wsS=data$WSS.S, wsthetaW=data$WSS.thetaW, wsD=data$WSS.D, wsDa=data$WSS.Da, 
                  gsHe=data$GSS.He, gsDj=data$GSS.Dj, gsWCst=data$GSS.WCst, gsS=data$GSS.S, 
                  gsthetaW=data$GSS.thetaW, gsD=data$GSS.D, gsDa=data$GSS.Da);


#save.image(file = "results/workspaceFile.RData");

## Prediction

## Based on Regression ABC - Parameter estimation

pred.par <- predict(object = regabc,
                       obs = obsgrf,
                  training = grf,
                     paral = TRUE);

pred.par$med
pred.par$quantiles

pdf(file = "results/figures/run1/regression_abc_posterior_orig.pdf")
densityPlot(object    = regabc,
            obs       = obsgrf,
            training  = grf,
            main      = expression("N"[e]),
            paral     = TRUE)

lines(density(priors$ne), col="grey")
abline(v=500, col = "red")
dev.off();

# Based on Classification ABC - Identification advantage mutation

predict(object = claabc$model.rf,  data = obslrf)$predictions

rowSums(predict(object=claabc$model.rf, predict.all = TRUE,data=obslrf)$predictions==1)

selected_locus <- predict(       object  = claabc,
                                    obs  = obslrf,
                                training = lrf,
                                   ntree = 500,
                                   paral = TRUE,
                          paral.predict  = TRUE);

(selected_locus)

selected_locus$allocation[which(selected_locus$allocation == "A") ];

pdf(file = "results/figures/run1/classification_abc_posterior.pdf")
plot(selected_locus$vote, main = "Loci Classification",
     col = ifelse(selected_locus$allocation == "A", "red", "blue"),
     pch = ifelse(selected_locus$allocation == "A", 19, 1))
abline(v=0.95, col = "red", lty = 2)
dev.off();

post_probA <- ifelse(selected_locus$allocation == "A", 
                     selected_locus$post.prob, 
                     (1-selected_locus$post.prob));

chr_pos <- as.character(data$ID);
chr_splitted <- strsplit(chr_pos, "\\:"); chr_splitted <- do.call(rbind, chr_splitted);


genescan <- data.frame(   snp = 1:length(chr_splitted[,2]),
                          chr = ifelse(chr_splitted[,1]=="chr1", 1, 2), 
                          pos = as.numeric(chr_splitted[,2]),
                       classy = selected_locus$allocation,
                         prob = post_probA, 
                          fst = obslrf$lsWCst,
                           he = obslrf$lsHe);

genescan[which(genescan$classy == "A"), ]

pdf(file = "results/figures/run1/classification_abc_manhattan.pdf")

p <- ggplot(genescan, aes(y = prob, x = pos)) +
     geom_point(aes(colour = factor(chr))) +
     scale_color_manual(values=c("black", "grey")) + 
     facet_wrap("chr") +
     xlab("Position") + 
     ylab("Posterior probability")
p
dev.off();

pdf(file = "results/figures/run1/classification_abc_manhattan_fst.pdf")

p <- ggplot(genescan, aes(y = fst, x = pos)) +
  geom_point(aes(colour = factor(chr))) +
  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
  facet_wrap("chr") +
  xlab("Position") + 
  ylab(expression("F"[ST])) +
  geom_hline(yintercept = 0.211504, color = "#fb8072")
p
dev.off();

pdf(file = "results/figures/run1/classification_abc_fstout.pdf")

p <- ggplot(genescan, aes(y = fst, x = he)) +
  geom_point(aes(colour = factor(chr))) +
  scale_color_manual(values=c("#e41a1c", "#377eb8")) + 
  facet_wrap("chr") +
  xlab(expression("H"[E])) + 
  ylab(expression("F"[ST]))
p
dev.off();

# Plot snp along the chromosome

chrstruct <- structure(list(chromosome = 1:2, size = c(50000L, 50000L)), 
                 .Names = c("chr", "size"), class = "data.frame", row.names = c(NA, -2L))

classycol <- ifelse(genescan$classy == "N", "#74a9cf", "#e34a33")

pdf(file = "results/figures/run1/classification_abc_chr.pdf")
bp <- barplot(chrstruct$size, border=NA, col="grey80", names.arg = c("chr 1", "chr 2"))
with(genescan,segments(bp[chr,]-0.5, pos, bp[chr,]+0.5, pos, col=classycol, lwd=2, lend=1))
dev.off();

save.image(file = "results/workspaceFile.RData");
