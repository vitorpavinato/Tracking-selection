#################################################
## Tracking Selection in Time-Series Data with ## 
##          Forward-time Simulations           ##
#################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("abc","abcrf"));
#library(abc);
#library(abcrf);

#sessionInfo()$otherPkgs$abc$Version;
#sessionInfo()$otherPkgs$abcrf$Version;

recordSessionInfo <- sessionInfo();
working_dir <- "/Users/vitorpavinato/Documents/My_repositories/Tracking-selection/SLiM_in_R"

if (is.na(match(working_dir,getwd()))){
  setwd(working_dir)
} else {
  print("Working directory already set up")
}

save(recordSessionInfo,file="results/sessionInfo.RData");
save.image(file = "results/workspaceFile.RData");

#################################################
##       Runing Forward-time Simulations       ##
#################################################

## Creating SLiM infiles

source("src/create_slim_infile.r");

# Simulating N from prior;
N_size <- as.integer(runif(10, 500, 1000));

for (n in 1:length(N_size)){
cat(paste("You simulated sample",n,"of size N =", N_size[n],"\n"));
};

ifelse(file.exists(dir("infiles/", pattern = ".slim$", full.names = TRUE)), 
       file.remove(dir("infiles/", pattern = ".slim$", full.names = TRUE)), NULL)
      
model1(      ne = N_size,
       filename = "infile_slim", 
         folder = "infiles/");

## Running SLiM simulations 

source("src/slim.r");

# Random seeds;
random_seeds <- as.integer(runif(length(N_size), 100000, 900000));

ifelse(file.exists(dir('data/', pattern = "output_", full.names = TRUE)), 
       file.remove(dir('data/', pattern = "output_", full.names = TRUE)), NULL)

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

slim(     parm = N_size, 
          seed = random_seeds,
      filename = 'outfile_slim',
     folderout = 'data/',
        infile = 'infile_slim',
      folderin = 'infiles/');

# system time: 10 simulation Ne < 1000
# user  system elapsed 
# 17.502   0.146  17.690

# SLiM2 runs v2:
# It runs SLiM2 and saves only the output;

slimclean(    parm = N_size, 
              seed = random_seeds,
            infile = 'infile_slim',      
          folderin = 'infiles/');

# system time: 10 simulation Ne < 1000
# user  system elapsed 
# 18.895   0.155  19.081

## Uploading SLiM simulations outputs

infot1 <- read.table("data/output_slim_N_520_mutInfo_t1.txt", header = TRUE)
dim(infot1)

genot1 <- read.table("data/output_slim_N_520_genotypes_t1.txt", header = FALSE)
genot1 <- t(genot1)

index = seq(from=1, to=length(genot1[1,]), by=1)
colnames(genot1) <- paste0("indiv",index,"@pop1", "")

cbind(infot1,genot1)
