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

# Simulating from prior;
fromprior <- function(nsims){
                   ne <- as.integer(runif(nsims, 500, 1000));
                 mud1 <- runif(nsims, 0.5, 1.0);
                 mud2 <- runif(nsims, 0.5, 1.0);
                mufv2 <- runif(nsims, 0, 0.5);
                  ts1 <- as.integer(runif(nsims, 1, 1000));
                int12 <- as.integer(runif(nsims, 1, 1000));
                  
                result <- cbind(ne, mud1, mud2, mufv2, ts1, int12);
                return(result)
                }

parm <- fromprior(10)

for (i in 1:nsim){
  cat(paste("You simulated sample",i,"of size N =", parm[,1][i],"\n"));
  cat(paste("You simulated selcoeff ",i," =", parm[,4][i],"\n"));
};

ifelse(file.exists(dir("infiles/", pattern = ".slim$", full.names = TRUE)), 
       file.remove(dir("infiles/", pattern = ".slim$", full.names = TRUE)), NULL)
      
model1(   nsims = 1,
       filename = "infile_slim", 
         folder = "infiles/");

## Running SLiM simulations 

source("src/slim.r");

# Random seeds;
random_seeds <- as.integer(runif(1, 100000, 900000));

ifelse(file.exists(dir('data/', pattern = "output_", full.names = TRUE)), 
       file.remove(dir('data/', pattern = "output_", full.names = TRUE)), NULL)

# SLiM2 runs v1:
# It runs SLiM2 and saves the parameters of each simulation in an outfile;
# Thanks to Eidos code each simulation now generate 2 files for each sampling;

slim(     nsims = 1, 
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

slimclean(   nsims = 1, 
              seed = random_seeds,
            infile = 'infile_slim',      
          folderin = 'infiles/');

# system time: 10 simulation Ne < 1000
# user  system elapsed 
# 18.895   0.155  19.081

## Uploading SLiM simulations outputs
library(plyr)

infot1 <- read.table("data/output_slim_1_mutInfo_t1.txt", header = TRUE)
head(infot1)
dim(infot1)[1]

genot1 <- read.table("data/output_slim_1_genotypes_t1.txt", header = FALSE)
head(genot1)
genot1 <- t(genot1)
dim(genot1)[2]

mindex = seq(from=1, to=dim(genot1)[1], by=1)
rownames(genot1) <- paste(mindex)
iindex = seq(from=1, to=dim(genot1)[2], by=1)
colnames(genot1) <- paste0("indiv",iindex,"@pop1", "")

sp1 <- cbind(infot1,genot1)

infot2 <- read.table("data/output_slim_1_mutInfo_t2.txt", header = TRUE)
head(infot2)
dim(infot2)[1]

genot2 <- read.table("data/output_slim_1_genotypes_t2.txt", header = FALSE)
head(genot2)
genot2 <- t(genot2)
dim(genot2)[2]

mindex = seq(from=1, to=dim(genot2)[1], by=1)
rownames(genot2) <- paste(mindex)
iindex = seq(from=1, to=dim(genot2)[2], by=1)
colnames(genot2) <- paste0("indiv",iindex,"@pop2", "")

sp2 <- cbind(infot2,genot2)

test <- merge(sp1, sp2, by.x= c(1,2,3,4,5), by.y = c(1,2,3,4,5), all=TRUE)
test[is.na(test)] <- 00
dim(test)

colnames(test)[3] <- "status"
test[,3]<- rep("NS", times=dim(test)[1])

colnames(test)[5] <- "alleles"
test[,5]<- rep("A,T", times=dim(test)[1])

write.table(test, file = "testmerge.txt", quote=FALSE, sep="\t", row.names = FALSE)

