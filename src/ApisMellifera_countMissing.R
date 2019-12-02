## Avalon population
##------------------

avalon_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/avalon_merged.txt", header = F)
# overall average - mean per sample
mean(apply(avalon_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) # 0.1082116

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(avalon_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(833626 * 7) # 0.1082116

## Humboldt population
##------------------

humboldt_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/humboldt_merged.txt", header = F)
# overall average - mean per sample
mean(apply(humboldt_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) # 0.003029194

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(humboldt_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(834683 * 12) # 0.003029194

## Davis population
##------------------

davis_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/davis_merged.txt", header = F)
# overall average - mean per sample
mean(apply(davis_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) #  0.004315563

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(davis_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(833528 * 7) # 0.004315563

## Stanislaus population
##------------------

stanislaus_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/stanislaus_merged.txt", header = F)
# overall average - mean per sample
mean(apply(stanislaus_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) # 0.00222511

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(stanislaus_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(834060 * 8) # 0.00222511

## Stebbins population
##------------------

stebbins_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/stebbins_merged.txt", header = F)
# overall average - mean per sample
mean(apply(stebbins_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) #  0.005241321

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(stebbins_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(834160 * 10) # 0.005241321

## Riverside population
##------------------

riverside_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/riverside_merged.txt", header = F)
# overall average - mean per sample
mean(apply(riverside_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) #  0.002839675

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(riverside_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(834532 * 10) # 0.002839675

## Placerita population
##------------------

placerita_d <- read.table(file="data/ApisMellifera/PRNJA385500_vcf_files/placerita_merged.txt", header = F)
# overall average - mean per sample
mean(apply(placerita_d[,-c(1:4)], 2, function(x){sum(x == "./.")/length(x)})) #  0.1012176

# overall proportion - overall sum of missing by total number of genotypes
sum(apply(placerita_d[,-c(1:4)], 2, function(x){sum(x == "./.")}))/(833970 * 11) # 0.1012176
