do_sim <- function(sim, nsim, slim_model, path_to_slim, slim_output_folder,
                   path_to_bgzip, path_to_tabix, path_to_bcftools,
                   egglib_input_folder, egglib_output_folder,
                   path_to_python, path_to_egglib_summstat, remove_files,  
                   mu_rate, mu_random, mu_min, mu_max, 
                   ne0_value, ne0_random, ne0_min, ne0_max,
                   ne1_value, ne1_random, ne1_min, ne1_max,  
                   gammaM_value, gammak_value, gammaM_gammak, gammaM_random, gammaM_min, gammaM_max, 
                   gammak_random, gammak_min, gammak_max, 
                   PrGWSel_value, PrGWSel_random, PrGWSel_min, PrGWSel_max, 
                   prbe_value, prbe_random, prbe_min, prbe_max, 
                   domN, domN_random, domN_min, domN_max, 
                   domB, domB_random, domB_min, domB_max, 
                   rr_rate, rr_random, rr_min, rr_max,
                   SS1, SS2, ts2, genomeS, fragS, chrN,
                   wss_wspan, sfs_bins,
                   random_simulations
                   ){
  
  # write progress of simulation on screen 
  if (sim==1 | sim%%10==0 | sim==nsim){
    cat(paste(sim,"of",nsim))  
  }
  
  #######################
  ##   SAMPLED VALUES  ##
  #######################
  ##-------------------
  
  # SLiM SEED
  sim_seed  <- as.integer(runif(n = 1, min = 100000000, max = 900000000));
  
  # MUTATION RATE
  if (mu_random){
    mu <- 10^runif(n = 1, min = log10(mu_min), max = log10(mu_max))
  } else {
    mu <- mu_rate
  }
  
  # THETA AND EFFECTIVE POPULATION SIZE 0
  if (ne0_random){
    theta_min = 4*ne0_min*mu
    theta_max = 4*ne0_max*mu
    theta <- 10^runif(n = 1, min = log10(theta_min), max = log10(theta_max))
    Ne0 <- as.integer(theta/(4*mu))
  } else {
    Ne0 <- ne0_value
    theta <- 4*Ne0*mu
  }
  
  # EFFECTIVE POPULATION SIZE 1
  if (ne1_random){
    Ne1 <- as.integer(10^runif(n = 1, min = log10(ne1_min), max = log10(ne1_max)))  
  } else {
    Ne1 <- ne1_value
  }
  
  # SAMPLING VALUES FOR SLiM GAMMA DISTRIBUTION OF FITNESS EFFECTS
  
  # gamma mean
  if (gammaM_random){
    gammaM <- round(10^runif(n = 1, min = log10(gammaM_min), max = log10(gammaM_max)), digits = 4)
  } else {
    gammaM <- gammaM_value
  }
  
  # gamma shape 
  if (gammak_random){
    if (gammaM_gammak){
      gammak = gammaM
    } else {
      gammak <- round(10^runif(n = 1, min = log10(gammak_min), max = log10(gammak_max)), digits = 4)
    }
  } else {
    gammak <- gammak_value
  }
  
  # GENOME-WIDE PROPORTION OF SELECTED GENOMIC ELEMENTS G2
  if (PrGWSel_random){
    PrGWSel <- round(runif(n = 1, min = PrGWSel_min, max = PrGWSel_max), digits = 6)
  } else {
    PrGWSel <-  PrGWSel_value
  }
  
  # PROPORTION OF GENOME-WIDE BENEFICIAL MUTATION IN GENOMIC ELEMENTS G2
  if (prbe_random){
    prbe <- round(runif(1, min = prbe_min, max = prbe_max), digits = 6)
  } else {
    prbe <- prbe_value
  }
  # DOMINANCE FOR GENOME-WIDE MUTATIONS
  
  # Neutral mutations - m1 and m2
  if (domN_random){
    dm1 = dm2 <- runif(n = 1, min = domN_min, max = domN_max)
  } else {
    dm1 = dm2 <- domN
  }
  
  # Beneficial mutations - m3
  if (domB_random){
    dm3 <- runif(n = 1, min = domB_min, max = domB_max)
  } else {
    dm3 <- domB
  }
  
  # RECOMBINATION RATE
  if (rr_random){
    rr <- 4.2 * (10^runif(1, min = log10(rr_min), max = log10(rr_max)))
  } else {
    rr <- rr_rate
  }
  
  ## RUM SLiM2
  ##-------------------
  
  # check if the folder exists
  if (!file_test("-d", slim_output_folder)){
    dir.create(file.path(slim_output_folder))
  }
  
  # generate text with slim command
  slim_output_fullpath <- paste0("'",getwd(), "/", slim_output_folder, "'")
  
  slim_simID     <- paste0("-d simID=", sim)
  slim_mu        <- paste0("-d mu=", mu)
  slim_theta     <- paste0("-d theta=", theta)
  slim_Ne1       <- paste0("-d Ne1=", Ne1)
  slim_gammaM    <- paste0("-d gammaM=", gammaM)
  slim_gammak    <- paste0("-d gammak=", gammak)
  slim_PrGWSel   <- paste0("-d PrGWSel=", PrGWSel)
  slim_prbe      <- paste0("-d prbe=", prbe)
  slim_dm1       <- paste0("-d dm1=", dm1)
  slim_dm2       <- paste0("-d dm2=", dm2)
  slim_dm3       <- paste0("-d dm3=", dm3)
  slim_rr        <- paste0("-d rr=", rr)
  slim_SS1       <- paste0("-d SS1=", SS1)
  slim_SS2       <- paste0("-d SS2=", SS2)
  slim_ts2       <- paste0("-d ts2=", ts2)
  slim_genomeS   <- paste0("-d genomeS=", genomeS)
  slim_fragS     <- paste0("-d fragS=", fragS)
  slim_chrN      <- paste0("-d chrN=", chrN)
  slim_output    <- paste0("-d ", '"',"outputpath=", slim_output_fullpath, '"')
  
  slim_run <- paste( path_to_slim,
                     "-s", sim_seed,                # seed  = simulation seed number
                     slim_simID,                    # simID = simulation id number
                     slim_mu,                       #    mu = mutation rate
                     slim_theta,                    # theta
                     slim_Ne1,                      #   ne1 = effective population size 1  
                     slim_gammaM,
                     slim_gammak,
                     slim_PrGWSel,
                     slim_prbe,
                     slim_dm1,
                     slim_dm2,
                     slim_dm3,
                     slim_rr,
                     slim_SS1,
                     slim_SS2,
                     slim_ts2,
                     slim_genomeS,
                     slim_fragS,
                     slim_chrN,
                     slim_output,
                     slim_model)                    # model
  
  # rum slim on system
  system(slim_run)
  
  ## HANDLE SLiM2 OUTPUT 1 - GENETIC DATA
  ##-----------------------------------------------------
  
  # sort vcf files
  slim_output_t1        <- paste0(slim_output_folder,"slim_output_t1_", sim, ".vcf")
  slim_output_t1_sorted <- paste0(slim_output_folder,"slim_output_t1_", sim, "_sorted" , ".vcf")
  
  slim_output_t2        <- paste0(slim_output_folder,"slim_output_t2_", sim, ".vcf")
  slim_output_t2_sorted <- paste0(slim_output_folder,"slim_output_t2_", sim, "_sorted" , ".vcf")
  
  sort_t1_vcf <- paste("grep '^#'", slim_output_t1, ">", slim_output_t1_sorted,
                       "&& grep -v '^#'", slim_output_t1, "| sort -k1,1 -k2,2n", ">>", slim_output_t1_sorted)
  
  sort_t2_vcf <- paste("grep '^#'", slim_output_t2, ">", slim_output_t2_sorted,
                       "&& grep -v '^#'", slim_output_t2, "| sort -k1,1 -k2,2n", ">>", slim_output_t2_sorted)
  
  system(sort_t1_vcf)
  system(sort_t2_vcf)
  
  # bgzip sorted vcf files
  system(paste(path_to_bgzip, slim_output_t1_sorted))
  system(paste(path_to_bgzip, slim_output_t2_sorted))
  
  # tabix bgziped files
  system(paste(path_to_tabix, paste0(slim_output_t1_sorted, ".gz")))
  system(paste(path_to_tabix, paste0(slim_output_t2_sorted, ".gz")))
  
  # merge and get the data
  slim_output_merged <- paste0(slim_output_folder,"slim_output_merged_", sim, ".txt")
  
  bcftools_query <- paste(path_to_bcftools, "merge --force-samples",
                          paste0(slim_output_t1_sorted, ".gz"),
                          paste0(slim_output_t2_sorted, ".gz"),
                          "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%S\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                          ">", slim_output_merged) 
  
  system(bcftools_query)
  
  ## READ SLiM OUTPUT 1 AND CONVERT TO EGGLIB INPUT
  ##-----------------------------------------------
  
  # assembly the header
  header_1       <- c("chrom", "position","alleles", "MID", "S", "DOM", "GO", "MT", "selection")
  sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
  sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
  
  full_header <- c(header_1, sample_names_1, sample_names_2)
  
  # imported the data
  slim_raw_data <- read.table(file = slim_output_merged, header = F, col.names = full_header, check.names = F, na.strings = "./.")
  
  # remove temporary SLiM outputs from the folder
  if (remove_files){
    file.remove(paste0(slim_output_t1))
    file.remove(paste0(slim_output_t1_sorted, ".gz"))
    file.remove(paste0(slim_output_t1_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_t2))
    file.remove(paste0(slim_output_t2_sorted, ".gz"))
    file.remove(paste0(slim_output_t2_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_merged))
  }
  
  # split the data
  slim_snp_geno <- slim_raw_data[, 10:ncol(slim_raw_data)]
  slim_snp_info <- slim_raw_data[, 1:length(header_1)]
  
  # change the genotype annotations
  slim_snp_geno <- as.matrix(slim_snp_geno)
  slim_snp_geno[is.na(slim_snp_geno)]   <- "11"
  slim_snp_geno[slim_snp_geno == "0|0"] <- "11"
  slim_snp_geno[slim_snp_geno =="1|1"]  <- "22"
  slim_snp_geno[slim_snp_geno == "0|1" | slim_snp_geno == "1|0"] <- "12"
  slim_snp_geno <- as.data.frame(slim_snp_geno)
  
  # mark monomophormic mutations (all 11 or 22)
  count_ref_geno    <- apply(slim_snp_geno, 1, function(x){sum(x == 11)})
  count_alt_geno    <- apply(slim_snp_geno, 1, function(x){sum(x == 22)})
  keep_snps     <- count_ref_geno < (SS1 + SS2) & count_alt_geno < (SS1 + SS2) # REMOVE FIXED MUTATIONS
  
  # re-assemble the data
  slim_data <- cbind(slim_snp_info, slim_snp_geno)
  
  # remove monomophormic mutations
  slim_data <- slim_data[keep_snps, ]
  
  # remove duplicated mutations
  # with new version of egglib summstats it supposedly possible to work with redundant positions
  # When it starts working comment this part
  slim_data <- slim_data[!duplicated(slim_data[ ,1:2]), ] # REMOVE DUPLICATED CHROM-POS
  
  # prepare egglib input data
  slim_to_egglib_data <- data.frame(chrom=slim_data$chrom, position=slim_data$position, status=slim_data$MT,
                                    selection=slim_data$selection, alleles=slim_data$alleles)
  
  # assembly final egglib input
  slim_to_egglib_data <- cbind(slim_to_egglib_data, slim_data[, (length(header_1)+1):ncol(slim_data)])
  
  # re-code the status column
  slim_to_egglib_data$status <- ifelse(slim_to_egglib_data$status == 1, "S", "NS")
  
  if (!file_test("-d", egglib_input_folder)){
    dir.create(file.path(egglib_input_folder))
  }
  
  # export egglib input file to the egglib input folder
  egglib_converted_file <- paste0("egglib_input", "_", sim, ".txt");
  write.table(slim_to_egglib_data, file = paste0(egglib_input_folder, egglib_converted_file), quote=FALSE, sep="\t", row.names = FALSE)
  
  # save only the information of the snps
  slim_to_egglib_snps <- data.frame(ID=paste0(slim_data$chrom, ":", slim_data$position), MID=slim_data$MID,
                                    MT=slim_data$MT, S=slim_data$S, DOM=slim_data$DOM, GO=slim_data$GO)
  
  ## HANDLE SLiM2 OUTPUT 2 - GENETIC LOAD
  ##-----------------------------------------------------
  
  # merge and get the data
  slim_output_geneticLoad <- paste0(slim_output_folder,"slim_output_load_", sim, ".txt")
  geneticLoad <- as.numeric(scan(file = slim_output_geneticLoad, quiet = T, na.strings = "NA", what = "character"))
  
  ## READ EGGLIB INPUT AND RUN EGGLIB SUMSTAT
  ##----------------------------------------------------- 
  
  # check if the folder exists
  if (!file_test("-d", egglib_output_folder)){
    dir.create(file.path(egglib_output_folder))
  }
  
  # generate text with slim command  
  egglib_run <- paste(path_to_python,
                      paste0(getwd(), "/", path_to_egglib_summstat),
                      paste0("input-file=", egglib_input_folder, egglib_converted_file),
                      paste0("output-file=", egglib_output_folder, "egglib_output", "_", sim, ".txt"),
                      paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
                      paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"), collapse = ",")),
                      paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "SFS"), collapse = ",")),
                      paste0("wspan=", wss_wspan),
                      paste0("SFS-bins=", sfs_bins),
                      paste0("select=", "all"));
  
  # rum egglib summstat on system
  system(egglib_run)
  
  ## READ EGGLIB OUTPUT AND CREATE THE REFERENCE TABLE
  ##-----------------------------------------------------
  
  # import egglib output
  egglib_summary_stats <- read.csv(file = paste0(egglib_output_folder,"egglib_output_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  
  # remove redundant summary statistics
  egglib_summary_stats <- egglib_summary_stats[, unique(names(egglib_summary_stats))]
  
  # rename the summary statistics
  colnames(egglib_summary_stats) <- gsub(":", "_", names(egglib_summary_stats))
  
  # GLOBAL SUMMARY STATISTICS
  #-------------------------------------------------------
  
  # take only calculated global statistics
  globla_GSS_stats <- egglib_summary_stats[1 , grepl("GSS" , unique(names(egglib_summary_stats)))]
  global_SFS_stats <- egglib_summary_stats[1 , grepl("SFS" , unique(names(egglib_summary_stats)))]
  
  # calculate additional GLOBAL summary statistics
  mean_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){mean(x, na.rm=T)})
  names(mean_LSS_stats) <- c("MEAN_LSS_He", "MEAN_LSS_Dj","MEAN_LSS_WCst", "MEAN_WSS_He", "MEAN_WSS_Dj", 
                             "MEAN_WSS_WCst", "MEAN_WSS_S", "MEAN_WSS_thetaW", "MEAN_WSS_D", "MEAN_WSS_Da", "MEAN_WSS_ZZ")
  
  var_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){var(x, na.rm=T)})
  names(var_LSS_stats) <- c("VAR_LSS_He", "VAR_LSS_Dj","VAR_LSS_WCst","VAR_WSS_He", "VAR_WSS_Dj", 
                               "VAR_WSS_WCst", "VAR_WSS_S", "VAR_WSS_thetaW", "VAR_WSS_D", "VAR_WSS_Da", "VAR_WSS_ZZ")
  
  kurt_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){kurtosis(x, na.rm=T)})
  names(kurt_LSS_stats) <- c("KURT_LSS_He", "KURT_LSS_Dj","KURT_LSS_WCst","KURT_WSS_He", "KURT_WSS_Dj",
                             "KURT_WSS_WCst", "KURT_WSS_S", "KURT_WSS_thetaW", "KURT_WSS_D", "KURT_WSS_Da", "KURT_WSS_ZZ")
  
  skew_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){skewness(x, na.rm=T)})
  names(skew_LSS_stats) <- c("SKEW_LSS_He", "SKEW_LSS_Dj","SKEW_LSS_WCst","SKEW_WSS_He", "SKEW_WSS_Dj",
                             "SKEW_WSS_WCst", "SKEW_WSS_S", "SKEW_WSS_thetaW", "SKEW_WSS_D", "SKEW_WSS_Da", "SKEW_WSS_ZZ")
  
  q05_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
  names(q05_LSS_stats) <- c("Q05_LSS_He", "Q05_LSS_Dj","Q05_LSS_WCst","Q05_WSS_He", "Q05_WSS_Dj",
                            "Q05_WSS_WCst", "Q05_WSS_S", "Q05_WSS_thetaW", "Q05_WSS_D", "Q05_WSS_Da", "Q05_WSS_ZZ")
  
  q95_LSS_stats <- apply(egglib_summary_stats[,-c(1, 13:29)], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
  names(q95_LSS_stats) <- gsub("Q05", "Q95", names(q05_LSS_stats))
  
  additional_global_summary_stats <-cbind(as.data.frame(t(mean_LSS_stats)), as.data.frame(t(var_LSS_stats)), as.data.frame(t(kurt_LSS_stats)), 
                                          as.data.frame(t(skew_LSS_stats)), as.data.frame(t(q05_LSS_stats)), as.data.frame(t(q95_LSS_stats)))
  
  # assemble GLOBAL summary statistics
  global_summary_stats <- cbind(globla_GSS_stats, global_SFS_stats, additional_global_summary_stats)
  
  # LOCUS-SPECIFIC SUMMARY STATISTICS
  #-------------------------------------------------------
  if (random_simulations){
    # sampling snps for the locus-specific reference table
    if (any(slim_to_egglib_snps$MT == 1)){
      m1 <- sample(which(slim_to_egglib_snps$MT == 1), size=1)
    } else {
      m1 <- 0
    }
        
    if (any(slim_to_egglib_snps$MT == 2)){
      m2 <- sample(which(slim_to_egglib_snps$MT == 2), size=1)
    } else {
      m2 <- 0
    }
  
    if (any(slim_to_egglib_snps$MT == 3)){
      m3 <- sample(which(slim_to_egglib_snps$MT == 3), size=1)
    } else {
      m3 <- 0
    }
  
    if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 2)){
      m1m2 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2), size=1)
    } else {
      m1m2 <- 0
    }
  
    if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 3)){
      m1m3 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 3), size=1)
    } else {
      m1m3 <- 0
    }
  
    if (any(slim_to_egglib_snps$MT == 2) | any(slim_to_egglib_snps$MT == 3)){
      m2m3 <- sample(which(slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
    } else {
      m2m3 <- 0
    }

    if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 2) | any(slim_to_egglib_snps$MT == 3)){
      m1m2m3 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
    } else {
      m1m2m3 <- 0
    }
  
    snps_in_reftable <- c(m1,m2,m3,m1m2,m1m3,m2m3,m1m2m3)
  
    # calculate sample minor allele frequency
    slim_to_egglib_genotypes <- slim_to_egglib_data[snps_in_reftable, (grep("alleles", names(slim_to_egglib_data)) + 1):(SS1 + SS2 + 5)]
  
    # MAF sample 1
    S1_genotypes <- slim_to_egglib_genotypes[, grepl(paste0("@pop", 1), names(slim_to_egglib_genotypes))]
    S1n11 <- apply(S1_genotypes==11, 1, sum, na.rm=T)
    S1n12 <- apply(S1_genotypes==12 | S1_genotypes==21, 1, sum, na.rm=T)
    S1n22 <- apply(S1_genotypes==22, 1, sum, na.rm=T)
    MAF1 <- (2*(S1n22) + S1n12)/((2*(S1n11) + S1n12)+(2*(S1n22) + S1n12))
    
    # MAF sample 2
    S2_genotypes <- slim_to_egglib_genotypes[, grepl(paste0("@pop", 2), names(slim_to_egglib_genotypes))]
    S2n11 <- apply(S2_genotypes==11, 1, sum, na.rm=T)
    S2n12 <- apply(S2_genotypes==12 | S2_genotypes==21, 1, sum, na.rm=T)
    S2n22 <- apply(S2_genotypes==22, 1, sum, na.rm=T)
    MAF2 <- (2*(S2n22) + S2n12)/((2*(S2n11) + S2n12)+(2*(S2n22) + S2n12))
  
    # assemble LOCUS-SPECIFIC summary statistics
    locus_lss_info <- slim_to_egglib_snps[snps_in_reftable, ]
    locus_lss_stats <- egglib_summary_stats[snps_in_reftable , grepl("LSS" , unique(names(egglib_summary_stats)))]
    locus_wss_stats <- egglib_summary_stats[snps_in_reftable , grepl("WSS" , unique(names(egglib_summary_stats)))]
  
    locus_summary_stats <- cbind(locus_lss_info, MAF1, MAF2, locus_lss_stats, locus_wss_stats)
    locus_summary_stats <- locus_summary_stats[ !duplicated(locus_summary_stats[ , 3]), ] 
  
  }
  
  ## RAW REFERENCE TABLE
  if (random_simulations){
    raw_reftable  <- suppressWarnings(cbind(sim=sim, theta=theta, mu=mu, rr=rr, Ne0=Ne0, Ne1=Ne1,
                                            gammaMean=gammaM, gammak=gammak,
                                            PrGWSel=PrGWSel, PropMSel=prbe, GeneticLoad=geneticLoad,
                                            locus_summary_stats, global_summary_stats))
  } else {
    raw_reftable  <- suppressWarnings(cbind(sim=sim, theta=theta, mu=mu, rr=rr, Ne0=Ne0, Ne1=Ne1,
                                            gammaMean=gammaM, gammak=gammak,
                                            PrGWSel=PrGWSel, PropMSel=prbe, GeneticLoad=geneticLoad,
                                            global_summary_stats))
  }
  
  # remove all intermediate files 
  if (remove_files){
    if (!sim %% 10 == 0){
      file.remove(paste0(egglib_input_folder, egglib_converted_file))
    }
    file.remove(paste0(egglib_output_folder,"egglib_output_", sim, ".txt"))
    file.remove(paste0(slim_output_geneticLoad))
  }
  
  ## OUTPUT RAW REFERENCE TABLE
  return(raw_reftable)
  
}
