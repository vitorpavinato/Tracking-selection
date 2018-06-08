do_sim <- function(sim, nsim, model,
                   mu_rate, mu_min, mu_max, 
                   ne0_min, ne0_max, ne1_min, ne1_max, pge2_min, pge2_max, mpb_min, mpb_max, 
                   gammaM_gammak, gammaM_min, gammaM_max, gammak_min, gammak_max, 
                   rr_rate, rr_min, rr_max, SS1, SS2, ts2,
                   slim_output_folder, egglib_input, save_extra_data, extra_output,  
                   python_path, egglib_summstat, wss_wspan, sfs_bins,
                   egglib_output, 
                   remove_files){
  
  # write progress of simulation on screen 
  if (sim==1 | sim%%10==0 | sim==nsim){
    cat(paste(sim,"of",nsim))  
  }
  
  #######################
  ##   SAMPLED VALUES  ##
  #######################
  ##-------------------
  
  # SLiM SEED
  sim_seed  <- as.integer(runif(1, 100000000, 900000000));
  
  # MUTATION RATE
  if (mu_rate == 0){
    mu <- 1e-7
  } else {
    #mu <- rlunif(1, min = mu_min, max = mu_max, base = exp(10)) 
    mu <- 10^runif(1, min = log10(mu_min), max = log10(mu_max)) # with works if you're unable to install KScorrect package
  }
  
  # THETA
  theta_min = 4*ne0_min*mu
  theta_max = 4*ne0_max*mu
  #theta <- rlunif(1, min = theta_min, max = theta_max, base = exp(10))
  theta <- 10^runif(1, min = log10(theta_min), max = log10(theta_max))
  
  # EFFECTIVE POPULATION SIZE 0 - Ne0
  ne0 = as.integer(theta/(4*mu))
  
  # EFFECTIVE POPULATION SIZE 1 - Ne1
  #ne1 <- as.integer(rlunif(1, min = ne1_min, max = ne1_max, base = exp(10)))
  ne1 <- as.integer(10^runif(1, min = log10(ne1_min), max = log10(ne1_max)))
  
  # SAMPLING VALUES FOR SLiM GAMMA DISTRIBUTION OF FITNESS EFFECTS
  # GAMMA MEAN
  #gammaM = round(rlunif(1, min = gammaM_min, max = gammaM_max, base = exp(10)), digits = 4)
  gammaM = round(10^runif(1, min = log10(gammaM_min), max = log10(gammaM_max)), digits = 4)
  
  # GAMMA SHAPE 
  if (gammaM_gammak){
    gammak = gammaM
  } else {
    #gammak = round(rlunif(1, min = gammak_min, max = gammak_max, base = exp(10)), digits = 4)
    gammak = round(10^runif(1, min = log10(gammak_min), max = log10(gammak_max)), digits = 4)
  }
  
  # PROPORTION OF GENOMIC ELEMENTS G2
  #pge2 = round(rlunif(1, min = pge2_min, max = pge2_max, base = exp(10)), digits = 6)
  pge2 = round(10^runif(1, min = log10(pge2_min), max = log10(pge2_max)), digits = 6)
  
  # PROPORTION OF BENEFICIAL MUTATION IN GENOMIC ELEMENTS G2
  #mpb = round(rlunif(1, min = mpb_min, max = mpb_max, base = exp(10)), digits = 6)
  mpb = round(10^runif(1, min = log10(mpb_min), max = log10(mpb_max)), digits = 6)
  
  # RECOMBINATION RATE
  if (rr_rate == 0){
    rr <- 4.2 * 1e-7
  } else {
    #rr <- 4.2 * (rlunif(1, min = rr_min, max = rr_max, base=exp(10)))
    rr <- 4.2 * (10^runif(1, min = log10(rr_min), max = log10(rr_max)))
  }
  
  ## RUM SLiM2
  ##-------------------
  
  # check if the folder exists
  if (!file_test("-d", slim_output_folder)){
    dir.create(file.path(slim_output_folder))
  }
  
  # generate text with slim command
  slim_output_fullpath <- paste0("'",getwd(), "/", slim_output_folder, "'")
  
  slim_simID    <- paste0("-d simID=", sim)
  slim_theta    <- paste0("-d theta=", theta)
  slim_mu       <- paste0("-d mu=", mu)
  slim_Ne1      <- paste0("-d Ne1=", ne1)
  slim_gamma2_M <- paste0("-d gamma2_M=", gammaM)
  slim_gamma2_k <- paste0("-d gamma2_k=", gammak)
  slim_pge2     <- paste0("-d pge2=", sprintf("%.6f", pge2))
  slim_mpb      <- paste0("-d mpb=", sprintf("%.6f", mpb))
  slim_rr       <- paste0("-d rr=", rr)
  slim_SS1      <- paste0("-d SS1=", SS1)
  slim_SS2      <- paste0("-d SS2=", SS2)
  slim_ts2      <- paste0("-d ts2=", ts2)
  slim_output <- paste0("-d ", '"',"outputpath=", slim_output_fullpath, '"')
  slim_model <- paste0(model, ".slim")
  
  slim_run <- paste( "/usr/local/bin/slim",
                     "-s", sim_seed,                # seed  = simulation seed number
                     slim_simID,                    # simID = simulation id number
                     slim_theta,                    # theta
                     slim_mu,                       #    mu = mutation rate
                     slim_Ne1,                      #   ne1 = effective population size 1  
                     slim_gamma2_M,
                     slim_gamma2_k,
                     slim_pge2,
                     slim_mpb,
                     slim_rr,
                     slim_SS1,
                     slim_SS2,
                     slim_ts2,
                     slim_output,
                     slim_model)                    # model
  
  # rum slim on system
  system(slim_run)
  
  ## HANDLING SLiM2 OUTPUT 1 & 2 - SWEEP, SELECTION & LOAD DATA
  ##-----------------------------------------------------------
  
  ######################
  # "SEGREGATION" DATA #
  ######################
  
  # LOAD T=1 DATA
  slim_segr_output_t1        <- paste0(slim_output_folder,"slim_segr_output_t1_", sim, ".txt")
  slim_segr_data_t1 <- read.table(file = slim_segr_output_t1, header = T, check.names = F)
  
  # LOAD T=2 DATA
  slim_segr_output_t2        <- paste0(slim_output_folder,"slim_segr_output_t2_", sim, ".txt")
  slim_segr_data_t2 <- read.table(file = slim_segr_output_t2, header = T, check.names = F)
  
  # Merge segregation data T=1 & T=2
  segr_merged_data  <- merge.data.frame(slim_segr_data_t1, by.x = c(1,2,3,4,5,6,7,8,9),
                                        slim_segr_data_t2, by.y = c(1,2,3,4,5,6,7,8,9), 
                                        all = T, sort = T);
  
  # Check if there is at leas one mutation - else NULL
  if (sum(is.na(segr_merged_data)) == (length(segr_merged_data)-2)){
    segr_merged_data = NULL
  } else {
    segr_merged_data <- segr_merged_data[complete.cases(segr_merged_data[,-c(10:11)]), ]
  }
  
  ######################
  # "SELECTION" DATA #
  ######################
  
  # LOAD T=2 DATA
  slim_sel_output_t2        <- paste0(slim_output_folder,"slim_sel_output_t2_", sim, ".txt")
  slim_sel_data_t2 <- read.table(file = slim_sel_output_t2, header = T, check.names = F)
  
  # Save GENETIC LOAD DATA
  genetic_load = slim_sel_data_t2[1,9]
  
  # Remove genetic load data  
  slim_sel_data_t2 <- slim_sel_data_t2[, -c(9)]
  
  # Check if there is at leas one selected site - else NULL
  if (sum(is.na(slim_sel_data_t2)) == (length(slim_sel_data_t2)-2)){
    slim_sel_data_t2 = NULL
  }
  
  # Take the mutation closest its midpoint
  if (!is.null(segr_merged_data)){
    midpoint_segr <- segr_merged_data[which.min(segr_merged_data$d2midpoint), ]
  } else {
    midpoint_segr <- NULL
  }

  # Calculate the ALPHA and the PROBABILITY OF SCAPING A SELECTIVE SWEEP (Pr_es)
  if (!is.null(midpoint_segr)){
    if (!is.null(slim_sel_data_t2)){
      alpha_sel = (rr*log(2*ne1))/slim_sel_data_t2$S
      dist_sel  = sqrt((midpoint_segr$position - slim_sel_data_t2$position)^2)
      PPr_es = prod(1 - exp(-alpha_sel*dist_sel)) 
      midpoint_segr$PPr_es <- PPr_es
     } else {
        PPr_es = 1
        midpoint_segr$PPr_es <- PPr_es
     }
  }  
  
  # Format the the output to produce the reference table
  if (!is.null(midpoint_segr)){
    midpoint_segr <- midpoint_segr[,-c(2,4,5)]
    midpoint_segr <- cbind(ID=paste0(midpoint_segr$chrom, ":", midpoint_segr$position), midpoint_segr[, -c(1,2)])
  } 
  
  ## HANDLE SLiM2 OUTPUT 3 - GENETIC DATA
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
  system(paste("bgzip", slim_output_t1_sorted))
  system(paste("bgzip", slim_output_t2_sorted))
  
  # tabix bgziped files
  system(paste("tabix", paste0(slim_output_t1_sorted, ".gz")))
  system(paste("tabix", paste0(slim_output_t2_sorted, ".gz")))
  
  # merge and get the data
  slim_output_merged <- paste0(slim_output_folder,"slim_output_merged_", sim, ".txt")
  
  bcftools_query <- paste("bcftools merge --force-samples",
                          paste0(slim_output_t1_sorted, ".gz"),
                          paste0(slim_output_t2_sorted, ".gz"),
                          "| bcftools query -f 'chr%CHROM\t%POS\t%MID\t%MT\tY\t%REF,%ALT[\t%GT]\n'",
                          ">", slim_output_merged) 
  
  system(bcftools_query)
  
  
  
  ## READ SLiM OUTPUT 1 AND CONVERT TO EGGLIB INPUT
  ##-----------------------------------------------
  
  # assembly the header
  header_1       <- c("chrom", "position","MID", "status", "selection","alleles")
  sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
  sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
  
  full_header <- c(header_1, sample_names_1, sample_names_2)
  
  # imported the data
  slim_data_raw <- read.table(file = slim_output_merged, header = F, col.names = full_header, check.names = F, na.strings = "./.")
  
  # remove temporary SLiM outputs from the folder
  if (remove_files){
    file.remove(paste0(slim_segr_output_t1))
    file.remove(paste0(slim_segr_output_t2))
    file.remove(paste0(slim_sel_output_t2))
    file.remove(paste0(slim_output_t1))
    file.remove(paste0(slim_output_t1_sorted))
    file.remove(paste0(slim_output_t1_sorted, ".gz"))
    file.remove(paste0(slim_output_t1_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_t2))
    file.remove(paste0(slim_output_t2_sorted))
    file.remove(paste0(slim_output_t2_sorted, ".gz"))
    file.remove(paste0(slim_output_t2_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_merged))
  }
  
  # split the data
  slim_data_gen <- slim_data_raw[, 7:(SS1+SS2+6)]
  slim_data_snp <- slim_data_raw[, 1:6]
  
  # change the genotypes annotation
  slim_data_gen <- as.matrix(slim_data_gen)
  slim_data_gen[is.na(slim_data_gen)]   <- "11"
  slim_data_gen[slim_data_gen == "0|0"] <- "11"
  slim_data_gen[slim_data_gen =="1|1"]  <- "22"
  slim_data_gen[slim_data_gen == "0|1" | slim_data_gen=="1|0"] <- "12"
  
  slim_data_gen        <- as.data.frame(slim_data_gen)
  
  # remove monomophormic mutations (all 11 or 22)
  geno_count_ref    <- apply(slim_data_gen, 1, function(x){sum(x == 11)})
  geno_count_alt    <- apply(slim_data_gen, 1, function(x){sum(x == 22)})
  
  kept_snps     <- geno_count_ref < (SS1 + SS2) & geno_count_alt < (SS1 + SS2) # REMOVE FIXED MUTATIONS
  
  slim_data <- cbind(slim_data_snp, slim_data_gen)
  slim_data <- slim_data[, -c(3)]
  slim_egglib <- slim_data[kept_snps, ]
  
  # remove duplicated mutations
  # with new version of egglib summstats it supposedly possible to work with redundant positions
  # When it starts working comment this part
  slim_egglib <- slim_egglib[!duplicated(slim_egglib[ ,1:2]), ] # REMOVE DUPLICATED CHROM-POS
                                                                
  # re-code the status column
  slim_egglib$status <- ifelse(slim_egglib$status == 1, "S", "NS")
  
  if (!file_test("-d", egglib_input)){
    dir.create(file.path(egglib_input))
  }
  
  # export egglib input file
  conv_file <- paste0("egglib_input", "_", sim, ".txt");
  write.table(slim_egglib, file = paste0(egglib_input, conv_file), quote=FALSE, sep="\t", row.names = FALSE)
  
  ## READ EGGLIB INPUT AND RUN EGGLIB SUMSTAT
  ##----------------------------------------------------- 
  
  # check if the folder exists
  if (!file_test("-d", egglib_output)){
    dir.create(file.path(egglib_output))
  }
  
  # generate text with slim command  
  egglib_run <- paste(python_path,
                      paste0(getwd(), "/", egglib_summstat),
                      paste0("input-file=", egglib_input, conv_file),
                      paste0("output-file=", egglib_output, "egglib_output", "_", sim, ".txt"),
                      paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
                      paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da"), collapse = ",")),
                      paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da"), collapse = ",")), # summstats_1.0.py add "SFS"
                      paste0("wspan=", wss_wspan),
                      #paste0("SFS-bins=", sfs_bins),
                      paste0("select=", "all"));
  
  # rum egglib summstat on system
  if(.Platform$OS.type == "unix") {
    slim_run <- paste( ".", egglib_run, sep="")
    system(egglib_run)
  }
  
  ## READ EGGLIB OUTPUT AND CREATE THE REFERENCE TABLE
  ##-----------------------------------------------------
  
  # Import egglib summary statistics output
  summ_stats <- read.csv(file = paste0(egglib_output,"egglib_output_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  summ_stats <- summ_stats[, unique(names(summ_stats))]
  colnames(summ_stats) <- gsub(":", "_", names(summ_stats))
  
  # Match midpoint mutation with locus-specific summary statistics 
  if (!is.null(midpoint_segr)){
    ref_locus <- summ_stats[summ_stats$ID==paste(midpoint_segr$ID), -1]
    if (dim(ref_locus)[1] == 1){
    ref_locus_stats <- cbind(midpoint_segr, ref_locus)
    } else {
      ref_locus_stats <- as.data.frame(matrix(NA, nrow = 1, ncol = 25))
      names(ref_locus_stats) <- c("ID","MID","S","DOM","GO","MF1","MF2","PPr_es",
                               "LSS_He","LSS_Dj","LSS_WCst","WSS_He","WSS_Dj","WSS_WCst","WSS_S","WSS_thetaW","WSS_D","WSS_Da",
                               "GSS_He","GSS_Dj","GSS_WCst","GSS_S","GSS_thetaW","GSS_D","GSS_Da")
    }
  } else {
    ref_locus_stats <- as.data.frame(matrix(NA, nrow = 1, ncol = 25))
    names(ref_locus_stats) <- c("ID","MID","S","DOM","GO","MF1","MF2","PPr_es",
                          "LSS_He","LSS_Dj","LSS_WCst","WSS_He","WSS_Dj","WSS_WCst","WSS_S","WSS_thetaW","WSS_D","WSS_Da",
                          "GSS_He","GSS_Dj","GSS_WCst","GSS_S","GSS_thetaW","GSS_D","GSS_Da")
  }
  
  # Take global summary statistics
  global_stats <- summ_stats[1 , grepl( "GSS" , unique(names(summ_stats)) ) ]
  
  # Calculate additional summary statistics
  mean_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){mean(x, na.rm=T)})
  names(mean_global_stats) <- c("MEAN_LSS_He", "MEAN_LSS_Dj","MEAN_LSS_WCst",
                                "MEAN_WSS_He", "MEAN_WSS_Dj", "MEAN_WSS_WCst", "MEAN_WSS_S", "MEAN_WSS_thetaW", "MEAN_WSS_D", "MEAN_WSS_Da")
  
  var_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){var(x, na.rm=T)})
  names(var_global_stats) <- c("VAR_LSS_He", "VAR_LSS_Dj","VAR_LSS_WCst",
                               "VAR_WSS_He", "VAR_WSS_Dj", "VAR_WSS_WCst", "VAR_WSS_S", "VAR_WSS_thetaW", "VAR_WSS_D", "VAR_WSS_Da")
  
  kurt_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){kurtosis(x, na.rm=T)})
  names(kurt_global_stats) <- c("KURT_LSS_He", "KURT_LSS_Dj","KURT_LSS_WCst",
                                "KURT_WSS_He", "KURT_WSS_Dj","KURT_WSS_WCst", "KURT_WSS_S", "KURT_WSS_thetaW", "KURT_WSS_D", "KURT_WSS_Da")
  
  skew_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){skewness(x, na.rm=T)})
  names(skew_global_stats) <- c("SKEW_LSS_He", "SKEW_LSS_Dj","SKEW_LSS_WCst",
                                "SKEW_WSS_He", "SKEW_WSS_Dj","SKEW_WSS_WCst", "SKEW_WSS_S", "SKEW_WSS_thetaW", "SKEW_WSS_D", "SKEW_WSS_Da")
  
  q05_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
  names(q05_global_stats) <- c("Q05_LSS_He", "Q05_LSS_Dj","Q05_LSS_WCst",
                               "Q05_WSS_He", "Q05_WSS_Dj","Q05_WSS_WCst", "Q05_WSS_S", "Q05_WSS_thetaW", "Q05_WSS_D", "Q05_WSS_Da")
  
  q95_global_stats <- apply(summ_stats[,-c(1, 12:18)], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
  names(q95_global_stats) <- gsub("Q05", "Q95", names(q05_global_stats))
  
  add_global_stats <-cbind(as.data.frame(t(mean_global_stats)), as.data.frame(t(var_global_stats)), as.data.frame(t(kurt_global_stats)), 
                           as.data.frame(t(skew_global_stats)), as.data.frame(t(q05_global_stats)), as.data.frame(t(q95_global_stats)))
  
  ## REFERENCE TABLE
  raw_reftable  <- cbind(sim=sim, seed=sim_seed, theta=theta, mu=mu, Ne0=ne0, Ne1=ne1, gammaMean=gammaM, gammak=gammak,
                         PropGSel=sprintf("%.4f", pge2), PropMSel=sprintf("%.4f", mpb), rr=sprintf("%1.e", rr),
                         genetic_load, global_stats, add_global_stats,
                         ref_locus_stats)
  
  
  # remove all intermediate files 
  if (remove_files){
    file.remove(paste0(egglib_input, conv_file))
    file.remove(paste0(egglib_output,"egglib_output_", sim, ".txt"))
  }
  
  ## OUTPUT RAW REFERENCE TABLE
  return(raw_reftable)
  
}
