do_sim <- function(sim, nsim, model,
                   mu_rate, mu_min, mu_max, 
                   ne0_min, ne0_max, ne1_min, ne1_max, pge2_min, pge2_max, mpb_min, mpb_max, 
                   gamma_min, gamma_max, rr_rate, rr_min, rr_max, SS1, SS2, ts2,
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
    mu <- rlunif(1, min = mu_min, max = mu_max, base=exp(10))
  }
  
  # THETA
  theta_min = 4*ne0_min*mu
  theta_max = 4*ne0_max*mu
  theta <- rlunif(1, min = theta_min, max = theta_max, base = exp(10))
  
  # EFFECTIVE POPULATION SIZE 0 - Ne0
  ne0 = as.integer(theta/(4*mu))
  
  # EFFECTIVE POPULATION SIZE 1 - Ne1
  ne1 <- as.integer(rlunif(1, min = ne1_min, max = ne1_max, base = exp(10)))
  
  # SAMPLING VALUES FOR SLiM GAMMA DISTRIBUTION OF FITNESS EFFECTS
  gamma_values = sort(runif(2, min = gamma_min, max = gamma_max))
  gamma2_M = gamma_values[1]
  gamma2_k = gamma_values[2]
  
  # PROPORTION OF GENOMIC ELEMENTS G2 
  pge2 = runif(1, min = pge2_min, max = pge2_max)
  
  # PROPORTION OF BENEFICIAL MUTATION IN GENOMIC ELEMENTS G2
  mpb = runif(1, min = mpb_min, max = mbp_max)
  
  # RECOMBINATION RATE
  if (rr_rate == 0){
    rr <- 4.2 * 1e-7
  } else {
    mu <- rlunif(1, min = rr_min, max = rr_max, base=exp(10))
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
  slim_gamma2_M <- paste0("-d gamma2_M=", gamma2_M)
  slim_gamma2_k <- paste0("-d gamma2_k=", gamma2_k)
  slim_pge2     <- paste0("-d pge2=", pge2)
  slim_mpb      <- paste0("-d mpb=", mpb)
  slim_rr       <- paste0("-d rr=", rr)
  slim_SS1      <- paste0("-d SS1", SS1)
  slim_SS2      <- paste0("-d SS1", SS2)
  slim_ts2      <- paste0("-d ts2", ts2)
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
                          "| bcftools query -f 'chr%CHROM\t%POS\t%MT\tY\t%REF,%ALT[\t%GT]\t%MID\t%DOM\t%S\t%GO\n'",
                          ">", slim_output_merged) 
  
  system(bcftools_query)
  
  ## READ SLiM2 OUTPUT 2 - GENETIC LOAD AND SWEEP DATA
  ##----------------------------------------------------------
  
  slim_Pres_output_t1        <- paste0(slim_output_folder,"slim_Pres_output_t1_", sim, ".txt")
  slim_Pres_output_t2        <- paste0(slim_output_folder,"slim_Pres_output_t2_", sim, ".txt")
  
  # imported the data
  slim_Pres_data_t1 <- read.table(file = slim_Pres_output_t1, header = T, check.names = F)
  slim_Pres_data_t1 <- slim_Pres_data_t1[, -c(12)]
  
  slim_Pres_data_t2 <- read.table(file = slim_Pres_output_t2, header = T, check.names = F)
  slim_Pres_data_t2 <- slim_Pres_data_t2[, -c(12)]
  
  Pres_data_merged  <- merge.data.frame(slim_Pres_data_t1, by.x = c(1,2,3,4,5,6,7,8),
                                        slim_Pres_data_t2, by.y = c(1,2,3,4,5,6,7,8), 
                                        all = T, sort = T);
  
  ## READ SLiM OUTPUT 1 AND CONVERT TO EGGLIB INPUT
  ##-----------------------------------------------
  
  # assembly the header
  header_1       <- c("chrom", "position", "status", "selection","alleles")
  sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
  sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
  header_2       <- c("MID", "DOM", "S", "GO")
  
  full_header <- c(header_1, sample_names_1, sample_names_2, header_2)
  
  # imported the data
  slim_data_raw <- read.table(file = slim_output_merged, header = F, col.names = full_header, check.names = F, na.strings = "./.")
  
  # remove temporary SLiM outputs from the folder
  if (remove_files){
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
  slim_data_gen <- slim_data_raw[, 6:(SS1+SS2+5)]
  slim_data_snp <- slim_data_raw[, 1:5]
  slim_data_ext <- slim_data_raw[, c(1:2, (SS1+SS2+6):(SS1+SS2+9))]
  
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
  slim_egglib <- slim_data[kept_snps, ]
  
  # remove duplicated mutations
  # with new version of egglib summstats it supposedly possible to work with redundant positions
  # When it starts working comment this part
  slim_egglib <- slim_egglib[!duplicated(slim_egglib[ ,1:2]), ] # REMOVE DUPLICATED CHROM-POS
                                                                
  # re-code the status column
  slim_egglib$status <- ifelse(slim_egglib$status == 1, "S", "NS")
  
  #### ORIGINAL_starts
  #ts_1 <- read.csv(file = paste0(slim_output_folder,"slim_output_t1_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  #ts_2 <- read.csv(file = paste0(slim_output_folder,"slim_output_t2_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  
  #rwm <- merge.data.frame(ts_1, ts_2, by.x = c(1,2,3,4,5,106,107,108,109), 
  #                        by.y = c(1,2,3,4,5,106,107,108,109), 
  #                        all = T, sort = T); 
  
  #prm <- rwm[, -c(6,7,8,9,110,211)]
  #prm[is.na(prm)] <- "11"
  
  #tempm <- prm[,-c(1,2,3,4,5)]
  #filter1 <- apply(tempm, 1, function(x){sum(x == 11)})
  #rm <- filter1 < dim(tempm)[2] & filter1 !=0
  #m <- prm[rm, ]
  #m <- m[!duplicated(m[ ,1:2]), ]
  #### ORIGINAL_ends
  
  if (!file_test("-d", egglib_input)){
    dir.create(file.path(egglib_input))
  }
  
  # export egglib input file
  conv_file <- paste0("egglib_input", "_", sim, ".txt");
  write.table(slim_egglib, file = paste0(egglib_input, conv_file), quote=FALSE, sep="\t", row.names = FALSE)
  
  # remove the information of monomorphic and duplicated snps
  # with new version of egglib summstats it supposedly possible to work with redundant positions
  # When it starts working comment this part
  slim_data_ext <- slim_data_ext[kept_snps, ]
  slim_data_ext <- slim_data_ext[!duplicated(slim_data_ext[, 1:2]), ] 
                                                                      
  
  #### ORIGINAL_starts
  #pre_extradata <- rwm[, c(1,2,6,7,8,9,110,211)]
  #extradata <- pre_extradata[rm, ]
  #extradata <- extradata[!duplicated(extradata[, 1:2]), ]
  #### ORIGINAL_ends
  
  if (save_extra_data) {
    
    if (!file_test("-d", extra_output)){
      dir.create(file.path(extra_output))
    }
    
    # export extra information
    extra_file <- paste0("extra_output", "_", sim, ".txt");
    write.table(slim_data_ext, file = paste0(extra_output, extra_file), quote=FALSE, sep="\t", row.names = FALSE)
  }
  
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
                      paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"), collapse = ",")),
                      paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "D", "Da", "ZZ"), collapse = ",")), # summstats_1.0.py add "SFS"
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
  
  summstats <- read.csv(file = paste0(egglib_output,"egglib_output_", sim, ".txt"), header = T, sep = "\t", check.names = F)
  summstats <- summstats[, unique(names(summstats))]
  colnames(summstats) <- gsub(":", "_", names(summstats))
  
  # recode ID column to match egglib output format
  slim_data_ext <- cbind(ID=paste0(slim_data_ext$chrom, ":", slim_data_ext$position), slim_data_ext[, -c(1,2)])
  
  # take only the beneficial mutation and calculate the proportion of it
  beneficial_muts <- slim_data_ext[which(slim_data_ext$muSel != 0.0), ]
  pro_beneficial_muts <- length(beneficial_muts$ID)/length(slim_data_ext$ID)
  
  # take only the neutral mutations
  neutral_muts <- slim_data_ext[which(slim_data_ext$muSel == 0.0), ] 
  
  ## GLOBAL REFERENCE TABLE
  glb_summstats <- summstats[1 , grepl( "GSS" , unique(names(summstats)) ) ]
  glb_reftable  <- cbind(sim=sim, seed=sim_seed, theta=theta, mu=mu, Ne0=ne0, Ne1=ne1, prop_beneficial=pro_beneficial_muts, glb_summstats)
  
  ## LOCI-SPECIFIC REFERENCE TABLE
  # I need to test how this part behave with redundant mutations
  
  # Sample one beneficial mutation
  if (is.data.frame(beneficial_muts) && nrow(beneficial_muts)==0){ 
    sampled_beneficial_ID <- NA
  } else {
    sampled_beneficial <- sample_n(beneficial_muts, size = 1)
    sampled_beneficial_ID <- sampled_beneficial$ID
  } 
  
  # Sample one neutral mutation
  if (is.data.frame(neutral_muts) && nrow(neutral_muts)==0){
    sampled_neutral_ID <- NA
  } else {
    sampled_neutral <- sample_n(neutral_muts, size = 1)
    sampled_neutral_ID <- sampled_neutral$ID
  } 
  
  if (is.na(sampled_beneficial_ID)){
    neutral  <- cbind(lc_classif="N", 
                      slim_data_ext[which(slim_data_ext$ID == sampled_neutral_ID), -1],
                      summstats[which(summstats$ID == sampled_neutral_ID), -1])
    loc_reftable <- neutral;
  } else {
    adaptive <- cbind(lc_classif="A",
                      slim_data_ext[which(slim_data_ext$ID == sampled_beneficial_ID), -1],
                      summstats[which(summstats$ID == sampled_beneficial_ID), -1])
    neutral  <- cbind(lc_classif="N", 
                      slim_data_ext[which(slim_data_ext$ID == sampled_neutral_ID), -1],
                      summstats[which(summstats$ID == sampled_neutral_ID), -1])
    loc_reftable <- rbind(adaptive,neutral);
  }
  
  # remove all intermediate files but the reference table 
  if (remove_files){
    file.remove(paste0(egglib_input, conv_file))
    file.remove(paste0(egglib_output,"egglib_output_", sim, ".txt"))
  }
  
  ## combine the GLOBAL and LOCUS-SPECIFIC REFERENCE TABLES
  res <- suppressWarnings(cbind(glb_reftable, loc_reftable))
  return(res)
  
}
