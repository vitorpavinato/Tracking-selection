## MODEL SELECTION
##------------------------------------------------------------------------------------------------

model_title <- switch(model_type, "DN", "BS", "SV")

## FIXED VALUES DEFINING GENOMIC STRUCTURE
##------------------------------------------------------------------------------------------------

# DEFINE CHROMOSOME SIZE AND RECOMBINATION LIMITS
if (chrN == 1){
  chrS = as.integer(0.25 * genomeS) 
  
  rr_limits = c((genomeS-1), genomeS, ((genomeS+chrS)-1))
} else {
  chrS = as.integer(genomeS/chrN)
  
  chrs_uppers = NULL
  chrs_lowers = NULL
  rr_limits = NULL
  for(i in seq(from = 1, to = (chrN-1))){
    uppers = ((i*chrS)-1)
    chrs_uppers = c(chrs_uppers, uppers)
    
    lowers = (uppers+1)
    chrs_lowers = c(chrs_lowers, lowers)
    
    limits = c(uppers, lowers)
    rr_limits = c(rr_limits, limits)
  }
}	

## IF IT IS A RADseq DATA
##----------------------------------------------------------------------------------------------
if (data_type == 2){
  
  # Set the total number of the RADseq reads
  radseq_readN = as.integer((genomeS/radseq_readL))
  
  # Set the STARTS and the ENDS for the RADseq reads
  rad_starts = NULL
  for(i in seq(from = 0, to = radseq_readN-1)){
    starts    = i*radseq_readL
    rad_starts = c(rad_starts, starts)
  }
  
  # Sample the starts/ends pairs of RADseq reads
  radseq_sampled <- sort(sample(rad_starts, as.integer(radseq_cov*radseq_readN), replace = FALSE))
  
  # Generate the vector of positionS for the sampled RADseq reads
  rad_interval = NULL
  for(i in seq_along(radseq_sampled)){
    interval   = radseq_sampled[i]:(radseq_sampled[i]+(radseq_readL-1))
    rad_interval = c(rad_interval, interval)
  }
}

## ADDITIONAL FUNCTIONS
##---------------------------------------------------------------------------------------------
# function to plug into the apply function to re-code the chromosome name
chromtagging <- function(x){
  
  for (i in seq(from = 1, to = (length(chrs_lowers)-1))){
    
    if (x < chrs_lowers[1]){ chrom_idd = paste0("chr", 1)}
    
    else if (x > chrs_lowers[length(chrs_lowers)]){ chrom_idd = paste0("chr", (length(chrs_lowers)+1))}
    
    else if (x > chrs_lowers[i] & x < chrs_lowers[i+1]){ chrom_idd = paste0("chr", (i+1))}
    
  }
  return(chrom_idd)
}

## SIMULATIONS
##---------------------------------------------------------------------------------------------
do_sim <- function(sim, nsim, 
                   path_to_slim_model, slim_model_prefix, model_type, model_title,
                   path_to_slim, slim_output_folder,
                   path_to_bgzip, path_to_tabix, path_to_bcftools,
                   egglib_input_folder, egglib_output_folder,
                   path_to_python, path_to_egglib_summstat, 
                   SS1, SS2, tau, genomeS, fragS, chrN, chrS, chrTAG, chromtagging,  
                   rr_limits, data_type, radseq_readL, radseq_readN, radseq_cov, missing_data, haplotype,
                   mu_rate, mu_random, mu_min, mu_max, 
                   neq_value, neq_random, neq_min, neq_max,
                   n_value, n_random, n_min, n_max,  
                   gammaM_value, gammak_value, gammaM_gammak, 
                   gammaM_random, gammaM_min, gammaM_max, 
                   gammak_random, gammak_min, gammak_max, 
                   PrGWSel_value, PrGWSel_random, PrGWSel_min, PrGWSel_max, 
                   prbe_value, prbe_random, prbe_min, prbe_max, 
                   domN, domN_random, domN_min, domN_max, 
                   domB, domB_random, domB_min, domB_max, 
                   rr_rate, rr_random, rr_min, rr_max,
                   selfing_rate, selfing_random, selfing_min, selfing_max,
                   wss_wspan_run, sfs_bins_run,
                   add_WSSwspan_SFSbins_1, add_wss_wspan_1, add_sfs_bins_1,
                   add_WSSwspan_SFSbins_2, add_wss_wspan_2, add_sfs_bins_2,
                   remove_files
                   ){
  
  # write progress of simulation on screen 
  if (sim==1 | sim%%10==0 | sim==nsim){
    cat(paste(sim,"of",nsim))  
  }
  
  ## SAMPLED VALUES FOR SIMULATION's PARAMETERS
  ##----------------------------------------------------------------------------------
  
  # SIMULATION SEED
  sim_seed  <- as.integer(runif(n = 1, min = 100000000, max = 900000000));
  
  # MUTATION RATE
  if (mu_random){
    mu <- 10^runif(n = 1, min = log10(mu_min), max = log10(mu_max))
  } else {
    mu <- mu_rate
  }
  
  # POPULATION SIZE Neq
  if (neq_random){
    Neq <- as.integer(10^runif(n = 1, min = log10(neq_min), max = log10(neq_max)))  
  } else {
    Neq <- neq_value
  }
  
  # POPULATION SIZE N
  if (n_random){
    N <- as.integer(10^runif(n = 1, min = log10(n_min), max = log10(n_max)))  
  } else {
    N <- n_value
  }
  
  # GENOME-WIDE DFE FOR BENEFICIAL MUTATIONS
  # gamma mean
  if (gammaM_random){
    gammaM <- 10^runif(n = 1, min = log10(gammaM_min), max = log10(gammaM_max))
  } else {
    gammaM <- gammaM_value
  }
  
  if (model_type == 2){
    gammaM = (-1)*gammaM
  }
    
  # gamma shape 
  if (gammak_random){
    if (gammaM_gammak){
      if(model_type != 2){
       gammak = gammaM
      }else {
        gammak = (-1)*gammaM
      }
    } else {
      gammak <- 10^runif(n = 1, min = log10(gammak_min), max = log10(gammak_max))
    }
  } else {
    gammak <- gammak_value
  }
  
  # PROPORTION OF THE GENOME THAT CONTAINS BENEFICIAL MUTATIONS - G2 ELEMENTS
  if (PrGWSel_random){
    PrGWSel <- runif(n = 1, min = PrGWSel_min, max = PrGWSel_max)
  } else {
    PrGWSel <-  PrGWSel_value
  }
  
  # PROPORTION OF GENOME-WIDE BENEFICIAL MUTATION IN G2 ELEMENTS
  if (prbe_random){
    prbe <- runif(1, min = prbe_min, max = prbe_max)
  } else {
    prbe <- prbe_value
  }
  
  # DOMINANCE FOR GENOME-WIDE MUTATIONS
  # GENOME Neutral mutations (m1 and m2) and EXTRA CHROMOSOME Neutral mutation (m4)
  if (domN_random){
    dm1 = dm2 <- runif(n = 1, min = domN_min, max = domN_max)
  } else {
    dm1 = dm2 = dm4 <- domN
  }
  
  # GENOME Beneficial mutations (m3)
  if (domB_random){
    dm3 <- runif(n = 1, min = domB_min, max = domB_max)
  } else {
    dm3 <- domB
  }
  
  # GENOME-WIDE RECOMBINATION RATE
  if (rr_random){
    rr <- 4.2 * (10^runif(1, min = log10(rr_min), max = log10(rr_max)))
  } else {
    rr <- rr_rate
  }
  
  # the distribution of rr across chromosome limits
  if (chrN == 1){
    rr_rates = c(rr, 0.5, rr)
  } else {
    rr_rates = rep(c(rr, 0.5), as.integer(length(rr_limits)/2))

    # updated rr_limits and rr_rates
    rr_limits = c(rr_limits, (genomeS-1), genomeS, ((genomeS+chrS)-1))
    rr_rates = c(rr_rates, rr, 0.5, rr)
  }	
  
  # SELFING RATE
  if (selfing_random){
    selfing <- (10^runif(1, min = log10(selfing_min), max = log10(selfing_max)))
  } else {
    selfing <- selfing_rate
  }
  
  # SELECTION ON STANDING VARIATION
  if (model_type == 3){
    tc = as.integer(runif(1, min = 0, max = tau))
  } else {
    tc = 0
  }
  
  ## RANDOM VALUES DEFINING GENOMIC ELEMENTS  
  ##-------------------------------------------------------------------------------
  
  # GENOME'S GENOME ELEMENT TYPE
  
  # Set the number of the GenomicElementType
  genomicElementN = as.integer((genomeS/fragS))
  
  # Set the STARTS and the ENDS for the initializeGenomicElement
  e_starts = NULL
  e_ends   = NULL
  for(i in seq(from = 0, to = genomicElementN-1)){
    starts    = i*fragS
    e_starts = c(e_starts, starts)
  
    ends = ((i+1)*fragS)-1
    e_ends = c(e_ends, ends)
  }
  
  # Sample the starts/ends pairs to assign to GenomicElementType G2
  indexes = seq(from = 0, to = (genomicElementN-1))
  
  g2_idx = sort(sample(indexes, as.integer(PrGWSel*genomicElementN), replace = FALSE))
  
  # The difference are the GenomicElementType G1
  g1_idx = setdiff(indexes, g2_idx)
  
  ## RUNNING SLiM
  ##-------------------------------------------------------------------------------
  
  # check if the folder exists
  if (!file_test("-d", slim_output_folder)){
    dir.create(file.path(slim_output_folder))
  }
  
  # generate text with slim command
  slim_model_p1 <- paste0(path_to_slim_model,slim_model_prefix,"_", model_title,"_COAL",".slim")
  slim_model_p2 <- paste0(path_to_slim_model,slim_model_prefix,"_", model_title,"_posCOAL",".slim")
  
  slim_output_fullpath <- paste0("'",getwd(), "/", slim_output_folder, "'")
  
  slim_output    <- paste0("-d ", '"',"outputpath=", slim_output_fullpath, '"')
  slim_simID     <- paste0("-d simID=", sim)
  slim_mu        <- paste0("-d mu=", mu)
  slim_Neq       <- paste0("-d Neq=", Neq)
  slim_N         <- paste0("-d N=", N)
  slim_gammaM    <- paste0("-d gammaM=", gammaM)          # p1 and p2 - DN and BS; p2 - SV
  slim_gammak    <- paste0("-d gammak=", gammak)          # p1 and p2 - DN and BS; p2 - SV
  slim_prbe      <- paste0("-d prbe=", prbe)
  slim_dm1       <- paste0("-d dm1=", dm1)
  slim_dm2       <- paste0("-d dm2=", dm2)
  slim_dm3       <- paste0("-d dm3=", dm3)
  slim_dm4       <- paste0("-d dm4=", dm4)
  slim_rr        <- paste0("-d rr=", rr)
  slim_selfing   <- paste0("-d selfing=", selfing)        # p1 only
  slim_tc        <- paste0("-d tc=", tc)                  # p2 - SV
  slim_genomeS   <- paste0("-d genomeS=", as.integer(genomeS))
  slim_chrS      <- paste0("-d chrS=", chrS)
  slim_rr_rates  <- paste0("-d rr_rates=", "'c(", paste(rr_rates, collapse = ","), ")'")
  slim_rr_limits <- paste0("-d rr_limits=", "'c(", paste(rr_limits, collapse = ","), ")'")
  slim_ge_starts <- paste0("-d e_starts=", "'c(", paste(e_starts, collapse = ","), ")'")
  slim_ge_ends   <- paste0("-d e_ends=", "'c(", paste(e_ends, collapse = ","), ")'")
  slim_g2_idx    <- paste0("-d g2_idx=", "'c(", paste(g2_idx, collapse = ","), ")'")
  slim_g1_idx    <- paste0("-d g1_idx=", "'c(", paste(g1_idx, collapse = ","), ")'")
  slim_SS1       <- paste0("-d SS1=", SS1)
  slim_SS2       <- paste0("-d SS2=", SS2)
  slim_tau       <- paste0("-d tau=", tau)
  
  # slim coalescent check run parameters
  slim_run_p1 <- paste(path_to_slim,
                       "-s", sim_seed,                # seed  = simulation seed number
                       slim_simID,                    # simID = simulation id number
                       slim_mu,                       
                       slim_Neq,                    
                       slim_gammaM,                   # p1 and p2 - DN and BS; p2 - SV
                       slim_gammak,                   # p1 and p2 - DN and BS; p2 - SV
                       slim_prbe,
                       slim_dm1,
                       slim_dm2,
                       slim_dm3,
                       slim_dm4,
                       slim_rr,
                       slim_selfing,                  # p1 only 
                       slim_genomeS,
                       slim_chrS,
                       slim_rr_rates,
                       slim_rr_limits,
                       slim_ge_starts,
                       slim_ge_ends,
                       slim_g2_idx,
                       slim_g1_idx,
                       slim_output,
                       slim_model_p1)                 # model part 1 - coalescent check
  
  # run slim coalescent check on system
  system(slim_run_p1)
  
  slim_run_p2 <- paste(path_to_slim,
                      "-s", sim_seed,                # seed  = simulation seed number
                      slim_simID,                    # simID = simulation id number
                      slim_mu,                       
                      slim_N,
                      slim_gammaM,                   # p1 and p2 - DN and BS; p2 - SV
                      slim_gammak,                   # p1 and p2 - DN and BS; p2 - SV
                      slim_prbe,
                      slim_dm1,
                      slim_dm2,
                      slim_dm3,
                      slim_dm4,
                      slim_rr,
                      slim_tc,                       # p2 - SV
                      slim_genomeS,
                      slim_chrS,
                      slim_rr_rates,
                      slim_rr_limits,
                      slim_ge_starts,
                      slim_ge_ends,
                      slim_g2_idx,
                      slim_g1_idx,
                      slim_SS1,
                      slim_SS2,
                      slim_tau,
                      slim_output,
                      slim_model_p2)                 # model part 2 - sampling
  
  # run slim on system
  system(slim_run_p2)
  
  ## HANDLING SLiM2 OUTPUT t1:t2 - GENETIC LOAD
  ##------------------------------------------------------------------------------
  
  # load the data
  slim_output_geneticLoad <- paste0(slim_output_folder,"slim_output_load_", sim, ".txt")
  
  if (file.exists(slim_output_geneticLoad)){
    info_GenLoad_file = file.info(slim_output_geneticLoad)
    
    if (info_GenLoad_file$size != 0){
      geneticLoad <- read.csv(file = slim_output_geneticLoad, sep = "\t", header = F, col.names = c("gen", "time", "genetiload"), na.strings = "NA")
      
      if (nrow(geneticLoad) != 0){
        
        # check if there is any -Inf/Inf and remove the rows containing it
        rows_to_keep <- !is.infinite(geneticLoad$genetiload)
        geneticLoad <- geneticLoad[rows_to_keep, ]
        
        # get the mean value
        averageGenLoad <- mean(geneticLoad[, "genetiload"], na.rm = TRUE) # all NA's averageGenLoad = NaN;
        
        # get the last value
        lastGenLoad <- geneticLoad[nrow(geneticLoad), "genetiload"]
        
      } else {
        averageGenLoad <- as.numeric(NA)
        lastGenLoad <- as.numeric(NA)
      }
      
    } else {
      averageGenLoad <- as.numeric(NA)
      lastGenLoad <- as.numeric(NA)
    }
  
  } else {
    averageGenLoad <- as.numeric(NA)
    lastGenLoad <- as.numeric(NA)
  }
  
  # remove temporary genetic load output from slim_output folder
  if (remove_files){
    if (file.exists(slim_output_geneticLoad)){
      file.remove(slim_output_geneticLoad)
    }
  }
  
  ## HANDLING SLiM2 OUTPUT t1:t2 - PEDIGREE NE
  ##------------------------------------------------------------------------------
  
  # load the data
  slim_output_pedigreene <- paste0(slim_output_folder,"slim_output_pedigreeNe_", sim, ".txt")
  
  if (file.exists(slim_output_pedigreene)){
    info_pedigreene_file = file.info(slim_output_pedigreene)
    
    if (info_pedigreene_file$size != 0){
      pedigreene <- read.csv(file = slim_output_pedigreene, sep = "\t", header = F, col.names = c("gen", "time", "ne"), na.strings = "NA")
      
      if (nrow(pedigreene) != 0){
        
        # check if there is any -Inf/Inf and susbstitute it for NA
        pedigreene[which(is.infinite(pedigreene$ne)), "ne"] <- NA
        
        # take the Pedigree Ne for each interval
        pedigreeNetimes <- pedigreene$ne
        pedigreeNetimes <- as.data.frame(t(pedigreene$ne))
        names(pedigreeNetimes) <- paste0("PedigreeNe",seq(from=0,to=(tau)))
        
        # get the Pedigree Ne for the whole period - harmonic mean
        pedigreeNetotal <- 1/mean(1/pedigreene$ne, na.rm = TRUE) # all NA's pedigreeNetotal = NA;
        
      } else {
        pedigreeNetimes <- as.data.frame(t(rep(NA, tau+1)))
        names(pedigreeNetimes) <- paste0("PedigreeNe",seq(from=0,to=(tau)))
        pedigreeNetotal <- as.numeric(NA)
      }
      
    } else {
      pedigreeNetimes <- as.data.frame(t(rep(NA, tau+1)))
      names(pedigreeNetimes) <- paste0("PedigreeNe",seq(from=0,to=(tau)))
      pedigreeNetotal <- as.numeric(NA)
    }
    
  } else {
    pedigreeNetimes <- as.data.frame(t(rep(NA, tau+1)))
    names(pedigreeNetimes) <- paste0("PedigreeNe",seq(from=0,to=(tau)))
    pedigreeNetotal <- as.numeric(NA)
  }
  
  # remove temporary genetic load output from slim_output folder
  if (remove_files){
    if (file.exists(slim_output_pedigreene)){
      file.remove(slim_output_pedigreene)
    }
  }
  
  ## HANDLING SLiM2 OUTPUT t1:t2 - POPULATION ALLELE FREQUENCIES
  ##-------------------------------------------------------------------------------
  
  ### GENOME
  if (tau >= 10){
    filenames_genome_1 <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_t","[0-9]_", sim), full.names=FALSE)
    filenames_genome_2 <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_t","[1-9][0-9]_", sim), full.names=FALSE)
    filenames_genome <- c(filenames_genome_1,filenames_genome_2)
    
    datalist_genome <- lapply(filenames_genome, function(x){read.table(file=paste0(slim_output_folder, x), header=T, na.strings = "NA")})
  
  } else {
    filenames_genome <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_t","[0-9]_", sim), full.names=FALSE)
    datalist_genome <- lapply(filenames_genome, function(x){read.table(file=paste0(slim_output_folder, x), header=T, na.strings = "NA")})
  }
  
  if (model_type == 3){
    all_merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome)
  } else {
    all_merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome)
  }
  
  if(!is.null(all_merged_genome)){
    
    # last generation file
    slim_output_lastgen <- paste0(slim_output_folder,"slim_output_lastgen_", sim, ".txt")
    
    # load file with the last generation 
    lastgen <- scan(file = slim_output_lastgen, what = integer(), na.strings = "NA", quiet = TRUE)
    
    ## IBD and VAR NE CALCULATION
    
    # remove new mutations for Ne calculation
    merged_genome <- all_merged_genome[which(all_merged_genome$GO <= lastgen), ]
    
    ## ALL MUTATION TYPES
    if (nrow(merged_genome) != 0){
      
      he_merged_genome <- merged_genome[ , grepl("^HE" , names(merged_genome))]
      he_merged_genome[is.na(he_merged_genome)] <- 0
      
      mean_he_merged_genome <- colMeans(he_merged_genome, na.rm = TRUE) 
      
      paaf_merged_genome <- merged_genome[ , grepl("^PAAF" , names(merged_genome))]
      paaf_merged_genome[is.na(paaf_merged_genome)] <- 0
      
      # IBD NE Genome-wide for each time interval
      IBDNeGWtimes <- -(1/(2*log(mean_he_merged_genome[2:length(mean_he_merged_genome)]/mean_he_merged_genome[1:(length(mean_he_merged_genome)-1)])))
      IBDNeGWtimes <- as.data.frame(t(IBDNeGWtimes))
      colnames(IBDNeGWtimes) <- paste0("IBDNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      # IBD NE Genome-wide total
      #IBDNeGWtotal <- -(tau/(2*log(mean_he_merged_genome[length(mean_he_merged_genome)]/mean_he_merged_genome[1])))
      IBDNeGWtotal <- -(tau/(2*log(mean_he_merged_genome[paste0("HE",tau)]/mean_he_merged_genome[paste0("HE0")])))
      names(IBDNeGWtotal) <- "IBDNeGWtotal"
      
      # VAR NE Genome-wide for each time interval
      VARNeGWtimes <- mean_he_merged_genome[1:tau]/(2*colMeans((paaf_merged_genome[,1:tau] - paaf_merged_genome[,2:(tau+1)])^2))
      VARNeGWtimes <- as.data.frame(t(VARNeGWtimes))
      colnames(VARNeGWtimes) <- paste0("VARNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      # VAR NE Genome-wide total
      #VARNeGWtotal <- tau*mean_he_merged_genome[1]/(2*mean((paaf_merged_genome[,1] - paaf_merged_genome[,(tau+1)])^2))
      VARNeGWtotal <- tau*mean_he_merged_genome[paste0("HE0")]/(2*mean((paaf_merged_genome[,paste0("PAAF0")] - paaf_merged_genome[,paste0("PAAF",tau)])^2))
      names(VARNeGWtotal) <- "VARNeGWtotal"
      
    } else {
      
      IBDNeGWtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(IBDNeGWtimes) <- paste0("IBDNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      IBDNeGWtotal <- as.numeric(NA)
      names(IBDNeGWtotal) <- "IBDNeGWtotal"
      
      VARNeGWtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(VARNeGWtimes) <- paste0("VARNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      VARNeGWtotal <- as.numeric(NA)
      names(VARNeGWtotal) <- "VARNeGWtotal"
    }
    
    ## ONLY NEUTRAL MUTATIONS
    if (any(merged_genome$MT == 1 | merged_genome$MT == 2 )){
      
      m1_and_m2 <- c(1, 2)
      merged_genome_neutral <- merged_genome[merged_genome$MT %in% m1_and_m2,]
      
      he_merged_genome_neutral <- merged_genome_neutral[ , grepl("HE" , names(merged_genome_neutral))]
      he_merged_genome_neutral[is.na(he_merged_genome_neutral)] <- 0
      
      mean_he_merged_genome_neutral <- colMeans(he_merged_genome_neutral, na.rm = TRUE)
      
      paaf_merged_genome_neutral <- merged_genome_neutral[ , grepl("^PAAF" , names(merged_genome_neutral))]
      paaf_merged_genome_neutral[is.na(paaf_merged_genome_neutral)] <- 0
      
      # IBD NE Genome-wide neutral for each time interval
      IBDNeGWNtimes <- -(1/(2*log(mean_he_merged_genome_neutral[2:length(mean_he_merged_genome_neutral)]/mean_he_merged_genome_neutral[1:(length(mean_he_merged_genome_neutral)-1)])))
      IBDNeGWNtimes <- as.data.frame(t(IBDNeGWNtimes))
      colnames(IBDNeGWNtimes) <- paste0("IBDNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      # IBD NE Genome-wide neutral total
      #IBDNeGWNtotal <- -(tau/(2*log(mean_he_merged_genome_neutral[length(mean_he_merged_genome_neutral)]/mean_he_merged_genome_neutral[1])))
      IBDNeGWNtotal <- -(tau/(2*log(mean_he_merged_genome_neutral[paste0("HE",tau)]/mean_he_merged_genome_neutral[paste0("HE0")])))
      names(IBDNeGWNtotal) <- "IBDNeGWNtotal"
      
      # VAR NE Genome-wide neutral for each time interval
      VARNeGWNtimes <- mean_he_merged_genome_neutral[1:tau]/(2*colMeans((paaf_merged_genome_neutral[,1:tau] - paaf_merged_genome_neutral[,2:(tau+1)])^2))
      VARNeGWNtimes <- as.data.frame(t(VARNeGWNtimes))
      colnames(VARNeGWNtimes) <- paste0("VARNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      # VAR NE Genome-wide neutral total
      #VARNeGWNtotal <- tau*mean_he_merged_genome_neutral[1]/(2*mean((paaf_merged_genome_neutral[,1] - paaf_merged_genome_neutral[,(tau+1)])^2))
      VARNeGWNtotal <- tau*mean_he_merged_genome_neutral[paste0("HE0")]/(2*mean((paaf_merged_genome_neutral[,paste0("PAAF0")] - paaf_merged_genome_neutral[,paste0("PAAF",tau)])^2))
      names(VARNeGWNtotal) <- "VARNeGWNtotal"
      
    } else {
      
      IBDNeGWNtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(IBDNeGWNtimes) <- paste0("IBDNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      IBDNeGWNtotal <- as.numeric(NA)
      names(IBDNeGWNtotal) <- "IBDNeGWNtotal"
      
      VARNeGWNtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(VARNeGWNtimes) <- paste0("VARNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
      
      VARNeGWNtotal <- as.numeric(NA)
      names(VARNeGWNtotal) <- "VARNeGWNtotal"
    }
    
    if (model_type == 3){
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
      if (any(all_merged_genome[, paste0("S", tc)] != 0, na.rm = TRUE)){
        
        actual_pop_prbe <- sum(all_merged_genome[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
        
        if (any(N*all_merged_genome[, paste0("S",tc)] > 1, na.rm = TRUE)){
          strong_positive_prbe <- sum(N*all_merged_genome[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
        } else {
          strong_positive_prbe <- as.numeric(0)
        }
        
        if (any(N*all_merged_genome[, paste0("S",tc)] < -1, na.rm = TRUE)){
          strong_negative_prbe <- sum(N*all_merged_genome[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
        } else {
          strong_negative_prbe <- as.numeric(0)
        }
        
        strong_pop_prbe <- strong_positive_prbe + strong_negative_prbe 
        
      } else {
        actual_pop_prbe <- as.numeric(0)
        strong_pop_prbe <- as.numeric(0)
      }
      
      ## IBD AND VAR NE ONLY SELECTED MUTATIONS
      if (any(merged_genome[, paste0("S",tc)] != 0, na.rm = TRUE)){
        
        merged_genome_selection <- merged_genome[which(merged_genome[, paste0("S",tc)] != 0), ]
      
        if (!is.null(merged_genome_selection)){
          
          he_merged_genome_selection <- merged_genome_selection[ , grepl("HE" , names(merged_genome_selection))]
          he_merged_genome_selection[is.na(he_merged_genome_selection)] <- 0
          
          mean_he_merged_genome_selection <- colMeans(he_merged_genome_selection, na.rm = TRUE)
          
          paaf_merged_genome_selection <- merged_genome_selection[ , grepl("^PAAF" , names(merged_genome_selection))]
          paaf_merged_genome_selection[is.na(paaf_merged_genome_selection)] <- 0
          
          # IBD NE Genome-wide selection for each time interval
          IBDNeGWStimes <- -(1/(2*log(mean_he_merged_genome_selection[2:length(mean_he_merged_genome_selection)]/mean_he_merged_genome_selection[1:(length(mean_he_merged_genome_selection)-1)])))
          IBDNeGWStimes <- as.data.frame(t(IBDNeGWStimes))
          colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
          
          # IBD NE Genome-wide selection total
          #IBDNeGWStotal <- -(tau/(2*log(mean_he_merged_genome_selection[length(mean_he_merged_genome_selection)]/mean_he_merged_genome_selection[1])))
          IBDNeGWStotal <- -(tau/(2*log(mean_he_merged_genome_selection[paste0("HE",tau)]/mean_he_merged_genome_selection[paste0("HE0")])))
          names(IBDNeGWStotal) <- "IBDNeGWStotal"
          
          # VAR NE Genome-wide selection for each time interval
          VARNeGWStimes <- mean_he_merged_genome_selection[1:tau]/(2*colMeans((paaf_merged_genome_selection[,1:tau] - paaf_merged_genome_selection[,2:(tau+1)])^2))
          VARNeGWStimes <- as.data.frame(t(VARNeGWStimes))
          colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
          
          # VAR NE Genome-wide neutral total
          #VARNeGWStotal <- tau*mean_he_merged_genome_selection[1]/(2*mean((paaf_merged_genome_selection[,1] - paaf_merged_genome_selection[,(tau+1)])^2))
          VARNeGWStotal <- tau*mean_he_merged_genome_selection[paste0("HE0")]/(2*mean((paaf_merged_genome_selection[,paste0("PAAF0")] - paaf_merged_genome_selection[,paste0("PAAF",tau)])^2))
          names(VARNeGWStotal) <- "VARNeGWStotal"
          
        } else {
          
          IBDNeGWStimes <- as.data.frame(t(rep(NA, tau)))
          colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
          
          IBDNeGWStotal <- as.numeric(NA)
          names(IBDNeGWStotal) <- "IBDNeGWStotal"
          
          VARNeGWStimes <- as.data.frame(t(rep(NA, tau)))
          colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
          
          VARNeGWStotal <- as.numeric(NA)
          names(VARNeGWStotal) <- "VARNeGWStotal"
        }
        
      } else {
        
        IBDNeGWStimes <- as.data.frame(t(rep(NA, tau)))
        colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        IBDNeGWStotal <- as.numeric(NA)
        names(IBDNeGWStotal) <- "IBDNeGWStotal"
        
        VARNeGWStimes <- as.data.frame(t(rep(NA, tau)))
        colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        VARNeGWStotal <- as.numeric(NA)
        names(VARNeGWStotal) <- "VARNeGWStotal"
      }
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
      if (any(all_merged_genome[, "S"] != 0, na.rm = TRUE)) {
        
        actual_pop_prbe <- sum(all_merged_genome[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
        
        if (any(N*all_merged_genome[, "S"] > 1, na.rm = TRUE)){
          strong_positive_prbe <- sum(N*all_merged_genome[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
        } else {
          strong_positive_prbe <- as.numeric(0)
        }
        
        if (any(N*all_merged_genome[, "S"] < -1, na.rm = TRUE)){
          strong_negative_prbe <- sum(N*all_merged_genome[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
        } else {
          strong_negative_prbe <- as.numeric(0)
        }
        
        strong_pop_prbe <- strong_positive_prbe + strong_negative_prbe 
        
      } else {
        actual_pop_prbe <- as.numeric(0)
        strong_pop_prbe <- as.numeric(0)
      }
      
      ## IBD AND VAR NE ONLY SELECTED MUTATIONS
      if (any(merged_genome[, "S"] != 0, na.rm = TRUE)){
        
        merged_genome_selection <- merged_genome[merged_genome$S != 0, ]
        
        he_merged_genome_selection <- merged_genome_selection[ , grepl("HE" , names(merged_genome_selection))]
        he_merged_genome_selection[is.na(he_merged_genome_selection)] <- 0
        
        mean_he_merged_genome_selection <- colMeans(he_merged_genome_selection, na.rm = TRUE)
        
        paaf_merged_genome_selection <- merged_genome_selection[ , grepl("^PAAF" , names(merged_genome_selection))]
        paaf_merged_genome_selection[is.na(paaf_merged_genome_selection)] <- 0
        
        # IBD NE Genome-wide selection for each time interval
        IBDNeGWStimes <- -(1/(2*log(mean_he_merged_genome_selection[2:length(mean_he_merged_genome_selection)]/mean_he_merged_genome_selection[1:(length(mean_he_merged_genome_selection)-1)])))
        IBDNeGWStimes <- as.data.frame(t(IBDNeGWStimes))
        colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        # IBD NE Genome-wide selection total
        #IBDNeGWStotal <- -(tau/(2*log(mean_he_merged_genome_selection[length(mean_he_merged_genome_selection)]/mean_he_merged_genome_selection[1])))
        IBDNeGWStotal <- -(tau/(2*log(mean_he_merged_genome_selection[paste0("HE",tau)]/mean_he_merged_genome_selection[paste0("HE0")])))
        names(IBDNeGWStotal) <- "IBDNeGWStotal"
        
        # VAR NE Genome-wide selection for each time interval
        VARNeGWStimes <- mean_he_merged_genome_selection[1:tau]/(2*colMeans((paaf_merged_genome_selection[,1:tau] - paaf_merged_genome_selection[,2:(tau+1)])^2))
        VARNeGWStimes <- as.data.frame(t(VARNeGWStimes))
        colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        # VAR NE Genome-wide neutral total
        #VARNeGWStotal <- tau*mean_he_merged_genome_selection[1]/(2*mean((paaf_merged_genome_selection[,1] - paaf_merged_genome_selection[,(tau+1)])^2))
        VARNeGWStotal <- tau*mean_he_merged_genome_selection[paste0("HE0")]/(2*mean((paaf_merged_genome_selection[,paste0("PAAF0")] - paaf_merged_genome_selection[,paste0("PAAF",tau)])^2))
        names(VARNeGWStotal) <- "VARNeGWStotal"
        
      } else {
        
        IBDNeGWStimes <- as.data.frame(t(rep(NA, tau)))
        colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        IBDNeGWStotal <- as.numeric(NA)
        names(IBDNeGWStotal) <- "IBDNeGWStotal"
        
        VARNeGWStimes <- as.data.frame(t(rep(NA, tau)))
        colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
        
        VARNeGWStotal <- as.numeric(NA)
        names(VARNeGWStotal) <- "VARNeGWStotal"
      }
    }
  } else {
    
    ## IBD and VAR NE CALCULATION
    
    ## ALL MUTATION TYPES
    IBDNeGWtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(IBDNeGWtimes) <- paste0("IBDNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    IBDNeGWtotal <- as.numeric(NA)
    names(IBDNeGWtotal) <- "IBDNeGWtotal"
    
    VARNeGWtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(VARNeGWtimes) <- paste0("VARNeGW",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    VARNeGWtotal <- as.numeric(NA)
    names(VARNeGWtotal) <- "VARNeGWtotal"
    
    ## ONLY NEUTRAL MUTATIONS
    IBDNeGWNtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(IBDNeGWNtimes) <- paste0("IBDNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    IBDNeGWNtotal <- as.numeric(NA)
    names(IBDNeGWNtotal) <- "IBDNeGWNtotal"
    
    VARNeGWNtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(VARNeGWNtimes) <- paste0("VARNeGWN",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    VARNeGWNtotal <- as.numeric(NA)
    names(VARNeGWNtotal) <- "VARNeGWNtotal"
    
    ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
    actual_pop_prbe <- as.numeric(0)
    strong_pop_prbe <- as.numeric(0)
    
    ## IBD AND VAR NE ONLY SELECTED MUTATIONS
    IBDNeGWStimes <- as.data.frame(t(rep(NA, tau)))
    colnames(IBDNeGWStimes) <- paste0("IBDNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    IBDNeGWStotal <- as.numeric(NA)
    names(IBDNeGWStotal) <- "IBDNeGWStotal"
    
    VARNeGWStimes <- as.data.frame(t(rep(NA, tau)))
    colnames(VARNeGWStimes) <- paste0("VARNeGWS",seq(from=0,to=(tau-1)),"_",seq(from=1,to=tau))
    
    VARNeGWStotal <- as.numeric(NA)
    names(VARNeGWStotal) <- "VARNeGWStotal"
    
  }
  
  # remove temporary PAAF file from slim_output folder
  if (remove_files){
    if(any(file.exists(paste0(slim_output_folder,filenames_genome)))){
      file.remove(paste0(slim_output_folder,filenames_genome))
    }
  }
  
  ## NEUTRAL MUTATION IN EXTRA CHROMOSOME
  if (tau >= 10){
    filenames_extraChr_1 <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_extraChr_t","[0-9]_", sim), full.names=FALSE)
    filenames_extraChr_2 <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_extraChr_t","[1-9][0-9]_", sim), full.names=FALSE)
    filenames_extraChr <- c(filenames_extraChr_1,filenames_extraChr_2)
    
    datalist_extraChr <- lapply(filenames_extraChr, function(x){read.table(file=paste0(slim_output_folder, x), header=T, na.strings = "NA")})
    
  } else {
    filenames_extraChr <- list.files(path=slim_output_folder, pattern = paste0("slim_output_paaf_extraChr_t","[0-9]_", sim), full.names=FALSE)
    datalist_extraChr <- lapply(filenames_extraChr, function(x){read.table(file=paste0(slim_output_folder, x), header=T, na.strings = "NA")})
  }
  
  merged_extraChr <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_extraChr)
  
  if(!is.null(merged_extraChr)){
    
    # last generation file
    slim_output_lastgen <- paste0(slim_output_folder,"slim_output_lastgen_", sim, ".txt")
    
    # load file with the last generation 
    lastgen <- scan(file = slim_output_lastgen, what = integer(), na.strings = "NA", quiet = TRUE)
    
    # remove new mutations
    merged_extraChr <- merged_extraChr[which(merged_extraChr$GO <= lastgen), ]
    
    if (nrow(merged_extraChr) != 0){
      
      he_merged_extraChr <- merged_extraChr[ , grepl("HE" , names(merged_extraChr))]
      he_merged_extraChr[is.na(he_merged_extraChr)] <- 0
      
      mean_he_merged_extraChr <- colMeans(he_merged_extraChr, na.rm = TRUE)
      
      paaf_merged_extraChr <- merged_extraChr[ , grepl("^PAAF" , names(merged_extraChr))]
      paaf_merged_extraChr[is.na(paaf_merged_extraChr)] <- 0
      
      # IBD NE extra chromosome for each time interval
      IBDNeChrtimes <- -(1/(2*log(mean_he_merged_extraChr[2:length(mean_he_merged_extraChr)]/mean_he_merged_extraChr[1:(length(mean_he_merged_extraChr)-1)])))
      IBDNeChrtimes <- as.data.frame(t(IBDNeChrtimes))
      colnames(IBDNeChrtimes) <- paste0("IBDNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
      
      # IBD NE extra chromosome total
      #IBDNeChrtotal <- -(tau/(2*log(mean_he_merged_extraChr[length(mean_he_merged_extraChr)]/mean_he_merged_extraChr[1])))
      IBDNeChrtotal <- -(tau/(2*log(mean_he_merged_extraChr[paste0("HE",tau)]/mean_he_merged_extraChr[paste0("HE0")])))
      names(IBDNeChrtotal) <- "IBDNeChrtotal"
      
      # VAR NE extra chromosome for each time interval
      VARNeChrtimes <- mean_he_merged_extraChr[1:tau]/(2*colMeans((paaf_merged_extraChr[,1:tau] - paaf_merged_extraChr[,2:(tau+1)])^2))
      VARNeChrtimes <- as.data.frame(t(VARNeChrtimes))
      colnames(VARNeChrtimes) <- paste0("VARNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
      
      # VAR NE extra chromosome total
      #VARNeChrtotal <- tau*mean_he_merged_extraChr[1]/(2*mean((paaf_merged_extraChr[,1] - paaf_merged_extraChr[,(tau+1)])^2))
      VARNeChrtotal <- tau*mean_he_merged_extraChr[paste0("HE0")]/(2*mean((paaf_merged_extraChr[,paste0("PAAF0")] - paaf_merged_extraChr[,paste0("PAAF",tau)])^2))
      names(VARNeChrtotal) <- "VARNeChrtotal"
      
    } else {
      
      IBDNeChrtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(IBDNeChrtimes) <- paste0("IBDNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
      
      IBDNeChrtotal <- as.numeric(NA)
      names(IBDNeChrtotal) <- "IBDNeChrtotal"
      
      VARNeChrtimes <- as.data.frame(t(rep(NA, tau)))
      colnames(VARNeChrtimes) <- paste0("VARNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
      
      VARNeChrtotal <- as.numeric(NA)
      names(VARNeChrtotal) <- "VARNeChrtotal"
    } 
  } else {
    
    IBDNeChrtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(IBDNeChrtimes) <- paste0("IBDNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
    
    IBDNeChrtotal <- as.numeric(NA)
    names(IBDNeChrtotal) <- "IBDNeChrtotal"
    
    VARNeChrtimes <- as.data.frame(t(rep(NA, tau)))
    colnames(VARNeChrtimes) <- paste0("VARNeChr",seq(from=0,to=(tau-1)),"_",seq(from=1, to=tau))
    
    VARNeChrtotal <- as.numeric(NA)
    names(VARNeChrtotal) <- "VARNeChrtotal"
  }
  
  # remove temporary PAAF file from slim_output folder
  if (remove_files){
    if(any(file.exists(paste0(slim_output_folder,filenames_extraChr)))){
      file.remove(paste0(slim_output_folder,filenames_extraChr))
    }
  }
  
  ## HANDLINDING SLiM2 OUTPUT - SAMPLED INDIVIDUALS
  ##---------------------------------------------------------------------------------
  
  # sort vcf files
  slim_output_sample_t1        <- paste0(slim_output_folder,"slim_output_sample_t1_", sim, ".vcf")
  slim_output_sample_t1_sorted <- paste0(slim_output_folder,"slim_output_sample_t1_", sim, "_sorted" , ".vcf")
  
  slim_output_sample_t2        <- paste0(slim_output_folder,"slim_output_sample_t2_", sim, ".vcf")
  slim_output_sample_t2_sorted <- paste0(slim_output_folder,"slim_output_sample_t2_", sim, "_sorted" , ".vcf")
  
  sort_sample_t1_vcf <- paste("grep '^#'", slim_output_sample_t1, ">", slim_output_sample_t1_sorted,
                              "&& grep -v '^#'", slim_output_sample_t1, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_t1_sorted)
  
  sort_sample_t2_vcf <- paste("grep '^#'", slim_output_sample_t2, ">", slim_output_sample_t2_sorted,
                             "&& grep -v '^#'", slim_output_sample_t2, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_t2_sorted)
  
  system(sort_sample_t1_vcf)
  system(sort_sample_t2_vcf)
  
  # bgzip sorted vcf files
  system(paste(path_to_bgzip, "-f", slim_output_sample_t1_sorted))
  system(paste(path_to_bgzip, "-f", slim_output_sample_t2_sorted))
  
  # tabix bgziped files
  system(paste(path_to_tabix, paste0(slim_output_sample_t1_sorted, ".gz")))
  system(paste(path_to_tabix, paste0(slim_output_sample_t2_sorted, ".gz")))
  
  # merge and get the data
  slim_output_sample_merged <- paste0(slim_output_folder,"slim_output_sample_merged_", sim, ".txt")
  
  bcftools_query <- paste(path_to_bcftools, "merge --force-samples",
                          paste0(slim_output_sample_t1_sorted, ".gz"),
                          paste0(slim_output_sample_t2_sorted, ".gz"),
                          "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                          "-e 'MT=4'",
                          ">", slim_output_sample_merged) 
  
  system(bcftools_query)
  
  ## READ SLiM OUTPUT 1 AND CONVERT TO EGGLIB INPUT
  
  # assembly the header
  header_1       <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
  sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
  sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
  
  full_header <- c(header_1, sample_names_1, sample_names_2)
  
  # imported the data
  slim_raw_data <- read.table(file = slim_output_sample_merged, header = F, col.names = full_header, check.names = F, na.strings = "./.")
  
  # remove temporary vcf and post processed vcf file from slim_output folder
  if (remove_files){
    file.remove(paste0(slim_output_sample_t1))
    file.remove(paste0(slim_output_sample_t1_sorted, ".gz"))
    file.remove(paste0(slim_output_sample_t1_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_sample_t2))
    file.remove(paste0(slim_output_sample_t2_sorted, ".gz"))
    file.remove(paste0(slim_output_sample_t2_sorted, ".gz.tbi"))
    file.remove(paste0(slim_output_sample_merged))
    file.remove(paste0(slim_output_folder,"slim_coalesced_",model_title,"_", sim, ".tree"))
    file.remove(paste0(slim_output_folder,"slim_output_lastgen_", sim, ".txt"))
  }
  
  # if it is a RADseq data
  if (data_type == 2){
    slim_raw_data <- slim_raw_data[which(slim_raw_data$position %in% rad_interval), ]
  }
  
  if (nrow(slim_raw_data) != 0){
    
    # split the data
    slim_snp_geno <- slim_raw_data[, 9:ncol(slim_raw_data)]
    slim_snp_info <- slim_raw_data[, 1:length(header_1)]
    
    # change the genotype annotations
    slim_snp_geno <- as.matrix(slim_snp_geno)
    slim_snp_geno[is.na(slim_snp_geno)]   <- "11"
    slim_snp_geno[slim_snp_geno == "0|0"] <- "11"
    slim_snp_geno[slim_snp_geno =="1|1"]  <- "22"
    
    if (haplotype){
      ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
      if (ref_or_alt == 0){
        slim_snp_geno[slim_snp_geno == "0|1" | slim_snp_geno == "1|0"] <- "11"
      } else {
        slim_snp_geno[slim_snp_geno == "0|1" | slim_snp_geno == "1|0"] <- "22"
      }
    } else {
      slim_snp_geno[slim_snp_geno == "0|1"] <- "12"
      slim_snp_geno[slim_snp_geno == "1|0"] <- "21"
    }
    
    slim_snp_geno <- as.data.frame(slim_snp_geno)
    
    # mark monomophormic mutations (all 11 or 22)
    count_ref_geno    <- apply(slim_snp_geno, 1, function(x){sum(x == 11)})
    count_alt_geno    <- apply(slim_snp_geno, 1, function(x){sum(x == 22)})
    keep_snps         <- count_ref_geno < (SS1 + SS2) & count_alt_geno < (SS1 + SS2) # MARK MONOMORPHIC MUTATIONS
    
    # adding missing data randomly - BEFORE OR AFTER REMOVE MONOMORPHIC?
    slim_snp_geno <- as.matrix(as.data.frame(lapply(slim_snp_geno, function(mi) mi[sample(c(TRUE, NA), prob = c((1-missing_data), missing_data), size = length(mi), replace = TRUE)])))
    slim_snp_geno[is.na(slim_snp_geno)] <- "00"
    slim_snp_geno <- as.data.frame(slim_snp_geno)
    
    colnames(slim_snp_geno) <- c(sample_names_1,sample_names_2) # I re-do it here because the missing data part is messing with the header
    
    # re-assemble the data
    slim_data <- cbind(slim_snp_info, slim_snp_geno)
    
    # remove monomorphic mutations
    slim_data <- slim_data[keep_snps, ]
    
    # remove duplicated mutations
    slim_data <- slim_data[!duplicated(slim_data[ ,1:2]), ]
    
    # re-code the chromosome name
    #if (chrTAG){
    #  if(chrN > 1){
    #    slim_data$chrom <- sapply(slim_data$position, chromtagging)
    #  }
    #}
    
    # re-code the chromosome name
    if(chrN > 1){
      if (chrTAG){
        slim_data$chrom <- sapply(slim_data$position, chromtagging)
      }
    }
    
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
    egglib_converted_file <- paste0("egglib_input_sample", "_", sim, ".txt");
    write.table(slim_to_egglib_data, file = paste0(egglib_input_folder, egglib_converted_file), quote=FALSE, sep="\t", row.names = FALSE)
    
    # save only the information of the snps
    if (model_type == 3){
      selcoeff_snps <- all_merged_genome[all_merged_genome$MID %in% slim_data$MID, paste0("S",tc)]
      } else {
      selcoeff_snps <- all_merged_genome[all_merged_genome$MID %in% slim_data$MID, "S"]
    }
    
    slim_to_egglib_snps <- cbind(ID=paste0(slim_data$chrom, ":", slim_data$position), 
                                 MID=slim_data$MID,
                                 MT=slim_data$MT, 
                                 S=selcoeff_snps,
                                 DOM=slim_data$DOM, 
                                 GO=slim_data$GO,
                                 slim_data[, (length(header_1)+1):ncol(slim_data)])
    
    ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
    ##-----------------------------------------------------------------------------------
    
    # check if the folder exists
    if (!file_test("-d", egglib_output_folder)){
      dir.create(file.path(egglib_output_folder))
    }
    
    # generate text with slim command  
    egglib_run <- paste(path_to_python,
                        paste0(getwd(), "/", path_to_egglib_summstat),
                        paste0("input-file=", egglib_input_folder, egglib_converted_file),
                        paste0("output-file=", egglib_output_folder, "egglib_output_sample", "_", sim, ".txt"),
                        paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
                        paste0("LSSp=", paste0(c("He", "Dj"), collapse = ",")),
                        paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "Pi", "D", "Da", "ZZ", "ZnS"), collapse = ",")),
                        paste0("WSSp=", paste0(c("He", "Dj", "S", "thetaW", "Pi", "D", "ZZ", "ZnS"), collapse = ",")),
                        paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW" , "Pi", "D", "Da", "SFS"), collapse = ",")),
                        paste0("GSSp=",paste0(c("He", "S", "thetaW", "Pi", "D", "Da"), collapse = ",")),
                        paste0("wspan=", wss_wspan_run),
                        paste0("SFS-bins=", sfs_bins_run),
                        paste0("select=", "all"));
    
    # run egglib summstat on system
    system(egglib_run)
    
    if (add_WSSwspan_SFSbins_1){
      egglib_run_1 <- paste(path_to_python,
                            paste0(getwd(), "/", path_to_egglib_summstat),
                            paste0("input-file=", egglib_input_folder, egglib_converted_file),
                            paste0("output-file=", egglib_output_folder, "egglib_output_sample_addwspan1", "_", sim, ".txt"),
                            paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "Pi", "D", "Da", "ZZ", "ZnS"), collapse = ",")),
                            paste0("WSSp=", paste0(c("He", "Dj", "S", "thetaW", "Pi", "D", "ZZ", "ZnS"), collapse = ",")),
                            paste0("GSS=", paste0(c("SFS"), collapse = ",")),
                            paste0("wspan=", add_wss_wspan_1),
                            paste0("SFS-bins=", add_sfs_bins_1),
                            paste0("select=", "all"));
      
      # run egglib summstat on system
      system(egglib_run_1)
    }
    
    if (add_WSSwspan_SFSbins_2){
      egglib_run_2 <- paste(path_to_python,
                            paste0(getwd(), "/", path_to_egglib_summstat),
                            paste0("input-file=", egglib_input_folder, egglib_converted_file),
                            paste0("output-file=", egglib_output_folder, "egglib_output_sample_addwspan2", "_", sim, ".txt"),
                            paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "Pi", "D", "Da", "ZZ", "ZnS"), collapse = ",")),
                            paste0("WSSp=", paste0(c("He", "Dj", "S", "thetaW", "Pi", "D", "ZZ", "ZnS"), collapse = ",")),
                            paste0("GSS=", paste0(c("SFS"), collapse = ",")),
                            paste0("wspan=", add_wss_wspan_2),
                            paste0("SFS-bins=", add_sfs_bins_2),
                            paste0("select=", "all"));
      
      # run egglib summstat on system
      system(egglib_run_2)
    }
    
    ## READ EGGLIB OUTPUT AND MAKE THE REFERENCE TABLE
    ##----------------------------------------------------------------------------------
    
    # import egglib output
    egglib_summary_stats <- read.csv(file = paste0(egglib_output_folder,"egglib_output_sample", "_", sim, ".txt"), header = T, sep = "\t", check.names = F)
    
    ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
    if(any(slim_to_egglib_snps$S != 0, na.rm = TRUE)){
      actual_sample_prbe <- sum(slim_to_egglib_snps$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
    } else {
      actual_sample_prbe <- as.numeric(0)
    }
    
    ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
    if(any(N*slim_to_egglib_snps$S > 1, na.rm = TRUE)){
      positive_sample_prbe <- sum(N*slim_to_egglib_snps$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
    } else {
      positive_sample_prbe <- as.numeric(0)
    }
    
    ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
    if(any(N*slim_to_egglib_snps$S < -1, na.rm = TRUE)){
      negative_sample_prbe <- sum(N*slim_to_egglib_snps$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
    } else {
      negative_sample_prbe <- as.numeric(0)
    }
    
    strong_sample_prbe = positive_sample_prbe + negative_sample_prbe
    
    ## GLOBAL SUMMARY STATISTICS
    # remove redundant summary statistics
    egglib_summary_stats <- egglib_summary_stats[, unique(names(egglib_summary_stats))]
    
    # rename the summary statistics
    colnames(egglib_summary_stats) <- gsub(":", "_", names(egglib_summary_stats))
    
    # egglib calculated GLOBAL statistics
    global_stats <- egglib_summary_stats[1 , grepl("^GSS" , unique(names(egglib_summary_stats)))]
    
    global_SFS   <- egglib_summary_stats[1 , grepl("^SFS" , unique(names(egglib_summary_stats)))]
    
    # calculate additional GLOBAL summary statistics
    mean_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){mean(x, na.rm=T)})
    
    var_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){var(x, na.rm=T)})
    
    kurt_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){kurtosis(x, na.rm=T)})
    
    skew_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){skewness(x, na.rm=T)})
    
    q05_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
    
    q95_locus_stats <- apply(egglib_summary_stats[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats))))[1]:length(egglib_summary_stats))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
    
    # assemble additional GLOBAL summary statistics
    add_global_stats <-cbind(as.data.frame(t(mean_locus_stats)), as.data.frame(t(var_locus_stats)), as.data.frame(t(kurt_locus_stats)), 
                             as.data.frame(t(skew_locus_stats)), as.data.frame(t(q05_locus_stats)), as.data.frame(t(q95_locus_stats)))
    
    # ASSEMBLY default GLOBAL summary statistics
    global_summary_stats <- cbind(global_stats, global_SFS, add_global_stats)
    colnames(global_summary_stats) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats)[2]))
    
    ## ADDITIONAL SUMMARY STATS IF DIFFERENTE WSPAN AND SFS BINS WERE SET UP
    if (add_WSSwspan_SFSbins_1){
      
      # import egglib output
      egglib_summary_stats_add1 <- read.csv(file = paste0(egglib_output_folder,"egglib_output_sample_addwspan1", "_", sim, ".txt"), header = T, sep = "\t", check.names = F)
      
      # remove redundant summary statistics
      egglib_summary_stats_add1 <- egglib_summary_stats_add1[, unique(names(egglib_summary_stats_add1))]
      
      # rename the summary statistics
      colnames(egglib_summary_stats_add1) <- gsub(":", "_", names(egglib_summary_stats_add1))
      
      # egglib calculated GLOBAL statistics
      global_SFS_bins1 <- egglib_summary_stats_add1[1 , grepl("^SFS" , unique(names(egglib_summary_stats_add1)))]
      
      # calculate additional GLOBAL summary statistics
      mean_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){mean(x, na.rm=T)})
      
      var_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){var(x, na.rm=T)})
      
      kurt_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){kurtosis(x, na.rm=T)})
      
      skew_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){skewness(x, na.rm=T)})
      
      q05_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
      
      q95_wss_wspan1 <- apply(egglib_summary_stats_add1[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add1)))))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
      
      # assemble additional GLOBAL summary statistics
      add_global_stats_add1 <-cbind(global_SFS_bins1,
                                    as.data.frame(t(mean_wss_wspan1)),as.data.frame(t(var_wss_wspan1)),as.data.frame(t(kurt_wss_wspan1)), 
                                    as.data.frame(t(skew_wss_wspan1)),as.data.frame(t(q05_wss_wspan1)),as.data.frame(t(q95_wss_wspan1))
      )
      
      colnames(add_global_stats_add1) <- paste0("ADD_1_GSS","_",seq(from=1,to=dim(add_global_stats_add1)[2]))
      
    }
    
    if (add_WSSwspan_SFSbins_2){
      
      # import egglib output
      egglib_summary_stats_add2 <- read.csv(file = paste0(egglib_output_folder,"egglib_output_sample_addwspan2", "_", sim, ".txt"), header = T, sep = "\t", check.names = F)
      
      # remove redundant summary statistics
      egglib_summary_stats_add2 <- egglib_summary_stats_add2[, unique(names(egglib_summary_stats_add2))]
      
      # rename the summary statistics
      colnames(egglib_summary_stats_add2) <- gsub(":", "_", names(egglib_summary_stats_add2))
      
      # egglib calculated GLOBAL statistics
      global_SFS_bins2 <- egglib_summary_stats_add2[1 , grepl("^SFS" , unique(names(egglib_summary_stats_add2)))]
      
      # calculate additional GLOBAL summary statistics
      mean_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){mean(x, na.rm=T)})
      
      var_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){var(x, na.rm=T)})
      
      kurt_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){kurtosis(x, na.rm=T)})
      
      skew_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){skewness(x, na.rm=T)})
      
      q05_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
      
      q95_wss_wspan2 <- apply(egglib_summary_stats_add2[,-c(1, which(grepl("^SFS_" , unique(names(egglib_summary_stats_add2)))))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
      
      # assemble additional GLOBAL summary statistics
      add_global_stats_add2 <-cbind(global_SFS_bins2,
                                    as.data.frame(t(mean_wss_wspan2)),as.data.frame(t(var_wss_wspan2)),as.data.frame(t(kurt_wss_wspan2)), 
                                    as.data.frame(t(skew_wss_wspan2)),as.data.frame(t(q05_wss_wspan2)),as.data.frame(t(q95_wss_wspan2))
      )
      
      colnames(add_global_stats_add2) <- paste0("ADD_2_GSS","_",seq(from=1,to=dim(add_global_stats_add2)[2]))
      
    }
    
    # ASSEMBLY FINAL GLOBAL summary statistics
    if (exists("add_global_stats_add1") & exists("add_global_stats_add2")){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add1, add_global_stats_add2)
    } else if (exists("add_global_stats_add1")){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add1)
    } else if (exists("add_global_stats_add2")){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add2)
    }
    
    # LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
    slim_to_egglib_snps$S[is.na(slim_to_egglib_snps$S)] <- 0
    
    raw_locusFST_table <- data.frame(Ns=N*slim_to_egglib_snps$S, LSS_WCst=egglib_summary_stats[, "LSS_WCst"])
    
    raw_locusFST_table$labelsNs <- ifelse(raw_locusFST_table$Ns > 1, 1, 0)
    
    ord_locusFST_table <- raw_locusFST_table[order(-raw_locusFST_table$LSS_WCst), ]
    
    if (any(ord_locusFST_table$Ns > 1)){
      if(!all(ord_locusFST_table$Ns > 1)){
        
        pred_locusFSTNS <- prediction(predictions = ord_locusFST_table$LSS_WCst, labels = ord_locusFST_table$labelsNs)
        perf_locusFSTNS <- performance(pred_locusFSTNS, "fpr", "prec")
        
        perf_table <- data.frame(fprNS=perf_locusFSTNS@y.values[[1]], precNS=perf_locusFSTNS@x.values[[1]],
                                 fdrNS=1-perf_locusFSTNS@x.values[[1]])
        
        perf_locusFST_table <- data.frame(ord_locusFST_table[1:dim(perf_table)[1], ], perf_table)
        
        FSTfdr <- as.data.frame(cbind(FSTfdrNS05 = perf_locusFST_table[which.min(perf_locusFST_table$fdrNS <= 0.05) ,"LSS_WCst"],
                                      FSTfdrNS10 = perf_locusFST_table[which.min(perf_locusFST_table$fdrNS <= 0.10) ,"LSS_WCst"],
                                      FSTfdrNS25 = perf_locusFST_table[which.min(perf_locusFST_table$fdrNS <= 0.25) ,"LSS_WCst"]))
        
        FSTprec <- as.data.frame(cbind(FSTprecNS75 = perf_locusFST_table[which.min(perf_locusFST_table$precNS >= 0.75) ,"LSS_WCst"],
                                       FSTprecNS95 = perf_locusFST_table[which.min(perf_locusFST_table$precNS >= 0.95) ,"LSS_WCst"]))
        
      } else {
        
        FSTfdr <- as.data.frame(cbind(FSTfdrNS05 = 0,
                                      FSTfdrNS10 = 0,
                                      FSTfdrNS25 = 0))
        
        FSTprec <- as.data.frame(cbind(FSTprecNS75 = 0,
                                       FSTprecNS95 = 0))
      }
    } else {
      
      FSTfdr <- as.data.frame(cbind(FSTfdrNS05 = 1,
                                    FSTfdrNS10 = 1,
                                    FSTfdrNS25 = 1))
      
      FSTprec <- as.data.frame(cbind(FSTprecNS75 = 1,
                                     FSTprecNS95 = 1))
      
    }
    
    ## LOCUS-SPECIFIC SUMMARY STATISTICS
    # sampling snps for the locus-specific reference table
    #if (any(slim_to_egglib_snps$MT == 1)){
    #  m1 <- sample(which(slim_to_egglib_snps$MT == 1), size=1)
    #} else {
    #  m1 <- 0
    #}
    #
    #if (any(slim_to_egglib_snps$MT == 2)){
    #  m2 <- sample(which(slim_to_egglib_snps$MT == 2), size=1)
    #} else {
    #  m2 <- 0
    #}
    #
    #if (any(slim_to_egglib_snps$MT == 3)){
    #  m3 <- sample(which(slim_to_egglib_snps$MT == 3), size=1)
    #} else {
    #  m3 <- 0
    #}
    
    #if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 2)){
    #  m1m2 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2), size=1)
    #} else {
    #  m1m2 <- 0
    #}
    
    #if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 3)){
    #  m1m3 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 3), size=1)
    #} else {
    #  m1m3 <- 0
    #}
    
    #if (any(slim_to_egglib_snps$MT == 2) | any(slim_to_egglib_snps$MT == 3)){
    #  m2m3 <- sample(which(slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
    #} else {
    #  m2m3 <- 0
    #}
    
    #if (any(slim_to_egglib_snps$MT == 1) | any(slim_to_egglib_snps$MT == 2) | any(slim_to_egglib_snps$MT == 3)){
    #  m1m2m3 <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
    #} else {
    #  m1m2m3 <- 0
    #}
    
    # vector of indexes
    #snps_in_reftable <- c(m1,m2,m3,m1m2,m1m3,m2m3,m1m2m3)
    
    # sample one RANDOM mutation
    snps_in_reftable <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
    sampled_snp <- slim_to_egglib_snps[snps_in_reftable, ]
    
    # calculate sample minor allele frequency
    sampled_snp_genotypes <- sampled_snp[, (grep("GO", names(sampled_snp)) + 1):(SS1 + SS2 + 5)]
    
    # sample alternative allele frequency - SAAF1
    S1_genotypes <- sampled_snp_genotypes[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes))]
    S1n11 <- apply(S1_genotypes==11, 1, sum, na.rm=T)
    S1n12 <- apply(S1_genotypes==12 | S1_genotypes==21, 1, sum, na.rm=T)
    S1n22 <- apply(S1_genotypes==22, 1, sum, na.rm=T)
    SAAF1 <- (2*(S1n22) + S1n12)/((2*(S1n11) + S1n12)+(2*(S1n22) + S1n12))
    
    # sample alternative allele frequency - SAAF2
    S2_genotypes <- sampled_snp_genotypes[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes))]
    S2n11 <- apply(S2_genotypes==11, 1, sum, na.rm=T)
    S2n12 <- apply(S2_genotypes==12 | S2_genotypes==21, 1, sum, na.rm=T)
    S2n22 <- apply(S2_genotypes==22, 1, sum, na.rm=T)
    SAAF2 <- (2*(S2n22) + S2n12)/((2*(S2n11) + S2n12)+(2*(S2n22) + S2n12))
    
    # assemble LOCUS-SPECIFIC summary statistics
    locus_lss_info  <- sampled_snp[, which(grepl("^ID" , unique(names(sampled_snp)))):which(grepl("^GO" , unique(names(sampled_snp))))]
    locus_lss_stats <- egglib_summary_stats[egglib_summary_stats$ID %in% sampled_snp$ID, grepl("^LSS" , unique(names(egglib_summary_stats)))]
    locus_wss_stats <- egglib_summary_stats[egglib_summary_stats$ID %in% sampled_snp$ID , grepl("^WSS" , unique(names(egglib_summary_stats)))]
    
    # ASSEMBLY default LOCUS-SPECIFIC summary statistics
    locus_summary_stats <- cbind(locus_lss_info, SAAF1, SAAF2, locus_lss_stats, locus_wss_stats)
    #locus_summary_stats <- locus_summary_stats[ !duplicated(locus_summary_stats[ , 3]), ]
    colnames(locus_summary_stats) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats)[2]))
    
    if (add_WSSwspan_SFSbins_1){
      locus_wss_stats_add1 <- egglib_summary_stats_add1[egglib_summary_stats_add1$ID %in% sampled_snp$ID , grepl("^WSS" , unique(names(egglib_summary_stats_add1)))]
      colnames(locus_wss_stats_add1) <- paste0("ADD_1_LSS","_",seq(from=1,to=dim(locus_wss_stats_add1)[2]))
    }
    
    if (add_WSSwspan_SFSbins_2){
      locus_wss_stats_add2 <- egglib_summary_stats_add2[egglib_summary_stats_add2$ID %in% sampled_snp$ID , grepl("^WSS" , unique(names(egglib_summary_stats_add2)))]
      colnames(locus_wss_stats_add2) <- paste0("ADD_2_LSS","_",seq(from=1,to=dim(locus_wss_stats_add2)[2]))
    }
    
    if (exists("locus_wss_stats_add1") & exists("locus_wss_stats_add2")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1, locus_wss_stats_add2)
    } else if (exists("locus_wss_stats_add1")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1)
    } else if (exists("locus_wss_stats_add2")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add2)
    }
  
  } else {
    actual_sample_prbe <- as.numeric(NA)
    strong_sample_prbe <- as.numeric(NA)
    global_stats <- as.data.frame(t(rep(NA, 20)))
    global_SFS <- as.data.frame(t(rep(NA, sfs_bins_run)))
    add_global_stats <- as.data.frame(t(rep(NA, 198)))
    global_summary_stats <- cbind(global_stats, global_SFS, add_global_stats)
    colnames(global_summary_stats) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats)[2]))
    
    if (add_WSSwspan_SFSbins_1){
      global_SFS_bins1 <- as.data.frame(t(rep(NA, add_sfs_bins_1)))
      add_global_stats_add1 <- cbind(global_SFS_bins1, as.data.frame(t(rep(NA, 156))))
      colnames(add_global_stats_add1) <- paste0("ADD_1_GSS","_",seq(from=1,to=dim(add_global_stats_add1)[2]))
    }
    
    if (add_WSSwspan_SFSbins_2){
      global_SFS_bins2 <- as.data.frame(t(rep(NA, add_sfs_bins_2)))
      add_global_stats_add2 <- cbind(global_SFS_bins2, as.data.frame(t(rep(NA, 156))))
      colnames(add_global_stats_add2) <- paste0("ADD_2_GSS","_",seq(from=1,to=dim(add_global_stats_add2)[2]))
    }
    
    if (add_WSSwspan_SFSbins_1 & add_WSSwspan_SFSbins_2){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add1, add_global_stats_add2)
    } else if (add_WSSwspan_SFSbins_1){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add1)
    } else if (add_WSSwspan_SFSbins_2){
      global_summary_stats <- cbind(global_summary_stats, add_global_stats_add2)
    }
    
    FSTfdr <- as.data.frame(cbind(FSTfdrNS05 = NA,
                                  FSTfdrNS10 = NA,
                                  FSTfdrNS25 = NA))
    
    FSTprec <- as.data.frame(cbind(FSTprecNS75 = NA,
                                   FSTprecNS95 = NA))
    
    locus_summary_stats <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats)[2]))
    
    if (add_WSSwspan_SFSbins_1){
      locus_wss_stats_add1 <- as.data.frame(t(rep(NA, 26)))
      colnames(locus_wss_stats_add1) <- paste0("ADD_1_LSS","_",seq(from=1,to=dim(locus_wss_stats_add1)[2]))
    }
    
    if (add_WSSwspan_SFSbins_2){
      locus_wss_stats_add2 <- as.data.frame(t(rep(NA, 26)))
      colnames(locus_wss_stats_add2) <- paste0("ADD_2_LSS","_",seq(from=1,to=dim(locus_wss_stats_add2)[2]))
    }
    
    if (exists("locus_wss_stats_add1") & exists("locus_wss_stats_add2")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1, locus_wss_stats_add2)
    } else if (exists("locus_wss_stats_add1")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1)
    } else if (exists("locus_wss_stats_add2")){
      locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add2)
    }
  }
  
  ## RAW REFERENCE TABLE
  ##------------------------------------------------------------------------------------------------
  
  raw_reftable  <- suppressWarnings(cbind(as.factor(model_type), 
                                          sim_seed, sim, 
                                          mu, rr, selfing, Neq, N,
                                          gammaM, gammak, tc,
                                          PrGWSel, prbe, 
                                          averageGenLoad, lastGenLoad,
                                          actual_pop_prbe, strong_pop_prbe,
                                          actual_sample_prbe, strong_sample_prbe,
                                          pedigreeNetimes, pedigreeNetotal,
                                          IBDNeGWtimes, IBDNeGWtotal, IBDNeGWNtimes, IBDNeGWNtotal, 
                                          IBDNeGWStimes, IBDNeGWStotal, IBDNeChrtimes, IBDNeChrtotal,
                                          VARNeGWtimes, VARNeGWtotal, VARNeGWNtimes, VARNeGWNtotal, 
                                          VARNeGWStimes, VARNeGWStotal, VARNeChrtimes, VARNeChrtotal,
                                          FSTfdr, FSTprec, locus_summary_stats, global_summary_stats))
  rownames(raw_reftable) <- sim
    
  # remove egglib file from egglib_input and egglib_output folders 
  if (remove_files){
    file.remove(paste0(egglib_input_folder, egglib_converted_file))
    file.remove(paste0(egglib_output_folder,"egglib_output_sample", "_", sim, ".txt"))
    file.remove(paste0(egglib_output_folder,"egglib_output_sample_addwspan1", "_", sim, ".txt"))
    file.remove(paste0(egglib_output_folder,"egglib_output_sample_addwspan2", "_", sim, ".txt"))
  }
  
  ## OUTPUT RAW REFERENCE TABLE
  return(raw_reftable)
  
}
