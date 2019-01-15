## MODEL SELECTION
##------------------------------------------------------------------------------------------------

model_title <- switch(model_type, "DN", "BS", "SV")

## FIXED VALUES DEFINING GENOMIC STRUCTURE
##------------------------------------------------------------------------------------------------

# DEFINE CHROMOSOME SIZE AND RECOMBINATION LIMITS
if (chrN == 1){
  rr_limits = c((genomeS-1))
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

## DEFINE RADSEQ LOCI - IF IT IS A RADseq DATA
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
  
  # remove list of RADseq loci starting position after use it 
  rm(rad_starts)
}

## ADDITIONAL FUNCTIONS
##---------------------------------------------------------------------------------------------

# function to re-code the chromosome name - use it with apply()
chromtagging <- function(x, chrs_lowers){
  
  for (i in seq(from = 1, to = (length(chrs_lowers)-1))){
    
    if (x < chrs_lowers[1]){ chrom_idd = paste0("chr", 1)}
    
    else if (x > chrs_lowers[length(chrs_lowers)]){ chrom_idd = paste0("chr", (length(chrs_lowers)+1))}
    
    else if (x > chrs_lowers[i] & x < chrs_lowers[i+1]){ chrom_idd = paste0("chr", (i+1))}
    
  }
  return(chrom_idd)
}

# function to prepare the lines of one mutation for the wfabc input format
countgen4wfabc <- function(input, t_points = 2){
  
  d = input
  g = vector("list", length = t_points)
  n_chrom = NULL
  n_Aalleles = NULL
  for (t in seq(t_points)){
    g[[t]] <- d[grepl(paste0("@pop", t), names(d))]
    n_chromT <- 2 * sum(g[[t]] != 00)
    n_AallelesT <- sum(g[[t]] == 12 | g[[t]] == 21) + 2 * sum(g[[t]] == 22)
    
    n_chrom <- cbind(n_chrom, n_chromT)
    n_Aalleles <- cbind(n_Aalleles, n_AallelesT)
  }
  return(rbind(n_chrom, n_Aalleles))
}

## SIMULATIONS
##---------------------------------------------------------------------------------------------
do_sim <- function(sim, nsim, 
                   path_to_slim_model, slim_model_prefix,
                   path_to_slim, slim_output_folder,
                   path_to_bgzip, path_to_tabix, path_to_bcftools,
                   egglib_input_folder, egglib_output_folder,
                   path_to_python, path_to_egglib_summstat,
                   remove_files, debug_sim, debug_output_folder,
                   wfabc_input_file, wfabc_input_folder,
                   model_type, model_title, genomeS, fragS, chrN, chrS,
                   chrTAG, chromtagging, data_type, 
                   radseq_readL, radseq_readN, radseq_cov,
                   missing_data, haplotype, SS1, SS2, tau,
                   mu_rate, mu_random, mu_min, mu_max, 
                   neq_value, neq_random, neq_min, neq_max,
                   ncs_value, ncs_random, ncs_min, ncs_max,  
                   gammaM_value, gammak_value, gammaM_gammak, 
                   gammaM_random, gammaM_min, gammaM_max, 
                   gammak_random, gammak_min, gammak_max, 
                   PrGWSel_value, PrGWSel_random, PrGWSel_min, PrGWSel_max, 
                   prbe_value, prbe_random, prbe_min, prbe_max, 
                   domN, domN_random, domN_min, domN_max, 
                   domB, domB_random, domB_min, domB_max, 
                   rr_rate, rr_random, rr_min, rr_max, rr_limits,
                   selfing_rate, selfing_random, selfing_min, selfing_max,
                   tc_random, tc_value, wss_wspan_run, sfs_bins_run,
                   add_WSSwspan_SFSbins_1, add_wss_wspan_1, add_sfs_bins_1,
                   add_WSSwspan_SFSbins_2, add_wss_wspan_2, add_sfs_bins_2
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
  
  # POPULATION SIZE EQUILIBRIUM PHASE - Neq
  if (neq_random){
    Neq <- as.integer(10^runif(n = 1, min = log10(neq_min), max = log10(neq_max)))  
  } else {
    Neq <- neq_value
  }
  
  # POPULATION CENSUS SIZE - Ncs
  if (ncs_random){
    Ncs <- as.integer(10^runif(n = 1, min = log10(ncs_min), max = log10(ncs_max)))  
  } else {
    Ncs <- ncs_value
  }
  
  # GENOME-WIDE DFE FOR BENEFICIAL MUTATIONS
  # gamma mean
  if (gammaM_random){
    gammaM <- 10^runif(n = 1, min = log10(gammaM_min), max = log10(gammaM_max))
  } else {
    gammaM <- gammaM_value
  }
  
  # gamma shape 
  if (gammak_random){
    if (gammaM_gammak){
       gammak = gammaM
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
  
  # PROPORTION OF BENEFICIAL MUTATION IN G2 ELEMENTS
  if (prbe_random){
    prbe <- 10^runif(1, min = log10(prbe_min), max = log10(prbe_max))
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
  
  # Beneficial mutations (m3)
  if (domB_random){
    dm3 <- runif(n = 1, min = domB_min, max = domB_max)
  } else {
    dm3 <- domB
  }
  
  # PER BASE RECOMBINATION RATE
  if (rr_random){
    rr <- 10^runif(1, min = log10(rr_min), max = log10(rr_max))
  } else {
    rr <- rr_rate
  }
  
  # the distribution of rr across chromosome limits
  if (chrN == 1){
    rr_rates = c(rr)
  } else {
    rr_rates = rep(c(rr, 0.5), as.integer(length(rr_limits)/2))

    # updated rr_limits and rr_rates
    rr_limits = c(rr_limits, (genomeS-1))
    rr_rates = c(rr_rates, rr)
  }	
  
  # SELFING RATE
  if (selfing_random){
    selfing <- (10^runif(1, min = log10(selfing_min), max = log10(selfing_max)))
  } else {
    selfing <- selfing_rate
  }
  
  # SELECTION ON STANDING VARIATION
  if (model_type == 3){
    if(tc_random){
      tc <- as.integer(runif(1, min = 0, max = tau))
    } else {
      tc <- tc_value
    }
  } else {
    tc <- 0
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
  slim_Ncs       <- paste0("-d Ncs=", Ncs)
  slim_gammaM    <- paste0("-d gammaM=", gammaM)          # p1 and p2 - DN and BS; p2 - SV
  slim_gammak    <- paste0("-d gammak=", gammak)          # p1 and p2 - DN and BS; p2 - SV
  slim_prbe      <- paste0("-d prbe=", prbe)
  slim_dm1       <- paste0("-d dm1=", dm1)
  slim_dm2       <- paste0("-d dm2=", dm2)
  slim_dm3       <- paste0("-d dm3=", dm3)
  slim_rr        <- paste0("-d rr=", rr)
  slim_selfing   <- paste0("-d selfing=", selfing)        # p1 only
  slim_tc        <- paste0("-d tc=", tc)                  # p2 - SV
  slim_genomeS   <- paste0("-d genomeS=", as.integer(genomeS))
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
                       slim_rr,
                       slim_selfing,                  # p1 only 
                       slim_genomeS,
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
  
  ## ADD conditional here to check if .tree exists
  if (file.exists(paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"))){
    
    slim_run_p2 <- paste(path_to_slim,
                         "-s", sim_seed,                # seed  = simulation seed number
                         slim_simID,                    # simID = simulation id number
                         slim_mu,                       
                         slim_Ncs,
                         slim_gammaM,                   # p1 and p2 - DN and BS; p2 - SV
                         slim_gammak,                   # p1 and p2 - DN and BS; p2 - SV
                         slim_prbe,
                         slim_dm1,
                         slim_dm2,
                         slim_dm3,
                         slim_rr,
                         slim_tc,                       # p2 - SV
                         slim_genomeS,
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
    
    ## HANDLING SLiM OUTPUT PEDIGREE NE EQUILIBRIUM PHASE
    ##------------------------------------------------------------------------------
    
    # load the data
    slim_output_ne1 <- paste0(slim_output_folder,"slim_output_ne1_", sim, ".txt")
    
    if (file.exists(slim_output_ne1)){
      info_ne1_file = file.info(slim_output_ne1)
      
      if (info_ne1_file$size != 0){
        ne1 <- read.csv(file = slim_output_ne1, sep = "\t", header = F, col.names = c("gen", "ne"), na.strings = "NA")
        
        if (nrow(ne1) != 0){
          
          # check if there is any -Inf/Inf and susbstitute it for NA
          ne1[which(is.infinite(ne1$ne)), "ne"] <- NA
          
          # get the Pedigree Ne for the whole period - harmonic mean
          meanNe1 <- 1/mean(1/ne1$ne, na.rm = TRUE) # all NA's pedigreeNetotal = NA;
          
        } else {
          meanNe1 <- as.numeric(NA)
        }
        
        # remove pedigreeNe data after use it
        rm(ne1)
        
      } else {
        meanNe1 <- as.numeric(NA)
      }
      
    } else {
      
      meanNe1 <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        debug_message <- "no slim_output_ne1 file found"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## HANDLING SLiM2 OUTPUT t1:t2 - PEDIGREE NE SAMPLING PHASE
    ##------------------------------------------------------------------------------
    
    # load the data
    slim_output_ne2 <- paste0(slim_output_folder,"slim_output_ne2_", sim, ".txt")
    
    if (file.exists(slim_output_ne2)){
      info_ne2_file = file.info(slim_output_ne2)
      
      if (info_ne2_file$size != 0){
        ne2 <- read.csv(file = slim_output_ne2, sep = "\t", header = F, col.names = c("gen", "time", "ne"), na.strings = "NA")
        
        if (nrow(ne2) != 0){
          
          # check if there is any -Inf/Inf and susbstitute it for NA
          ne2[which(is.infinite(ne2$ne)), "ne"] <- NA
          
          # take the Pedigree Ne for each interval
          timesNe2 <- as.data.frame(t(ne2$ne))
          names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau)))
          
          # get the Pedigree Ne for the whole period - harmonic mean
          meanNe2 <- 1/mean(1/ne2$ne, na.rm = TRUE) # all NA's pedigreeNetotal = NA;
          
        } else {
          timesNe2 <- as.data.frame(t(rep(NA, tau+1)))
          names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau)))
          meanNe2 <- as.numeric(NA)
        }
        
        # remove pedigreeNe data after use it
        rm(ne2)
        
      } else {
        timesNe2 <- as.data.frame(t(rep(NA, tau+1)))
        names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau)))
        meanNe2 <- as.numeric(NA)
      }
      
    } else {
      timesNe2 <- as.data.frame(t(rep(NA, tau+1)))
      names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau)))
      meanNe2 <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        debug_message <- "no slim_output_ne2 file found"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
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
        
        # remove genetic load data after use it
        rm(geneticLoad)
        
      } else {
        averageGenLoad <- as.numeric(NA)
        lastGenLoad <- as.numeric(NA)
      }
      
    } else {
      
      averageGenLoad <- as.numeric(NA)
      lastGenLoad <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        debug_message <- "no slim_output_load file found"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## HANDLING SLiM2 OUTPUT t1:t2 - POPULATION MUTATIONS DATA
    ##-------------------------------------------------------------------------------
    
    filenames_genome <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=0, to=tau, by=1), "_", sim, ".txt")
      
    if(all(file.exists(c(filenames_genome)))){
      
      ### GENOME
      datalist_genome <- lapply(filenames_genome, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome)
        
        if(!is.null(all_merged_genome)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe    <- sum(all_merged_genome[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
            actual_pop_SelMean <- mean(all_merged_genome[all_merged_genome[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd   <- sd(all_merged_genome[all_merged_genome[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2)){
              
              if (any(meanNe2*all_merged_genome[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe    <- sum(meanNe2*all_merged_genome[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
                strong_positive_SelMean <- mean(all_merged_genome[meanNe2*all_merged_genome[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd   <- sd(all_merged_genome[meanNe2*all_merged_genome[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe    <- as.numeric(NA)
                strong_positive_SelMean <- as.numeric(NA)
                strong_positive_SelSd   <- as.numeric(NA)
              }
              
              if (any(meanNe2*all_merged_genome[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe    <- sum(meanNe2*all_merged_genome[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome[, paste0("S",tc)][!is.na(all_merged_genome[, paste0("S",tc)])])
                strong_negative_SelMean <- mean(all_merged_genome[meanNe2*all_merged_genome[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd   <- sd(all_merged_genome[meanNe2*all_merged_genome[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe    <- as.numeric(NA)
                strong_negative_SelMean <- as.numeric(NA)
                strong_negative_SelSd   <- as.numeric(NA)
              }
              
              if(is.na(strong_positive_prbe)  & is.na(strong_negative_prbe)){
                
                strong_pop_prbe    <- as.numeric(NA)
                strong_pop_SelMean <- as.numeric(NA)
                strong_pop_SelSd   <- as.numeric(NA)
                 
              } else {
                strong_pop_prbe    <- sum(strong_positive_prbe, strong_negative_prbe, na.rm = TRUE)
                strong_pop_SelMean <- sum(strong_positive_SelMean, strong_negative_SelMean, na.rm = TRUE)
                strong_pop_SelSd   <- sum(strong_positive_SelSd, strong_negative_SelSd, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe    <- as.numeric(NA)
              strong_pop_SelMean <- as.numeric(NA)
              strong_pop_SelSd   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe    <- as.numeric(NA)
            actual_pop_SelMean <- as.numeric(NA)
            actual_pop_SelSd   <- as.numeric(NA)
            
            strong_pop_prbe    <- as.numeric(NA)
            strong_pop_SelMean <- as.numeric(NA)
            strong_pop_SelSd   <- as.numeric(NA)
          }
          
        } else {
          actual_pop_prbe    <- as.numeric(NA)
          actual_pop_SelMean <- as.numeric(NA)
          actual_pop_SelSd   <- as.numeric(NA)
          
          strong_pop_prbe    <- as.numeric(NA)
          strong_pop_SelMean <- as.numeric(NA)
          strong_pop_SelSd   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome)
        
        if(!is.null(all_merged_genome)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe    <- sum(all_merged_genome[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
            actual_pop_SelMean <- mean(all_merged_genome[all_merged_genome[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd   <- sd(all_merged_genome[all_merged_genome[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2)){
              if (any(meanNe2*all_merged_genome[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe    <- sum(meanNe2*all_merged_genome[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
                strong_positive_SelMean <- mean(all_merged_genome[meanNe2*all_merged_genome[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd   <- sd(all_merged_genome[meanNe2*all_merged_genome[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe    <- as.numeric(NA)
                strong_positive_SelMean <- as.numeric(NA)
                strong_positive_SelSd   <- as.numeric(NA)
              }
              
              if (any(meanNe2*all_merged_genome[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe    <- sum(meanNe2*all_merged_genome[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome[, "S"][!is.na(all_merged_genome[, "S"])])
                strong_negative_SelMean <- mean(all_merged_genome[meanNe2*all_merged_genome[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd   <- sd(all_merged_genome[meanNe2*all_merged_genome[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe    <- as.numeric(NA)
                strong_negative_SelMean <- as.numeric(NA)
                strong_negative_SelSd   <- as.numeric(NA)
              }
              
              if(is.na(strong_positive_prbe)  & is.na(strong_negative_prbe)){
                
                strong_pop_prbe    <- as.numeric(NA)
                strong_pop_SelMean <- as.numeric(NA)
                strong_pop_SelSd   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe    <- sum(strong_positive_prbe, strong_negative_prbe, na.rm = TRUE)
                strong_pop_SelMean <- sum(strong_positive_SelMean, strong_negative_SelMean, na.rm = TRUE)
                strong_pop_SelSd   <- sum(strong_positive_SelSd, strong_negative_SelSd, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe    <- as.numeric(NA)
              strong_pop_SelMean <- as.numeric(NA)
              strong_pop_SelSd   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe    <- as.numeric(NA)
            actual_pop_SelMean <- as.numeric(NA)
            actual_pop_SelSd   <- as.numeric(NA)
            
            strong_pop_prbe    <- as.numeric(NA)
            strong_pop_SelMean <- as.numeric(NA)
            strong_pop_SelSd   <- as.numeric(NA)
          }
          
        } else {
          actual_pop_prbe    <- as.numeric(NA)
          actual_pop_SelMean <- as.numeric(NA)
          actual_pop_SelSd   <- as.numeric(NA)
          
          strong_pop_prbe    <- as.numeric(NA)
          strong_pop_SelMean <- as.numeric(NA)
          strong_pop_SelSd   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe    <- as.numeric(NA)
      actual_pop_SelMean <- as.numeric(NA)
      actual_pop_SelSd   <- as.numeric(NA)
      
      strong_pop_prbe    <- as.numeric(NA)
      strong_pop_SelMean <- as.numeric(NA)
      strong_pop_SelSd   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        file.copy(from = filenames_genome, to = debug_output_folder)
        
        debug_message <- "no pmuts file found"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## HANDLINDING SLiM2 OUTPUT - SAMPLED INDIVIDUALS
    ##---------------------------------------------------------------------------------
    
    # sort vcf files
    slim_output_sample_t1        <- paste0(slim_output_folder,"slim_output_sample_t1_", sim, ".vcf")
    slim_output_sample_t2        <- paste0(slim_output_folder,"slim_output_sample_t2_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_t1, slim_output_sample_t1)))){
      
      slim_output_sample_t1_sorted <- paste0(slim_output_folder,"slim_output_sample_t1_", sim, "_sorted" , ".vcf")
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
                              ">", slim_output_sample_merged) 
      
      system(bcftools_query)
      
      if(file.exists(slim_output_sample_merged)){
        
        # assembly the header
        header_1       <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1 <- paste0("indiv", seq(from=1, to=SS1, by=1), "@pop1", "")
        sample_names_2 <- paste0("indiv", seq(from=1, to=SS2, by=1), "@pop2", "")
        
        full_header <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data <- read.table(file = slim_output_sample_merged, header = F, col.names = full_header, check.names = F, na.strings = "./.")
        
        # if it is a RADseq data
        if (data_type == 2){
          slim_raw_data <- slim_raw_data[which(slim_raw_data$position %in% rad_interval), ]
        }
        
        if (nrow(slim_raw_data) != 0){
          
          # split the data
          slim_snp_geno <- slim_raw_data[, 9:ncol(slim_raw_data)]
          slim_snp_info <- slim_raw_data[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data)
          
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
          
          # remove raw snp data information after use it
          rm(slim_snp_geno)
          rm(slim_snp_info)
          
          # remove monomorphic mutations
          slim_data <- slim_data[keep_snps, ]
          
          # remove duplicated mutations
          slim_data <- slim_data[!duplicated(slim_data[ ,1:2]), ]
          
          # make WFABC input file
          if (wfabc_input_file){
            
            slim2wfabc <- slim_data[, -c(1:8)]
            slim2wfabc <- as.data.frame(t(slim2wfabc))
            
            wfabc_data  <- do.call(rbind, sapply(slim2wfabc, countgen4wfabc, t_points=2, simplify = F))
            
            if (!file_test("-d", wfabc_input_folder)){
              dir.create(file.path(wfabc_input_folder))
            }
            
            write(paste(dim(slim2wfabc)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_", sim,".txt")) 
            write(paste(0, (tau), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_", sim,".txt"), append = TRUE) 
            write.table(wfabc_data, file=paste0(wfabc_input_folder, "wfabc_input_sample_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE) 
          }
          
          # re-code the chromosome name
          if(chrN > 1){
            if (chrTAG){
              slim_data$chrom <- sapply(slim_data$position, chromtagging, chrs_lowers=chrs_lowers)
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
          egglib_converted_file <- paste0("egglib_input_sample", "_", sim, ".txt")
          write.table(slim_to_egglib_data, file = paste0(egglib_input_folder,egglib_converted_file), quote=FALSE, sep="\t", row.names = FALSE)
          
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
          
          # remove snp datasets after use it
          rm(slim_data)
          rm(slim_to_egglib_data)
          
          ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
          ##-----------------------------------------------------------------------------------
          
          if(file.exists(paste0(egglib_input_folder,egglib_converted_file))){
            
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
            
            # import egglib output
            egglib_output_summstats <- paste0(egglib_output_folder,"egglib_output_sample", "_", sim, ".txt") 
            
            if(file.exists(egglib_output_summstats)){
              
              egglib_summary_stats <- read.csv(file = paste0(egglib_output_folder,"egglib_output_sample", "_", sim, ".txt"), header = T, sep = "\t", check.names = F)
              
              ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
              if(any(slim_to_egglib_snps$S != 0, na.rm = TRUE)){
                actual_sample_prbe    <- sum(slim_to_egglib_snps$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
                actual_sample_SelMean <- mean(slim_to_egglib_snps[slim_to_egglib_snps[, "S"] != 0, "S"], na.rm = TRUE)
                actual_sample_SelSd   <- sd(slim_to_egglib_snps[slim_to_egglib_snps[, "S"] != 0, "S"], na.rm = TRUE)
                
              } else {
                actual_sample_prbe    <- as.numeric(NA)
                actual_sample_SelMean <- as.numeric(NA)
                actual_sample_SelSd   <- as.numeric(NA)
              }
              
             if (!is.na(meanNe2)){
                ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                if(any(meanNe2*slim_to_egglib_snps$S > 1, na.rm = TRUE)){
                  positive_sample_prbe    <- sum(meanNe2*slim_to_egglib_snps$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
                  positive_sample_SelMean <- mean(slim_to_egglib_snps[meanNe2*slim_to_egglib_snps[, "S"] > 1, "S"], na.rm = TRUE)
                  positive_sample_SelSd   <- sd(slim_to_egglib_snps[meanNe2*slim_to_egglib_snps[, "S"] > 1, "S"], na.rm = TRUE)
                  
                } else {
                  positive_sample_prbe    <- as.numeric(NA)
                  positive_sample_SelMean <- as.numeric(NA)
                  positive_sample_SelSd   <- as.numeric(NA)
                }
                
                ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                if(any(meanNe2*slim_to_egglib_snps$S < -1, na.rm = TRUE)){
                  negative_sample_prbe    <- sum(meanNe2*slim_to_egglib_snps$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps$S[!is.na(slim_to_egglib_snps$S)])
                  negative_sample_SelMean <- mean(slim_to_egglib_snps[meanNe2*slim_to_egglib_snps[, "S"] < -1, "S"], na.rm = TRUE)
                  negative_sample_SelSd   <- sd(slim_to_egglib_snps[meanNe2*slim_to_egglib_snps[, "S"] < -1, "S"], na.rm = TRUE)
                  
                } else {
                  negative_sample_prbe    <- as.numeric(NA)
                  negative_sample_SelMean <- as.numeric(NA)
                  negative_sample_SelSd   <- as.numeric(NA)
                }
                
               if(is.na(positive_sample_prbe)  & is.na(negative_sample_prbe)){
                 
                 strong_sample_prbe    <- as.numeric(NA)
                 strong_sample_SelMean <- as.numeric(NA)
                 strong_sample_SelSd   <- as.numeric(NA)
                 
               } else {
                 strong_sample_prbe    <- sum(positive_sample_prbe, negative_sample_prbe, na.rm = TRUE)
                 strong_sample_SelMean <- sum(positive_sample_SelMean, negative_sample_SelMean, na.rm = TRUE)
                 strong_sample_SelSd   <- sum(positive_sample_SelSd, negative_sample_SelSd, na.rm = TRUE)
               }
               
              } else {
                strong_sample_prbe    <- as.numeric(NA)
                strong_sample_SelMean <- as.numeric(NA)
                strong_sample_SelSd   <- as.numeric(NA)
              }
              
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
              
              locusFST_test_table <- data.frame(Ns=meanNe2*slim_to_egglib_snps$S, LSS_WCst=egglib_summary_stats[, "LSS_WCst"], Ns_test=ifelse(meanNe2*slim_to_egglib_snps$S > 1, 1, 0))
        
              locusFST_test_table <- locusFST_test_table[order(-locusFST_test_table$LSS_WCst), ]

                            
              if (any(locusFST_test_table$Ns > 1)){
                if(!all(locusFST_test_table$Ns > 1)){
                  
                  pred_locusFSTNS <- prediction(predictions = locusFST_test_table$LSS_WCst, labels = locusFST_test_table$Ns_test)
                  perf_locusFSTNS <- performance(pred_locusFSTNS, "ppv", "fpr")
                  
                  perf_table <- data.frame(ppvNS=perf_locusFSTNS@y.values[[1]], fdrNS=1-perf_locusFSTNS@y.values[[1]])
                  
                  perf_locusFST_table <- data.frame(locusFST_test_table[1:dim(perf_table)[1], ], perf_table)
                  
                  FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = min(perf_locusFST_table[which(perf_locusFST_table$fdrNS <= 0.005),"LSS_WCst"]),
                                                FSTfdrNS01  = min(perf_locusFST_table[which(perf_locusFST_table$fdrNS <= 0.01) ,"LSS_WCst"]),
                                                FSTfdrNS02  = min(perf_locusFST_table[which(perf_locusFST_table$fdrNS <= 0.02) ,"LSS_WCst"]),
                                                FSTfdrNS05  = min(perf_locusFST_table[which(perf_locusFST_table$fdrNS <= 0.05) ,"LSS_WCst"]),
                                                FSTfdrNS10  = min(perf_locusFST_table[which(perf_locusFST_table$fdrNS <= 0.10) ,"LSS_WCst"])))
                  
                } else {
                  
                  # all strongly selected mutations
                  FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = min(locusFST_test_table$LSS_WCst),
                                                FSTfdrNS01  = min(locusFST_test_table$LSS_WCst),
                                                FSTfdrNS02  = min(locusFST_test_table$LSS_WCst),
                                                FSTfdrNS05  = min(locusFST_test_table$LSS_WCst),
                                                FSTfdrNS10  = min(locusFST_test_table$LSS_WCst)))
                
                }
              } else {
                
                # all neutral mutations
                FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = max(locusFST_test_table$LSS_WCst),
                                              FSTfdrNS01 = max(locusFST_test_table$LSS_WCst),
                                              FSTfdrNS02 = max(locusFST_test_table$LSS_WCst),
                                              FSTfdrNS05 = max(locusFST_test_table$LSS_WCst),
                                              FSTfdrNS10 = max(locusFST_test_table$LSS_WCst)))
                
              }
              
              ## LOCUS-SPECIFIC SUMMARY STATISTICS
              # sampling ONE RANDOM mutation for the locus-specific reference table
              snps_in_reftable <- sample(which(slim_to_egglib_snps$MT == 1 | slim_to_egglib_snps$MT == 2 | slim_to_egglib_snps$MT == 3), size=1)
              sampled_snp <- slim_to_egglib_snps[snps_in_reftable, ]
              
              # remove complete snp table after use it
              rm(slim_to_egglib_snps)
              
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
              
              # remove summary statistics data after use it
              rm(egglib_summary_stats)
              
              # ASSEMBLY default LOCUS-SPECIFIC summary statistics
              locus_summary_stats <- cbind(locus_lss_info, SAAF1, SAAF2, locus_lss_stats, locus_wss_stats)
              #locus_summary_stats <- locus_summary_stats[ !duplicated(locus_summary_stats[ , 3]), ]
              colnames(locus_summary_stats) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats)[2]))
              
              if (add_WSSwspan_SFSbins_1){
                locus_wss_stats_add1 <- egglib_summary_stats_add1[egglib_summary_stats_add1$ID %in% sampled_snp$ID , grepl("^WSS" , unique(names(egglib_summary_stats_add1)))]
                colnames(locus_wss_stats_add1) <- paste0("ADD_1_LSS","_",seq(from=1,to=dim(locus_wss_stats_add1)[2]))
                
                # remove additional summary statistics 1 data after use it
                rm(egglib_summary_stats_add1)
              }
              
              if (add_WSSwspan_SFSbins_2){
                locus_wss_stats_add2 <- egglib_summary_stats_add2[egglib_summary_stats_add2$ID %in% sampled_snp$ID , grepl("^WSS" , unique(names(egglib_summary_stats_add2)))]
                colnames(locus_wss_stats_add2) <- paste0("ADD_2_LSS","_",seq(from=1,to=dim(locus_wss_stats_add2)[2]))
                
                # remove additional summary statistics 2 data after use it
                rm(egglib_summary_stats_add2)
              }
              
              if (exists("locus_wss_stats_add1") & exists("locus_wss_stats_add2")){
                locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1, locus_wss_stats_add2)
              } else if (exists("locus_wss_stats_add1")){
                locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add1)
              } else if (exists("locus_wss_stats_add2")){
                locus_summary_stats <- cbind(locus_summary_stats, locus_wss_stats_add2)
              }
              
            } else {
              
              actual_sample_prbe    <- as.numeric(NA)
              actual_sample_SelMean <- as.numeric(NA)
              actual_sample_SelSd   <- as.numeric(NA)
              
              strong_sample_prbe    <- as.numeric(NA)
              strong_sample_SelMean <- as.numeric(NA)
              strong_sample_SelSd   <- as.numeric(NA)
              
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
              
              FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                            FSTfdrNS01  = NA,
                                            FSTfdrNS02  = NA,
                                            FSTfdrNS05  = NA,
                                            FSTfdrNS10  = NA))
              
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
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                file.copy(from = c(slim_output_sample_t1, slim_output_sample_t2), to = debug_output_folder)
                file.copy(from = slim_output_sample_merged, to = debug_output_folder)
                file.copy(from = paste0(egglib_input_folder, egglib_converted_file), to = debug_output_folder)
                
                debug_message <- "egglib output file not found"
                
                debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                                      sim_seed, sim, 
                                                      mu, rr, selfing, Neq, Ncs,
                                                      gammaM, gammak, tc,
                                                      PrGWSel, prbe, debug_message))
                rownames(debug_dump) <- sim
                write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
              }
              
            }
            
          } else {
            
            actual_sample_prbe    <- as.numeric(NA)
            actual_sample_SelMean <- as.numeric(NA)
            actual_sample_SelSd   <- as.numeric(NA)
            
            strong_sample_prbe    <- as.numeric(NA)
            strong_sample_SelMean <- as.numeric(NA)
            strong_sample_SelSd   <- as.numeric(NA)
            
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
            
            FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                          FSTfdrNS01  = NA,
                                          FSTfdrNS02  = NA,
                                          FSTfdrNS05  = NA,
                                          FSTfdrNS10  = NA))
            
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
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              file.copy(from = c(slim_output_sample_t1, slim_output_sample_t2), to = debug_output_folder)
              file.copy(from = slim_output_sample_merged, to = debug_output_folder)
              
              debug_message <- "egglib input file not found"
              
              debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                                    sim_seed, sim, 
                                                    mu, rr, selfing, Neq, Ncs,
                                                    gammaM, gammak, tc,
                                                    PrGWSel, prbe, debug_message))
              rownames(debug_dump) <- sim
              write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
            }
          }
          
        } else {
          
          actual_sample_prbe    <- as.numeric(NA)
          actual_sample_SelMean <- as.numeric(NA)
          actual_sample_SelSd   <- as.numeric(NA)
          
          strong_sample_prbe    <- as.numeric(NA)
          strong_sample_SelMean <- as.numeric(NA)
          strong_sample_SelSd   <- as.numeric(NA)
          
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
          
          FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                        FSTfdrNS01  = NA,
                                        FSTfdrNS02  = NA,
                                        FSTfdrNS05  = NA,
                                        FSTfdrNS10  = NA))
          
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
        
      } else {
        
        actual_sample_prbe    <- as.numeric(NA)
        actual_sample_SelMean <- as.numeric(NA)
        actual_sample_SelSd   <- as.numeric(NA)
        
        strong_sample_prbe    <- as.numeric(NA)
        strong_sample_SelMean <- as.numeric(NA)
        strong_sample_SelSd   <- as.numeric(NA)
        
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
        
        FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                      FSTfdrNS01  = NA,
                                      FSTfdrNS02  = NA,
                                      FSTfdrNS05  = NA,
                                      FSTfdrNS10  = NA))
        
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
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          file.copy(from = c(slim_output_sample_t1, slim_output_sample_t2), to = debug_output_folder)
          
          debug_message <- "merged vcf files not found"
          
          debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                                sim_seed, sim, 
                                                mu, rr, selfing, Neq, Ncs,
                                                gammaM, gammak, tc,
                                                PrGWSel, prbe, debug_message))
          rownames(debug_dump) <- sim
          write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
        }
        
      }
      
    } else {
      
      actual_sample_prbe    <- as.numeric(NA)
      actual_sample_SelMean <- as.numeric(NA)
      actual_sample_SelSd   <- as.numeric(NA)
      
      strong_sample_prbe    <- as.numeric(NA)
      strong_sample_SelMean <- as.numeric(NA)
      strong_sample_SelSd   <- as.numeric(NA)
      
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
      
      FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                    FSTfdrNS01  = NA,
                                    FSTfdrNS02  = NA,
                                    FSTfdrNS05  = NA,
                                    FSTfdrNS10  = NA))
      
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
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        file.copy(from = c(slim_output_sample_t1, slim_output_sample_t2), to = debug_output_folder)
        
        debug_message <- "one or both vcf files were not found"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }

    if (remove_files){
      file.remove(paste0(slim_output_folder,"slim_coalesced_",model_title,"_", sim, ".tree"))
      file.remove(slim_output_geneticLoad)
      file.remove(slim_output_ne1)
      file.remove(slim_output_ne2)
      file.remove(filenames_genome)
      file.remove(paste0(slim_output_sample_t1))
      file.remove(paste0(slim_output_sample_t1_sorted, ".gz"))
      file.remove(paste0(slim_output_sample_t1_sorted, ".gz.tbi"))
      file.remove(paste0(slim_output_sample_t2))
      file.remove(paste0(slim_output_sample_t2_sorted, ".gz"))
      file.remove(paste0(slim_output_sample_t2_sorted, ".gz.tbi"))
      file.remove(paste0(slim_output_sample_merged))
      file.remove(paste0(egglib_input_folder,egglib_converted_file))
      file.remove(paste0(egglib_output_folder,"egglib_output_sample", "_", sim, ".txt"))
      if (add_WSSwspan_SFSbins_1){
        file.remove(paste0(egglib_output_folder,"egglib_output_sample_addwspan1", "_", sim, ".txt"))
      }
      if (add_WSSwspan_SFSbins_2){
        file.remove(paste0(egglib_output_folder,"egglib_output_sample_addwspan2", "_", sim, ".txt"))
      }
    }
    
  } else {
    
    if (debug_sim){
      
      # check if the folder exists
      if (!file_test("-d", debug_output_folder)){
        dir.create(file.path(debug_output_folder))
      }
      
      debug_message <- "no .tree file found"
      debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                            sim_seed, sim, 
                                            mu, rr, selfing, Neq, Ncs,
                                            gammaM, gammak, tc,
                                            PrGWSel, prbe, debug_message))
      rownames(debug_dump) <- sim
      write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
    }
    
    meanNe1 <- as.numeric(NA)
    meanNe2 <- as.numeric(NA)
    timesNe2 <- as.data.frame(t(rep(NA, tau+1)))
    names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau)))
    averageGenLoad <- as.numeric(NA)
    lastGenLoad <- as.numeric(NA)
    
    ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
    actual_pop_prbe    <- as.numeric(NA)
    actual_pop_SelMean <- as.numeric(NA)
    actual_pop_SelSd   <- as.numeric(NA)
    
    strong_pop_prbe    <- as.numeric(NA)
    strong_pop_SelMean <- as.numeric(NA)
    strong_pop_SelSd   <- as.numeric(NA)
    
    actual_sample_prbe    <- as.numeric(NA)
    actual_sample_SelMean <- as.numeric(NA)
    actual_sample_SelSd   <- as.numeric(NA)
    
    strong_sample_prbe    <- as.numeric(NA)
    strong_sample_SelMean <- as.numeric(NA)
    strong_sample_SelSd   <- as.numeric(NA)
    
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
    
    FSTfdr <- as.data.frame(cbind(FSTfdrNS005 = NA,
                                  FSTfdrNS01  = NA,
                                  FSTfdrNS02  = NA,
                                  FSTfdrNS05  = NA,
                                  FSTfdrNS10  = NA))
    
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
                                          mu, rr, selfing, Neq, Ncs,
                                          gammaM, gammak, tc,
                                          PrGWSel, prbe, meanNe1, meanNe2, timesNe2,
                                          averageGenLoad, lastGenLoad,
                                          actual_pop_prbe, actual_pop_SelMean, actual_pop_SelSd, 
                                          strong_pop_prbe, strong_pop_SelMean, strong_pop_SelSd,
                                          actual_sample_prbe, actual_sample_SelMean, actual_sample_SelSd,
                                          strong_sample_prbe, strong_sample_SelMean, strong_sample_SelSd,
                                          FSTfdr, locus_summary_stats, 
                                          global_summary_stats))
  rownames(raw_reftable) <- sim
  
  ## OUTPUT RAW REFERENCE TABLE
  return(raw_reftable)
  
}
