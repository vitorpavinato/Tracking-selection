## MODEL SELECTION
##------------------------------------------------------------------------------------------------

model_title <- switch(model_type, "DN", "BS", "SV")

## FIXED VALUES DEFINING GENOMIC STRUCTURE
##------------------------------------------------------------------------------------------------

# DEFINE CHROMOSOME SIZE AND RECOMBINATION LIMITS
if (chrN == 1){
  rr_limits = c((genomeS-1))
} else {
  chrS = round(genomeS/chrN)
  
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
  
  # set the total number of RADseq loci expected
  radseq_lociN = round((genomeS/radseq_readL))
  
  radseq_starts = c(1, 1:(radseq_lociN-1) * radseq_readL + 1)
  
  # Sample a number of RADseq loci STARTs based on the sequencing effort
  radseq_sampled <- sort(sample(radseq_starts, round(radseq_cov*radseq_lociN), replace = FALSE))
  
  # Generate the vector of positionS for the sampled RADseq reads
  radseq_interval = NULL
  for(i in seq_along(radseq_sampled)){
    interval   = radseq_sampled[i]:(radseq_sampled[i]+(radseq_readL-1))
    radseq_interval = c(radseq_interval, interval)
  }
  
  # remove list of RADseq loci starting position after use it 
  rm(radseq_starts)
}

## ADDITIONAL FUNCTIONS
##---------------------------------------------------------------------------------------------

# function to re-code the chromosome name - use it with apply()
chromtagging <- function(x, chrsLimits){
  
  for (i in seq(from = 1, to = (length(chrsLimits)-1))){
    
    if (x < chrsLimits[1]){ chrom_idd = paste0("chr", 1)}
    
    else if (x > chrsLimits[length(chrsLimits)]){ chrom_idd = paste0("chr", (length(chrsLimits)+1))}
    
    else if (x > chrsLimits[i] & x < chrsLimits[i+1]){ chrom_idd = paste0("chr", (i+1))}
    
  }
  return(chrom_idd)
}

# function to re-code the RADseq loci name - use it with apply()
radseqtagging <- function(x, tagSampled, readLength){
  
  radtag_idd = NULL
  for (i in seq_along(tagSampled)){
    
    if (  (x >= tagSampled[i] & x < tagSampled[i] + (readLength-1)) | (x == tagSampled[i] + (readLength-1)) ){ radtag_idd = paste0("rad", i)}
  }
  return(radtag_idd)
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
                   wfabc_input_file, wfabc_input_folder, egglib_input_selcoeff,               # add egglib_input_selcoeff
                   model_type, model_title, genomeS, fragS, 
                   chrN, chrTAG, chromtagging, chrs_lowers,                                   # removeer chrS
                   data_type, radseq_readL, radseq_interval, radseq_sampled, one_snp_radseq,  # add chrs_lowers  #radseq_readN to radseq_lociN ## remove radseq_readN ## remove radseq_cov
                   missing_data, haplotype, SSs, tau,                                         # add radseq_sampled and radseq_interval and one_snp_radseq; SS1 and SS2 to SSs
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
                   tc_random, tc_value, wss_wspan_run, sfs_bins_run
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
    Neq <- round(10^runif(n = 1, min = log10(neq_min), max = log10(neq_max)))  
  } else {
    Neq <- neq_value
  }
  
  # POPULATION CENSUS SIZE - Ncs
  if (ncs_random){
    Ncs <- round(10^runif(n = 1, min = log10(ncs_min), max = log10(ncs_max)))  
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
    rr_rates = rep(c(rr, 0.5), round(length(rr_limits)/2))

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
      tc <- round(runif(1, min = 0, max = tau))
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
  genomicElementN = round((genomeS/fragS)) # round
  
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
  
  g2_idx = sort(sample(indexes, round(PrGWSel*genomicElementN), replace = FALSE)) #round
  
  # The difference are the GenomicElementType G1
  g1_idx = setdiff(indexes, g2_idx)
  
  ## RUNNING SLiM
  ##-------------------------------------------------------------------------------
  
  # check if the folder exists
  if (!file_test("-d", slim_output_folder)){
    dir.create(file.path(slim_output_folder))
  }
  
  # generate text with slim command
  slim_model_p1 <- paste0(path_to_slim_model,slim_model_prefix,"_", model_title,"_COAL_bees",".slim")
  slim_model_p2 <- paste0(path_to_slim_model,slim_model_prefix,"_", model_title,"_posCOAL_bees",".slim")
  
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
  slim_genomeS   <- paste0("-d genomeS=", genomeS) # remove as.integer
  slim_rr_rates  <- paste0("-d rr_rates=", "'c(", paste(rr_rates, collapse = ","), ")'")
  slim_rr_limits <- paste0("-d rr_limits=", "'c(", paste(rr_limits, collapse = ","), ")'")
  slim_ge_starts <- paste0("-d e_starts=", "'c(", paste(e_starts, collapse = ","), ")'")
  slim_ge_ends   <- paste0("-d e_ends=", "'c(", paste(e_ends, collapse = ","), ")'")
  slim_g2_idx    <- paste0("-d g2_idx=", "'c(", paste(g2_idx, collapse = ","), ")'")
  slim_g1_idx    <- paste0("-d g1_idx=", "'c(", paste(g1_idx, collapse = ","), ")'")
  slim_SSs       <- paste0("-d SSs=", "'c(", paste(SSs, collapse = ","), ")'")
  slim_tau       <- paste0("-d tau=", "'c(", paste(tau, collapse = ","), ")'")
  
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
                         slim_SSs,
                         slim_tau,
                         slim_output,
                         slim_model_p2)                 # model part 2 - sampling
    
    # run slim on system
    system(slim_run_p2)
    
    ## HANDLING SLiM OUTPUT NE EQUILIBRIUM PERIOD - NE(1)
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
          
          # get the Ne for the whole period - harmonic mean
          meanNe1 <- 1/mean(1/ne1$ne, na.rm = TRUE)
          
        } else {
          meanNe1 <- as.numeric(NA)
        }
        
        # remove Ne(1) data after use it
        rm(ne1)
        
      } else {
        meanNe1 <- as.numeric(NA)
      }
      
      if (remove_files){
        file.remove(slim_output_ne1)
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
    
    ## HANDLING SLiM2 OUTPUT NE SAMPLING PERIOD - NE(2)
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
          
          # get the Ne for each generation in the period
          timesNe2 <- as.data.frame(t(ne2$ne))
          names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau[6])))
          
          # get the Ne for the whole interval of each pair of samples - harmonic mean
          meanNe2_Ava <- 1/mean(1/ne2$ne, na.rm = TRUE)
          meanNe2_Hum <- 1/mean(1/ne2$ne[(tau[1]+1):(tau[6]+1)], na.rm = TRUE)
          meanNe2_Dav <- 1/mean(1/ne2$ne[(tau[2]+1):(tau[6]+1)], na.rm = TRUE)
          meanNe2_Sta <- 1/mean(1/ne2$ne[(tau[3]+1):(tau[6]+1)], na.rm = TRUE)
          meanNe2_Ste <- 1/mean(1/ne2$ne[(tau[4]+1):(tau[6]+1)], na.rm = TRUE)
          meanNe2_Riv <- 1/mean(1/ne2$ne[(tau[5]+1):(tau[6]+1)], na.rm = TRUE)
          meanNe2_Pla = meanNe2_Riv
          
        } else {
          timesNe2 <- as.data.frame(t(rep(NA, (tau[6]+1))))
          names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau[6])))
          meanNe2_Ava <- as.numeric(NA)
          meanNe2_Hum <- as.numeric(NA)
          meanNe2_Dav <- as.numeric(NA)
          meanNe2_Sta <- as.numeric(NA)
          meanNe2_Ste <- as.numeric(NA)
          meanNe2_Riv = meanNe2_Pla <- as.numeric(NA)
          
        }
        
        # remove pedigreeNe data after use it
        rm(ne2)
        
      } else {
        timesNe2 <- as.data.frame(t(rep(NA, (tau[6]+1))))
        names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau[6])))
        meanNe2_Ava <- as.numeric(NA)
        meanNe2_Hum <- as.numeric(NA)
        meanNe2_Dav <- as.numeric(NA)
        meanNe2_Sta <- as.numeric(NA)
        meanNe2_Ste <- as.numeric(NA)
        meanNe2_Riv = meanNe2_Pla <- as.numeric(NA)
      }
      
      if (remove_files){
        file.remove(slim_output_ne2)
      }
      
    } else {
      timesNe2 <- as.data.frame(t(rep(NA, (tau[6]+1))))
      names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau[6])))
      meanNe2_Ava <- as.numeric(NA)
      meanNe2_Hum <- as.numeric(NA)
      meanNe2_Dav <- as.numeric(NA)
      meanNe2_Sta <- as.numeric(NA)
      meanNe2_Ste <- as.numeric(NA)
      meanNe2_Riv = meanNe2_Pla <- as.numeric(NA)
      
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
    
    ## HANDLING SLiM2 OUTPUT GENETIC LOAD OF SAMPLING PERIOD
    ##------------------------------------------------------------------------------
    
    # load the data
    slim_output_geneticLoad <- paste0(slim_output_folder,"slim_output_load_", sim, ".txt")
    
    if (file.exists(slim_output_geneticLoad)){
      info_geneticLoad_file = file.info(slim_output_geneticLoad)
      
      if (info_geneticLoad_file$size != 0){
        geneticLoad <- read.csv(file = slim_output_geneticLoad, sep = "\t", header = F, col.names = c("gen", "time", "genetciload"), na.strings = "NA")
        
        if (nrow(geneticLoad) != 0){
          
          # subset the genetic load for each pair of samples
          # check if there is any -Inf/Inf and remove the rows containing it
          geneticLoad_Ava <- geneticLoad$genetciload
          averageGenLoad_Ava <- mean(geneticLoad_Ava[!is.infinite(geneticLoad_Ava)], na.rm = TRUE)
          
          geneticLoad_Hum <- geneticLoad$genetciload[(tau[1]+1):(tau[6]+1)]
          averageGenLoad_Hum <- mean(geneticLoad_Hum[!is.infinite(geneticLoad_Hum)], na.rm = TRUE)
          
          geneticLoad_Dav <- geneticLoad$genetciload[(tau[2]+1):(tau[6]+1)]
          averageGenLoad_Dav <- mean(geneticLoad_Dav[!is.infinite(geneticLoad_Dav)], na.rm = TRUE)
          
          geneticLoad_Sta <- geneticLoad$genetciload[(tau[3]+1):(tau[6]+1)]
          averageGenLoad_Sta <- mean(geneticLoad_Sta[!is.infinite(geneticLoad_Sta)], na.rm = TRUE)
          
          geneticLoad_Ste <- geneticLoad$genetciload[(tau[4]+1):(tau[6]+1)]
          averageGenLoad_Ste <- mean(geneticLoad_Ste[!is.infinite(geneticLoad_Ste)], na.rm = TRUE)
          
          geneticLoad_Riv <- geneticLoad$genetciload[(tau[5]+1):(tau[6]+1)]
          averageGenLoad_Riv <- mean(geneticLoad_Riv[!is.infinite(geneticLoad_Riv)], na.rm = TRUE)
          
          averageGenLoad_Pla = averageGenLoad_Riv
          
        } else {
          averageGenLoad_Ava <- as.numeric(NA)
          averageGenLoad_Hum <- as.numeric(NA)
          averageGenLoad_Dav <- as.numeric(NA)
          averageGenLoad_Sta <- as.numeric(NA)
          averageGenLoad_Ste <- as.numeric(NA)
          averageGenLoad_Riv = averageGenLoad_Pla <- as.numeric(NA)
        }
        
        # remove genetic load data after use it
        rm(geneticLoad)
        
      } else {
        averageGenLoad_Ava <- as.numeric(NA)
        averageGenLoad_Hum <- as.numeric(NA)
        averageGenLoad_Dav <- as.numeric(NA)
        averageGenLoad_Sta <- as.numeric(NA)
        averageGenLoad_Ste <- as.numeric(NA)
        averageGenLoad_Riv = averageGenLoad_Pla <- as.numeric(NA)
      }
      
      if (remove_files){
        file.remove(slim_output_geneticLoad)
      }
      
    } else {
      
      averageGenLoad_Ava <- as.numeric(NA)
      averageGenLoad_Hum <- as.numeric(NA)
      averageGenLoad_Dav <- as.numeric(NA)
      averageGenLoad_Sta <- as.numeric(NA)
      averageGenLoad_Ste <- as.numeric(NA)
      averageGenLoad_Riv = averageGenLoad_Pla <- as.numeric(NA)
      
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
    
    ## HANDLING SLiM2 OUTPUT POPULATION MUTATIONS DATA OF SAMPLING PERIOD
    ##-------------------------------------------------------------------------------
    
    ## Avalon population
    ##-----------------
    
    filenames_genome_Ava <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=0, to=tau[6], by=1), "_", sim, ".txt")
      
    if(all(file.exists(c(filenames_genome_Ava)))){
      
      ### GENOME
      datalist_genome_Ava <- lapply(filenames_genome_Ava, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Ava <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Ava)
        
        if(!is.null(all_merged_genome_Ava)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Ava[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Ava    <- sum(all_merged_genome_Ava[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Ava[, paste0("S",tc)][!is.na(all_merged_genome_Ava[, paste0("S",tc)])])
            actual_pop_SelMean_Ava <- mean(all_merged_genome_Ava[all_merged_genome_Ava[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Ava   <- sd(all_merged_genome_Ava[all_merged_genome_Ava[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Ava)){
              
              if (any(meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Ava    <- sum(meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Ava[, paste0("S",tc)][!is.na(all_merged_genome_Ava[, paste0("S",tc)])])
                strong_positive_SelMean_Ava <- mean(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Ava   <- sd(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Ava    <- as.numeric(0)
                strong_positive_SelMean_Ava <- as.numeric(0)
                strong_positive_SelSd_Ava   <- as.numeric(0)
              }
              
              if (any(meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Ava    <- sum(meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Ava[, paste0("S",tc)][!is.na(all_merged_genome_Ava[, paste0("S",tc)])])
                strong_negative_SelMean_Ava <- mean(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Ava   <- sd(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Ava    <- as.numeric(0)
                strong_negative_SelMean_Ava <- as.numeric(0)
                strong_negative_SelSd_Ava   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Ava)  & is.na(strong_negative_prbe_Ava)){
                
                strong_pop_prbe_Ava    <- as.numeric(NA)
                strong_pop_SelMean_Ava <- as.numeric(NA)
                strong_pop_SelSd_Ava   <- as.numeric(NA)
                 
              } else {
                strong_pop_prbe_Ava    <- sum(strong_positive_prbe_Ava, strong_negative_prbe_Ava, na.rm = TRUE)
                strong_pop_SelMean_Ava <- sum(strong_positive_SelMean_Ava, strong_negative_SelMean_Ava, na.rm = TRUE)
                strong_pop_SelSd_Ava   <- sum(strong_positive_SelSd_Ava, strong_negative_SelSd_Ava, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Ava    <- as.numeric(NA)
              strong_pop_SelMean_Ava <- as.numeric(NA)
              strong_pop_SelSd_Ava   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Ava    <- as.numeric(0)
            actual_pop_SelMean_Ava <- as.numeric(0)
            actual_pop_SelSd_Ava   <- as.numeric(0)
            
            strong_pop_prbe_Ava    <- as.numeric(0)
            strong_pop_SelMean_Ava <- as.numeric(0)
            strong_pop_SelSd_Ava   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Ava    <- as.numeric(NA)
          actual_pop_SelMean_Ava <- as.numeric(NA)
          actual_pop_SelSd_Ava   <- as.numeric(NA)
          
          strong_pop_prbe_Ava    <- as.numeric(NA)
          strong_pop_SelMean_Ava <- as.numeric(NA)
          strong_pop_SelSd_Ava   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Ava <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Ava)
        
        if(!is.null(all_merged_genome_Ava)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Ava[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Ava    <- sum(all_merged_genome_Ava[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Ava[, "S"][!is.na(all_merged_genome_Ava[, "S"])])
            actual_pop_SelMean_Ava <- mean(all_merged_genome_Ava[all_merged_genome_Ava[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Ava   <- sd(all_merged_genome_Ava[all_merged_genome_Ava[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Ava)){
              if (any(meanNe2_Ava*all_merged_genome_Ava[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Ava    <- sum(meanNe2_Ava*all_merged_genome_Ava[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Ava[, "S"][!is.na(all_merged_genome_Ava[, "S"])])
                strong_positive_SelMean_Ava <- mean(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Ava   <- sd(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Ava    <- as.numeric(0)
                strong_positive_SelMean_Ava <- as.numeric(0)
                strong_positive_SelSd_Ava   <- as.numeric(0)
              }
              
              if (any(meanNe2_Ava*all_merged_genome_Ava[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Ava    <- sum(meanNe2_Ava*all_merged_genome_Ava[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Ava[, "S"][!is.na(all_merged_genome_Ava[, "S"])])
                strong_negative_SelMean_Ava <- mean(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Ava   <- sd(all_merged_genome_Ava[meanNe2_Ava*all_merged_genome_Ava[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Ava    <- as.numeric(0)
                strong_negative_SelMean_Ava <- as.numeric(0)
                strong_negative_SelSd_Ava   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Ava)  & is.na(strong_negative_prbe_Ava)){
                
                strong_pop_prbe_Ava    <- as.numeric(NA)
                strong_pop_SelMean_Ava <- as.numeric(NA)
                strong_pop_SelSd_Ava   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Ava    <- sum(strong_positive_prbe_Ava, strong_negative_prbe_Ava, na.rm = TRUE)
                strong_pop_SelMean_Ava <- sum(strong_positive_SelMean_Ava, strong_negative_SelMean_Ava, na.rm = TRUE)
                strong_pop_SelSd_Ava   <- sum(strong_positive_SelSd_Ava, strong_negative_SelSd_Ava, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Ava    <- as.numeric(NA)
              strong_pop_SelMean_Ava <- as.numeric(NA)
              strong_pop_SelSd_Ava   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Ava    <- as.numeric(0)
            actual_pop_SelMean_Ava <- as.numeric(0)
            actual_pop_SelSd_Ava   <- as.numeric(0)
            
            strong_pop_prbe_Ava    <- as.numeric(0)
            strong_pop_SelMean_Ava <- as.numeric(0)
            strong_pop_SelSd_Ava   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Ava    <- as.numeric(NA)
          actual_pop_SelMean_Ava <- as.numeric(NA)
          actual_pop_SelSd_Ava   <- as.numeric(NA)
          
          strong_pop_prbe_Ava    <- as.numeric(NA)
          strong_pop_SelMean_Ava <- as.numeric(NA)
          strong_pop_SelSd_Ava   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Ava)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Ava    <- as.numeric(NA)
      actual_pop_SelMean_Ava <- as.numeric(NA)
      actual_pop_SelSd_Ava   <- as.numeric(NA)
      
      strong_pop_prbe_Ava    <- as.numeric(NA)
      strong_pop_SelMean_Ava <- as.numeric(NA)
      strong_pop_SelSd_Ava   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Ava)))){
          file.copy(from = filenames_genome_Ava, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Avalon population"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## Humboldt population
    ##-----------------
    
    filenames_genome_Hum <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=tau[1], to=tau[6], by=1), "_", sim, ".txt")
    
    if(all(file.exists(c(filenames_genome_Hum)))){
      
      ### GENOME
      datalist_genome_Hum <- lapply(filenames_genome_Hum, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Hum <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Hum)
        
        if(!is.null(all_merged_genome_Hum)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Hum[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Hum    <- sum(all_merged_genome_Hum[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Hum[, paste0("S",tc)][!is.na(all_merged_genome_Hum[, paste0("S",tc)])])
            actual_pop_SelMean_Hum <- mean(all_merged_genome_Hum[all_merged_genome_Hum[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Hum   <- sd(all_merged_genome_Hum[all_merged_genome_Hum[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Hum)){
              
              if (any(meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Hum    <- sum(meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Hum[, paste0("S",tc)][!is.na(all_merged_genome_Hum[, paste0("S",tc)])])
                strong_positive_SelMean_Hum <- mean(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Hum   <- sd(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Hum    <- as.numeric(0)
                strong_positive_SelMean_Hum <- as.numeric(0)
                strong_positive_SelSd_Hum   <- as.numeric(0)
              }
              
              if (any(meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Hum    <- sum(meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Hum[, paste0("S",tc)][!is.na(all_merged_genome_Hum[, paste0("S",tc)])])
                strong_negative_SelMean_Hum <- mean(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Hum   <- sd(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Hum    <- as.numeric(0)
                strong_negative_SelMean_Hum <- as.numeric(0)
                strong_negative_SelSd_Hum   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Hum)  & is.na(strong_negative_prbe_Hum)){
                
                strong_pop_prbe_Hum    <- as.numeric(NA)
                strong_pop_SelMean_Hum <- as.numeric(NA)
                strong_pop_SelSd_Hum   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Hum    <- sum(strong_positive_prbe_Hum, strong_negative_prbe_Hum, na.rm = TRUE)
                strong_pop_SelMean_Hum <- sum(strong_positive_SelMean_Hum, strong_negative_SelMean_Hum, na.rm = TRUE)
                strong_pop_SelSd_Hum   <- sum(strong_positive_SelSd_Hum, strong_negative_SelSd_Hum, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Hum    <- as.numeric(NA)
              strong_pop_SelMean_Hum <- as.numeric(NA)
              strong_pop_SelSd_Hum   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Hum    <- as.numeric(0)
            actual_pop_SelMean_Hum <- as.numeric(0)
            actual_pop_SelSd_Hum   <- as.numeric(0)
            
            strong_pop_prbe_Hum    <- as.numeric(0)
            strong_pop_SelMean_Hum <- as.numeric(0)
            strong_pop_SelSd_Hum   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Hum    <- as.numeric(NA)
          actual_pop_SelMean_Hum <- as.numeric(NA)
          actual_pop_SelSd_Hum   <- as.numeric(NA)
          
          strong_pop_prbe_Hum    <- as.numeric(NA)
          strong_pop_SelMean_Hum <- as.numeric(NA)
          strong_pop_SelSd_Hum   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Hum <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Hum)
        
        if(!is.null(all_merged_genome_Hum)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Hum[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Hum    <- sum(all_merged_genome_Hum[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Hum[, "S"][!is.na(all_merged_genome_Hum[, "S"])])
            actual_pop_SelMean_Hum <- mean(all_merged_genome_Hum[all_merged_genome_Hum[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Hum   <- sd(all_merged_genome_Hum[all_merged_genome_Hum[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Hum)){
              if (any(meanNe2_Hum*all_merged_genome_Hum[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Hum    <- sum(meanNe2_Hum*all_merged_genome_Hum[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Hum[, "S"][!is.na(all_merged_genome_Hum[, "S"])])
                strong_positive_SelMean_Hum <- mean(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Hum   <- sd(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Hum    <- as.numeric(0)
                strong_positive_SelMean_Hum <- as.numeric(0)
                strong_positive_SelSd_Hum   <- as.numeric(0)
              }
              
              if (any(meanNe2_Hum*all_merged_genome_Hum[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Hum    <- sum(meanNe2_Hum*all_merged_genome_Hum[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Hum[, "S"][!is.na(all_merged_genome_Hum[, "S"])])
                strong_negative_SelMean_Hum <- mean(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Hum   <- sd(all_merged_genome_Hum[meanNe2_Hum*all_merged_genome_Hum[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Hum    <- as.numeric(0)
                strong_negative_SelMean_Hum <- as.numeric(0)
                strong_negative_SelSd_Hum   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Hum)  & is.na(strong_negative_prbe_Hum)){
                
                strong_pop_prbe_Hum    <- as.numeric(NA)
                strong_pop_SelMean_Hum <- as.numeric(NA)
                strong_pop_SelSd_Hum   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Hum    <- sum(strong_positive_prbe_Hum, strong_negative_prbe_Hum, na.rm = TRUE)
                strong_pop_SelMean_Hum <- sum(strong_positive_SelMean_Hum, strong_negative_SelMean_Hum, na.rm = TRUE)
                strong_pop_SelSd_Hum   <- sum(strong_positive_SelSd_Hum, strong_negative_SelSd_Hum, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Hum    <- as.numeric(NA)
              strong_pop_SelMean_Hum <- as.numeric(NA)
              strong_pop_SelSd_Hum   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Hum    <- as.numeric(0)
            actual_pop_SelMean_Hum <- as.numeric(0)
            actual_pop_SelSd_Hum   <- as.numeric(0)
            
            strong_pop_prbe_Hum    <- as.numeric(0)
            strong_pop_SelMean_Hum <- as.numeric(0)
            strong_pop_SelSd_Hum   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Hum    <- as.numeric(NA)
          actual_pop_SelMean_Hum <- as.numeric(NA)
          actual_pop_SelSd_Hum   <- as.numeric(NA)
          
          strong_pop_prbe_Hum    <- as.numeric(NA)
          strong_pop_SelMean_Hum <- as.numeric(NA)
          strong_pop_SelSd_Hum   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Hum)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Hum    <- as.numeric(NA)
      actual_pop_SelMean_Hum <- as.numeric(NA)
      actual_pop_SelSd_Hum   <- as.numeric(NA)
      
      strong_pop_prbe_Hum    <- as.numeric(NA)
      strong_pop_SelMean_Hum <- as.numeric(NA)
      strong_pop_SelSd_Hum   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Hum)))){
          file.copy(from = filenames_genome_Hum, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Humboldt population"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## Davis population
    ##-----------------
    
    filenames_genome_Dav <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=tau[2], to=tau[6], by=1), "_", sim, ".txt")
    
    if(all(file.exists(c(filenames_genome_Dav)))){
      
      ### GENOME
      datalist_genome_Dav <- lapply(filenames_genome_Dav, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Dav <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Dav)
        
        if(!is.null(all_merged_genome_Dav)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Dav[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Dav    <- sum(all_merged_genome_Dav[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Dav[, paste0("S",tc)][!is.na(all_merged_genome_Dav[, paste0("S",tc)])])
            actual_pop_SelMean_Dav <- mean(all_merged_genome_Dav[all_merged_genome_Dav[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Dav   <- sd(all_merged_genome_Dav[all_merged_genome_Dav[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Dav)){
              
              if (any(meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Dav    <- sum(meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Dav[, paste0("S",tc)][!is.na(all_merged_genome_Dav[, paste0("S",tc)])])
                strong_positive_SelMean_Dav <- mean(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Dav   <- sd(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Dav    <- as.numeric(0)
                strong_positive_SelMean_Dav <- as.numeric(0)
                strong_positive_SelSd_Dav   <- as.numeric(0)
              }
              
              if (any(meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Dav    <- sum(meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Dav[, paste0("S",tc)][!is.na(all_merged_genome_Dav[, paste0("S",tc)])])
                strong_negative_SelMean_Dav <- mean(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Dav   <- sd(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Dav    <- as.numeric(0)
                strong_negative_SelMean_Dav <- as.numeric(0)
                strong_negative_SelSd_Dav   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Dav)  & is.na(strong_negative_prbe_Dav)){
                
                strong_pop_prbe_Dav    <- as.numeric(NA)
                strong_pop_SelMean_Dav <- as.numeric(NA)
                strong_pop_SelSd_Dav   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Dav    <- sum(strong_positive_prbe_Dav, strong_negative_prbe_Dav, na.rm = TRUE)
                strong_pop_SelMean_Dav <- sum(strong_positive_SelMean_Dav, strong_negative_SelMean_Dav, na.rm = TRUE)
                strong_pop_SelSd_Dav   <- sum(strong_positive_SelSd_Dav, strong_negative_SelSd_Dav, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Dav    <- as.numeric(NA)
              strong_pop_SelMean_Dav <- as.numeric(NA)
              strong_pop_SelSd_Dav   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Dav    <- as.numeric(0)
            actual_pop_SelMean_Dav <- as.numeric(0)
            actual_pop_SelSd_Dav   <- as.numeric(0)
            
            strong_pop_prbe_Dav    <- as.numeric(0)
            strong_pop_SelMean_Dav <- as.numeric(0)
            strong_pop_SelSd_Dav   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Dav    <- as.numeric(NA)
          actual_pop_SelMean_Dav <- as.numeric(NA)
          actual_pop_SelSd_Dav   <- as.numeric(NA)
          
          strong_pop_prbe_Dav    <- as.numeric(NA)
          strong_pop_SelMean_Dav <- as.numeric(NA)
          strong_pop_SelSd_Dav   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Dav <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Dav)
        
        if(!is.null(all_merged_genome_Dav)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Dav[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Dav    <- sum(all_merged_genome_Dav[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Dav[, "S"][!is.na(all_merged_genome_Dav[, "S"])])
            actual_pop_SelMean_Dav <- mean(all_merged_genome_Dav[all_merged_genome_Dav[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Dav   <- sd(all_merged_genome_Dav[all_merged_genome_Dav[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Dav)){
              if (any(meanNe2_Dav*all_merged_genome_Dav[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Dav    <- sum(meanNe2_Dav*all_merged_genome_Dav[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Dav[, "S"][!is.na(all_merged_genome_Dav[, "S"])])
                strong_positive_SelMean_Dav <- mean(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Dav   <- sd(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Dav    <- as.numeric(0)
                strong_positive_SelMean_Dav <- as.numeric(0)
                strong_positive_SelSd_Dav   <- as.numeric(0)
              }
              
              if (any(meanNe2_Dav*all_merged_genome_Dav[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Dav    <- sum(meanNe2_Dav*all_merged_genome_Dav[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Dav[, "S"][!is.na(all_merged_genome_Dav[, "S"])])
                strong_negative_SelMean_Dav <- mean(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Dav   <- sd(all_merged_genome_Dav[meanNe2_Dav*all_merged_genome_Dav[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Dav    <- as.numeric(0)
                strong_negative_SelMean_Dav <- as.numeric(0)
                strong_negative_SelSd_Dav   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Dav)  & is.na(strong_negative_prbe_Dav)){
                
                strong_pop_prbe_Dav    <- as.numeric(NA)
                strong_pop_SelMean_Dav <- as.numeric(NA)
                strong_pop_SelSd_Dav   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Dav    <- sum(strong_positive_prbe_Dav, strong_negative_prbe_Dav, na.rm = TRUE)
                strong_pop_SelMean_Dav <- sum(strong_positive_SelMean_Dav, strong_negative_SelMean_Dav, na.rm = TRUE)
                strong_pop_SelSd_Dav   <- sum(strong_positive_SelSd_Dav, strong_negative_SelSd_Dav, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Dav    <- as.numeric(NA)
              strong_pop_SelMean_Dav <- as.numeric(NA)
              strong_pop_SelSd_Dav   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Dav    <- as.numeric(0)
            actual_pop_SelMean_Dav <- as.numeric(0)
            actual_pop_SelSd_Dav   <- as.numeric(0)
            
            strong_pop_prbe_Dav    <- as.numeric(0)
            strong_pop_SelMean_Dav <- as.numeric(0)
            strong_pop_SelSd_Dav   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Dav    <- as.numeric(NA)
          actual_pop_SelMean_Dav <- as.numeric(NA)
          actual_pop_SelSd_Dav   <- as.numeric(NA)
          
          strong_pop_prbe_Dav    <- as.numeric(NA)
          strong_pop_SelMean_Dav <- as.numeric(NA)
          strong_pop_SelSd_Dav   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Dav)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Dav    <- as.numeric(NA)
      actual_pop_SelMean_Dav <- as.numeric(NA)
      actual_pop_SelSd_Dav   <- as.numeric(NA)
      
      strong_pop_prbe_Dav    <- as.numeric(NA)
      strong_pop_SelMean_Dav <- as.numeric(NA)
      strong_pop_SelSd_Dav   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Dav)))){
          file.copy(from = filenames_genome_Dav, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Davis population"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## Stanislaus population
    ##-----------------
    
    filenames_genome_Sta <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=tau[3], to=tau[6], by=1), "_", sim, ".txt")
    
    if(all(file.exists(c(filenames_genome_Sta)))){
      
      ### GENOME
      datalist_genome_Sta <- lapply(filenames_genome_Sta, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Sta <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Sta)
        
        if(!is.null(all_merged_genome_Sta)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Sta[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Sta    <- sum(all_merged_genome_Sta[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Sta[, paste0("S",tc)][!is.na(all_merged_genome_Sta[, paste0("S",tc)])])
            actual_pop_SelMean_Sta <- mean(all_merged_genome_Sta[all_merged_genome_Sta[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Sta   <- sd(all_merged_genome_Sta[all_merged_genome_Sta[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Sta)){
              
              if (any(meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Sta    <- sum(meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Sta[, paste0("S",tc)][!is.na(all_merged_genome_Sta[, paste0("S",tc)])])
                strong_positive_SelMean_Sta <- mean(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Sta   <- sd(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Sta    <- as.numeric(0)
                strong_positive_SelMean_Sta <- as.numeric(0)
                strong_positive_SelSd_Sta   <- as.numeric(0)
              }
              
              if (any(meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Sta    <- sum(meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Sta[, paste0("S",tc)][!is.na(all_merged_genome_Sta[, paste0("S",tc)])])
                strong_negative_SelMean_Sta <- mean(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Sta   <- sd(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Sta    <- as.numeric(0)
                strong_negative_SelMean_Sta <- as.numeric(0)
                strong_negative_SelSd_Sta   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Sta)  & is.na(strong_negative_prbe_Sta)){
                
                strong_pop_prbe_Sta    <- as.numeric(NA)
                strong_pop_SelMean_Sta <- as.numeric(NA)
                strong_pop_SelSd_Sta   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Sta    <- sum(strong_positive_prbe_Sta, strong_negative_prbe_Sta, na.rm = TRUE)
                strong_pop_SelMean_Sta <- sum(strong_positive_SelMean_Sta, strong_negative_SelMean_Sta, na.rm = TRUE)
                strong_pop_SelSd_Sta   <- sum(strong_positive_SelSd_Sta, strong_negative_SelSd_Sta, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Sta    <- as.numeric(NA)
              strong_pop_SelMean_Sta <- as.numeric(NA)
              strong_pop_SelSd_Sta   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Sta    <- as.numeric(0)
            actual_pop_SelMean_Sta <- as.numeric(0)
            actual_pop_SelSd_Sta   <- as.numeric(0)
            
            strong_pop_prbe_Sta    <- as.numeric(0)
            strong_pop_SelMean_Sta <- as.numeric(0)
            strong_pop_SelSd_Sta   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Sta    <- as.numeric(NA)
          actual_pop_SelMean_Sta <- as.numeric(NA)
          actual_pop_SelSd_Sta   <- as.numeric(NA)
          
          strong_pop_prbe_Sta    <- as.numeric(NA)
          strong_pop_SelMean_Sta <- as.numeric(NA)
          strong_pop_SelSd_Sta   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Sta <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Sta)
        
        if(!is.null(all_merged_genome_Sta)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Sta[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Sta    <- sum(all_merged_genome_Sta[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Sta[, "S"][!is.na(all_merged_genome_Sta[, "S"])])
            actual_pop_SelMean_Sta <- mean(all_merged_genome_Sta[all_merged_genome_Sta[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Sta   <- sd(all_merged_genome_Sta[all_merged_genome_Sta[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Sta)){
              if (any(meanNe2_Sta*all_merged_genome_Sta[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Sta    <- sum(meanNe2_Sta*all_merged_genome_Sta[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Sta[, "S"][!is.na(all_merged_genome_Sta[, "S"])])
                strong_positive_SelMean_Sta <- mean(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Sta   <- sd(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Sta    <- as.numeric(0)
                strong_positive_SelMean_Sta <- as.numeric(0)
                strong_positive_SelSd_Sta   <- as.numeric(0)
              }
              
              if (any(meanNe2_Sta*all_merged_genome_Sta[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Sta    <- sum(meanNe2_Sta*all_merged_genome_Sta[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Sta[, "S"][!is.na(all_merged_genome_Sta[, "S"])])
                strong_negative_SelMean_Sta <- mean(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Sta   <- sd(all_merged_genome_Sta[meanNe2_Sta*all_merged_genome_Sta[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Sta    <- as.numeric(0)
                strong_negative_SelMean_Sta <- as.numeric(0)
                strong_negative_SelSd_Sta   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Sta)  & is.na(strong_negative_prbe_Sta)){
                
                strong_pop_prbe_Sta    <- as.numeric(NA)
                strong_pop_SelMean_Sta <- as.numeric(NA)
                strong_pop_SelSd_Sta   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Sta    <- sum(strong_positive_prbe_Sta, strong_negative_prbe_Sta, na.rm = TRUE)
                strong_pop_SelMean_Sta <- sum(strong_positive_SelMean_Sta, strong_negative_SelMean_Sta, na.rm = TRUE)
                strong_pop_SelSd_Sta   <- sum(strong_positive_SelSd_Sta, strong_negative_SelSd_Sta, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Sta    <- as.numeric(NA)
              strong_pop_SelMean_Sta <- as.numeric(NA)
              strong_pop_SelSd_Sta   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Sta    <- as.numeric(0)
            actual_pop_SelMean_Sta <- as.numeric(0)
            actual_pop_SelSd_Sta   <- as.numeric(0)
            
            strong_pop_prbe_Sta    <- as.numeric(0)
            strong_pop_SelMean_Sta <- as.numeric(0)
            strong_pop_SelSd_Sta   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Sta    <- as.numeric(NA)
          actual_pop_SelMean_Sta <- as.numeric(NA)
          actual_pop_SelSd_Sta   <- as.numeric(NA)
          
          strong_pop_prbe_Sta    <- as.numeric(NA)
          strong_pop_SelMean_Sta <- as.numeric(NA)
          strong_pop_SelSd_Sta   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Sta)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Sta    <- as.numeric(NA)
      actual_pop_SelMean_Sta <- as.numeric(NA)
      actual_pop_SelSd_Sta   <- as.numeric(NA)
      
      strong_pop_prbe_Sta    <- as.numeric(NA)
      strong_pop_SelMean_Sta <- as.numeric(NA)
      strong_pop_SelSd_Sta   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Sta)))){
          file.copy(from = filenames_genome_Sta, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Stanislaus population"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## Stebbins population
    ##-----------------
    
    filenames_genome_Ste <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=tau[4], to=tau[6], by=1), "_", sim, ".txt")
    
    if(all(file.exists(c(filenames_genome_Ste)))){
      
      ### GENOME
      datalist_genome_Ste <- lapply(filenames_genome_Ste, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Ste <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Ste)
        
        if(!is.null(all_merged_genome_Ste)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Ste[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Ste    <- sum(all_merged_genome_Ste[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Ste[, paste0("S",tc)][!is.na(all_merged_genome_Ste[, paste0("S",tc)])])
            actual_pop_SelMean_Ste <- mean(all_merged_genome_Ste[all_merged_genome_Ste[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Ste   <- sd(all_merged_genome_Ste[all_merged_genome_Ste[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Ste)){
              
              if (any(meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Ste    <- sum(meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Ste[, paste0("S",tc)][!is.na(all_merged_genome_Ste[, paste0("S",tc)])])
                strong_positive_SelMean_Ste <- mean(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Ste   <- sd(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Ste    <- as.numeric(0)
                strong_positive_SelMean_Ste <- as.numeric(0)
                strong_positive_SelSd_Ste   <- as.numeric(0)
              }
              
              if (any(meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Ste    <- sum(meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Ste[, paste0("S",tc)][!is.na(all_merged_genome_Ste[, paste0("S",tc)])])
                strong_negative_SelMean_Ste <- mean(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Ste   <- sd(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Ste    <- as.numeric(0)
                strong_negative_SelMean_Ste <- as.numeric(0)
                strong_negative_SelSd_Ste   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Ste)  & is.na(strong_negative_prbe_Ste)){
                
                strong_pop_prbe_Ste    <- as.numeric(NA)
                strong_pop_SelMean_Ste <- as.numeric(NA)
                strong_pop_SelSd_Ste   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Ste    <- sum(strong_positive_prbe_Ste, strong_negative_prbe_Ste, na.rm = TRUE)
                strong_pop_SelMean_Ste <- sum(strong_positive_SelMean_Ste, strong_negative_SelMean_Ste, na.rm = TRUE)
                strong_pop_SelSd_Ste   <- sum(strong_positive_SelSd_Ste, strong_negative_SelSd_Ste, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Ste    <- as.numeric(NA)
              strong_pop_SelMean_Ste <- as.numeric(NA)
              strong_pop_SelSd_Ste   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Ste    <- as.numeric(0)
            actual_pop_SelMean_Ste <- as.numeric(0)
            actual_pop_SelSd_Ste   <- as.numeric(0)
            
            strong_pop_prbe_Ste    <- as.numeric(0)
            strong_pop_SelMean_Ste <- as.numeric(0)
            strong_pop_SelSd_Ste   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Ste    <- as.numeric(NA)
          actual_pop_SelMean_Ste <- as.numeric(NA)
          actual_pop_SelSd_Ste   <- as.numeric(NA)
          
          strong_pop_prbe_Ste    <- as.numeric(NA)
          strong_pop_SelMean_Ste <- as.numeric(NA)
          strong_pop_SelSd_Ste   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Ste <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Ste)
        
        if(!is.null(all_merged_genome_Ste)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Ste[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Ste    <- sum(all_merged_genome_Ste[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Ste[, "S"][!is.na(all_merged_genome_Ste[, "S"])])
            actual_pop_SelMean_Ste <- mean(all_merged_genome_Ste[all_merged_genome_Ste[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Ste   <- sd(all_merged_genome_Ste[all_merged_genome_Ste[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Ste)){
              if (any(meanNe2_Ste*all_merged_genome_Ste[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Ste    <- sum(meanNe2_Ste*all_merged_genome_Ste[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Ste[, "S"][!is.na(all_merged_genome_Ste[, "S"])])
                strong_positive_SelMean_Ste <- mean(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Ste   <- sd(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Ste    <- as.numeric(0)
                strong_positive_SelMean_Ste <- as.numeric(0)
                strong_positive_SelSd_Ste   <- as.numeric(0)
              }
              
              if (any(meanNe2_Ste*all_merged_genome_Ste[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Ste    <- sum(meanNe2_Ste*all_merged_genome_Ste[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Ste[, "S"][!is.na(all_merged_genome_Ste[, "S"])])
                strong_negative_SelMean_Ste <- mean(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Ste   <- sd(all_merged_genome_Ste[meanNe2_Ste*all_merged_genome_Ste[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Ste    <- as.numeric(0)
                strong_negative_SelMean_Ste <- as.numeric(0)
                strong_negative_SelSd_Ste   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Ste)  & is.na(strong_negative_prbe_Ste)){
                
                strong_pop_prbe_Ste    <- as.numeric(NA)
                strong_pop_SelMean_Ste <- as.numeric(NA)
                strong_pop_SelSd_Ste   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Ste    <- sum(strong_positive_prbe_Ste, strong_negative_prbe_Ste, na.rm = TRUE)
                strong_pop_SelMean_Ste <- sum(strong_positive_SelMean_Ste, strong_negative_SelMean_Ste, na.rm = TRUE)
                strong_pop_SelSd_Ste   <- sum(strong_positive_SelSd_Ste, strong_negative_SelSd_Ste, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Ste    <- as.numeric(NA)
              strong_pop_SelMean_Ste <- as.numeric(NA)
              strong_pop_SelSd_Ste   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Ste    <- as.numeric(0)
            actual_pop_SelMean_Ste <- as.numeric(0)
            actual_pop_SelSd_Ste   <- as.numeric(0)
            
            strong_pop_prbe_Ste    <- as.numeric(0)
            strong_pop_SelMean_Ste <- as.numeric(0)
            strong_pop_SelSd_Ste   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Ste    <- as.numeric(NA)
          actual_pop_SelMean_Ste <- as.numeric(NA)
          actual_pop_SelSd_Ste   <- as.numeric(NA)
          
          strong_pop_prbe_Ste    <- as.numeric(NA)
          strong_pop_SelMean_Ste <- as.numeric(NA)
          strong_pop_SelSd_Ste   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Ste)
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Ste    <- as.numeric(NA)
      actual_pop_SelMean_Ste <- as.numeric(NA)
      actual_pop_SelSd_Ste   <- as.numeric(NA)
      
      strong_pop_prbe_Ste    <- as.numeric(NA)
      strong_pop_SelMean_Ste <- as.numeric(NA)
      strong_pop_SelSd_Ste   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Ste)))){
          file.copy(from = filenames_genome_Ste, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Stebbins population"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
    }
    
    ## Riverside and Placerita population
    ##-----------------
    
    filenames_genome_Riv <- paste0(slim_output_folder, "slim_output_pmuts_t", seq(from=tau[5], to=tau[6], by=1), "_", sim, ".txt")
    
    if(all(file.exists(c(filenames_genome_Riv)))){
      
      ### GENOME
      datalist_genome_Riv <- lapply(filenames_genome_Riv, function(x){read.table(file= x, header=T, na.strings = "NA")})
      
      if (model_type == 3){
        all_merged_genome_Riv <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4), all = TRUE)}, datalist_genome_Riv)
        
        if(!is.null(all_merged_genome_Riv)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Riv[, paste0("S", tc)] != 0, na.rm = TRUE)){
            
            actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- sum(all_merged_genome_Riv[, paste0("S",tc)] != 0, na.rm = TRUE)/length(all_merged_genome_Riv[, paste0("S",tc)][!is.na(all_merged_genome_Riv[, paste0("S",tc)])])
            actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- mean(all_merged_genome_Riv[all_merged_genome_Riv[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- sd(all_merged_genome_Riv[all_merged_genome_Riv[, paste0("S",tc)] != 0, paste0("S",tc)], na.rm = TRUE)
            
            if (!is.na(meanNe2_Riv)){
              
              if (any(meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] > 1, na.rm = TRUE)){
                strong_positive_prbe_Riv    = strong_positive_prbe_Pla    <- sum(meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] > 1, na.rm = TRUE)/length(all_merged_genome_Riv[, paste0("S",tc)][!is.na(all_merged_genome_Riv[, paste0("S",tc)])])
                strong_positive_SelMean_Riv = strong_positive_SelMean_Pla <- mean(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                strong_positive_SelSd_Riv   = strong_positive_SelSd_Pla   <- sd(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] > 1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Riv    = strong_positive_prbe_Pla    <- as.numeric(0)
                strong_positive_SelMean_Riv = strong_positive_SelMean_Pla <- as.numeric(0)
                strong_positive_SelSd_Riv   = strong_positive_SelSd_Pla   <- as.numeric(0)
              }
              
              if (any(meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] < -1, na.rm = TRUE)){
                strong_negative_prbe_Riv    = strong_negative_prbe_Pla    <- sum(meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] < -1, na.rm = TRUE)/length(all_merged_genome_Riv[, paste0("S",tc)][!is.na(all_merged_genome_Riv[, paste0("S",tc)])])
                strong_negative_SelMean_Riv = strong_negative_SelMean_Pla <- mean(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                strong_negative_SelSd_Riv   = strong_negative_SelSd_Pla   <- sd(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, paste0("S",tc)] < -1, paste0("S",tc)], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Riv    = strong_negative_prbe_Pla    <- as.numeric(0)
                strong_negative_SelMean_Riv = strong_negative_SelMean_Pla <- as.numeric(0)
                strong_negative_SelSd_Riv   = strong_negative_SelSd_Pla   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Riv)  & is.na(strong_negative_prbe_Riv)){
                
                strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
                strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
                strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- sum(strong_positive_prbe_Riv, strong_negative_prbe_Riv, na.rm = TRUE)
                strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- sum(strong_positive_SelMean_Riv, strong_negative_SelMean_Riv, na.rm = TRUE)
                strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- sum(strong_positive_SelSd_Riv, strong_negative_SelSd_Riv, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
              strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
              strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(0)
            actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(0)
            actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(0)
            
            strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(0)
            strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(0)
            strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(NA)
          actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(NA)
          actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(NA)
          
          strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
          strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
          strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
        }
        
      } else {
        all_merged_genome_Riv <- Reduce(function(x,y) {merge(x,y, by=c(1,2,3,4,5), all = TRUE)}, datalist_genome_Riv)
        
        if(!is.null(all_merged_genome_Riv)){
          
          ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION 
          if (any(all_merged_genome_Riv[, "S"] != 0, na.rm = TRUE)) {
            
            actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- sum(all_merged_genome_Riv[, "S"] != 0, na.rm = TRUE)/length(all_merged_genome_Riv[, "S"][!is.na(all_merged_genome_Riv[, "S"])])
            actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- mean(all_merged_genome_Riv[all_merged_genome_Riv[, "S"] != 0, "S"], na.rm = TRUE)
            actual_pop_SelSd_Riv   =  actual_pop_SelSd_Pla  <- sd(all_merged_genome_Riv[all_merged_genome_Riv[, "S"] != 0, "S"], na.rm = TRUE)
            
            if (!is.na(meanNe2_Riv)){
              if (any(meanNe2_Riv*all_merged_genome_Riv[, "S"] > 1, na.rm = TRUE)){
                strong_positive_prbe_Riv    = strong_positive_prbe_Pla    <- sum(meanNe2_Riv*all_merged_genome_Riv[, "S"] > 1, na.rm = TRUE)/length(all_merged_genome_Riv[, "S"][!is.na(all_merged_genome_Riv[, "S"])])
                strong_positive_SelMean_Riv = strong_positive_SelMean_Pla <- mean(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, "S"] > 1, "S"], na.rm = TRUE)
                strong_positive_SelSd_Riv   = strong_positive_SelSd_Pla   <- sd(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, "S"] > 1, "S"], na.rm = TRUE)
                
              } else {
                strong_positive_prbe_Riv    = strong_positive_prbe_Pla    <- as.numeric(0)
                strong_positive_SelMean_Riv = strong_positive_SelMean_Pla <- as.numeric(0)
                strong_positive_SelSd_Riv   = strong_positive_SelSd_Pla   <- as.numeric(0)
              }
              
              if (any(meanNe2_Riv*all_merged_genome_Riv[, "S"] < -1, na.rm = TRUE)){
                strong_negative_prbe_Riv    = strong_negative_prbe_Pla    <- sum(meanNe2_Riv*all_merged_genome_Riv[, "S"] < -1, na.rm = TRUE)/length(all_merged_genome_Riv[, "S"][!is.na(all_merged_genome_Riv[, "S"])])
                strong_negative_SelMean_Riv = strong_negative_SelMean_Pla <- mean(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, "S"] < -1, "S"], na.rm = TRUE)
                strong_negative_SelSd_Riv   = strong_negative_SelSd_Pla   <- sd(all_merged_genome_Riv[meanNe2_Riv*all_merged_genome_Riv[, "S"] < -1, "S"], na.rm = TRUE)
                
              } else {
                strong_negative_prbe_Riv    = strong_negative_prbe_Pla    <- as.numeric(0)
                strong_negative_SelMean_Riv = strong_negative_SelMean_Pla <- as.numeric(0)
                strong_negative_SelSd_Riv   = strong_negative_SelSd_Pla   <- as.numeric(0)
              }
              
              if(is.na(strong_positive_prbe_Riv)  & is.na(strong_negative_prbe_Riv)){
                
                strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
                strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
                strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
                
              } else {
                strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- sum(strong_positive_prbe_Riv, strong_negative_prbe_Riv, na.rm = TRUE)
                strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- sum(strong_positive_SelMean_Riv, strong_negative_SelMean_Riv, na.rm = TRUE)
                strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- sum(strong_positive_SelSd_Riv, strong_negative_SelSd_Riv, na.rm = TRUE)
              }
              
            } else {
              strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
              strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
              strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
            }
            
          } else {
            actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(0)
            actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(0)
            actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(0)
            
            strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(0)
            strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(0)
            strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(0)
          }
          
        } else {
          actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(NA)
          actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(NA)
          actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(NA)
          
          strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
          strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
          strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
        }
      }
      
      # remove population allele frequency data after use it
      rm(datalist_genome_Riv)
      
      if (remove_files){
        file.remove(filenames_genome_Ava)
      }
      
    } else {
      
      ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
      actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(NA)
      actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(NA)
      actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(NA)
      
      strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
      strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
      strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(any(file.exists(c(filenames_genome_Riv)))){
          file.copy(from = filenames_genome_Riv, to = debug_output_folder)
        }
        
        debug_message <- "no pmuts file found for Riverside and Placerita populations"
        
        debug_dump  <- suppressWarnings(cbind(as.factor(model_type), 
                                              sim_seed, sim, 
                                              mu, rr, selfing, Neq, Ncs,
                                              gammaM, gammak, tc,
                                              PrGWSel, prbe, debug_message))
        rownames(debug_dump) <- sim
        write.table(debug_dump,file=paste0(debug_output_folder, "debug_",sim,".txt"), row.names = F, quote = F)
      }
      
      if (remove_files){
        if (any(file.exists(filenames_genome_Ava))){
          file.remove(filenames_genome_Ava)
        }
      }
    }
    
    ## HANDLINDING SLiM2 OUTPUT - SAMPLED INDIVIDUALS
    ##---------------------------------------------------------------------------------
    
    ## Avalon population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts1        <- paste0(slim_output_folder,"slim_output_sample_ts1_", sim, ".vcf")
    slim_output_sample_ts7        <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts1, slim_output_sample_ts7)))){
      
      slim_output_sample_ts1_sorted <- paste0(slim_output_folder,"slim_output_sample_ts1_", sim, "_sorted" , ".vcf")
      slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts1_vcf <- paste("grep '^#'", slim_output_sample_ts1, ">", slim_output_sample_ts1_sorted,
                                  "&& grep -v '^#'", slim_output_sample_ts1, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts1_sorted)
      
      sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                  "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
      
      system(sort_sample_ts1_vcf)
      system(sort_sample_ts7_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts1_sorted))
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts1_sorted, ".gz"))) 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts1_sorted, ".gz")))
        system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        
      }
      
      # merge and get the data
      slim_output_sample_merged_Ava <- paste0(slim_output_folder,"slim_output_sample_merged_Ava_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts1_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Ava) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Ava)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[1], by=1), "@pop1", "")
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Ava, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(8:10)]
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Ava_", sim,".txt")) 
              write(paste(0, (tau[6]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Ava_", sim,".txt"), append = TRUE) 
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Ava_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Ava <- paste0("egglib_input_sample_Ava", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Ava), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Ava[all_merged_genome_Ava$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Ava[all_merged_genome_Ava$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Ava <- paste0("egglib_input_sample_selcoeff_Ava", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder, egglib_selcoeff_file_Ava), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Ava))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                             paste0(getwd(), "/", path_to_egglib_summstat),
                             paste0("input-file=", egglib_input_folder, egglib_converted_file_Ava),
                             paste0("output-file=", egglib_output_folder, "egglib_output_sample_Ava", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstats_Ava <- paste0(egglib_output_folder,"egglib_output_sample_Ava", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstats_Ava)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstats_Ava, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Ava    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Ava <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Ava   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Ava    <- as.numeric(0)
                  actual_sample_SelMean_Ava <- as.numeric(0)
                  actual_sample_SelSd_Ava   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Ava)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Ava*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Ava    <- sum(meanNe2_Ava*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Ava <- mean(slim_to_egglib_snps_[meanNe2_Ava*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Ava   <- sd(slim_to_egglib_snps_[meanNe2_Ava*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Ava    <- as.numeric(0)
                    positive_sample_SelMean_Ava <- as.numeric(0)
                    positive_sample_SelSd_Ava   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Ava*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Ava    <- sum(meanNe2_Ava*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Ava <- mean(slim_to_egglib_snps_[meanNe2_Ava*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Ava   <- sd(slim_to_egglib_snps_[meanNe2_Ava*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Ava    <- as.numeric(0)
                    negative_sample_SelMean_Ava <- as.numeric(0)
                    negative_sample_SelSd_Ava   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Ava)  & is.na(negative_sample_prbe_Ava)){
                    
                    strong_sample_prbe_Ava    <- as.numeric(NA)
                    strong_sample_SelMean_Ava <- as.numeric(NA)
                    strong_sample_SelSd_Ava   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Ava    <- sum(positive_sample_prbe_Ava, negative_sample_prbe_Ava, na.rm = TRUE)
                    strong_sample_SelMean_Ava <- sum(positive_sample_SelMean_Ava, negative_sample_SelMean_Ava, na.rm = TRUE)
                    strong_sample_SelSd_Ava   <- sum(positive_sample_SelSd_Ava, negative_sample_SelSd_Ava, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Ava    <- as.numeric(NA)
                  strong_sample_SelMean_Ava <- as.numeric(NA)
                  strong_sample_SelSd_Ava   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Ava*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Ava*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Ava = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Ava = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Ava = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Ava = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Ava = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Ava = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Ava = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Ava = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Ava = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Ava = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava,
                                                      FSTfdrNS01_Ava,
                                                      FSTfdrNS02_Ava,
                                                      FSTfdrNS05_Ava,
                                                      FSTfdrNS10_Ava))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Ava  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Ava  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Ava  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Ava  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Ava = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Ava = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Ava = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Ava = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[1] + (SSs[7]-3) + 6)]
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Ava <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Ava", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Ava    <- as.numeric(NA)
                actual_sample_SelMean_Ava <- as.numeric(NA)
                actual_sample_SelSd_Ava   <- as.numeric(NA)
                
                strong_sample_prbe_Ava    <- as.numeric(NA)
                strong_sample_SelMean_Ava <- as.numeric(NA)
                strong_sample_SelSd_Ava   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
                
                FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                                  FSTfdrNS01_Ava  = NA,
                                                  FSTfdrNS02_Ava  = NA,
                                                  FSTfdrNS05_Ava  = NA,
                                                  FSTfdrNS10_Ava  = NA))
                
                locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts1)){
                    file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Ava)){
                    file.copy(from = slim_output_sample_merged_Ava, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Ava))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Ava), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Ava))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Ava), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Ava    <- as.numeric(NA)
              actual_sample_SelMean_Ava <- as.numeric(NA)
              actual_sample_SelSd_Ava   <- as.numeric(NA)
              
              strong_sample_prbe_Ava    <- as.numeric(NA)
              strong_sample_SelMean_Ava <- as.numeric(NA)
              strong_sample_SelSd_Ava   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
              add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
              global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
              
              FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                                FSTfdrNS01_Ava  = NA,
                                                FSTfdrNS02_Ava  = NA,
                                                FSTfdrNS05_Ava  = NA,
                                                FSTfdrNS10_Ava  = NA))
              
              locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts1)){
                  file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Ava)){
                  file.copy(from = slim_output_sample_merged_Ava, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Ava))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Ava))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Ava))){ #### ADD SNP INFO HERE
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Ava))
              }
            }
            
          } else {
            
            actual_sample_prbe_Ava    <- as.numeric(NA)
            actual_sample_SelMean_Ava <- as.numeric(NA)
            actual_sample_SelSd_Ava   <- as.numeric(NA)
            
            strong_sample_prbe_Ava    <- as.numeric(NA)
            strong_sample_SelMean_Ava <- as.numeric(NA)
            strong_sample_SelSd_Ava   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
            add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
            global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
            
            FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                              FSTfdrNS01_Ava  = NA,
                                              FSTfdrNS02_Ava  = NA,
                                              FSTfdrNS05_Ava  = NA,
                                              FSTfdrNS10_Ava  = NA))
            
            locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts1)){
                file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Ava)){
                file.copy(from = slim_output_sample_merged_Ava, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Ava    <- as.numeric(NA)
          actual_sample_SelMean_Ava <- as.numeric(NA)
          actual_sample_SelSd_Ava   <- as.numeric(NA)
          
          strong_sample_prbe_Ava    <- as.numeric(NA)
          strong_sample_SelMean_Ava <- as.numeric(NA)
          strong_sample_SelSd_Ava   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
          add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
          global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
          
          FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                            FSTfdrNS01_Ava  = NA,
                                            FSTfdrNS02_Ava  = NA,
                                            FSTfdrNS05_Ava  = NA,
                                            FSTfdrNS10_Ava  = NA))
          
          locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
        
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts1)){
              file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Ava)){
              file.copy(from = slim_output_sample_merged_Ava, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Ava    <- as.numeric(NA)
        actual_sample_SelMean_Ava <- as.numeric(NA)
        actual_sample_SelSd_Ava   <- as.numeric(NA)
        
        strong_sample_prbe_Ava    <- as.numeric(NA)
        strong_sample_SelMean_Ava <- as.numeric(NA)
        strong_sample_SelSd_Ava   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
        add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
        global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
        
        FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                          FSTfdrNS01_Ava  = NA,
                                          FSTfdrNS02_Ava  = NA,
                                          FSTfdrNS05_Ava  = NA,
                                          FSTfdrNS10_Ava  = NA))
        
        locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts1)){
            file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts1))
        if (file.exists(paste0(slim_output_sample_ts1_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts1_sorted, ".gz"))
        }
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts1_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts1_sorted, ".gz.csi"))
          }
        } else {
          if (file.exists(paste0(slim_output_sample_ts1_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts1_sorted, ".gz.tbi"))
          }
        }
        file.remove(paste0(slim_output_sample_merged_Ava))
      }
      
    } else {
      
      actual_sample_prbe_Ava    <- as.numeric(NA)
      actual_sample_SelMean_Ava <- as.numeric(NA)
      actual_sample_SelSd_Ava   <- as.numeric(NA)
      
      strong_sample_prbe_Ava    <- as.numeric(NA)
      strong_sample_SelMean_Ava <- as.numeric(NA)
      strong_sample_SelSd_Ava   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
      add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
      global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
      
      FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                        FSTfdrNS01_Ava  = NA,
                                        FSTfdrNS02_Ava  = NA,
                                        FSTfdrNS05_Ava  = NA,
                                        FSTfdrNS10_Ava  = NA))
      
      locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts1)){
          file.copy(from = slim_output_sample_ts1, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    
    ## Humboldt population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts2        <- paste0(slim_output_folder,"slim_output_sample_ts2_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts2, slim_output_sample_ts7)))){
      
      slim_output_sample_ts2_sorted <- paste0(slim_output_folder,"slim_output_sample_ts2_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts2_vcf <- paste("grep '^#'", slim_output_sample_ts2, ">", slim_output_sample_ts2_sorted,
                                   "&& grep -v '^#'", slim_output_sample_ts2, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts2_sorted)
      
      system(sort_sample_ts2_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts2_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts2_sorted, ".gz"))) 
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts2_sorted, ".gz")))
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      # merge and get the data
      slim_output_sample_merged_Hum <- paste0(slim_output_folder,"slim_output_sample_merged_Hum_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts2_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Hum) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Hum)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[2], by=1), "@pop1", "") ##HERE
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Hum, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(13:14)] ## remove extra inds here
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Hum_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[1]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Hum_", sim,".txt"), append = TRUE) ##HERE
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Hum_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Hum <- paste0("egglib_input_sample_Hum", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Hum), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Hum[all_merged_genome_Hum$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Hum[all_merged_genome_Hum$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Hum <- paste0("egglib_input_sample_selcoeff_Hum", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Hum), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Hum))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Hum),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Hum", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstats_Hum <- paste0(egglib_output_folder,"egglib_output_sample_Hum", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstats_Hum)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstats_Hum, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Hum    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Hum <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Hum   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Hum    <- as.numeric(0)
                  actual_sample_SelMean_Hum <- as.numeric(0)
                  actual_sample_SelSd_Hum   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Hum)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Hum*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Hum    <- sum(meanNe2_Hum*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Hum <- mean(slim_to_egglib_snps_[meanNe2_Hum*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Hum   <- sd(slim_to_egglib_snps_[meanNe2_Hum*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Hum    <- as.numeric(0)
                    positive_sample_SelMean_Hum <- as.numeric(0)
                    positive_sample_SelSd_Hum   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Hum*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Hum    <- sum(meanNe2_Hum*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Hum <- mean(slim_to_egglib_snps_[meanNe2_Hum*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Hum   <- sd(slim_to_egglib_snps_[meanNe2_Hum*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Hum    <- as.numeric(0)
                    negative_sample_SelMean_Hum <- as.numeric(0)
                    negative_sample_SelSd_Hum   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Hum)  & is.na(negative_sample_prbe_Hum)){
                    
                    strong_sample_prbe_Hum    <- as.numeric(NA)
                    strong_sample_SelMean_Hum <- as.numeric(NA)
                    strong_sample_SelSd_Hum   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Hum    <- sum(positive_sample_prbe_Hum, negative_sample_prbe_Hum, na.rm = TRUE)
                    strong_sample_SelMean_Hum <- sum(positive_sample_SelMean_Hum, negative_sample_SelMean_Hum, na.rm = TRUE)
                    strong_sample_SelSd_Hum   <- sum(positive_sample_SelSd_Hum, negative_sample_SelSd_Hum, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Hum    <- as.numeric(NA)
                  strong_sample_SelMean_Hum <- as.numeric(NA)
                  strong_sample_SelSd_Hum   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Hum*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Hum*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Hum = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Hum = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Hum = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Hum = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Hum = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Hum = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Hum = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Hum = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Hum = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Hum = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum,
                                                      FSTfdrNS01_Hum,
                                                      FSTfdrNS02_Hum,
                                                      FSTfdrNS05_Hum,
                                                      FSTfdrNS10_Hum))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Hum  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Hum  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Hum  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Hum  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Hum = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Hum = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Hum = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Hum = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[2] + (SSs[7]-2) + 6)] #HERE
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Hum <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Hum", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Hum    <- as.numeric(NA)
                actual_sample_SelMean_Hum <- as.numeric(NA)
                actual_sample_SelSd_Hum   <- as.numeric(NA)
                
                strong_sample_prbe_Hum    <- as.numeric(NA)
                strong_sample_SelMean_Hum <- as.numeric(NA)
                strong_sample_SelSd_Hum   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
                
                FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                                  FSTfdrNS01_Hum  = NA,
                                                  FSTfdrNS02_Hum  = NA,
                                                  FSTfdrNS05_Hum  = NA,
                                                  FSTfdrNS10_Hum  = NA))
                
                locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts2)){
                    file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Hum)){
                    file.copy(from = slim_output_sample_merged_Hum, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Hum))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Hum), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Hum))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Hum), to = debug_output_folder)
                  }
                
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
              
              actual_sample_prbe_Hum    <- as.numeric(NA)
              actual_sample_SelMean_Hum <- as.numeric(NA)
              actual_sample_SelSd_Hum   <- as.numeric(NA)
              
              strong_sample_prbe_Hum    <- as.numeric(NA)
              strong_sample_SelMean_Hum <- as.numeric(NA)
              strong_sample_SelSd_Hum   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
              add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
              global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
              
              FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                                FSTfdrNS01_Hum  = NA,
                                                FSTfdrNS02_Hum  = NA,
                                                FSTfdrNS05_Hum  = NA,
                                                FSTfdrNS10_Hum  = NA))
              
              locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts2)){
                  file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Hum)){
                  file.copy(from = slim_output_sample_merged_Hum, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Hum))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Hum))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Hum))){ #### ADD SNP INFO HERE
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Hum))
              }
            }
            
          } else {
            
            actual_sample_prbe_Hum    <- as.numeric(NA)
            actual_sample_SelMean_Hum <- as.numeric(NA)
            actual_sample_SelSd_Hum   <- as.numeric(NA)
            
            strong_sample_prbe_Hum    <- as.numeric(NA)
            strong_sample_SelMean_Hum <- as.numeric(NA)
            strong_sample_SelSd_Hum   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
            add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
            global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
            
            FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                              FSTfdrNS01_Hum  = NA,
                                              FSTfdrNS02_Hum  = NA,
                                              FSTfdrNS05_Hum  = NA,
                                              FSTfdrNS10_Hum  = NA))
            
            locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts2)){
                file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Hum)){
                file.copy(from = slim_output_sample_merged_Hum, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Hum    <- as.numeric(NA)
          actual_sample_SelMean_Hum <- as.numeric(NA)
          actual_sample_SelSd_Hum   <- as.numeric(NA)
          
          strong_sample_prbe_Hum    <- as.numeric(NA)
          strong_sample_SelMean_Hum <- as.numeric(NA)
          strong_sample_SelSd_Hum   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
          add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
          global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
          
          FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                            FSTfdrNS01_Hum  = NA,
                                            FSTfdrNS02_Hum  = NA,
                                            FSTfdrNS05_Hum  = NA,
                                            FSTfdrNS10_Hum  = NA))
          
          locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts2)){
              file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Hum)){
              file.copy(from = slim_output_sample_merged_Hum, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Hum    <- as.numeric(NA)
        actual_sample_SelMean_Hum <- as.numeric(NA)
        actual_sample_SelSd_Hum   <- as.numeric(NA)
        
        strong_sample_prbe_Hum    <- as.numeric(NA)
        strong_sample_SelMean_Hum <- as.numeric(NA)
        strong_sample_SelSd_Hum   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
        add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
        global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
        
        FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                          FSTfdrNS01_Hum  = NA,
                                          FSTfdrNS02_Hum  = NA,
                                          FSTfdrNS05_Hum  = NA,
                                          FSTfdrNS10_Hum  = NA))
        
        locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts2)){
            file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts2))
        if (file.exists(paste0(slim_output_sample_ts2_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts2_sorted, ".gz"))
        }
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts2_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts2_sorted, ".gz.csi"))
          }
        } else {
          if (file.exists(paste0(slim_output_sample_ts2_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts2_sorted, ".gz.tbi"))
          }
        }
        file.remove(paste0(slim_output_sample_merged_Hum))
      }
      
    } else {
      
      actual_sample_prbe_Hum    <- as.numeric(NA)
      actual_sample_SelMean_Hum <- as.numeric(NA)
      actual_sample_SelSd_Hum   <- as.numeric(NA)
      
      strong_sample_prbe_Hum    <- as.numeric(NA)
      strong_sample_SelMean_Hum <- as.numeric(NA)
      strong_sample_SelSd_Hum   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
      add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
      global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
      
      FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                        FSTfdrNS01_Hum  = NA,
                                        FSTfdrNS02_Hum  = NA,
                                        FSTfdrNS05_Hum  = NA,
                                        FSTfdrNS10_Hum  = NA))
      
      locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts2)){
          file.copy(from = slim_output_sample_ts2, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    ## Davis population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts3        <- paste0(slim_output_folder,"slim_output_sample_ts3_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts3, slim_output_sample_ts7)))){
      
      slim_output_sample_ts3_sorted <- paste0(slim_output_folder,"slim_output_sample_ts3_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts3_vcf <- paste("grep '^#'", slim_output_sample_ts3, ">", slim_output_sample_ts3_sorted,
                                   "&& grep -v '^#'", slim_output_sample_ts3, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts3_sorted)
      
      system(sort_sample_ts3_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts3_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts3_sorted, ".gz"))) 
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts3_sorted, ".gz")))
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      # merge and get the data
      slim_output_sample_merged_Dav <- paste0(slim_output_folder,"slim_output_sample_merged_Dav_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts3_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Dav) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Dav)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[3], by=1), "@pop1", "") ##HERE
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Dav, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(9:10)] ## remove extra inds here
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Dav_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[2]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Dav_", sim,".txt"), append = TRUE) ##HERE
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Dav_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Dav <- paste0("egglib_input_sample_Dav", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Dav), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Dav[all_merged_genome_Dav$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Dav[all_merged_genome_Dav$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Dav <- paste0("egglib_input_sample_selcoeff_Dav", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Dav), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Dav))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Dav),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Dav", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstats_Dav <- paste0(egglib_output_folder,"egglib_output_sample_Dav", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstats_Dav)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstats_Dav, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Dav    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Dav <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Dav   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Dav    <- as.numeric(0)
                  actual_sample_SelMean_Dav <- as.numeric(0)
                  actual_sample_SelSd_Dav   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Dav)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Dav*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Dav    <- sum(meanNe2_Dav*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Dav <- mean(slim_to_egglib_snps_[meanNe2_Dav*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Dav   <- sd(slim_to_egglib_snps_[meanNe2_Dav*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Dav    <- as.numeric(0)
                    positive_sample_SelMean_Dav <- as.numeric(0)
                    positive_sample_SelSd_Dav   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Dav*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Dav    <- sum(meanNe2_Dav*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Dav <- mean(slim_to_egglib_snps_[meanNe2_Dav*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Dav   <- sd(slim_to_egglib_snps_[meanNe2_Dav*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Dav    <- as.numeric(0)
                    negative_sample_SelMean_Dav <- as.numeric(0)
                    negative_sample_SelSd_Dav   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Dav)  & is.na(negative_sample_prbe_Dav)){
                    
                    strong_sample_prbe_Dav    <- as.numeric(NA)
                    strong_sample_SelMean_Dav <- as.numeric(NA)
                    strong_sample_SelSd_Dav   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Dav    <- sum(positive_sample_prbe_Dav, negative_sample_prbe_Dav, na.rm = TRUE)
                    strong_sample_SelMean_Dav <- sum(positive_sample_SelMean_Dav, negative_sample_SelMean_Dav, na.rm = TRUE)
                    strong_sample_SelSd_Dav   <- sum(positive_sample_SelSd_Dav, negative_sample_SelSd_Dav, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Dav    <- as.numeric(NA)
                  strong_sample_SelMean_Dav <- as.numeric(NA)
                  strong_sample_SelSd_Dav   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Dav*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Dav*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Dav = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Dav = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Dav = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Dav = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Dav = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Dav = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Dav = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Dav = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Dav = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Dav = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav,
                                                      FSTfdrNS01_Dav,
                                                      FSTfdrNS02_Dav,
                                                      FSTfdrNS05_Dav,
                                                      FSTfdrNS10_Dav))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Dav  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Dav  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Dav  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Dav  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Dav = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Dav = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Dav = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Dav = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[3] + (SSs[7]-2) + 6)] #HERE
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Dav <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Dav", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Dav    <- as.numeric(NA)
                actual_sample_SelMean_Dav <- as.numeric(NA)
                actual_sample_SelSd_Dav   <- as.numeric(NA)
                
                strong_sample_prbe_Dav    <- as.numeric(NA)
                strong_sample_SelMean_Dav <- as.numeric(NA)
                strong_sample_SelSd_Dav   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
                
                FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                                  FSTfdrNS01_Dav  = NA,
                                                  FSTfdrNS02_Dav  = NA,
                                                  FSTfdrNS05_Dav  = NA,
                                                  FSTfdrNS10_Dav  = NA))
                
                locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts3)){
                    file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Dav)){
                    file.copy(from = slim_output_sample_merged_Dav, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Dav))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Dav), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Dav))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Dav), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Dav    <- as.numeric(NA)
              actual_sample_SelMean_Dav <- as.numeric(NA)
              actual_sample_SelSd_Dav   <- as.numeric(NA)
              
              strong_sample_prbe_Dav    <- as.numeric(NA)
              strong_sample_SelMean_Dav <- as.numeric(NA)
              strong_sample_SelSd_Dav   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
              add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
              global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
              
              FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                                FSTfdrNS01_Dav  = NA,
                                                FSTfdrNS02_Dav  = NA,
                                                FSTfdrNS05_Dav  = NA,
                                                FSTfdrNS10_Dav  = NA))
              
              locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts3)){
                  file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Dav)){
                  file.copy(from = slim_output_sample_merged_Dav, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Dav))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Dav))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Dav))){
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Dav))
              }
            }
            
          } else {
            
            actual_sample_prbe_Dav    <- as.numeric(NA)
            actual_sample_SelMean_Dav <- as.numeric(NA)
            actual_sample_SelSd_Dav   <- as.numeric(NA)
            
            strong_sample_prbe_Dav    <- as.numeric(NA)
            strong_sample_SelMean_Dav <- as.numeric(NA)
            strong_sample_SelSd_Dav   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
            add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
            global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
            
            FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                              FSTfdrNS01_Dav  = NA,
                                              FSTfdrNS02_Dav  = NA,
                                              FSTfdrNS05_Dav  = NA,
                                              FSTfdrNS10_Dav  = NA))
            
            locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts3)){
                file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Dav)){
                file.copy(from = slim_output_sample_merged_Dav, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Dav    <- as.numeric(NA)
          actual_sample_SelMean_Dav <- as.numeric(NA)
          actual_sample_SelSd_Dav   <- as.numeric(NA)
          
          strong_sample_prbe_Dav    <- as.numeric(NA)
          strong_sample_SelMean_Dav <- as.numeric(NA)
          strong_sample_SelSd_Dav   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
          add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
          global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
          
          FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                            FSTfdrNS01_Dav  = NA,
                                            FSTfdrNS02_Dav  = NA,
                                            FSTfdrNS05_Dav  = NA,
                                            FSTfdrNS10_Dav  = NA))
          
          locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts3)){
              file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Dav)){
              file.copy(from = slim_output_sample_merged_Dav, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Dav    <- as.numeric(NA)
        actual_sample_SelMean_Dav <- as.numeric(NA)
        actual_sample_SelSd_Dav   <- as.numeric(NA)
        
        strong_sample_prbe_Dav    <- as.numeric(NA)
        strong_sample_SelMean_Dav <- as.numeric(NA)
        strong_sample_SelSd_Dav   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
        add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
        global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
        
        FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                          FSTfdrNS01_Dav  = NA,
                                          FSTfdrNS02_Dav  = NA,
                                          FSTfdrNS05_Dav  = NA,
                                          FSTfdrNS10_Dav  = NA))
        
        locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts3)){
            file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts3))
        if (file.exists(paste0(slim_output_sample_ts3_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts3_sorted, ".gz"))
        }
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts3_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts3_sorted, ".gz.csi"))
          }
        } else {
          if (file.exists(paste0(slim_output_sample_ts3_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts3_sorted, ".gz.tbi"))
          }
        }
        file.remove(paste0(slim_output_sample_merged_Dav))
      }
      
    } else {
      
      actual_sample_prbe_Dav    <- as.numeric(NA)
      actual_sample_SelMean_Dav <- as.numeric(NA)
      actual_sample_SelSd_Dav   <- as.numeric(NA)
      
      strong_sample_prbe_Dav    <- as.numeric(NA)
      strong_sample_SelMean_Dav <- as.numeric(NA)
      strong_sample_SelSd_Dav   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
      add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
      global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
      
      FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                        FSTfdrNS01_Dav  = NA,
                                        FSTfdrNS02_Dav  = NA,
                                        FSTfdrNS05_Dav  = NA,
                                        FSTfdrNS10_Dav  = NA))
      
      locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts3)){
          file.copy(from = slim_output_sample_ts3, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    ## Stanislaus population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts4        <- paste0(slim_output_folder,"slim_output_sample_ts4_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts4, slim_output_sample_ts7)))){
      
      slim_output_sample_ts4_sorted <- paste0(slim_output_folder,"slim_output_sample_ts4_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts4_vcf <- paste("grep '^#'", slim_output_sample_ts4, ">", slim_output_sample_ts4_sorted,
                                   "&& grep -v '^#'", slim_output_sample_ts4, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts4_sorted)
      
      system(sort_sample_ts4_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts4_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts4_sorted, ".gz"))) 
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts4_sorted, ".gz")))
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      # merge and get the data
      slim_output_sample_merged_Sta <- paste0(slim_output_folder,"slim_output_sample_merged_Sta_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts4_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Sta) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Sta)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[4], by=1), "@pop1", "") ##HERE
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Sta, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(9:10)] ## remove extra inds here
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Sta_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[3]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Sta_", sim,".txt"), append = TRUE) ##HERE
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Sta_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Sta <- paste0("egglib_input_sample_Sta", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Sta), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Sta[all_merged_genome_Sta$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Sta[all_merged_genome_Sta$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Sta <- paste0("egglib_input_sample_selcoeff_Sta", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Sta), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Sta))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Sta),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Sta", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstats_Sta <- paste0(egglib_output_folder,"egglib_output_sample_Sta", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstats_Sta)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstats_Sta, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Sta    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Sta <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Sta   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Sta    <- as.numeric(0)
                  actual_sample_SelMean_Sta <- as.numeric(0)
                  actual_sample_SelSd_Sta   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Sta)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Sta*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Sta    <- sum(meanNe2_Sta*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Sta <- mean(slim_to_egglib_snps_[meanNe2_Sta*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Sta   <- sd(slim_to_egglib_snps_[meanNe2_Sta*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Sta    <- as.numeric(0)
                    positive_sample_SelMean_Sta <- as.numeric(0)
                    positive_sample_SelSd_Sta   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Sta*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Sta    <- sum(meanNe2_Sta*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Sta <- mean(slim_to_egglib_snps_[meanNe2_Sta*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Sta   <- sd(slim_to_egglib_snps_[meanNe2_Sta*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Sta    <- as.numeric(0)
                    negative_sample_SelMean_Sta <- as.numeric(0)
                    negative_sample_SelSd_Sta   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Sta)  & is.na(negative_sample_prbe_Sta)){
                    
                    strong_sample_prbe_Sta    <- as.numeric(NA)
                    strong_sample_SelMean_Sta <- as.numeric(NA)
                    strong_sample_SelSd_Sta   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Sta    <- sum(positive_sample_prbe_Sta, negative_sample_prbe_Sta, na.rm = TRUE)
                    strong_sample_SelMean_Sta <- sum(positive_sample_SelMean_Sta, negative_sample_SelMean_Sta, na.rm = TRUE)
                    strong_sample_SelSd_Sta   <- sum(positive_sample_SelSd_Sta, negative_sample_SelSd_Sta, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Sta    <- as.numeric(NA)
                  strong_sample_SelMean_Sta <- as.numeric(NA)
                  strong_sample_SelSd_Sta   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Sta*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Sta*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Sta = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Sta = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Sta = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Sta = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Sta = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Sta = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Sta = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Sta = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Sta = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Sta = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta,
                                                      FSTfdrNS01_Sta,
                                                      FSTfdrNS02_Sta,
                                                      FSTfdrNS05_Sta,
                                                      FSTfdrNS10_Sta))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Sta  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Sta  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Sta  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Sta  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Sta = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Sta = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Sta = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Sta = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[4] + (SSs[7]-2) + 6)] #HERE
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Sta <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Sta", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Sta    <- as.numeric(NA)
                actual_sample_SelMean_Sta <- as.numeric(NA)
                actual_sample_SelSd_Sta   <- as.numeric(NA)
                
                strong_sample_prbe_Sta    <- as.numeric(NA)
                strong_sample_SelMean_Sta <- as.numeric(NA)
                strong_sample_SelSd_Sta   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
                
                FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                                  FSTfdrNS01_Sta  = NA,
                                                  FSTfdrNS02_Sta  = NA,
                                                  FSTfdrNS05_Sta  = NA,
                                                  FSTfdrNS10_Sta  = NA))
                
                locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts4)){
                    file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Sta)){
                    file.copy(from = slim_output_sample_merged_Sta, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Sta))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Sta), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Sta))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Sta), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Sta    <- as.numeric(NA)
              actual_sample_SelMean_Sta <- as.numeric(NA)
              actual_sample_SelSd_Sta   <- as.numeric(NA)
              
              strong_sample_prbe_Sta    <- as.numeric(NA)
              strong_sample_SelMean_Sta <- as.numeric(NA)
              strong_sample_SelSd_Sta   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
              add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
              global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
              
              FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                                FSTfdrNS01_Sta  = NA,
                                                FSTfdrNS02_Sta  = NA,
                                                FSTfdrNS05_Sta  = NA,
                                                FSTfdrNS10_Sta  = NA))
              
              locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts4)){
                  file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Sta)){
                  file.copy(from = slim_output_sample_merged_Sta, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Sta))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Sta))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Sta))){
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Sta))
              }
            }
            
          } else {
            
            actual_sample_prbe_Sta    <- as.numeric(NA)
            actual_sample_SelMean_Sta <- as.numeric(NA)
            actual_sample_SelSd_Sta   <- as.numeric(NA)
            
            strong_sample_prbe_Sta    <- as.numeric(NA)
            strong_sample_SelMean_Sta <- as.numeric(NA)
            strong_sample_SelSd_Sta   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
            add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
            global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
            
            FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                              FSTfdrNS01_Sta  = NA,
                                              FSTfdrNS02_Sta  = NA,
                                              FSTfdrNS05_Sta  = NA,
                                              FSTfdrNS10_Sta  = NA))
            
            locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts4)){
                file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Sta)){
                file.copy(from = slim_output_sample_merged_Sta, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Sta    <- as.numeric(NA)
          actual_sample_SelMean_Sta <- as.numeric(NA)
          actual_sample_SelSd_Sta   <- as.numeric(NA)
          
          strong_sample_prbe_Sta    <- as.numeric(NA)
          strong_sample_SelMean_Sta <- as.numeric(NA)
          strong_sample_SelSd_Sta   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
          add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
          global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
          
          FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                            FSTfdrNS01_Sta  = NA,
                                            FSTfdrNS02_Sta  = NA,
                                            FSTfdrNS05_Sta  = NA,
                                            FSTfdrNS10_Sta  = NA))
          
          locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts4)){
              file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Sta)){
              file.copy(from = slim_output_sample_merged_Sta, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Sta    <- as.numeric(NA)
        actual_sample_SelMean_Sta <- as.numeric(NA)
        actual_sample_SelSd_Sta   <- as.numeric(NA)
        
        strong_sample_prbe_Sta    <- as.numeric(NA)
        strong_sample_SelMean_Sta <- as.numeric(NA)
        strong_sample_SelSd_Sta   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
        add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
        global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
        
        FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                          FSTfdrNS01_Sta  = NA,
                                          FSTfdrNS02_Sta  = NA,
                                          FSTfdrNS05_Sta  = NA,
                                          FSTfdrNS10_Sta  = NA))
        
        locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts4)){
            file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts4))
        
        if (file.exists(paste0(slim_output_sample_ts4_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts4_sorted, ".gz"))
        }
        
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts4_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts4_sorted, ".gz.csi"))
          }
          
        } else {
          if (file.exists(paste0(slim_output_sample_ts4_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts4_sorted, ".gz.tbi"))
          }
        }
        file.remove(paste0(slim_output_sample_merged_Sta))
      }
      
    } else {
      
      actual_sample_prbe_Sta    <- as.numeric(NA)
      actual_sample_SelMean_Sta <- as.numeric(NA)
      actual_sample_SelSd_Sta   <- as.numeric(NA)
      
      strong_sample_prbe_Sta    <- as.numeric(NA)
      strong_sample_SelMean_Sta <- as.numeric(NA)
      strong_sample_SelSd_Sta   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
      add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
      global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
      
      FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                        FSTfdrNS01_Sta  = NA,
                                        FSTfdrNS02_Sta  = NA,
                                        FSTfdrNS05_Sta  = NA,
                                        FSTfdrNS10_Sta  = NA))
      
      locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts4)){
          file.copy(from = slim_output_sample_ts4, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    ## Stebbins population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts5        <- paste0(slim_output_folder,"slim_output_sample_ts5_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts5, slim_output_sample_ts7)))){
      
      slim_output_sample_ts5_sorted <- paste0(slim_output_folder,"slim_output_sample_ts5_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts5_vcf <- paste("grep '^#'", slim_output_sample_ts5, ">", slim_output_sample_ts5_sorted,
                                   "&& grep -v '^#'", slim_output_sample_ts5, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts5_sorted)
      
      system(sort_sample_ts5_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts5_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts5_sorted, ".gz"))) 
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts5_sorted, ".gz")))
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      # merge and get the data
      slim_output_sample_merged_Ste <- paste0(slim_output_folder,"slim_output_sample_merged_Ste_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts5_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Ste) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Ste)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[5], by=1), "@pop1", "") ##HERE
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Ste, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(11:13)] ## remove extra inds here
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Ste_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[4]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Ste_", sim,".txt"), append = TRUE) ##HERE
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Ste_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Ste <- paste0("egglib_input_sample_Ste", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Ste), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Ste[all_merged_genome_Ste$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Ste[all_merged_genome_Ste$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Ste <- paste0("egglib_input_sample_selcoeff_Ste", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Ste), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Ste))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Ste),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Ste", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstats_Ste <- paste0(egglib_output_folder,"egglib_output_sample_Ste", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstats_Ste)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstats_Ste, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Ste    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Ste <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Ste   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Ste    <- as.numeric(0)
                  actual_sample_SelMean_Ste <- as.numeric(0)
                  actual_sample_SelSd_Ste   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Ste)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Ste*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Ste    <- sum(meanNe2_Ste*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Ste <- mean(slim_to_egglib_snps_[meanNe2_Ste*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Ste   <- sd(slim_to_egglib_snps_[meanNe2_Ste*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Ste    <- as.numeric(0)
                    positive_sample_SelMean_Ste <- as.numeric(0)
                    positive_sample_SelSd_Ste   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Ste*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Ste    <- sum(meanNe2_Ste*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Ste <- mean(slim_to_egglib_snps_[meanNe2_Ste*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Ste   <- sd(slim_to_egglib_snps_[meanNe2_Ste*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Ste    <- as.numeric(0)
                    negative_sample_SelMean_Ste <- as.numeric(0)
                    negative_sample_SelSd_Ste   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Ste)  & is.na(negative_sample_prbe_Ste)){
                    
                    strong_sample_prbe_Ste    <- as.numeric(NA)
                    strong_sample_SelMean_Ste <- as.numeric(NA)
                    strong_sample_SelSd_Ste   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Ste    <- sum(positive_sample_prbe_Ste, negative_sample_prbe_Ste, na.rm = TRUE)
                    strong_sample_SelMean_Ste <- sum(positive_sample_SelMean_Ste, negative_sample_SelMean_Ste, na.rm = TRUE)
                    strong_sample_SelSd_Ste   <- sum(positive_sample_SelSd_Ste, negative_sample_SelSd_Ste, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Ste    <- as.numeric(NA)
                  strong_sample_SelMean_Ste <- as.numeric(NA)
                  strong_sample_SelSd_Ste   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Ste*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Ste*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Ste = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Ste = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Ste = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Ste = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Ste = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Ste = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Ste = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Ste = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Ste = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Ste = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste,
                                                      FSTfdrNS01_Ste,
                                                      FSTfdrNS02_Ste,
                                                      FSTfdrNS05_Ste,
                                                      FSTfdrNS10_Ste))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Ste  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Ste  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Ste  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Ste  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Ste = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Ste = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Ste = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Ste = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[5] + (SSs[7]-3) + 6)] #HERE
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Ste <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Ste", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Ste    <- as.numeric(NA)
                actual_sample_SelMean_Ste <- as.numeric(NA)
                actual_sample_SelSd_Ste   <- as.numeric(NA)
                
                strong_sample_prbe_Ste    <- as.numeric(NA)
                strong_sample_SelMean_Ste <- as.numeric(NA)
                strong_sample_SelSd_Ste   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
                
                FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                                  FSTfdrNS01_Ste  = NA,
                                                  FSTfdrNS02_Ste  = NA,
                                                  FSTfdrNS05_Ste  = NA,
                                                  FSTfdrNS10_Ste  = NA))
                
                locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts5)){
                    file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Ste)){
                    file.copy(from = slim_output_sample_merged_Ste, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Ste))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Ste), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Ste))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Ste), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Ste    <- as.numeric(NA)
              actual_sample_SelMean_Ste <- as.numeric(NA)
              actual_sample_SelSd_Ste   <- as.numeric(NA)
              
              strong_sample_prbe_Ste    <- as.numeric(NA)
              strong_sample_SelMean_Ste <- as.numeric(NA)
              strong_sample_SelSd_Ste   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
              add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
              global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
              
              FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                                FSTfdrNS01_Ste  = NA,
                                                FSTfdrNS02_Ste  = NA,
                                                FSTfdrNS05_Ste  = NA,
                                                FSTfdrNS10_Ste  = NA))
              
              locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts5)){
                  file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Ste)){
                  file.copy(from = slim_output_sample_merged_Ste, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Ste))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Ste))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Ste))){
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Ste))
              }
            }
            
          } else {
            
            actual_sample_prbe_Ste    <- as.numeric(NA)
            actual_sample_SelMean_Ste <- as.numeric(NA)
            actual_sample_SelSd_Ste   <- as.numeric(NA)
            
            strong_sample_prbe_Ste    <- as.numeric(NA)
            strong_sample_SelMean_Ste <- as.numeric(NA)
            strong_sample_SelSd_Ste   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
            add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
            global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
            
            FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                              FSTfdrNS01_Ste  = NA,
                                              FSTfdrNS02_Ste  = NA,
                                              FSTfdrNS05_Ste  = NA,
                                              FSTfdrNS10_Ste  = NA))
            
            locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts5)){
                file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Ste)){
                file.copy(from = slim_output_sample_merged_Ste, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Ste    <- as.numeric(NA)
          actual_sample_SelMean_Ste <- as.numeric(NA)
          actual_sample_SelSd_Ste   <- as.numeric(NA)
          
          strong_sample_prbe_Ste    <- as.numeric(NA)
          strong_sample_SelMean_Ste <- as.numeric(NA)
          strong_sample_SelSd_Ste   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
          add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
          global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
          
          FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                            FSTfdrNS01_Ste  = NA,
                                            FSTfdrNS02_Ste  = NA,
                                            FSTfdrNS05_Ste  = NA,
                                            FSTfdrNS10_Ste  = NA))
          
          locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts5)){
              file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Ste)){
              file.copy(from = slim_output_sample_merged_Ste, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Ste    <- as.numeric(NA)
        actual_sample_SelMean_Ste <- as.numeric(NA)
        actual_sample_SelSd_Ste   <- as.numeric(NA)
        
        strong_sample_prbe_Ste    <- as.numeric(NA)
        strong_sample_SelMean_Ste <- as.numeric(NA)
        strong_sample_SelSd_Ste   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
        add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
        global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
        
        FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                          FSTfdrNS01_Ste  = NA,
                                          FSTfdrNS02_Ste  = NA,
                                          FSTfdrNS05_Ste  = NA,
                                          FSTfdrNS10_Ste  = NA))
        
        locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts5)){
            file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts5))
        
        if (file.exists(paste0(slim_output_sample_ts5_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts5_sorted, ".gz"))
        }
        
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts5_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts5_sorted, ".gz.csi"))
          }
          
        } else {
          if (file.exists(paste0(slim_output_sample_ts5_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts5_sorted, ".gz.tbi"))
          }
        }
        file.remove(paste0(slim_output_sample_merged_Ste))
      }
      
    } else {
      
      actual_sample_prbe_Ste    <- as.numeric(NA)
      actual_sample_SelMean_Ste <- as.numeric(NA)
      actual_sample_SelSd_Ste   <- as.numeric(NA)
      
      strong_sample_prbe_Ste    <- as.numeric(NA)
      strong_sample_SelMean_Ste <- as.numeric(NA)
      strong_sample_SelSd_Ste   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
      add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
      global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
      
      FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                        FSTfdrNS01_Ste  = NA,
                                        FSTfdrNS02_Ste  = NA,
                                        FSTfdrNS05_Ste  = NA,
                                        FSTfdrNS10_Ste  = NA))
      
      locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts5)){
          file.copy(from = slim_output_sample_ts5, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    ## Riverside population
    ##-----------------
    
    # sort vcf files
    slim_output_sample_ts6        <- paste0(slim_output_folder,"slim_output_sample_ts6_", sim, ".vcf")
    
    if(all(file.exists(c(slim_output_sample_ts6, slim_output_sample_ts7)))){
      
      slim_output_sample_ts6_sorted <- paste0(slim_output_folder,"slim_output_sample_ts6_", sim, "_sorted" , ".vcf")
      
      sort_sample_ts6_vcf <- paste("grep '^#'", slim_output_sample_ts6, ">", slim_output_sample_ts6_sorted,
                                   "&& grep -v '^#'", slim_output_sample_ts6, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts6_sorted)
      
      system(sort_sample_ts6_vcf)
      
      # bgzip sorted vcf files
      system(paste(path_to_bgzip, "-f", slim_output_sample_ts6_sorted))
      
      # tabix bgziped files
      if (genomeS > 2^29){
        # csi instead of tbi for large chromosome 
        system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts6_sorted, ".gz"))) 
      } else {
        system(paste(path_to_tabix, paste0(slim_output_sample_ts6_sorted, ".gz")))
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      # merge and get the data
      slim_output_sample_merged_Riv <- paste0(slim_output_folder,"slim_output_sample_merged_Riv_", sim, ".txt")
      
      bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                               paste0(slim_output_sample_ts6_sorted, ".gz"),
                               paste0(slim_output_sample_ts7_sorted, ".gz"),
                               "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                               ">", slim_output_sample_merged_Riv) 
      
      system(bcftools_query_)
      
      if(file.exists(slim_output_sample_merged_Riv)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[6], by=1), "@pop1", "")
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Riv, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(3:5)]
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Riv_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[5]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Riv_", sim,".txt"), append = TRUE)
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Riv_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Riv <- paste0("egglib_input_sample_Riv", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Riv), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Riv[all_merged_genome_Riv$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Riv[all_merged_genome_Riv$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Riv <- paste0("egglib_input_sample_selcoeff_Riv", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Riv), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Riv))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Riv),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Riv", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstat_Riv <- paste0(egglib_output_folder,"egglib_output_sample_Riv", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstat_Riv)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstat_Riv, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Riv    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Riv <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Riv   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Riv    <- as.numeric(0)
                  actual_sample_SelMean_Riv <- as.numeric(0)
                  actual_sample_SelSd_Riv   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Riv)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Riv*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Riv    <- sum(meanNe2_Riv*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Riv <- mean(slim_to_egglib_snps_[meanNe2_Riv*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Riv   <- sd(slim_to_egglib_snps_[meanNe2_Riv*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Riv    <- as.numeric(0)
                    positive_sample_SelMean_Riv <- as.numeric(0)
                    positive_sample_SelSd_Riv   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Riv*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Riv    <- sum(meanNe2_Riv*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Riv <- mean(slim_to_egglib_snps_[meanNe2_Riv*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Riv   <- sd(slim_to_egglib_snps_[meanNe2_Riv*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Riv    <- as.numeric(0)
                    negative_sample_SelMean_Riv <- as.numeric(0)
                    negative_sample_SelSd_Riv   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Riv)  & is.na(negative_sample_prbe_Riv)){
                    
                    strong_sample_prbe_Riv    <- as.numeric(NA)
                    strong_sample_SelMean_Riv <- as.numeric(NA)
                    strong_sample_SelSd_Riv   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Riv    <- sum(positive_sample_prbe_Riv, negative_sample_prbe_Riv, na.rm = TRUE)
                    strong_sample_SelMean_Riv <- sum(positive_sample_SelMean_Riv, negative_sample_SelMean_Riv, na.rm = TRUE)
                    strong_sample_SelSd_Riv   <- sum(positive_sample_SelSd_Riv, negative_sample_SelSd_Riv, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Riv    <- as.numeric(NA)
                  strong_sample_SelMean_Riv <- as.numeric(NA)
                  strong_sample_SelSd_Riv   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Riv*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Riv*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Riv = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Riv = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Riv = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Riv = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Riv = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Riv = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Riv = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Riv = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Riv = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Riv = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv,
                                                      FSTfdrNS01_Riv,
                                                      FSTfdrNS02_Riv,
                                                      FSTfdrNS05_Riv,
                                                      FSTfdrNS10_Riv))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Riv  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Riv  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Riv  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Riv  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Riv = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Riv = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Riv = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Riv = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):((SSs[6]-3) + SSs[7] + 6)]
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Riv <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Riv", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Riv    <- as.numeric(NA)
                actual_sample_SelMean_Riv <- as.numeric(NA)
                actual_sample_SelSd_Riv   <- as.numeric(NA)
                
                strong_sample_prbe_Riv    <- as.numeric(NA)
                strong_sample_SelMean_Riv <- as.numeric(NA)
                strong_sample_SelSd_Riv   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
                
                FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                                  FSTfdrNS01_Riv  = NA,
                                                  FSTfdrNS02_Riv  = NA,
                                                  FSTfdrNS05_Riv  = NA,
                                                  FSTfdrNS10_Riv  = NA))
                
                locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts6)){
                    file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Riv)){
                    file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Riv))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Riv), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Riv))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Riv), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Riv    <- as.numeric(NA)
              actual_sample_SelMean_Riv <- as.numeric(NA)
              actual_sample_SelSd_Riv   <- as.numeric(NA)
              
              strong_sample_prbe_Riv    <- as.numeric(NA)
              strong_sample_SelMean_Riv <- as.numeric(NA)
              strong_sample_SelSd_Riv   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
              add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
              global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
              
              FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                                FSTfdrNS01_Riv  = NA,
                                                FSTfdrNS02_Riv  = NA,
                                                FSTfdrNS05_Riv  = NA,
                                                FSTfdrNS10_Riv  = NA))
              
              locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts6)){
                  file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Riv)){
                  file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Riv))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Riv))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Riv))){
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Riv))
              }
            }
            
          } else {
            
            actual_sample_prbe_Riv    <- as.numeric(NA)
            actual_sample_SelMean_Riv <- as.numeric(NA)
            actual_sample_SelSd_Riv   <- as.numeric(NA)
            
            strong_sample_prbe_Riv    <- as.numeric(NA)
            strong_sample_SelMean_Riv <- as.numeric(NA)
            strong_sample_SelSd_Riv   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
            add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
            global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
            
            FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                              FSTfdrNS01_Riv  = NA,
                                              FSTfdrNS02_Riv  = NA,
                                              FSTfdrNS05_Riv  = NA,
                                              FSTfdrNS10_Riv  = NA))
            
            locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts6)){
                file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Riv)){
                file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Riv    <- as.numeric(NA)
          actual_sample_SelMean_Riv <- as.numeric(NA)
          actual_sample_SelSd_Riv   <- as.numeric(NA)
          
          strong_sample_prbe_Riv    <- as.numeric(NA)
          strong_sample_SelMean_Riv <- as.numeric(NA)
          strong_sample_SelSd_Riv   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
          add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
          global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
          
          FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                            FSTfdrNS01_Riv  = NA,
                                            FSTfdrNS02_Riv  = NA,
                                            FSTfdrNS05_Riv  = NA,
                                            FSTfdrNS10_Riv  = NA))
          
          locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts6)){
              file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Riv)){
              file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Riv    <- as.numeric(NA)
        actual_sample_SelMean_Riv <- as.numeric(NA)
        actual_sample_SelSd_Riv   <- as.numeric(NA)
        
        strong_sample_prbe_Riv    <- as.numeric(NA)
        strong_sample_SelMean_Riv <- as.numeric(NA)
        strong_sample_SelSd_Riv   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
        add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
        global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
        
        FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                          FSTfdrNS01_Riv  = NA,
                                          FSTfdrNS02_Riv  = NA,
                                          FSTfdrNS05_Riv  = NA,
                                          FSTfdrNS10_Riv  = NA))
        
        locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts6)){
            file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      actual_sample_prbe_Riv    <- as.numeric(NA)
      actual_sample_SelMean_Riv <- as.numeric(NA)
      actual_sample_SelSd_Riv   <- as.numeric(NA)
      
      strong_sample_prbe_Riv    <- as.numeric(NA)
      strong_sample_SelMean_Riv <- as.numeric(NA)
      strong_sample_SelSd_Riv   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
      add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
      global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
      
      FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                        FSTfdrNS01_Riv  = NA,
                                        FSTfdrNS02_Riv  = NA,
                                        FSTfdrNS05_Riv  = NA,
                                        FSTfdrNS10_Riv  = NA))
      
      locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts6)){
          file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    
    ## Placerita population
    ##-----------------
    
    if(all(file.exists(c(slim_output_sample_ts6, slim_output_sample_ts7)))){
      
      if(length(slim_output_sample_ts6_sorted) == 0){
        slim_output_sample_ts6_sorted <- paste0(slim_output_folder,"slim_output_sample_ts6_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts6_vcf <- paste("grep '^#'", slim_output_sample_ts6, ">", slim_output_sample_ts6_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts6, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts6_sorted)
        
        system(sort_sample_ts6_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts6_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts6_sorted, ".gz"))) 
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts6_sorted, ".gz")))
        }
      }
      
      if(length(slim_output_sample_ts7_sorted) == 0){
        slim_output_sample_ts7_sorted <- paste0(slim_output_folder,"slim_output_sample_ts7_", sim, "_sorted" , ".vcf")
        
        sort_sample_ts7_vcf <- paste("grep '^#'", slim_output_sample_ts7, ">", slim_output_sample_ts7_sorted,
                                     "&& grep -v '^#'", slim_output_sample_ts7, "| sort -k1,1 -k2,2n", ">>", slim_output_sample_ts7_sorted)
        system(sort_sample_ts7_vcf)
        
        # bgzip sorted vcf files
        system(paste(path_to_bgzip, "-f", slim_output_sample_ts7_sorted))
        
        # tabix bgziped files
        if (genomeS > 2^29){
          # csi instead of tbi for large chromosome 
          system(paste(path_to_tabix, "-C", paste0(slim_output_sample_ts7_sorted, ".gz")))
        } else {
          system(paste(path_to_tabix, paste0(slim_output_sample_ts7_sorted, ".gz")))
        }
      }
      
      if(length(slim_output_sample_merged_Riv) == 0){
        
        # merge and get the data - use the same as Riverside population
        slim_output_sample_merged_Riv <- paste0(slim_output_folder,"slim_output_sample_merged_Riv_", sim, ".txt")
        
        bcftools_query_ <- paste(path_to_bcftools, "merge --force-samples",
                                 paste0(slim_output_sample_ts6_sorted, ".gz"),
                                 paste0(slim_output_sample_ts7_sorted, ".gz"),
                                 "|", path_to_bcftools, "query -f 'chr%CHROM\t%POS\t%REF,%ALT\t%MID\t%DOM\t%GO\t%MT\tY[\t%GT]\n'",
                                 ">", slim_output_sample_merged_Riv) 
        
        system(bcftools_query_)
      }
      
      if(file.exists(slim_output_sample_merged_Riv)){
        
        # assembly the header
        header_1           <- c("chrom", "position","alleles", "MID", "DOM", "GO", "MT", "selection")
        sample_names_1     <- paste0("indiv", seq(from=1, to=SSs[6], by=1), "@pop1", "")
        sample_names_2     <- paste0("indiv", seq(from=1, to=SSs[7], by=1), "@pop2", "")
        
        full_header_       <- c(header_1, sample_names_1, sample_names_2)
        
        # imported the data
        slim_raw_data_ <- read.table(file = slim_output_sample_merged_Riv, header = F, col.names = full_header_, check.names = F, na.strings = "./.")
        
        rm(sample_names_1)
        rm(full_header_)
        
        # if it is a RADseq data
        if (data_type == 2){
          
          slim_raw_data_ <- slim_raw_data_[which(slim_raw_data_$position %in% radseq_interval), ]
          
          slim_raw_data_$chrom <- sapply(slim_raw_data_$position, radseqtagging, tagSampled=radseq_sampled, readLength=radseq_readL)
          
          if (one_snp_radseq){
            slim_raw_data_ <- slim_raw_data_[!duplicated(slim_raw_data_[ ,1]), ]
          }
        }
        
        if (nrow(slim_raw_data_) != 0){
          
          # split the data
          slim_snp_geno_ <- slim_raw_data_[, 9:ncol(slim_raw_data_)]
          slim_snp_geno_ <- slim_snp_geno_[,-c(12:13)]
          slim_snp_info_ <- slim_raw_data_[, 1:length(header_1)]
          
          # remove raw snp data after use it
          rm(slim_raw_data_)
          
          # change the genotype annotations
          slim_snp_geno_ <- as.matrix(slim_snp_geno_)
          slim_snp_geno_[is.na(slim_snp_geno_)]   <- "11"
          slim_snp_geno_[slim_snp_geno_ == "0|0"] <- "11"
          slim_snp_geno_[slim_snp_geno_ =="1|1"]  <- "22"
          
          if (haplotype){
            ref_or_alt <- rbinom(n=1, size = 1, prob = 0.5)
            if (ref_or_alt == 0){
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "11"
            } else {
              slim_snp_geno_[slim_snp_geno_ == "0|1" | slim_snp_geno_ == "1|0"] <- "22"
            }
          } else {
            slim_snp_geno_[slim_snp_geno_ == "0|1"] <- "12"
            slim_snp_geno_[slim_snp_geno_ == "1|0"] <- "21"
          }
          
          # adding missing data randomly
          slim_snp_geno_[sample(1:length(slim_snp_geno_), size=round(length(slim_snp_geno_)*missing_data), replace = FALSE)] <- NA
          slim_snp_geno_[is.na(slim_snp_geno_)] <- "00"
          
          slim_snp_geno_ <- as.data.frame(slim_snp_geno_)
          
          # mark monomophormic mutations (all "11" + "00" or "22")
          count_ref_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "11" | x == "00")})
          count_alt_geno_  <- apply(slim_snp_geno_, 1, function(x){all(x == "22" | x == "00")})
          keep_snps_       <- !(count_ref_geno_ | count_alt_geno_) # MARK MONOMORPHIC MUTATIONS
          
          rm(count_ref_geno_)
          rm(count_alt_geno_)
          
          # re-assemble the data
          slim_data_ <- cbind(slim_snp_info_, slim_snp_geno_)
          
          # remove raw snp data information after use it
          rm(slim_snp_geno_)
          rm(slim_snp_info_)
          
          # remove monomorphic mutations
          slim_data_ <- slim_data_[keep_snps_, ]
          
          # remove vector of kept snps after use it
          rm(keep_snps_)
          
          # remove duplicated mutations
          slim_data_ <- slim_data_[!duplicated(slim_data_[ ,1:2]), ]
          
          if (nrow(slim_data_) != 0){
            
            # make WFABC input file
            if (wfabc_input_file){
              
              slim2wfabc_ <- slim_data_[, -c(1:8)]
              slim2wfabc_ <- as.data.frame(t(slim2wfabc_))
              
              wfabc_data_  <- do.call(rbind, sapply(slim2wfabc_, countgen4wfabc, t_points=2, simplify = F))
              
              if (!file_test("-d", wfabc_input_folder)){
                dir.create(file.path(wfabc_input_folder))
              }
              
              write(paste(dim(slim2wfabc_)[2], 2),file=paste0(wfabc_input_folder, "wfabc_input_sample_Pla_", sim,".txt")) 
              write(paste(0, (tau[6]-tau[5]), sep = ","),file=paste0(wfabc_input_folder, "wfabc_input_sample_Pla_", sim,".txt"), append = TRUE)
              write.table(wfabc_data_, file=paste0(wfabc_input_folder, "wfabc_input_sample_Pla_", sim,".txt"),append=TRUE, quote=FALSE, sep=",", row.names = FALSE, col.names = FALSE)
              
              rm(slim2wfabc_)
              rm(wfabc_data_)
              
            }
            
            # re-code the chromosome name
            if(chrN > 1){
              if (chrTAG){
                slim_data_$chrom <- sapply(slim_data_$position, chromtagging, chrsLimits=chrs_lowers)
              }
            } 
            
            # prepare egglib input data
            slim_to_egglib_data_ <- data.frame(chrom=slim_data_$chrom, 
                                               position=slim_data_$position, 
                                               status=slim_data_$MT,
                                               selection=slim_data_$selection, 
                                               alleles=slim_data_$alleles)
            
            # assembly final egglib input
            slim_to_egglib_data_ <- cbind(slim_to_egglib_data_, slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # re-code the status column
            slim_to_egglib_data_$status <- ifelse(slim_to_egglib_data_$status == 1, "S", "NS")
            
            if (!file_test("-d", egglib_input_folder)){
              dir.create(file.path(egglib_input_folder))
            }
            
            # export egglib input file to the egglib input folder
            egglib_converted_file_Pla <- paste0("egglib_input_sample_Pla", "_", sim, ".txt")
            write.table(slim_to_egglib_data_, file = paste0(egglib_input_folder,egglib_converted_file_Pla), quote=FALSE, sep="\t", row.names = FALSE)
            
            # save only the information of the snps
            if (model_type == 3){
              selcoeff_snps_ <- all_merged_genome_Riv[all_merged_genome_Riv$MID %in% slim_data_$MID, paste0("S",tc)]
            } else {
              selcoeff_snps_ <- all_merged_genome_Riv[all_merged_genome_Riv$MID %in% slim_data_$MID, "S"]
            }
            
            slim_to_egglib_snps_ <- cbind(ID  = paste0(slim_data_$chrom, ":", slim_data_$position), 
                                          MID = slim_data_$MID,
                                          MT  = slim_data_$MT, 
                                          S   = selcoeff_snps_,
                                          DOM = slim_data_$DOM, 
                                          GO  = slim_data_$GO,
                                          slim_data_[, (length(header_1)+1):ncol(slim_data_)])
            
            # expor the complete data of mutations present in egglib input file
            egglib_selcoeff_file_Pla <- paste0("egglib_input_sample_selcoeff_Pla", "_", sim, ".txt")
            
            if (egglib_input_selcoeff){
              write.table(slim_to_egglib_snps_, file = paste0(egglib_input_folder,egglib_selcoeff_file_Pla), quote=FALSE, sep="\t", row.names = FALSE)
            }
            
            # remove snp datasets after use it
            rm(slim_data_)
            rm(slim_to_egglib_data_)
            rm(selcoeff_snps_)
            
            ## RUNNUNG EGGLIB - CALCULATING SUMMARY STATISTICS
            ##-----------------------------------------------------------------------------------
            
            if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Pla))){
              
              # check if the folder exists
              if (!file_test("-d", egglib_output_folder)){
                dir.create(file.path(egglib_output_folder))
              }
              
              # generate text with egglib command  
              egglib_run_ <- paste(path_to_python,
                                   paste0(getwd(), "/", path_to_egglib_summstat),
                                   paste0("input-file=", egglib_input_folder, egglib_converted_file_Pla),
                                   paste0("output-file=", egglib_output_folder, "egglib_output_sample_Pla", "_", sim, ".txt"),
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
              system(egglib_run_)
              
              # import egglib output
              egglib_output_summstat_Pla <- paste0(egglib_output_folder,"egglib_output_sample_Pla", "_", sim, ".txt") 
              
              if(file.exists(egglib_output_summstat_Pla)){
                
                egglib_summary_stats_ <- read.csv(file = egglib_output_summstat_Pla, header = T, sep = "\t", check.names = F)
                
                ## ACTUAL Pr SELECTED MUTATIONS IN THE SAMPLE
                ##-------------------------------------------
                if(any(slim_to_egglib_snps_$S != 0, na.rm = TRUE)){
                  actual_sample_prbe_Pla    <- sum(slim_to_egglib_snps_$S != 0, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                  actual_sample_SelMean_Pla <- mean(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  actual_sample_SelSd_Pla   <- sd(slim_to_egglib_snps_[slim_to_egglib_snps_[, "S"] != 0, "S"], na.rm = TRUE)
                  
                } else {
                  actual_sample_prbe_Pla    <- as.numeric(0)
                  actual_sample_SelMean_Pla <- as.numeric(0)
                  actual_sample_SelSd_Pla   <- as.numeric(0)
                }
                
                if (!is.na(meanNe2_Pla)){
                  ## Pr MUTATIONS UNDER STRONG POSITIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Pla*slim_to_egglib_snps_$S > 1, na.rm = TRUE)){
                    positive_sample_prbe_Pla    <- sum(meanNe2_Pla*slim_to_egglib_snps_$S > 1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    positive_sample_SelMean_Pla <- mean(slim_to_egglib_snps_[meanNe2_Pla*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    positive_sample_SelSd_Pla   <- sd(slim_to_egglib_snps_[meanNe2_Pla*slim_to_egglib_snps_[, "S"] > 1, "S"], na.rm = TRUE)
                    
                  } else {
                    positive_sample_prbe_Pla    <- as.numeric(0)
                    positive_sample_SelMean_Pla <- as.numeric(0)
                    positive_sample_SelSd_Pla   <- as.numeric(0)
                  }
                  
                  ## Pr MUTATIONS UNDER STRONG NEGATIVE SELECTION --[Ns > 1]-- IN THE SAMPLE
                  ##------------------------------------------------------------------------
                  if(any(meanNe2_Pla*slim_to_egglib_snps_$S < -1, na.rm = TRUE)){
                    negative_sample_prbe_Pla    <- sum(meanNe2_Pla*slim_to_egglib_snps_$S < -1, na.rm = TRUE)/length(slim_to_egglib_snps_$S[!is.na(slim_to_egglib_snps_$S)])
                    negative_sample_SelMean_Pla <- mean(slim_to_egglib_snps_[meanNe2_Pla*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    negative_sample_SelSd_Pla   <- sd(slim_to_egglib_snps_[meanNe2_Pla*slim_to_egglib_snps_[, "S"] < -1, "S"], na.rm = TRUE)
                    
                  } else {
                    negative_sample_prbe_Pla    <- as.numeric(0)
                    negative_sample_SelMean_Pla <- as.numeric(0)
                    negative_sample_SelSd_Pla   <- as.numeric(0)
                  }
                  
                  if(is.na(positive_sample_prbe_Pla)  & is.na(negative_sample_prbe_Pla)){
                    
                    strong_sample_prbe_Pla    <- as.numeric(NA)
                    strong_sample_SelMean_Pla <- as.numeric(NA)
                    strong_sample_SelSd_Pla   <- as.numeric(NA)
                    
                  } else {
                    strong_sample_prbe_Pla    <- sum(positive_sample_prbe_Pla, negative_sample_prbe_Pla, na.rm = TRUE)
                    strong_sample_SelMean_Pla <- sum(positive_sample_SelMean_Pla, negative_sample_SelMean_Pla, na.rm = TRUE)
                    strong_sample_SelSd_Pla   <- sum(positive_sample_SelSd_Pla, negative_sample_SelSd_Pla, na.rm = TRUE)
                  }
                  
                } else {
                  strong_sample_prbe_Pla    <- as.numeric(NA)
                  strong_sample_SelMean_Pla <- as.numeric(NA)
                  strong_sample_SelSd_Pla   <- as.numeric(NA)
                }
                
                ## GLOBAL SUMMARY STATISTICS
                ##---------------------------
                
                # remove redundant summary statistics
                egglib_summary_stats_ <- egglib_summary_stats_[, unique(names(egglib_summary_stats_))]
                
                # rename the summary statistics
                colnames(egglib_summary_stats_) <- gsub(":", "_", names(egglib_summary_stats_))
                
                # egglib calculated GLOBAL statistics
                global_stats_ <- egglib_summary_stats_[1 , grepl("^GSS" , unique(names(egglib_summary_stats_)))]
                
                global_SFS_   <- egglib_summary_stats_[1 , grepl("^SFS" , unique(names(egglib_summary_stats_)))]
                
                # calculate additional GLOBAL summary statistics
                mean_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){mean(x, na.rm=T)})
                
                var_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){var(x, na.rm=T)})
                
                kurt_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){kurtosis(x, na.rm=T)})
                
                skew_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){skewness(x, na.rm=T)})
                
                q05_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.05, na.rm=T)})
                
                q95_locus_stats_ <- apply(egglib_summary_stats_[,-c(1, which(grepl("^GSS" , unique(names(egglib_summary_stats_))))[1]:length(egglib_summary_stats_))], 2, function(x){quantile(x, probs = 0.95, na.rm=T)})
                
                # assemble additional GLOBAL summary statistics
                add_global_stats_ <-cbind(as.data.frame(t(mean_locus_stats_)), as.data.frame(t(var_locus_stats_)), as.data.frame(t(kurt_locus_stats_)), 
                                          as.data.frame(t(skew_locus_stats_)), as.data.frame(t(q05_locus_stats_)), as.data.frame(t(q95_locus_stats_)))
                
                rm(mean_locus_stats_)
                rm(var_locus_stats_)
                rm(kurt_locus_stats_)
                rm(skew_locus_stats_)
                rm(q05_locus_stats_)
                rm(q95_locus_stats_)
                
                # ASSEMBLY default GLOBAL summary statistics
                global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
                
                rm(global_stats_)
                rm(global_SFS_)
                rm(add_global_stats_)
                
                ## LOCUS-SPECIFIC FST - ACCURACY, PRECISION, SENSITIVITY, FPR AND FDR
                ##-------------------------------------------------------------------
                slim_to_egglib_snps_$S[is.na(slim_to_egglib_snps_$S)] <- 0
                
                locusFST_test_table_ <- data.frame(Ns=meanNe2_Pla*slim_to_egglib_snps_$S, LSS_WCst=egglib_summary_stats_[, "LSS_WCst"], Ns_test=ifelse(meanNe2_Pla*slim_to_egglib_snps_$S > 1, 1, 0))
                
                locusFST_test_table_ <- locusFST_test_table_[order(-locusFST_test_table_$LSS_WCst), ]
                
                if (any(locusFST_test_table_$Ns > 1)){
                  if(!all(locusFST_test_table_$Ns > 1)){
                    
                    pred_locusFSTNS_ <- prediction(predictions = locusFST_test_table_$LSS_WCst, labels = locusFST_test_table_$Ns_test)
                    perf_locusFSTNS_ <- performance(pred_locusFSTNS_, "ppv", "fpr")
                    
                    perf_table_ <- data.frame(ppvNS=perf_locusFSTNS_@y.values[[1]], fdrNS=1-perf_locusFSTNS_@y.values[[1]])
                    
                    perf_locusFST_table_ <- data.frame(locusFST_test_table_[1:dim(perf_table_)[1], ], perf_table_)
                    
                    whichfdrNS005_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.005),"LSS_WCst"]
                    whichfdrNS01_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.01),"LSS_WCst"]
                    whichfdrNS02_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.02),"LSS_WCst"]
                    whichfdrNS05_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.05),"LSS_WCst"]
                    whichfdrNS10_ <- perf_locusFST_table_[which(perf_locusFST_table_$fdrNS <= 0.10),"LSS_WCst"]
                    
                    if (length(whichfdrNS005_) == 0){
                      FSTfdrNS005_Pla = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS005_Pla = min(whichfdrNS005_)
                    }
                    
                    if (length(whichfdrNS01_) == 0){
                      FSTfdrNS01_Pla = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS01_Pla = min(whichfdrNS01_)
                    }
                    
                    if (length(whichfdrNS02_) == 0){
                      FSTfdrNS02_Pla = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS02_Pla = min(whichfdrNS02_)
                    }
                    
                    if (length(whichfdrNS05_) == 0){
                      FSTfdrNS05_Pla = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS05_Pla = min(whichfdrNS05_)
                    }
                    
                    if (length(whichfdrNS10_) == 0){
                      FSTfdrNS10_Pla = max(locusFST_test_table_$LSS_WCst)
                    } else {
                      FSTfdrNS10_Pla = min(whichfdrNS10_)
                    }
                    
                    FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla,
                                                      FSTfdrNS01_Pla,
                                                      FSTfdrNS02_Pla,
                                                      FSTfdrNS05_Pla,
                                                      FSTfdrNS10_Pla))
                    
                    # remove files
                    rm(pred_locusFSTNS_)
                    rm(perf_locusFSTNS_)
                    rm(perf_table_)
                    rm(perf_locusFST_table_)
                    rm(whichfdrNS005_)
                    rm(whichfdrNS01_) 
                    rm(whichfdrNS02_) 
                    rm(whichfdrNS05_) 
                    rm(whichfdrNS10_)
                    
                  } else {
                    
                    # all strongly selected mutations
                    FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS01_Pla  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS02_Pla  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS05_Pla  = min(locusFST_test_table_$LSS_WCst),
                                                      FSTfdrNS10_Pla  = min(locusFST_test_table_$LSS_WCst)))
                    
                  }
                } else {
                  
                  # all neutral mutations
                  FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS01_Pla = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS02_Pla = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS05_Pla = max(locusFST_test_table_$LSS_WCst),
                                                    FSTfdrNS10_Pla = max(locusFST_test_table_$LSS_WCst)))
                  
                }
                
                rm(locusFST_test_table_)
                
                ## LOCUS-SPECIFIC SUMMARY STATISTICS
                ##----------------------------------
                
                # sampling ONE RANDOM mutation for the locus-specific reference table
                snps_in_reftable_ <- sample(which(slim_to_egglib_snps_$MT == 1 | slim_to_egglib_snps_$MT == 2 | slim_to_egglib_snps_$MT == 3), size=1)
                sampled_snp_ <- slim_to_egglib_snps_[snps_in_reftable_, ]
                
                # remove complete snp table after use it
                rm(slim_to_egglib_snps_)
                rm(snps_in_reftable_)
                
                # calculate sample minor allele frequency
                sampled_snp_genotypes_ <- sampled_snp_[, (grep("GO", names(sampled_snp_)) + 1):(SSs[6] + (SSs[7]-2) + 6)] #HERE
                
                # sample alternative allele frequency - SAAF1
                S1_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 1), names(sampled_snp_genotypes_))]
                S1n11_ <- apply(S1_genotypes_==11, 1, sum, na.rm=T)
                S1n12_ <- apply(S1_genotypes_==12 | S1_genotypes_==21, 1, sum, na.rm=T)
                S1n22_ <- apply(S1_genotypes_==22, 1, sum, na.rm=T)
                SAAF1_ <- (2*(S1n22_) + S1n12_)/((2*(S1n11_) + S1n12_)+(2*(S1n22_) + S1n12_))
                
                # sample alternative allele frequency - SAAF2
                S2_genotypes_ <- sampled_snp_genotypes_[, grepl(paste0("@pop", 2), names(sampled_snp_genotypes_))]
                S2n11_ <- apply(S2_genotypes_==11, 1, sum, na.rm=T)
                S2n12_ <- apply(S2_genotypes_==12 | S2_genotypes_==21, 1, sum, na.rm=T)
                S2n22_ <- apply(S2_genotypes_==22, 1, sum, na.rm=T)
                SAAF2_ <- (2*(S2n22_) + S2n12_)/((2*(S2n11_) + S2n12_)+(2*(S2n22_) + S2n12_))
                
                # assemble LOCUS-SPECIFIC summary statistics
                locus_lss_info_  <- sampled_snp_[, which(grepl("^ID" , unique(names(sampled_snp_)))):which(grepl("^GO" , unique(names(sampled_snp_))))]
                locus_lss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID, grepl("^LSS" , unique(names(egglib_summary_stats_)))]
                locus_wss_stats_ <- egglib_summary_stats_[egglib_summary_stats_$ID %in% sampled_snp_$ID , grepl("^WSS" , unique(names(egglib_summary_stats_)))]
                
                # remove summary statistics data after use it
                rm(egglib_summary_stats_)
                rm(sampled_snp_)
                rm(sampled_snp_genotypes_)
                
                # ASSEMBLY default LOCUS-SPECIFIC summary statistics
                locus_summary_stats_Pla <- cbind(locus_lss_info_, SAAF1_, SAAF2_, locus_lss_stats_, locus_wss_stats_)
                colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
                
                rm(S1_genotypes_)
                rm(S1n11_)
                rm(S1n12_)
                rm(S1n22_)
                rm(SAAF1_)
                rm(S2_genotypes_)
                rm(S2n11_)
                rm(S2n12_)
                rm(S2n22_)
                rm(SAAF2_)
                rm(locus_lss_info_)
                rm(locus_lss_stats_)
                rm(locus_wss_stats_)
                
                if (remove_files){
                  file.remove(paste0(egglib_output_folder,"egglib_output_sample_Pla", "_", sim, ".txt"))
                }
                
              } else {
                
                actual_sample_prbe_Pla    <- as.numeric(NA)
                actual_sample_SelMean_Pla <- as.numeric(NA)
                actual_sample_SelSd_Pla   <- as.numeric(NA)
                
                strong_sample_prbe_Pla    <- as.numeric(NA)
                strong_sample_SelMean_Pla <- as.numeric(NA)
                strong_sample_SelSd_Pla   <- as.numeric(NA)
                
                global_stats_ <- as.data.frame(t(rep(NA, 20)))
                global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
                add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
                global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
                colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
                
                FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                                  FSTfdrNS01_Pla  = NA,
                                                  FSTfdrNS02_Pla  = NA,
                                                  FSTfdrNS05_Pla  = NA,
                                                  FSTfdrNS10_Pla  = NA))
                
                locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
                colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
                
                if (debug_sim){
                  
                  # check if the folder exists
                  if (!file_test("-d", debug_output_folder)){
                    dir.create(file.path(debug_output_folder))
                  }
                  
                  file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                  if(file.exists(slim_output_sample_ts6)){
                    file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_ts7)){
                    file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                  }
                  
                  if(file.exists(slim_output_sample_merged_Riv)){
                    file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_converted_file_Pla))){
                    file.copy(from = paste0(egglib_input_folder, egglib_converted_file_Pla), to = debug_output_folder)
                  }
                  
                  if(file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Pla))){
                    file.copy(from = paste0(egglib_input_folder, egglib_selcoeff_file_Pla), to = debug_output_folder)
                  }
                  
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
              
              actual_sample_prbe_Pla    <- as.numeric(NA)
              actual_sample_SelMean_Pla <- as.numeric(NA)
              actual_sample_SelSd_Pla   <- as.numeric(NA)
              
              strong_sample_prbe_Pla    <- as.numeric(NA)
              strong_sample_SelMean_Pla <- as.numeric(NA)
              strong_sample_SelSd_Pla   <- as.numeric(NA)
              
              global_stats_ <- as.data.frame(t(rep(NA, 20)))
              global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
              add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
              global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
              colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
              
              FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                                FSTfdrNS01_Pla  = NA,
                                                FSTfdrNS02_Pla  = NA,
                                                FSTfdrNS05_Pla  = NA,
                                                FSTfdrNS10_Pla  = NA))
              
              locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
              colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
              
              if (debug_sim){
                
                # check if the folder exists
                if (!file_test("-d", debug_output_folder)){
                  dir.create(file.path(debug_output_folder))
                }
                
                file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
                if(file.exists(slim_output_sample_ts6)){
                  file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_ts7)){
                  file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
                }
                
                if(file.exists(slim_output_sample_merged_Riv)){
                  file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
                }
                
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
            
            if (remove_files){
              if (file.exists(paste0(egglib_input_folder, egglib_converted_file_Pla))){
                file.remove(paste0(egglib_input_folder, egglib_converted_file_Pla))
              }
              if (file.exists(paste0(egglib_input_folder, egglib_selcoeff_file_Pla))){
                file.remove(paste0(egglib_input_folder, egglib_selcoeff_file_Pla))
              }
            }
            
          } else {
            
            actual_sample_prbe_Pla    <- as.numeric(NA)
            actual_sample_SelMean_Pla <- as.numeric(NA)
            actual_sample_SelSd_Pla   <- as.numeric(NA)
            
            strong_sample_prbe_Pla    <- as.numeric(NA)
            strong_sample_SelMean_Pla <- as.numeric(NA)
            strong_sample_SelSd_Pla   <- as.numeric(NA)
            
            global_stats_ <- as.data.frame(t(rep(NA, 20)))
            global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
            add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
            global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
            colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
            
            FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                              FSTfdrNS01_Pla  = NA,
                                              FSTfdrNS02_Pla  = NA,
                                              FSTfdrNS05_Pla  = NA,
                                              FSTfdrNS10_Pla  = NA))
            
            locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
            colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
            
            if (debug_sim){
              
              # check if the folder exists
              if (!file_test("-d", debug_output_folder)){
                dir.create(file.path(debug_output_folder))
              }
              
              file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
              
              if(file.exists(slim_output_sample_ts6)){
                file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_ts7)){
                file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
              }
              
              if(file.exists(slim_output_sample_merged_Riv)){
                file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
              }
              
              debug_message <- "Not enough polymorphism"
              
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
          
          actual_sample_prbe_Pla    <- as.numeric(NA)
          actual_sample_SelMean_Pla <- as.numeric(NA)
          actual_sample_SelSd_Pla   <- as.numeric(NA)
          
          strong_sample_prbe_Pla    <- as.numeric(NA)
          strong_sample_SelMean_Pla <- as.numeric(NA)
          strong_sample_SelSd_Pla   <- as.numeric(NA)
          
          global_stats_ <- as.data.frame(t(rep(NA, 20)))
          global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
          add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
          global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
          colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
          
          FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                            FSTfdrNS01_Pla  = NA,
                                            FSTfdrNS02_Pla  = NA,
                                            FSTfdrNS05_Pla  = NA,
                                            FSTfdrNS10_Pla  = NA))
          
          locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
          colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
          
          if (debug_sim){
            
            # check if the folder exists
            if (!file_test("-d", debug_output_folder)){
              dir.create(file.path(debug_output_folder))
            }
            
            file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
            
            if(file.exists(slim_output_sample_ts6)){
              file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_ts7)){
              file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
            }
            
            if(file.exists(slim_output_sample_merged_Riv)){
              file.copy(from = slim_output_sample_merged_Riv, to = debug_output_folder)
            }
            
            debug_message <- "RADseq sampling removed all genotypic data"
            
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
        
        actual_sample_prbe_Pla    <- as.numeric(NA)
        actual_sample_SelMean_Pla <- as.numeric(NA)
        actual_sample_SelSd_Pla   <- as.numeric(NA)
        
        strong_sample_prbe_Pla    <- as.numeric(NA)
        strong_sample_SelMean_Pla <- as.numeric(NA)
        strong_sample_SelSd_Pla   <- as.numeric(NA)
        
        global_stats_ <- as.data.frame(t(rep(NA, 20)))
        global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
        add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
        global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
        colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
        
        FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                          FSTfdrNS01_Pla  = NA,
                                          FSTfdrNS02_Pla  = NA,
                                          FSTfdrNS05_Pla  = NA,
                                          FSTfdrNS10_Pla  = NA))
        
        locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
        colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
        
        if (debug_sim){
          
          # check if the folder exists
          if (!file_test("-d", debug_output_folder)){
            dir.create(file.path(debug_output_folder))
          }
          
          file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
          
          if(file.exists(slim_output_sample_ts6)){
            file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
          }
          
          if(file.exists(slim_output_sample_ts7)){
            file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
          }
          
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
      
      if (remove_files){
        file.remove(paste0(slim_output_sample_ts6))
        file.remove(paste0(slim_output_sample_ts7))
        
        if (file.exists(paste0(slim_output_sample_ts6_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts6_sorted, ".gz"))
        }
        if (file.exists(paste0(slim_output_sample_ts7_sorted, ".gz"))){
          file.remove(paste0(slim_output_sample_ts7_sorted, ".gz"))
        }
        if (genomeS > 2^29){
          if (file.exists(paste0(slim_output_sample_ts6_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts6_sorted, ".gz.csi"))
          }
          if (file.exists(paste0(slim_output_sample_ts7_sorted, ".gz.csi"))){
            file.remove(paste0(slim_output_sample_ts7_sorted, ".gz.csi"))
          }
        } else {
          if (file.exists(paste0(slim_output_sample_ts6_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts6_sorted, ".gz.tbi"))
          }
          if (file.exists(paste0(slim_output_sample_ts7_sorted, ".gz.tbi"))){
            file.remove(paste0(slim_output_sample_ts7_sorted, ".gz.tbi"))
          }
        }
        
        file.remove(paste0(slim_output_sample_merged_Riv))
      }
      
    } else {
      
      actual_sample_prbe_Pla    <- as.numeric(NA)
      actual_sample_SelMean_Pla <- as.numeric(NA)
      actual_sample_SelSd_Pla   <- as.numeric(NA)
      
      strong_sample_prbe_Pla    <- as.numeric(NA)
      strong_sample_SelMean_Pla <- as.numeric(NA)
      strong_sample_SelSd_Pla   <- as.numeric(NA)
      
      global_stats_ <- as.data.frame(t(rep(NA, 20)))
      global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run))) 
      add_global_stats_ <- as.data.frame(t(rep(NA, 198))) 
      global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
      colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
      
      FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                        FSTfdrNS01_Pla  = NA,
                                        FSTfdrNS02_Pla  = NA,
                                        FSTfdrNS05_Pla  = NA,
                                        FSTfdrNS10_Pla  = NA))
      
      locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
      colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
      
      if (debug_sim){
        
        # check if the folder exists
        if (!file_test("-d", debug_output_folder)){
          dir.create(file.path(debug_output_folder))
        }
        
        file.copy(from = paste0(slim_output_folder, "slim_coalesced_",model_title,"_", sim, ".tree"), to = debug_output_folder)
        if(file.exists(slim_output_sample_ts6)){
          file.copy(from = slim_output_sample_ts6, to = debug_output_folder)
        }
        
        if(file.exists(slim_output_sample_ts7)){
          file.copy(from = slim_output_sample_ts7, to = debug_output_folder)
        }
        
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
    }

  } else {
    
    ## NE EQUILIBRIUM PERIOD - NE(1)
    meanNe1 <- as.numeric(NA)
    
    ## NE SAMPLING PERIOD - NE(2)
    timesNe2 <- as.data.frame(t(rep(NA, (tau[6]+1))))
    names(timesNe2) <- paste0("timesNe2_",seq(from=0,to=(tau[6])))
    meanNe2_Ava <- as.numeric(NA)
    meanNe2_Hum <- as.numeric(NA)
    meanNe2_Dav <- as.numeric(NA)
    meanNe2_Sta <- as.numeric(NA)
    meanNe2_Ste <- as.numeric(NA)
    meanNe2_Riv = meanNe2_Pla <- as.numeric(NA)
    
    ## GENETIC LOAD OF SAMPLING PERIOD
    averageGenLoad_Ava <- as.numeric(NA)
    averageGenLoad_Hum <- as.numeric(NA)
    averageGenLoad_Dav <- as.numeric(NA)
    averageGenLoad_Sta <- as.numeric(NA)
    averageGenLoad_Ste <- as.numeric(NA)
    averageGenLoad_Riv = averageGenLoad_Pla <- as.numeric(NA)
    
    ## PROPORTION OF SELECTED MUTATIONS IN THE POPULATION
    actual_pop_prbe_Ava    <- as.numeric(NA)
    actual_pop_SelMean_Ava <- as.numeric(NA)
    actual_pop_SelSd_Ava   <- as.numeric(NA)
    
    strong_pop_prbe_Ava    <- as.numeric(NA)
    strong_pop_SelMean_Ava <- as.numeric(NA)
    strong_pop_SelSd_Ava   <- as.numeric(NA)
    
    actual_pop_prbe_Hum    <- as.numeric(NA)
    actual_pop_SelMean_Hum <- as.numeric(NA)
    actual_pop_SelSd_Hum   <- as.numeric(NA)
    
    strong_pop_prbe_Hum    <- as.numeric(NA)
    strong_pop_SelMean_Hum <- as.numeric(NA)
    strong_pop_SelSd_Hum   <- as.numeric(NA)
    
    actual_pop_prbe_Dav    <- as.numeric(NA)
    actual_pop_SelMean_Dav <- as.numeric(NA)
    actual_pop_SelSd_Dav   <- as.numeric(NA)
    
    strong_pop_prbe_Dav    <- as.numeric(NA)
    strong_pop_SelMean_Dav <- as.numeric(NA)
    strong_pop_SelSd_Dav   <- as.numeric(NA)
    
    actual_pop_prbe_Sta    <- as.numeric(NA)
    actual_pop_SelMean_Sta <- as.numeric(NA)
    actual_pop_SelSd_Sta   <- as.numeric(NA)
    
    strong_pop_prbe_Sta    <- as.numeric(NA)
    strong_pop_SelMean_Sta <- as.numeric(NA)
    strong_pop_SelSd_Sta   <- as.numeric(NA)
    
    actual_pop_prbe_Ste    <- as.numeric(NA)
    actual_pop_SelMean_Ste <- as.numeric(NA)
    actual_pop_SelSd_Ste   <- as.numeric(NA)
    
    strong_pop_prbe_Ste    <- as.numeric(NA)
    strong_pop_SelMean_Ste <- as.numeric(NA)
    strong_pop_SelSd_Ste   <- as.numeric(NA)
    
    actual_pop_prbe_Riv    = actual_pop_prbe_Pla    <- as.numeric(NA)
    actual_pop_SelMean_Riv = actual_pop_SelMean_Pla <- as.numeric(NA)
    actual_pop_SelSd_Riv   = actual_pop_SelSd_Pla   <- as.numeric(NA)
    
    strong_pop_prbe_Riv    = strong_pop_prbe_Pla    <- as.numeric(NA)
    strong_pop_SelMean_Riv = strong_pop_SelMean_Pla <- as.numeric(NA)
    strong_pop_SelSd_Riv   = strong_pop_SelSd_Pla   <- as.numeric(NA)
    
    ## SUMMARY STATISTICS
    
    global_stats_ <- as.data.frame(t(rep(NA, 20)))
    global_SFS_ <- as.data.frame(t(rep(NA, sfs_bins_run)))
    add_global_stats_ <- as.data.frame(t(rep(NA, 198)))
    
    ## Avalon population
    ##-----------------
    
    actual_sample_prbe_Ava    <- as.numeric(NA)
    actual_sample_SelMean_Ava <- as.numeric(NA)
    actual_sample_SelSd_Ava   <- as.numeric(NA)
    
    strong_sample_prbe_Ava    <- as.numeric(NA)
    strong_sample_SelMean_Ava <- as.numeric(NA)
    strong_sample_SelSd_Ava   <- as.numeric(NA)
    
    global_summary_stats_Ava <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Ava) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ava)[2]),"_","Ava")
    
    FSTfdr_Ava <- as.data.frame(cbind(FSTfdrNS005_Ava = NA,
                                      FSTfdrNS01_Ava  = NA,
                                      FSTfdrNS02_Ava  = NA,
                                      FSTfdrNS05_Ava  = NA,
                                      FSTfdrNS10_Ava  = NA))
    
    locus_summary_stats_Ava <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Ava) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ava)[2]),"_","Ava")
    
    ## Humboldt population
    ##-----------------
    
    actual_sample_prbe_Hum    <- as.numeric(NA)
    actual_sample_SelMean_Hum <- as.numeric(NA)
    actual_sample_SelSd_Hum   <- as.numeric(NA)
    
    strong_sample_prbe_Hum    <- as.numeric(NA)
    strong_sample_SelMean_Hum <- as.numeric(NA)
    strong_sample_SelSd_Hum   <- as.numeric(NA)
    
    global_summary_stats_Hum <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Hum) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Hum)[2]),"_","Hum")
    
    FSTfdr_Hum <- as.data.frame(cbind(FSTfdrNS005_Hum = NA,
                                      FSTfdrNS01_Hum  = NA,
                                      FSTfdrNS02_Hum  = NA,
                                      FSTfdrNS05_Hum  = NA,
                                      FSTfdrNS10_Hum  = NA))
    
    locus_summary_stats_Hum <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Hum) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Hum)[2]),"_","Hum")
    
    ## Davis population
    ##-----------------
    
    actual_sample_prbe_Dav    <- as.numeric(NA)
    actual_sample_SelMean_Dav <- as.numeric(NA)
    actual_sample_SelSd_Dav   <- as.numeric(NA)
    
    strong_sample_prbe_Dav    <- as.numeric(NA)
    strong_sample_SelMean_Dav <- as.numeric(NA)
    strong_sample_SelSd_Dav   <- as.numeric(NA)
    
    global_summary_stats_Dav <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Dav) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Dav)[2]),"_","Dav")
    
    FSTfdr_Dav <- as.data.frame(cbind(FSTfdrNS005_Dav = NA,
                                      FSTfdrNS01_Dav  = NA,
                                      FSTfdrNS02_Dav  = NA,
                                      FSTfdrNS05_Dav  = NA,
                                      FSTfdrNS10_Dav  = NA))
    
    locus_summary_stats_Dav <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Dav) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Dav)[2]),"_","Dav")
    
    ## Stanislaus population
    ##-----------------
    
    actual_sample_prbe_Sta    <- as.numeric(NA)
    actual_sample_SelMean_Sta <- as.numeric(NA)
    actual_sample_SelSd_Sta   <- as.numeric(NA)
    
    strong_sample_prbe_Sta    <- as.numeric(NA)
    strong_sample_SelMean_Sta <- as.numeric(NA)
    strong_sample_SelSd_Sta   <- as.numeric(NA)
    
    global_summary_stats_Sta <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Sta) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Sta)[2]),"_","Sta")
    
    FSTfdr_Sta <- as.data.frame(cbind(FSTfdrNS005_Sta = NA,
                                      FSTfdrNS01_Sta  = NA,
                                      FSTfdrNS02_Sta  = NA,
                                      FSTfdrNS05_Sta  = NA,
                                      FSTfdrNS10_Sta  = NA))
    
    locus_summary_stats_Sta <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Sta) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Sta)[2]),"_","Sta")
    
    ## Stebbins population
    ##-----------------
    
    actual_sample_prbe_Ste    <- as.numeric(NA)
    actual_sample_SelMean_Ste <- as.numeric(NA)
    actual_sample_SelSd_Ste   <- as.numeric(NA)
    
    strong_sample_prbe_Ste    <- as.numeric(NA)
    strong_sample_SelMean_Ste <- as.numeric(NA)
    strong_sample_SelSd_Ste   <- as.numeric(NA)
    
    global_summary_stats_Ste <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Ste) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Ste)[2]),"_","Ste")
    
    FSTfdr_Ste <- as.data.frame(cbind(FSTfdrNS005_Ste = NA,
                                      FSTfdrNS01_Ste  = NA,
                                      FSTfdrNS02_Ste  = NA,
                                      FSTfdrNS05_Ste  = NA,
                                      FSTfdrNS10_Ste  = NA))
    
    locus_summary_stats_Ste <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Ste) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Ste)[2]),"_","Ste")
    
    ## Riverside population
    ##-----------------
    
    actual_sample_prbe_Riv    <- as.numeric(NA)
    actual_sample_SelMean_Riv <- as.numeric(NA)
    actual_sample_SelSd_Riv   <- as.numeric(NA)
    
    strong_sample_prbe_Riv    <- as.numeric(NA)
    strong_sample_SelMean_Riv <- as.numeric(NA)
    strong_sample_SelSd_Riv   <- as.numeric(NA)
    
    global_summary_stats_Riv <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Riv) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Riv)[2]),"_","Riv")
    
    FSTfdr_Riv <- as.data.frame(cbind(FSTfdrNS005_Riv = NA,
                                      FSTfdrNS01_Riv  = NA,
                                      FSTfdrNS02_Riv  = NA,
                                      FSTfdrNS05_Riv  = NA,
                                      FSTfdrNS10_Riv  = NA))
    
    locus_summary_stats_Riv <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Riv) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Riv)[2]),"_","Riv")
    
    ## Placerita population
    ##-----------------
    
    actual_sample_prbe_Pla    <- as.numeric(NA)
    actual_sample_SelMean_Pla <- as.numeric(NA)
    actual_sample_SelSd_Pla   <- as.numeric(NA)
    
    strong_sample_prbe_Pla    <- as.numeric(NA)
    strong_sample_SelMean_Pla <- as.numeric(NA)
    strong_sample_SelSd_Pla   <- as.numeric(NA)
    
    global_summary_stats_Pla <- cbind(global_stats_, global_SFS_, add_global_stats_)
    colnames(global_summary_stats_Pla) <- paste0("GSS","_",seq(from=1,to=dim(global_summary_stats_Pla)[2]),"_","Pla")
    
    FSTfdr_Pla <- as.data.frame(cbind(FSTfdrNS005_Pla = NA,
                                      FSTfdrNS01_Pla  = NA,
                                      FSTfdrNS02_Pla  = NA,
                                      FSTfdrNS05_Pla  = NA,
                                      FSTfdrNS10_Pla  = NA))
    
    locus_summary_stats_Pla <- as.data.frame(t(rep(NA, 41)))
    colnames(locus_summary_stats_Pla) <- paste0("LSS","_",seq(from=1,to=dim(locus_summary_stats_Pla)[2]),"_","Pla")
    
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
    
  }
  
  ## RAW REFERENCE TABLE
  ##------------------------------------------------------------------------------------------------
  
  raw_reftable  <- suppressWarnings(cbind(as.factor(model_type), 
                                          sim_seed, sim, 
                                          mu, rr, selfing, Neq, Ncs,
                                          gammaM, gammak, tc,
                                          PrGWSel, prbe, meanNe1, timesNe2,
                                          meanNe2_Ava, averageGenLoad_Ava,
                                          actual_pop_prbe_Ava, actual_pop_SelMean_Ava, actual_pop_SelSd_Ava,  
                                          strong_pop_prbe_Ava, strong_pop_SelMean_Ava, strong_pop_SelSd_Ava,  
                                          actual_sample_prbe_Ava, actual_sample_SelMean_Ava, actual_sample_SelSd_Ava,  
                                          strong_sample_prbe_Ava, strong_sample_SelMean_Ava, strong_sample_SelSd_Ava,
                                          FSTfdr_Ava,locus_summary_stats_Ava, global_summary_stats_Ava,
                                          meanNe2_Hum, averageGenLoad_Hum,
                                          actual_pop_prbe_Hum, actual_pop_SelMean_Hum, actual_pop_SelSd_Hum,  
                                          strong_pop_prbe_Hum, strong_pop_SelMean_Hum, strong_pop_SelSd_Hum,  
                                          actual_sample_prbe_Hum, actual_sample_SelMean_Hum, actual_sample_SelSd_Hum,  
                                          strong_sample_prbe_Hum, strong_sample_SelMean_Hum, strong_sample_SelSd_Hum,
                                          FSTfdr_Hum,locus_summary_stats_Hum, global_summary_stats_Hum,
                                          meanNe2_Dav, averageGenLoad_Dav,
                                          actual_pop_prbe_Dav, actual_pop_SelMean_Dav, actual_pop_SelSd_Dav,  
                                          strong_pop_prbe_Dav, strong_pop_SelMean_Dav, strong_pop_SelSd_Dav,  
                                          actual_sample_prbe_Dav, actual_sample_SelMean_Dav, actual_sample_SelSd_Dav,  
                                          strong_sample_prbe_Dav, strong_sample_SelMean_Dav, strong_sample_SelSd_Dav,
                                          FSTfdr_Dav,locus_summary_stats_Dav, global_summary_stats_Dav,
                                          meanNe2_Sta, averageGenLoad_Sta,
                                          actual_pop_prbe_Sta, actual_pop_SelMean_Sta, actual_pop_SelSd_Sta,  
                                          strong_pop_prbe_Sta, strong_pop_SelMean_Sta, strong_pop_SelSd_Sta,  
                                          actual_sample_prbe_Sta, actual_sample_SelMean_Sta, actual_sample_SelSd_Sta,  
                                          strong_sample_prbe_Sta, strong_sample_SelMean_Sta, strong_sample_SelSd_Sta,
                                          FSTfdr_Sta,locus_summary_stats_Sta, global_summary_stats_Sta, 
                                          meanNe2_Ste, averageGenLoad_Ste,
                                          actual_pop_prbe_Ste, actual_pop_SelMean_Ste, actual_pop_SelSd_Ste,  
                                          strong_pop_prbe_Ste, strong_pop_SelMean_Ste, strong_pop_SelSd_Ste,  
                                          actual_sample_prbe_Ste, actual_sample_SelMean_Ste, actual_sample_SelSd_Ste,  
                                          strong_sample_prbe_Ste, strong_sample_SelMean_Ste, strong_sample_SelSd_Ste,
                                          FSTfdr_Ste,locus_summary_stats_Ste, global_summary_stats_Ste,
                                          meanNe2_Riv, averageGenLoad_Riv,
                                          actual_pop_prbe_Riv, actual_pop_SelMean_Riv, actual_pop_SelSd_Riv,  
                                          strong_pop_prbe_Riv, strong_pop_SelMean_Riv, strong_pop_SelSd_Riv,  
                                          actual_sample_prbe_Riv, actual_sample_SelMean_Riv, actual_sample_SelSd_Riv,  
                                          strong_sample_prbe_Riv, strong_sample_SelMean_Riv, strong_sample_SelSd_Riv,
                                          FSTfdr_Riv,locus_summary_stats_Riv, global_summary_stats_Riv,
                                          meanNe2_Pla, averageGenLoad_Pla,
                                          actual_pop_prbe_Pla, actual_pop_SelMean_Pla, actual_pop_SelSd_Pla,  
                                          strong_pop_prbe_Pla, strong_pop_SelMean_Pla, strong_pop_SelSd_Pla,  
                                          actual_sample_prbe_Pla, actual_sample_SelMean_Pla, actual_sample_SelSd_Pla,  
                                          strong_sample_prbe_Pla, strong_sample_SelMean_Pla, strong_sample_SelSd_Pla,
                                          FSTfdr_Pla,locus_summary_stats_Pla, global_summary_stats_Pla
                                          )
                                    )
  rownames(raw_reftable) <- sim
  
  ## OUTPUT RAW REFERENCE TABLE
  return(raw_reftable)
  
}
