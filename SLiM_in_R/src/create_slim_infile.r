##########################################################################################
# simuprior(nsim, ne = N_size, ..., int12_max = value);                                  #
#                                                                                        #
# This function creates a table with simulation parameters                               #
##########################################################################################                            

simuprior <- function(nsim, 
                      ne = 1000, ne_min = NULL, ne_max = NULL,
                      n = 100,   n_min = NULL,  n_man = NULL,
                      mud1 = 0.5, mud1_min = NULL, mud1_max = NULL,
                      mud2 = 0.5, mud2_min = NULL, mud2_max = NULL,
                      mufv2 = 0.1, mufv2_min = NULL, mufv2_max = NULL,
                      ts1 = 10, ts1_min = NULL, ts1_max = NULL,
                      int12 = 10, int12_min = NULL, int12_max = NULL
){
  
  simnumb <- seq(1,nsim,1);
  if (!is.null(ne_min)){
    ne <- as.integer(runif(nsim, ne_min, ne_max));
  } else {
    ne = rep(ne, nsim);
  }
  if (!is.null(n_min)){
    n <- as.integer(runif(nsim, n_min, n_max));
  } else {
    n = rep(n, nsim);
  }
  if (!is.null(mud1_min)){
    mud1 <- runif(nsim, mud1_min, mud1_max);
  } else {
    mud1 = rep(mud1, nsim);
  }
  if (!is.null(mud2_min)){
    mud2 <- runif(nsim, mud2_min, mud2_max);
  } else {
    mud2 = rep(mud2, nsim);
  }
  if (!is.null(mufv2_min)){
    mufv2 <- runif(nsim, mufv2_min, mufv2_max);
  } else {
    mufv2 = rep(mufv2, nsim);
  }
  if (!is.null(ts1_min)){
    ts1 <- as.integer(runif(nsim, ts1_min, ts1_max));
  } else {
    ts1 = rep(ts1, nsim);
  }
  if (!is.null(int12_min)){
    int12 <- as.integer(runif(nsim, int12_min_min, int12_max));
  } else {
    int12 = rep(int12, nsim);
  }
  
  parm <- cbind(simnumb, ne, n, mud1, mud2, mufv2, ts1, int12);
  return(parm)
}

##########################################################################################
# toymodel(nsim, ne = N_size, filename = 'infile_N_size', folder = 'data/');             #
#                                                                                        #
# This function creates the infile for a toymodel to run SLiM2 simulations.              #
# Toymodel helped the initial stages of R/SLiM implementation                            #                             
#																						                                             #
# ne = a vector containing population size N to be simulated;                            #
# ne_min and ne_max are the defined range values for random sampling taking from uniform #
# filename = prefix for a filename system;                                               #
# folder = path where to save the generated SLiM2 infiles.                               #
##########################################################################################

toymodel <- function(nsim,
                     ne = 1000, ne_min = NULL, ne_max = NULL,
                     filename, 
                     folder)
{
        simnumb <- seq(1,nsim,1);
        if (!is.null(ne_min)){
          ne <- as.integer(runif(nsim, ne_min, ne_max));
        } else {
          ne <- rep(ne, nsim);
        } 
        parm <- cbind(simnumb, ne);
          for (i in 1:nsim){
                sink(file = paste0(folder, filename, '_', simnumb[i], '.slim'), type = 'output');
    
                      cat(paste('initialize() {', '\n',
                                '\t', 'initializeMutationRate(1e-7);', '\n',
                                '\t', 'initializeMutationType("m1", 0.5, "f", 0.0);', '\n',
                                '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
                                '\t', 'initializeGenomicElement(g1, 0, 999999);', '\n',
                                '\t', 'initializeRecombinationRate(1e-8);', '\n',
                                '}', '\n\n',
                  
                                '1 { sim.addSubpop("p1",', ne[i], '); }', '\n\n',
                                ne[i]*10, 'late() { cat(sim.mutations.size() + "\\n");', '\n\n',
                  
                                '\t','sim.simulationFinished();', '\n\n',
                                
                                '}'));
    
                sink();
          }; 
        return(parm);
};

##########################################################################################
# MODEL 1 - HARD SWEEP                                         													 #
# 																						                                           #
# model1(nsim, ne = N_size, ..., filename = 'infile_slim', folder = 'infiles/');         #
#                                                                                        #
# This function creates the infile for model 1 to run SLiM2 simulations.                 #
##########################################################################################

model1 <- function(nsim, 
                   ne = 1000, ne_min = NULL, ne_max = NULL, 
                   n = 100,   n_min = NULL,  n_man = NULL,
                   mur = 2.5e-8, mur_min = NULL, mur_max = NULL,
                   mud1 = 0.5, mud1_min = NULL, mud1_max = NULL, 
                   mufd1 = "f", mufv1 = 0.0, mufv1_min = NULL, mufv1_max = NULL,
                   mud2 = 0.5, mud2_min = NULL, mud2_max = NULL, 
                   mufd2 = "f", mufv2 = 0.0, mufv2_min = NULL, mufv2_max = NULL,
                   ts1 = 10, ts1_min = NULL, ts1_max = NULL,  
                   int12 = 10, int12_min = NULL, int12_max = NULL, 
                   genomesize = 100000,
                   filename, folder)
{
  
        gb1 = 0; gb2 = ((genomesize/2)-1); gb3 = ((genomesize/2)+1); gb4 = (genomesize-1);
        
        simnumb <- seq(1,nsim,1);
        if (!is.null(ne_min)){
          ne <- as.integer(runif(nsim, ne_min, ne_max));
        } else {
          ne = as.integer(rep(ne, nsim));
        }
        if (!is.null(n_min)){
          n <- as.integer(runif(nsim, n_min, n_max));
        } else {
          n = as.integer(rep(n, nsim));
        }
        if (!is.null(mur_min)){
          mur <- runif(nsim, mur_min, mur_max);
        } else {
          mur = rep(mur, nsim);
        }
        if (!is.null(mud1_min)){
          mud1 <- runif(nsim, mud1_min, mud1_max);
        } else {
          mud1 = rep(mud1, nsim);
        }
        if (!is.null(mufv1_min)){
          mufv1 <- runif(nsim, mufv1_min, mufv1_max);
        } else {
          mufv1 = rep(mufv1, nsim);
        }
        if (!is.null(mud2_min)){
          mud2 <- runif(nsim, mud2_min, mud2_max);
        } else {
          mud2 = rep(mud2, nsim);
        }
        if (!is.null(mufv2_min)){
          mufv2 <- runif(nsim, mufv2_min, mufv2_max);
        } else {
          mufv2 = rep(mufv2, nsim);
        }
        if (!is.null(ts1_min)){
          ts1 <- as.integer(runif(nsim, ts1_min, ts1_max));
        } else {
          ts1 = as.integer(rep(ts1, nsim));
        }
        if (!is.null(int12_min)){
          int12 <- as.integer(runif(nsim, int12_min, int12_max));
        } else {
          int12 = as.integer(rep(int12, nsim));
        }
        ts2 = ts1+int12;
        
        parm <- cbind(simnumb, ne, n, sprintf("%1.1e",mur), sprintf("%.1f",mud1), sprintf("%.1f",mufv1), sprintf("%.1f",mud2), sprintf("%.1f", mufv2), ts1, int12, ts2);
        colnames(parm) <- c("Sim", "NE", "sampleS", "murate", "neuDom", "neuFitness", 
                            "benDom", "benFitness", "sampleT1", "INT12", "sampleT2");
        
          for(i in 1:nsim){
                sink(file = paste0(folder, filename, '_', simnumb[i], '.slim'), type = 'output');
        
                      cat(paste0('initialize() {', '\n',
                                '\t', 'initializeMutationRate(', sprintf("%1.1e", mur[i]), ');', '\n',
                                '\t', 'initializeMutationType("m1", ', sprintf("%.1f", mud1[i]), ', "', mufd1, '", ', sprintf("%.1f", mufv1[i]), ');', '\n',
                                '\t', 'initializeMutationType("m2", ', sprintf("%.1f", mud2[i]), ', "', mufd2, '", ', sprintf("%.1f", mufv2[i]), ');', '\n',
                                '\t', 'm1.mutationStackPolicy = "f";', '\n',
                                '\t', 'm2.mutationStackPolicy = "f";', '\n',
                                '\t', 'm1.convertToSubstitution = T;', '\n',
                                '\t', 'm2.convertToSubstitution = T;', '\n',
                                '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
                                '\t', 'initializeGenomicElementType("g2", m1, 1.0);', '\n',
                                '\t', 'initializeGenomicElement(g1, ', gb1, ', ', gb2, ');', '\n',
                                '\t', 'initializeGenomicElement(g2, ', gb3, ', ', gb4, ');', '\n',
                                '\t', 'initializeRecombinationRate(c(1e-8, 0.5, 1e-8), c(', gb2, ', ', (gb2+1), ', ', gb4, '));', '\n',
                                '}', '\n\n'));
                               
                      cat(paste0('1 { sim.addSubpop("p0",', ne[i], '); ', '}', '\n\n'));
                               
                      cat(paste0((ne[i]*10), ' late() {', '\n',
                                '\t', 'dnovo = sample(sim.subpopulations.genomes, 1);', '\n',
                                '\t', 'if (dnovo.countOfMutationsOfType(m2) == 0)', '\n',
                                '\t\t', ' dnovo.addNewDrawnMutation(m2, 1);', '\n\n', 
                                
                                '\t', 'm2.convertToSubstitution = F;', '\n',
                                '}', '\n\n'));
                               
                      cat(paste0((ne[i]*10)+(ts1[i]-(1)), ' late() {', '\n',
                                '\t', 'm1.convertToSubstitution = F;', '\n',
                                '}', '\n\n'));
                      
                      cat(paste0((ne[i]*10)+(ts1[i]), ' late() {', '\n',
                                '\t', 'u = sim.subpopulations.individuals;' , '\n',
                                '\t', 'g = u.genomes;', '\n',
                                '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                                '\t', 'mp = m.position;', '\n',
                                '\t', 'mi = m.id;', '\n',
                                '\t', 'mc = m.selectionCoeff;', '\n\n',
                                
                                '\t', 'chr = ifelse(mp > 50000, paste("chr2"), paste("chr1"));', '\n',
                                '\t', 'mps = apply(mp, "if (applyValue > 50000) paste(applyValue - 50000); else paste(applyValue);");', '\n',
                                '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                                
                                '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"ID"+"\\t"+"selection"+"\\t"+"selCoeff"+"\\n");', '\n',
                                '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+sel+"\\t"+format("%.2f",mc)+"\\n");', '\n',
                                '\t', 'mutInf = header + mutInf;', '\n',
                                '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', simnumb[i], '_mutInfo_t1.txt", mutInf);', '\n\n',
                                
                                '\t', 'su = sample(u, ', n[i], ', F);', '\n\n',
                                
                                '\t', 'genos = NULL;', '\n',
                                '\t', 'for (si in su)','\n',
                                '\t', '{', '\n',
                                '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                                '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                                
                                '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                                '\t\t', 'genos = c(genos, genoLine);', '\n',
                                '\t', '}', '\n\n',
                                
                                '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                                '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', simnumb[i], '_genotypes_t1.txt", genomatrix);', '\n\n',
                                
                                '}', '\n\n'));
                               
                      cat(paste0((ne[i]*10)+(ts1[i]+(1)), ' late() {', '\n',
                                '\t', 'm1.convertToSubstitution = T;', '\n',
                                '}', '\n\n'));
            
                      cat(paste0((ne[i]*10)+(ts2[i]-(1)), ' late() {', '\n',
                                '\t', 'm1.convertToSubstitution = F;', '\n',
                                '}', '\n\n'));
                      
                      cat(paste0((ne[i]*10)+(ts2[i]), ' late() {', '\n',
                                 '\t', 'u = sim.subpopulations.individuals;' , '\n',
                                 '\t', 'g = u.genomes;', '\n',
                                 '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                                 '\t', 'mp = m.position;', '\n',
                                 '\t', 'mi = m.id;', '\n',
                                 '\t', 'mc = m.selectionCoeff;', '\n\n',
                                 
                                 '\t', 'chr = ifelse(mp > 50000, paste("chr2"), paste("chr1"));', '\n',
                                 '\t', 'mps = apply(mp, "if (applyValue > 50000) paste(applyValue - 50000); else paste(applyValue);");', '\n',
                                 '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                                 
                                 '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"ID"+"\\t"+"selection"+"\\t"+"selCoeff"+"\\n");', '\n',
                                 '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+sel+"\\t"+format("%.2f",mc)+"\\n");', '\n',
                                 '\t', 'mutInf = header + mutInf;', '\n',
                                 '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', simnumb[i], '_mutInfo_t2.txt", mutInf);', '\n\n',
                                 
                                 '\t', 'su = sample(u, ', n[i], ', F);', '\n\n',
                                 
                                 '\t', 'genos = NULL;', '\n',
                                 '\t', 'for (si in su)','\n',
                                 '\t', '{', '\n',
                                 '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                                 '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                                 
                                 '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                                 '\t\t', 'genos = c(genos, genoLine);', '\n',
                                 '\t', '}', '\n\n',
                                 
                                 '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                                 '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', simnumb[i], '_genotypes_t2.txt", genomatrix);', '\n\n',
                               
                                 '\t','sim.simulationFinished();', '\n', 
                               
                                 '}', '\n\n'));
                sink();
          };
        return(parm);
};