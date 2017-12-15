##########################################################################################
# toymodel(ne = N_size, filename = 'infile_N_size', folder = 'data/');                   #
#                                                                                        #
# This function creates the infile for a toymodel to run SLiM2 simulations.              #
# Toymodel helped the R/SLiM implementation                                              #                             
#																						                                             #
# ne = a vector containing population size N to be simulated;                            #
# filename = prefix for a filename system;                                               #
# folder = path where to save the generated SLiM2 infiles.                               #
##########################################################################################

toymodel <- function(ne, filename, folder){
  for (i in 1:length(ne)){
    sink(file = paste0(folder, filename, '_', ne[i], '.slim'), type = 'output');
    
    cat(paste('initialize() {', '\n',
              '\t', 'initializeMutationRate(1e-7);', '\n',
              '\t', 'initializeMutationType("m1", 0.5, "f", 0.0);', '\n',
              '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
              '\t', 'initializeGenomicElement(g1, 0, 999999);', '\n',
              '\t', 'initializeRecombinationRate(1e-8);', '\n',
              '}', '\n\n',
              
              '1 { sim.addSubpop("p1",', ne[i], '); }', '\n\n',
              ne[i]*10, 'late() { cat(sim.mutations.size() + "\\n"); }', '\n\n',
              
              '\t','sim.simulationFinished();'));
    
    sink();
  };
};

##########################################################################################
# MODEL 1 - HARD SWEEP                                         													 #
# 																						                                           #
# model1(ne = N_size, filename = 'infile_N_size', folder = 'data/');                     #
#                                                                                        #
# This function creates the infile for model 1 to run SLiM2 simulations.                 #
##########################################################################################

model1 <- function(ne, n=100,
                   mur=2.5e-8, 
                   mud1=0.5, mufd1="f", mufv1=0.0, 
                   mud2=0.5, mufd2="f", mufv2=0.1,
                   ts1=10, int12=10, ts2=ts1+int12,
                   filename, folder){
  
  for (i in 1:length(ne)){
    sink(file = paste0(folder, filename, '_', ne[i], '.slim'), type = 'output');
    
      cat(paste0('initialize() {', '\n',
                '\t', 'initializeMutationRate(', mur, ');', '\n',
                '\t', 'initializeMutationType("m1", ', mud1, ', "', mufd1, '", ', mufv1, ');', '\n',
                '\t', 'initializeMutationType("m2", ', mud2, ', "', mufd2, '", ', mufv2, ');', '\n',
                '\t', 'm2.mutationStackPolicy = "l";', '\n',
                '\t', 'm1.convertToSubstitution = T;', '\n',
                '\t', 'm2.convertToSubstitution = T;', '\n',
                '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
                '\t', 'initializeGenomicElementType("g2", m1, 1.0);', '\n',
                '\t', 'initializeGenomicElement(g1, 0, 49999);', '\n',
                '\t', 'initializeGenomicElement(g2, 50001, 99999);', '\n',
                '\t', 'initializeRecombinationRate(c(1e-8, 0.5, 1e-8), c(49999, 50000, 99999));', '\n',
                '}', '\n\n'));
               
      cat(paste0('1 { sim.addSubpop("p0",', ne[i], '); ', '}', '\n\n'));
               
      cat(paste0((ne[i]*10), ' late() {', '\n',
                '\t', 'dnovo = sample(sim.subpopulations.genomes, 1);', '\n',
                '\t', 'if (dnovo.countOfMutationsOfType(m2) == 0)', '\n',
                '\t\t', 'dnovo.addNewDrawnMutation(m2, 0);', '\n\n', 
                
                '\t', 'm2.convertToSubstitution = F;', '\n',
                '}', '\n\n'));
               
      cat(paste0((ne[i]*10)+(ts1-(1)), ' late() {', '\n',
                '\t', 'm1.convertToSubstitution = F;', '\n',
                '}', '\n\n'));
    
      cat(paste0((ne[i]*10)+(ts1), ' late() {', '\n',
                '\t', 'u = sim.subpopulations.individuals;' , '\n',
                '\t', 'g = u.genomes;', '\n',
                '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                '\t', 'mp = m.position;', '\n',
                '\t', 'mi = m.id;', '\n',
                '\t', 'mc = m.selectionCoeff;', '\n\n',
                
                '\t', 'chr = ifelse(mp > 50000, paste("chr2"), paste("chr1"));', '\n',
                '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                
                '\t', 'cat(paste("chrom"+"\\t"+"pos"+"\\t"+"id"+"\\t"+"selection"+"\\n"));', '\n',
                '\t', 'cat(paste(chr+"\\t"+format("%5d",mp)+"\\t"+format("%d",mi)+"\\t"+sel+"\\n"));', '\n\n',
                
                '\t', 'su = sample(u, ', n, ', F);', '\n\n',
                
                '\t', 'for (si in su){', '\n',
                '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n',
                '\t\t', 'cat(paste(hp1 + "" + hp2) + "\\n");', '\n',
                '\t', '}', '\n',
                '}', '\n\n'));
               
      cat(paste0((ne[i]*10)+(ts2-(1)), ' late() {', '\n',
               '\t', 'm1.convertToSubstitution = F;', '\n',
               '}', '\n\n'));
    
      cat(paste0((ne[i]*10)+(ts2), ' late() {', '\n',
                 '\t', 'u = sim.subpopulations.individuals;' , '\n',
                 '\t', 'g = u.genomes;', '\n',
                 '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                 '\t', 'mp = m.position;', '\n',
                 '\t', 'mi = m.id;', '\n',
                 '\t', 'mc = m.selectionCoeff;', '\n\n',
                 
                 '\t', 'chr = ifelse(mp > 50000, paste("chr2"), paste("chr1"));', '\n',
                 '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                 
                 '\t', 'cat(paste("chrom"+"\\t"+"pos"+"\\t"+"id"+"\\t"+"selection"+"\\n"));', '\n',
                 '\t', 'cat(paste(chr+"\\t"+format("%5d",mp)+"\\t"+format("%d",mi)+"\\t"+sel+"\\n"));', '\n\n',
                 
                 '\t', 'su = sample(u, ', n, ', F);', '\n\n',
                 
                 '\t', 'for (si in su){', '\n',
                 '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                 '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n',
                 '\t\t', 'cat(paste(hp1 + "" + hp2) + "\\n");', '\n',
                 '\t', '}', '\n\n',
               
                 '\t','sim.simulationFinished();', '\n\n', 
               
                 '}', '\n\n'));
    sink();
  };
}; ## ate aqui esta ok;

## TO DO or FURTER IMPLEMENTATIONS
# Allow user definition for all parameters and set default values (Miguel 8/11);
# Sample all the genomes/individuals
# Fixed number of individuals that are samples or it can vary given the N?
# Time between sample 1 and 2 can vary os it would be fixed?
# Be aware of drift effect - it is ok for neutral, but not ok for the benefitial.
## If the benefitial arises in t1, it should be present in t2 as well;
