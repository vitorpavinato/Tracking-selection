##########################################################################################
# toymodel(parm = N_size, filename = 'infile_N_size', folder = 'data/');                 #
#                                                                                        #
# This function creates the infile for a toymodel to run SLiM2 simulations.              #
# Toymodel helped the R/SLiM implementation                                              #                             
#																						 #
# parm = a vector containing population size N to be simulated;                          #
# filename = prefix for a filename system;                                               #
# folder = path where to save the generated SLiM2 infiles.                               #
##########################################################################################

toymodel <- function(parm, filename, folder){
  for (i in 1:length(parm)){
    sink(file = paste0(folder, filename, '_', parm[i], '.slim'), type = 'output');
    
    cat(paste('initialize() {', '\n',
              '\t', 'initializeMutationRate(1e-7);', '\n',
              '\t', 'initializeMutationType("m1", 0.5, "f", 0.0);', '\n',
              '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
              '\t', 'initializeGenomicElement(g1, 0, 999999);', '\n',
              '\t', 'initializeRecombinationRate(1e-8);', '\n',
              '}', '\n\n',
              
              '1 { sim.addSubpop("p1",', parm[i], '); }', '\n\n',
              parm[i]*10, 'late() { cat(sim.mutations.size() + "\\n"); }', '\n\n',
              
              '\t','sim.simulationFinished();'));
    
    sink();
  };
};

##########################################################################################
# MODEL 1 - HARD SWEEP with NO LINKAGE 													 #
# 																						 #
# model1(parm = N_size, filename = 'infile_N_size', folder = 'data/');                   #
#                                                                                        #
# This function creates the infile for model 1 to run SLiM2 simulations.                 #
#                                                                                        #
# parm = a vector containing population size N to be simulated;                          #
# filename = prefix for a filename system;                                               #
# folder = path where to save the generated SLiM2 infiles.                               #
##########################################################################################

model1 <- function(parm, filename, folder){
  for (i in 1:length(parm)){
    sink(file = paste0(folder, filename, '_', parm[i], '.slim'), type = 'output');
    
    cat(paste0('initialize() {', '\n',
               '\t', 'initializeMutationRate(1e-7);', '\n',
               '\t', 'initializeMutationType("m1", 0.5, "f", 0.0);', '\n',
               '\t', 'initializeMutationType("m2", 0.5, "f", 0.8);', '\n',
               '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
               '\t', 'initializeGenomicElementType("g2", m1, 1.0);', '\n',
               '\t', 'initializeGenomicElement(g1, 0, 49999);', '\n',
               '\t', 'initializeGenomicElement(g2, 50001, 99999);', '\n',
               '\t', 'initializeRecombinationRate(c(1e-8, 0.5, 1e-8), c(49999, 50000, 99999));', '\n',
               '}', '\n\n',
               
               '1 { sim.addSubpop("p1",', parm[i], '); }', '\n\n',
               
               (parm[i]*10), 'late() {', '\n',
               '\t', 'mutation = sample(p1.genomes, 1);', '\n',
               '\t', 'mutation.addNewDrawnMutation(m2, 0); }', '\n\n',
               
               (parm[i]*10)+10, 'late() {', '\n',
               '\t', 'allIndividuals = sim.subpopulations.individuals;', '\n',
               '\t', 'sampledIndividuals = sample(allIndividuals, 100);', '\n',
               '\t', 'sampledIndividuals.genomes.outputVCF(filePath = NULL, outputMultiallelics = F); }', '\n\n',
               
               (parm[i]*10)+20, 'late() {', '\n',
               '\t', 'allIndividuals = sim.subpopulations.individuals;', '\n',
               '\t', 'sampledIndividuals = sample(allIndividuals, 100);', '\n',
               '\t', 'sampledIndividuals.genomes.outputVCF(filePath = NULL, outputMultiallelics = F); }', '\n\n',
               
               '\t','sim.simulationFinished();'));
    
    sink();
  };
};

## TO DO or FURTER IMPLEMENTATIONS
# Allow user definition for all parameters and set default values (Miguel 8/11);
# Sample all the genomes/individuals
# Fixed number of individuals that are samples or it can vary given the N?
# Time between sample 1 and 2 can vary os it would be fixed?
# Be aware of drift effect - it is ok for neutral, but not ok for the benefitial.
 ## If the benefitial arises in t1, it should be present in t2 as well.
