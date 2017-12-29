##########################################################################################
# toymodel(ne = N_size, filename = 'infile_N_size', folder = 'data/');                   #
#                                                                                        #
# This function creates the infile for a toymodel to run SLiM2 simulations.              #
# Toymodel helped the initial stages of R/SLiM implementation                                              #                             
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
# model1(ne = N_size, filename = 'infile_slim', folder = 'infiles/');                    #
#                                                                                        #
# This function creates the infile for model 1 to run SLiM2 simulations.                 #
##########################################################################################

model1 <- function(nsims=100,
                   ne=100, n=100,
                   mur=2.5e-8, 
                   mud1=0.5, mufd1="f", mufv1=0.0, 
                   mud2=0.5, mufd2="f", mufv2=0.1,
                   genome=100000,
                   ts1=10, int12=10, ts2=ts1+int12,
                   filename, folder){
  
          gb1 = 0;
          gb2 = ((genome/2)-1);
          gb3 = ((genome/2)+1);
          gb4 = (genome-1);
  
          for (i in 1:nsims){
              sink(file = paste0(folder, filename, '_', i, '.slim'), type = 'output');
    
                cat(paste0('initialize() {', '\n',
                          '\t', 'initializeMutationRate(', mur, ');', '\n',
                          '\t', 'initializeMutationType("m1", ', mud1, ', "', mufd1, '", ', mufv1, ');', '\n',
                          '\t', 'initializeMutationType("m2", ', mud2, ', "', mufd2, '", ', mufv2, ');', '\n',
                          '\t', 'm2.mutationStackPolicy = "l";', '\n',
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
                          '\t\t', ' dnovo.addNewDrawnMutation(m2, 0);', '\n\n', 
                          
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
                          '\t', 'mps = apply(mp, "if (applyValue > 50000) paste(applyValue - 50000); else paste(applyValue);");', '\n',
                          '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                          
                          '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"ID"+"\\t"+"selection"+"\\t"+"selCoeff"+"\\n");', '\n',
                          '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+sel+"\\t"+format("%.2f",mc)+"\\n");', '\n',
                          '\t', 'mutInf = header + mutInf;', '\n',
                          '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', i, '_mutInfo_t1.txt", mutInf);', '\n\n',
                          
                          '\t', 'su = sample(u, ', n, ', F);', '\n\n',
                          
                          '\t', 'genos = NULL;', '\n',
                          '\t', 'for (si in su)','\n',
                          '\t', '{', '\n',
                          '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                          '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                          
                          '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                          '\t\t', 'genos = c(genos, genoLine);', '\n',
                          '\t', '}', '\n\n',
                          
                          '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                          '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', i, '_genotypes_t1.txt", genomatrix);', '\n\n',
                          
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
                           '\t', 'mps = apply(mp, "if (applyValue > 50000) paste(applyValue - 50000); else paste(applyValue);");', '\n',
                           '\t', 'sel = ifelse(mc != 0, paste("Y"), paste("N"));', '\n\n',
                           
                           '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"ID"+"\\t"+"selection"+"\\t"+"selCoeff"+"\\n");', '\n',
                           '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+sel+"\\t"+format("%.2f",mc)+"\\n");', '\n',
                           '\t', 'mutInf = header + mutInf;', '\n',
                           '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', i, '_mutInfo_t2.txt", mutInf);', '\n\n',
                           
                           '\t', 'su = sample(u, ', n, ', F);', '\n\n',
                           
                           '\t', 'genos = NULL;', '\n',
                           '\t', 'for (si in su)','\n',
                           '\t', '{', '\n',
                           '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                           '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                           
                           '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                           '\t\t', 'genos = c(genos, genoLine);', '\n',
                           '\t', '}', '\n\n',
                           
                           '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                           '\t', 'writeFile("', getwd(), '/data/', 'output_slim_', i, '_genotypes_t2.txt", genomatrix);', '\n\n',
                         
                           '\t','sim.simulationFinished();', '\n', 
                         
                           '}', '\n\n'));
              sink();
          };
};
