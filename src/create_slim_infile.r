##########################################################################################
# toymodel(filename, folder);             #
#                                                                                        #
# This function creates the infile for a toymodel to run SLiM2 simulations.              #
# Toymodel helped the initial stages of R/SLiM implementation                            #                             
#																						                                             #
# ne = a vector containing population size N to be simulated;                            #
# ne_min and ne_max are the defined range values for random sampling taking from uniform #
# filename = prefix for a filename system;                                               #
# folder = path where to save the generated SLiM2 infiles.                               #
##########################################################################################

toymodel <- function(ne = 1000, mur = 1e-7, 
                     mud = 0.5, mufd = "f", mufv = 0.0, 
                     genomesize = 1000000, rrate = 1e-8, 
                     filename, folder)
{
        gb1 = 0; gb2 = (genomesize-1);
        
        parm <- as.data.frame(cbind(ne, mur, mud, mufv));
        
        for (i in 1:dim(parm)[1]){
         	sink(file = paste0(folder, filename, '_', i, '.slim'), type = 'output');
          
                cat(paste0('initialize() {', '\n',
                          '\t', 'initializeMutationRate(', sprintf("%1.1e", parm$mur[i]), ');', '\n',
                          '\t', 'initializeMutationType("m1", ', sprintf("%.1f", parm$mud[i]), ', "', mufd, '", ', sprintf("%.3f", parm$mufv[i]), ');', '\n',
                          '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
                          '\t', 'initializeGenomicElement(g1, ', gb1, ', ',  gb2, ');', '\n',
                          '\t', 'initializeRecombinationRate(', sprintf("%1.1e", rrate), ');', '\n',
                          '}', '\n\n',
            
                          '1 { sim.addSubpop("p1",', parm$ne[i], '); }', '\n\n',
                          parm$ne[i]*10, ' late() { cat(sim.mutations.size() + "\\n");', '\n\n',
            
                          '\t','sim.simulationFinished();', '\n\n',
                          
                          '}'));
            sink();
        };
        return(parm);
};

##########################################################################################
# MODEL 1 - SIMPLE HARD SWEEP TEMPORAL SAMPLES                                       		 #
# 																						                                           #
# model1(filename, folder, output, folderout);                                           #
#                                                                                        #
# This function creates the infile for model 1 to run SLiM2 simulations.                 #
##########################################################################################

model1 <- function(ne = 500, n = 100, mur = 1e-7,
                   mud1 = 0.5, mufd1 = "f", mufv1 = 0.0,
                   mud2 = 0.5, mufd2 = "f", mufv2 = 0.5,
                   ts1 = 25, int12 = 25, 
                   genomesize = 100000, rrate = 1e-8,
                   filename, folder, output, folderout)
{
        gb1 = 0; gb2 = ((genomesize/2)-1); gb3 = ((genomesize/2)+1); gb4 = (genomesize-1);
        
        parm <- as.data.frame(cbind(ne, n, mur, mud1, mufv1, mud2, mufv2, ts1, int12, genomesize, rrate));
        
        for(i in 1:dim(parm)[1]){
            sink(file = paste0(folder, filename, '_', i, '.slim'), type = 'output');
        
                cat(paste0('initialize() {', '\n',
                          '\t', 'initializeMutationRate(', sprintf("%1.1e", parm$mur[i]), ');', '\n',
                          '\t', 'initializeMutationType("m1", ', sprintf("%.1f", parm$mud1[i]), ', "', mufd1, '", ', sprintf("%.3f", parm$mufv1[i]), ');', '\n',
                          '\t', 'initializeMutationType("m2", ', sprintf("%.1f", parm$mud2[i]), ', "', mufd2, '", ', sprintf("%.3f", parm$mufv2[i]), ');', '\n',
                          '\t', 'm1.mutationStackPolicy = "s";', '\n',
                          '\t', 'm2.mutationStackPolicy = "s";', '\n',
                          '\t', 'm1.convertToSubstitution = T;', '\n',
                          '\t', 'm2.convertToSubstitution = T;', '\n',
                          '\t', 'initializeGenomicElementType("g1", m1, 1.0);', '\n',
                          '\t', 'initializeGenomicElementType("g2", m1, 1.0);', '\n',
                          '\t', 'initializeGenomicElement(g1, ', gb1, ', ', gb2, ');', '\n',
                          '\t', 'initializeGenomicElement(g2, ', gb3, ', ', gb4, ');', '\n',
                          '\t', 'initializeRecombinationRate(c(', sprintf("%1.1e", rrate), ', 0.5, ', sprintf("%1.1e", rrate), '), c(', gb2, ', ', (gb2+1), ', ', gb4, '));', '\n',
                          '}', '\n\n'));
                         
                cat(paste0('1 { sim.addSubpop("p0",', parm$ne[i], '); ', '}', '\n\n'));
                         
                cat(paste0((parm$ne[i]*10), ' late() {', '\n',
                          '\t', 'dnovo = sample(sim.subpopulations.genomes, 1);', '\n',
                          '\t', 'if (dnovo.countOfMutationsOfType(m2) == 0)', '\n',
                          '\t\t', ' dnovo.addNewDrawnMutation(m2, 1);', '\n\n', 
                          
                          '\t', 'm2.convertToSubstitution = F;', '\n',
                          '}', '\n\n'));
                         
                cat(paste0((parm$ne[i]*10)+(parm$ts1[i]-(1)), ' late() {', '\n',
                          '\t', 'm1.convertToSubstitution = F;', '\n',
                          '}', '\n\n'));
                
                cat(paste0((parm$ne[i]*10)+(parm$ts1[i]), ' late() {', '\n',
                          '\t', 'u = sim.subpopulations.individuals;' , '\n',
                          '\t', 'g = u.genomes;', '\n',
                          '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                          '\t', 'mp = m.position;', '\n',
                          '\t', 'mi = m.id;', '\n',
                          '\t', 'mf = sim.mutationFrequencies(subpops=NULL,mutations=m);', '\n',
                          '\t', 'mc = m.selectionCoeff;', '\n\n',
                          
                          '\t', 'chr = ifelse(mp > ', genomesize/2, ', paste("chr2"), paste("chr1"));', '\n',
                          '\t', 'mps = sapply(mp, "if (applyValue > ', genomesize/2, ') paste(applyValue - ', genomesize/2, '); else paste(applyValue);");', '\n',
                          '\t', 'status = repEach(paste("NS"), size(m));', '\n',
                          '\t', 'alleles = repEach(paste("A,T"), size(m));', '\n\n',
                          
                          '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"muID"+"\\t"+"mufreq"+"\\t"+"status"+"\\t"+"selCoeff"+"\\t"+"alleles"+"\\n");', '\n',
                          '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+format("%.3f",mf)+"\\t"+status+"\\t"+format("%.3f",mc)+"\\t"+alleles+"\\n");', '\n',
                          '\t', 'mutInf = header + mutInf;', '\n',
                          '\t', 'writeFile("', getwd(), '/', folderout, '/', output, '_', i, '_mutInfo_t1.txt", mutInf);', '\n\n',
                          
                          '\t', 'su = sample(u, ', parm$n[i], ', F);', '\n\n',
                          
                          '\t', 'genos = NULL;', '\n',
                          '\t', 'for (si in su)','\n',
                          '\t', '{', '\n',
                          '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                          '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                          
                          '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                          '\t\t', 'genos = c(genos, genoLine);', '\n',
                          '\t', '}', '\n\n',
                          
                          '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                          '\t', 'writeFile("', getwd(), '/', folderout, '/', output, '_', i, '_genotypes_t1.txt", genomatrix);', '\n\n',
                          
                          '}', '\n\n'));
                         
                cat(paste0((parm$ne[i]*10)+(parm$ts1[i]+(1)), ' late() {', '\n',
                          '\t', 'm1.convertToSubstitution = T;', '\n',
                          '}', '\n\n'));
            
                cat(paste0((parm$ne[i]*10)+(parm$ts1[i]+parm$int12[i]-(1)), ' late() {', '\n',
                          '\t', 'm1.convertToSubstitution = F;', '\n',
                          '}', '\n\n'));
                
                cat(paste0((parm$ne[i]*10)+(parm$ts1[i]+parm$int12[i]), ' late() {', '\n',
                           '\t', 'u = sim.subpopulations.individuals;' , '\n',
                           '\t', 'g = u.genomes;', '\n',
                           '\t', 'm = sortBy(unique(g.mutations), "position");' , '\n', 
                           '\t', 'mp = m.position;', '\n',
                           '\t', 'mi = m.id;', '\n',
                           '\t', 'mf = sim.mutationFrequencies(subpops=NULL,mutations=m);', '\n',
                           '\t', 'mc = m.selectionCoeff;', '\n\n',
                           
                           '\t', 'chr = ifelse(mp > ', genomesize/2, ', paste("chr2"), paste("chr1"));', '\n',
                           '\t', 'mps = sapply(mp, "if (applyValue > ', genomesize/2, ') paste(applyValue - ', genomesize/2, '); else paste(applyValue);");', '\n',
                           '\t', 'status = repEach(paste("NS"), size(m));', '\n',
                           '\t', 'alleles = repEach(paste("A,T"), size(m));', '\n\n',
                           
                           '\t', 'header = paste("chrom"+"\\t"+"position"+"\\t"+"muID"+"\\t"+"mufreq"+"\\t"+"status"+"\\t"+"selCoeff"+"\\t"+"alleles"+"\\n");', '\n',
                           '\t', 'mutInf = paste(chr+"\\t"+mps+"\\t"+format("%d",mi)+"\\t"+format("%.3f",mf)+"\\t"+status+"\\t"+format("%.3f",mc)+"\\t"+alleles+"\\n");', '\n',
                           '\t', 'mutInf = header + mutInf;', '\n',
                           '\t', 'writeFile("', getwd(), '/', folderout, '/', output, '_', i, '_mutInfo_t2.txt", mutInf);', '\n\n',
                           
                           '\t', 'su = sample(u, ', parm$n[i], ', F);', '\n\n',
                           
                           '\t', 'genos = NULL;', '\n',
                           '\t', 'for (si in su)','\n',
                           '\t', '{', '\n',
                           '\t\t', 'hp1 = (asInteger(match(m,sortBy(unique(si.genomes[0].mutations), "position")) >=0)+1);', '\n',
                           '\t\t', 'hp2 = (asInteger(match(m,sortBy(unique(si.genomes[1].mutations), "position")) >=0)+1);', '\n\n',
                           
                           '\t\t', 'genoLine = paste(hp1 + "" + hp2, "\\t");', '\n',
                           '\t\t', 'genos = c(genos, genoLine);', '\n',
                           '\t', '}', '\n\n',
                           
                           '\t', 'genomatrix = paste(genos, "\\n");', '\n',
                           '\t', 'writeFile("', getwd(), '/', folderout, '/', output, '_', i, '_genotypes_t2.txt", genomatrix);', '\n\n',
                         
                           '\t','sim.simulationFinished();', '\n', 
                         
                           '}', '\n\n'));
            sink();
        };
        return(parm);
};