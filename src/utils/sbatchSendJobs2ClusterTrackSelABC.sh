#!/bin/bash

for i in {2..10};
do
  echo '#!/bin/bash'                                       > trackingselsim.$i.sh
  echo "#SBATCH -J selsim.$i"                             >> trackingselsim.$i.sh;
  echo "#SBATCH -o /fs/scratch/PAS1715/TrackSel_analysis/Tracking-selection/slurmoutput/selsim.$i.o%A"      >> trackingselsim.$i.sh;
  echo "#SBATCH -e /fs/scratch/PAS1715/TrackSel_analysis/Tracking-selection/slurmoutput/selsim.$i.e%A"      >> trackingselsim.$i.sh;
  echo '#SBATCH -t 168:00:00'                              >> trackingselsim.$i.sh;
  echo '#SBATCH --mem=8G'                                 >> trackingselsim.$i.sh;
  echo '#SBATCH --account PAS1715'                         >> trackingselsim.$i.sh;
  echo "mkdir batch.$i"                                   >> trackingselsim.$i.sh;
  echo "cd batch.$i"                                      >> trackingselsim.$i.sh;
  echo 'workingDir=$PWD'                                  >> trackingselsim.$i.sh;
  echo 'mkdir bin/'                                       >> trackingselsim.$i.sh;
  echo 'mkdir src/'                                       >> trackingselsim.$i.sh;
  echo 'mkdir -p src/models/'                             >> trackingselsim.$i.sh;
  echo 'mkdir results/'                                   >> trackingselsim.$i.sh;
  echo 'module purge'                                     >> trackingselsim.$i.sh;
  echo 'module load R/3.6.3-gnu9.1'                       >> trackingselsim.$i.sh;  
  echo 'cd $SLURM_SUBMIT_DIR'                             >> trackingselsim.$i.sh;
  echo 'cp -rp bin/bcftools $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp bin/bgzip $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp bin/tabix $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp bin/python $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp bin/slim $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp bin/summstats.py $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp src/proof/*.R $workingDir/src/'                  >> trackingselsim.$i.sh;
  echo 'cp -rp src/models/proof/*.slim $workingDir/src/models/' >> trackingselsim.$i.sh;
  echo 'cp -rp renv/ $workingDir'                                  >> trackingselsim.$i.sh;
  echo 'cp -rp renv.lock $workingDir'                    >> trackingselsim.$i.sh;
  echo 'cd $workingDir'                                   >> trackingselsim.$i.sh;
  echo 'set -x'                                           >> trackingselsim.$i.sh;
  echo 'Rscript ./src/main.R $RANDOM'                     >> trackingselsim.$i.sh;
  echo 'set +x'                                           >> trackingselsim.$i.sh;
  echo 'cp -rp results/reference_table.RData ./'         >> trackingselsim.$i.sh; 
  echo 'rm -r bin/'                                      >> trackingselsim.$i.sh;
  echo 'rm -r renv/'	                                 >> trackingselsim.$i.sh; 
  echo 'rm -r results/'                                  >> trackingselsim.$i.sh;
  echo 'rm renv.lock'                                     >> trackingselsim.$i.sh;
  sbatch trackingselsim.$i.sh
done
