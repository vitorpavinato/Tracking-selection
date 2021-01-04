#!/bin/bash

for i in {1..1000};
do
  echo '#!/bin/bash'                                       > trackingselsim.$i.sh
  echo "#SBATCH -J selsim.$i"                             >> trackingselsim.$i.sh;
  echo "#SBATCH -o selsim.$i.o%A"                         >> trackingselsim.$i.sh;
  echo "#SBATCH -e selsim.$i.e%A"                         >> trackingselsim.$i.sh;
  echo '#SBATCH --partition=imag'                         >> trackingselsim.$i.sh;
  #echo '#SBATCH -t 00:00:00'                              >> trackingselsim.$i.sh;
  #echo '#SBATCH --mem=8G'                                 >> trackingselsim.$i.sh;
  #echo '#SBATCH -w muse003'                               >> trackingselsim.$i.sh;
  echo "mkdir batch.$i"                                   >> trackingselsim.$i.sh;
  echo "cd batch.$i"                                      >> trackingselsim.$i.sh;
  echo 'workingDir=$PWD'                                  >> trackingselsim.$i.sh;
  echo 'mkdir bin/'                                       >> trackingselsim.$i.sh;
  echo 'mkdir src/'                                       >> trackingselsim.$i.sh;
  echo 'mkdir src/models/'                                >> trackingselsim.$i.sh;
  echo 'mkdir results/'                                   >> trackingselsim.$i.sh;
  echo 'module purge'                                     >> trackingselsim.$i.sh;
  echo 'module load cv-standard R/3.4.3'                  >> trackingselsim.$i.sh; 
  echo 'module load cv-standard python/2.7.12'            >> trackingselsim.$i.sh;
  echo 'cd $SLURM_SUBMIT_DIR'                             >> trackingselsim.$i.sh;
  echo 'cp -rp bin/* $workingDir/bin/'                    >> trackingselsim.$i.sh;
  echo 'cp -rp src/*.R $workingDir/src/'                  >> trackingselsim.$i.sh;
  echo 'cp -rp src/models/*.slim $workingDir/src/models/' >> trackingselsim.$i.sh;
  echo 'cd $workingDir'                                   >> trackingselsim.$i.sh;
  echo 'set -x'                                           >> trackingselsim.$i.sh;
  echo 'Rscript ./src/main_bees.R $RANDOM'                >> trackingselsim.$i.sh;
  echo 'set +x'                                           >> trackingselsim.$i.sh;
  echo 'cp -rp $workingDir/results/* $workingDir'         >> trackingselsim.$i.sh; 
  echo 'rm -r $workingDir/bin/'                           >> trackingselsim.$i.sh;
  echo 'rm -r $workingDir/src/'	                          >> trackingselsim.$i.sh; 
  echo 'rm -r $workingDir/results/'                       >> trackingselsim.$i.sh;
  sbatch trackingselsim.$i.sh
done
