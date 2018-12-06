#!/bin/bash

for i in {1..2};
do
  echo '#!/bin/bash'                                      > trackingselsim.$i.sh
  echo "#SBATCH -J selsim.$i"                            >> trackingselsim.$i.sh;
  echo '#SBATCH -o output.out'                           >> trackingselsim.$i.sh;
  echo '#SBATCH -e error.out'                            >> trackingselsim.$i.sh;
  #echo '#SBATCH -t 01:00:00'                             >> trackingselsim.$i.sh;
  echo '#SBATCH --mem=4G'                                >> trackingselsim.$i.sh;
  echo '#SBATCH --mail-type=BEGIN,END,FAIL'              >> trackingselsim.$i.sh;
  echo 'module purge'                                    >> trackingselsim.$i.sh;
  echo 'module load system/R-3.5.1'                      >> trackingselsim.$i.sh;
  #echo 'module load bioinfo/SLiM-3.2'                    >> trackingselsim.$i.sh;
  #echo 'module load bioinfo/tabix-0.2.5'                 >> trackingselsim.$i.sh;
  #echo 'module load bioinfo/bcftools-1.6'                >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/bin/'                       	     >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/'                              >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/models/'                       >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/results/'                   	     >> trackingselsim.$i.sh;
  echo 'cp -rp bin/* $TMPDIR/bin/'             	  	     >> trackingselsim.$i.sh;
  echo 'cp -rp src/*.R $TMPDIR/src/'                     >> trackingselsim.$i.sh;
  echo 'cp -rp src/models/*.slim $TMPDIR/src/models/'    >> trackingselsim.$i.sh;
  echo "mkdir batch.$i"					                 >> trackingselsim.$i.sh;
  echo "cd batch.$i"					                 >> trackingselsim.$i.sh;
  echo 'workingDir=$PWD'                           	     >> trackingselsim.$i.sh;
  echo 'cd $TMPDIR'                               	     >> trackingselsim.$i.sh;
  echo 'set -x'                                   	     >> trackingselsim.$i.sh;
  echo 'Rscript ./src/main.R $RANDOM'             	     >> trackingselsim.$i.sh;
  echo 'set +x'                                   	     >> trackingselsim.$i.sh;
  echo 'cp -rp $TMPDIR/results/* $workingDir'   	     >> trackingselsim.$i.sh;
  sbatch trackingselsim.$i.sh
done
