#!/bin/bash

for i in {1..3};
do
  echo '#PBS -A PAS1554'                                 > trackingselsim.$i.sh;
  echo "#PBS -N selsim.$i"                               >> trackingselsim.$i.sh;
  echo '#PBS -l walltime=168:00:00'                      >> trackingselsim.$i.sh;
  echo '#PBS -l nodes=1:ppn=28'                          >> trackingselsim.$i.sh;
  echo '#PBS -j oe'                                      >> trackingselsim.$i.sh;
  echo 'module load intel/16.0.8'                        >> trackingselsim.$i.sh;
  echo 'module load R/3.4.2'                             >> trackingselsim.$i.sh;
  echo 'module load slim/3.2.1'                          >> trackingselsim.$i.sh;
  echo 'module load hstlib/1.9.2'                        >> trackingselsim.$i.sh;
  echo 'module load bcftools/1.9.2'                      >> trackingselsim.$i.sh;
  echo 'module load python/2.7'                          >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/bin/'                       	     >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/'                              >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/models/'                       >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/results/'                   	     >> trackingselsim.$i.sh;
  echo 'cd $PBS_O_WORKDIR'                               >> trackingselsim.$i.sh;
  echo 'cp -rp bin/* $TMPDIR/bin/'             	  	     >> trackingselsim.$i.sh;
  echo 'cp -rp src/*.R $TMPDIR/src/'                     >> trackingselsim.$i.sh;
  echo 'cp -rp src/models/*.slim $TMPDIR/src/models/'    >> trackingselsim.$i.sh;
  echo "mkdir batch.$i"					                 >> trackingselsim.$i.sh;
  echo "cd batch.$i"					                 >> trackingselsim.$i.sh;
  echo 'workingDir=$PWD'                           	     >> trackingselsim.$i.sh;
  echo 'cd $TMPDIR'                               	     >> trackingselsim.$i.sh;
  echo 'set -x'                                   	     >> trackingselsim.$i.sh;
  echo 'Rscript ./src/main_bee.R $RANDOM'             	 >> trackingselsim.$i.sh;
  echo 'set +x'                                   	     >> trackingselsim.$i.sh;
  echo 'cp -rp $TMPDIR/results/* $workingDir'   	     >> trackingselsim.$i.sh;
  qsub trackingselsim.$i.sh
done
