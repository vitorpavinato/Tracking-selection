#!/bin/bash

for i in {1..1000};
do
  echo '#!/bin/bash'                                      > trackingselsim.$i.sh
  echo '#$ -S /bin/bash'                          	     >> trackingselsim.$i.sh;
  echo "#$ -N selsim.$i"                          	     >> trackingselsim.$i.sh;
  echo '#$ -q long.q'                             	     >> trackingselsim.$i.sh;
  echo '#$ -cwd'                                  	     >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/bin/'                       	     >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/'                              >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/src/models/'                       >> trackingselsim.$i.sh;
  echo 'mkdir $TMPDIR/results/'                   	     >> trackingselsim.$i.sh;
  echo 'cp -rp bin/* $TMPDIR/bin/'             	  	     >> trackingselsim.$i.sh;
  echo 'cp -rp src/*.R $TMPDIR/src/'                     >> trackingselsim.$i.sh;
  echo 'cp -rp src/models/*.slim $TMPDIR/src/models/'    >> trackingselsim.$i.sh;
  echo "mkdir batch.$i"					                         >> trackingselsim.$i.sh;
  echo "cd batch.$i"					                           >> trackingselsim.$i.sh;
  echo 'workingDir=$PWD'                           	     >> trackingselsim.$i.sh;
  echo 'cd $TMPDIR'                               	     >> trackingselsim.$i.sh;
  echo 'set -x'                                   	     >> trackingselsim.$i.sh;
  echo 'Rscript ./src/main.R $RANDOM'             	     >> trackingselsim.$i.sh;
  echo 'set +x'                                   	     >> trackingselsim.$i.sh;
  echo 'cp -rp $TMPDIR/results/* $workingDir'   	       >> trackingselsim.$i.sh;
  qsub trackingselsim.$i.sh
done


