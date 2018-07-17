#!/bin/bash

for i in {1..1000};
do
  echo '#!/bin/bash'                              > trackingselsim2.$i.sh
  echo '#$ -S /bin/bash'                          >> trackingselsim2.$i.sh;
  echo "#$ -N selsim2.$i"                          >> trackingselsim2.$i.sh;
  echo '#$ -q long.q'                             >> trackingselsim2.$i.sh;
  echo '#$ -cwd'                                  >> trackingselsim2.$i.sh;
  echo 'mkdir $TMPDIR/bin/'                       >> trackingselsim2.$i.sh;
  echo 'mkdir $TMPDIR/src/'                       >> trackingselsim2.$i.sh;
  echo 'mkdir $TMPDIR/src/models/'                       >> trackingselsim2.$i.sh;
  echo 'mkdir $TMPDIR/results/'                   >> trackingselsim2.$i.sh;
  echo 'cp -rp bin/* $TMPDIR/bin/'             	  >> trackingselsim2.$i.sh;
  echo 'cp -rp src/*.R $TMPDIR/src/'              >> trackingselsim2.$i.sh;
  echo 'cp -rp src/models/*.slim $TMPDIR/src/models/'    >> trackingselsim2.$i.sh;
  echo "mkdir batch.$i"					>> trackingselsim2.$i.sh;
  echo "cd batch.$i"					>> trackingselsim2.$i.sh;
  echo 'workingDir=$PWD'                          >> trackingselsim2.$i.sh;
  echo 'cd $TMPDIR'                               >> trackingselsim2.$i.sh;
  echo 'set -x'                                   >> trackingselsim2.$i.sh;
  echo 'Rscript ./src/main.R $RANDOM'             >> trackingselsim2.$i.sh;
  echo 'set +x'                                   >> trackingselsim2.$i.sh;
  echo 'cp -rp $TMPDIR/results/* $workingDir'   >> trackingselsim2.$i.sh;
  qsub trackingselsim2.$i.sh
done


