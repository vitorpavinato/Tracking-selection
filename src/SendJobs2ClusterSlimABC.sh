#!/bin/bash

for i in {1..50};
do
   echo "#!/bin/bash"                              > trackingselsim.$i.sh
   echo "#$ -S /bin/bash"                          >> trackingselsim.$i.sh;
   echo "#$ -N selsim.$i"                          >> trackingselsim.$i.sh;
   echo "#$ -q mem.q"                              >> trackingselsim.$i.sh;
   echo "#$ -cwd"                                  >> trackingselsim.$i.sh;
   echo 'workingDir=$PWD'                          >> trackingselsim.$i.sh;
   echo 'mkdir $TMPDIR/results/'                   >> trackingselsim.$i.sh;
   echo 'mkdir $TMPDIR/src/models/'                >> trackingselsim.$i.sh;
   echo 'mkdir $TMPDIR/bin/'                       >> trackingselsim.$i.sh;
   echo 'cp -rp results/* $TMPDIR/results/'        >> trackingselsim.$i.sh;
   echo 'cp -rp src/* $TMPDIR/src/'                >> trackingselsim.$i.sh;
   echo 'cp -rp bin/*.py $TMPDIR/bin/'             >> trackingselsim.$i.sh;
   echo 'cp -rp Tracking-selection.Rproj $TMPDIR/' >> trackingselsim.$i.sh;
   echo 'cd $TMPDIR'                               >> trackingselsim.$i.sh;
   echo 'set -x'                                   >> trackingselsim.$i.sh;
   echo "R --no-save --args directory \$TMPDIR seed $i$p DIYABC_exe_name \"/home/bin/Diyabc/2.1.0/x64/bin/general\" run_in_cluster \"c(T)\" project \"P$p\" simulated_target_data \"c(T)\" num_of_sims 1000000 proportion_of_sims_kept 0.01 num_of_points 50 prior_TAU_max 4 number_of_replicates 1000 true_gsm 0.$p prior_GSM_min 0.0 prior_GSM_max 1.0 prior_SNI_min 0.0 prior_SNI_max 0.0 motif 1 scenarios_number $i < ABCskylineplot.R" >> ABCskysim.$p.$i.sh;
   echo 'set +x'                                   >> trackingselsim.$i.sh;
   echo 'cp -rp $TMPDIR/* $workingDir/'            >> trackingselsim.$i.sh;
   qsub trackingselsim.$i.sh
done
