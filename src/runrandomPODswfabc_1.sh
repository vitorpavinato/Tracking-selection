#!/bin/bash

wfabc_folder=/Users/vitorpavinato/Softwares/PopGen/WFABC_v1.1/binaries/Mac
INFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/Cluster_data/pooled_reftable_pods/copy_all_wfabc_inputs
OUTFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/Cluster_data/pooled_reftable_pods/results

for i in {1..2000}

#for i in $(ls -dl batch.*/ | cut -d"." -f2 | rev | cut -c 2- | rev | uniq) # used to produce the first version that contain the batch* in the row name**

do 
   #/home/pavinato/Softwares/WFABC_v1.1/binaries/Linux/wfabc_1 batch.${i}/wfabc_input/wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print "batch."i, $1}' >> wfabc_1_output.txt # same**
   ${wfabc_folder}/wfabc_1 -nboots 10000 -nthreads 8 ${INFOLDER}/wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print i, $1, $2}' >> ${OUTFOLDER}/wfabc_1_output_pods_new.txt
done
