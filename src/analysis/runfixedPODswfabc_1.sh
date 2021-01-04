#!/bin/bash

#wfabc_folder=/Users/vitorpavinato/Softwares/PopGen/WFABC_v1.1/binaries/Mac
#INFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/neutral/results/wfabc_input
#OUTFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/neutral/results

wfabc_folder=/Users/correapavinato.1/Softwares/PopGen/WFABC_v1.1/binaries/Mac
INFOLDER=/Users/correapavinato.1/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/mutation_unlimited/add_s0.01/wfabc_input
OUTFOLDER=/Users/correapavinato.1/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/mutation_unlimited/add_s0.01

for i in {1..100}

do 
   ${wfabc_folder}/wfabc_1 -nboots 10000 -nthreads 8 ${INFOLDER}/wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print i, $1, $2}' >> ${OUTFOLDER}/wfabc_1_output_mutation_unlimited_.txt
done
