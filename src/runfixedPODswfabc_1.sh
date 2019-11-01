#!/bin/bash

wfabc_folder=/Users/vitorpavinato/Softwares/PopGen/WFABC_v1.1/binaries/Mac
INFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods_2/mutation_unlimited/results/wfabc_input
OUTFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods_2/mutation_unlimited/results

for i in {1..100}

do 
   ${wfabc_folder}/wfabc_1 ${INFOLDER}/wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print i, $1}' >> ${OUTFOLDER}/wfabc_1_output.txt
done
