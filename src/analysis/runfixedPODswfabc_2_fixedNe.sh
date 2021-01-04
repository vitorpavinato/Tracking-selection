#!/bin/bash

#need to run ABC-RF to estimate Ne2 then create a vector with Ne estimates and save it as .txt

wfabc_folder=/Users/vitorpavinato/Softwares/PopGen/WFABC_v1.1/binaries/Mac
INFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/mutation_limited/results/wfabc_input
#NE_FILE=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/mutation_limited/results/wfabc_1_output_mutation_limited.txt

for i in {1..1}

do 
   ne=`awk '{print $2}' ${NE_FILE} | tail -n +${i} | head -n 1`   
   ${wfabc_folder}/wfabc_2 -ploidy 2 -fixed_N ${ne} -min_s 0.001 -max_s 1.0 ${INFOLDER}/wfabc_input_sample_${i}.txt
done
