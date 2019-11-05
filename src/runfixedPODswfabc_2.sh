#!/bin/bash

wfabc_folder=/Users/vitorpavinato/Softwares/PopGen/WFABC_v1.1/binaries/Mac
INFOLDER=/Users/vitorpavinato/My_repositories/Tracking-selection/results/pipeline_v5/fixedpods/neutral/results/wfabc_input

for i in {1..1}

do    
   ${wfabc_folder}/wfabc_2 -ploidy 2 -min_s 0.001 -max_s 1.0 ${INFOLDER}/wfabc_input_sample_${i}.txt
done
