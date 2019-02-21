#!/bin/bash

for i in {1..100}

do 
   /home/pavinato/Softwares/WFABC_v1.1/binaries/Linux/wfabc_1 wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print "batch."i, $1}' >> wfabc_1_output.txt
done
