#!/bin/bash

#for i in {1..5}

for i in $(ls -dl batch.*/ | cut -d"." -f2 | rev | cut -c 2- | rev | uniq) 

do 
   /home/pavinato/Softwares/WFABC_v1.1/binaries/Linux/wfabc_1 batch.${i}/wfabc_input/wfabc_input_sample_${i}.txt | awk -v i="$i" 'BEGIN{FS="\t"} {print "batch."i, $1}' >> wfabc_1_output.txt
done
