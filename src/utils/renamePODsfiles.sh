#!/bin/bash

for i in {1001..2000}

do
  rename "s/_1/_$i/" batch.$i/slim_output/*.tree
  rename "s/_1/_$i/" batch.$i/slim_output/*.vcf
  rename "s/_1/_$i/" batch.$i/slim_output/*_sorted.vcf.gz
  rename "s/_1/_$i/" batch.$i/slim_output/*_sorted.vcf.gz.tbi
  rename "s/_1/_$i/" batch.$i/slim_output/*.txt
  rename "s/_1/_$i/" batch.$i/egglib_input/*.txt
  rename "s/_1/_$i/" batch.$i/egglib_output/*.txt
  rename "s/_1/_$i/" batch.$i/wfabc_input/*.txt
done
