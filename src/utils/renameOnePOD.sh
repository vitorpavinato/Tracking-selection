#!/bin/bash

# usage: renameOnePOD.sh <original name> <new name>

rename "s/_$1/_$2/" batch.$2/slim_output/*.tree
rename "s/_$1/_$2/" batch.$2/slim_output/*.vcf
rename "s/_$1/_$2/" batch.$2/slim_output/*_sorted.vcf.gz
rename "s/_$1/_$2/" batch.$2/slim_output/*_sorted.vcf.gz.tbi
rename "s/_$1/_$2/" batch.$2/slim_output/*.txt
rename "s/_$1/_$2/" batch.$2/egglib_input/*.txt
rename "s/_$1/_$2/" batch.$2/egglib_output/*.txt
rename "s/_$1/_$2/" batch.$2/wfabc_input/*.txt
