#!/bin/bash

# usage mergeVcfs.sh <<list of vcfs>> <<outpu file prefix>>

# process the vcf files in the list
for ID in $(cat $1 ) ;

do
# sort the vcf file
grep '^#' ${ID}.vcf > ${ID}_sorted.vcf && grep -v '^#' ${ID}.vcf | sort -k1,1 -k2,2n >> ${ID}_sorted.vcf ;

# bgzip the vcf file
bgzip -f ${ID}_sorted.vcf ;

# tabix the vcf file
tabix ${ID}_sorted.vcf.gz ;
done

# set the list of sample names to merge
modified_names=`awk '{ printf(" %s_sorted.vcf.gz", $1); }' $1`

# merge and get only the necessary field of the merged vcf file
bcftools merge --force-samples ${modified_names} |  bcftools query -f 'chr%CHROM\t%POS\t%REF,%ALT\tY[\t%GT]\n' > ${2}_merged.txt ;

# remove intermediate files
rm *.gz
rm *.tbi
