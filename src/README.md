## **FUTURE IMPLEMENTATIONS**
### **Pipeline**

**Priorities**
- [X] Test how the pipeline handles/behaves SLiM vcf output with MULTIALLELIC=TRUE
- Implement this simulation values on the lss and gss reference table;
- Include a option to remove SNPs with MAF < 0.05 during data conversion;

**Thoughts about vcf MULTIALLELIC**
- Bioinformatics processing works OK (sort, bgzip, tabix and bcftools merge and get the data we need);
- The processing of the merged data in R also works fine - remember to comment the part that the pipeline remove redudancy;
- egglib summstats new version (the one that works with the new version of egglib: egglib3.0.0b21) still can't handle data with redundant mutations;
- I could not test the last part of the pipeline - prepare the locus-specific reference table;


**List of the last implementations**
- [X] Merge SLiM vcf outout with MULTIALLELIC=FALSE;


