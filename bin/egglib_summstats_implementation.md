## **egglib summstats.py**

**all v0 bugs were fixed in v1**

**List of bugs v0:**
- SFS calculation;
- select=list error;
- select=rand and select-num error;
- Some GSS are repeated before GSS columns in the output;
- Mutation with same position caused a problem in egglib;

	- Mutations in same position were likely with the output code for the first version of the pipeline; however one of the redundanct position was removed during the egglib file conversion in R (pipeline v1 and v2);
	- With the new code to produce the ouputs it is unlikely to have same IDs in the sample; but if it is the case, they can be removed with the R code;
	- If we are going to use the vcf-style slim output, if we use the argument MULTIALLELIC=F, BOTH mutation with same position are removed.