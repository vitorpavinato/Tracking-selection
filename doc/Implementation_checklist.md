## **FUTURE IMPLEMENTATIONS**

**Priorities**
- fix the time for first sample - Time of equilibrium teq = t1: [done]
- Implement the measurement of the hitchhiking influence on neutral markers:
- Implement the proportion of SNPs as being benefitial taken from the prior: [done]
- Implement the proportion of fragments (g1 and g2 elements) that have/not have benefitial mutations from prior:
- Implement this simulation values on the lss and gss reference table:
- re-implement output path in slim code: [done]
- remove SNPs with MAF < 0.05 during data conversion:
- re-think the utility of "convertToSubstitution" of neutral mutations: - model_0 [done]
																	    - model_1 [done]
- optimize the way slim output the mutation information - see item 15.3 SLiM manual:
	- strategy 1 - re-implement the custom output - replace for loops for apply;
	- strategy 2 - use vcf output for individuals - try to order mutations by position;

**List of the last implementations
- Sample theta and Ne from prior: sample theta and use it to define the population size for M-D equilibrium - Ne0, then change the population size to Ne1: done
- Work on the SLiM code to allow the simulation of Ne where Ne < ns: done
- Redifine genome architecture: done
