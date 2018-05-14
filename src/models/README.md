## **FUTURE IMPLEMENTATIONS**
### **SLiM model**

###**Priorities**
- Implement the measurement of the hitchhiking influence on neutral markers;
- [X] Optimize the way slim output the mutation information - see item 15.3 SLiM manual:
  - plan A - re-implement the custom output - replace for loops for apply;
  - [X] plan B - use vcf output for individuals - try to order mutations by position;

**List of the last implementations**
- [X] Vcf MULTIALLELIC test: (done - if it is FALSE it removes ALL mutations when nm >1);
- [X] Implement the proportion of fragments (g1 and g2 elements) that have/not have benefitial mutations from prior;
- [X] Change the fragment size for 40Kbp to have more fragments and reduce the genome size to 135e+5 (but real simulation will be with 135e+6bp);
- [X] Re-difine genome architecture;
- [X] Implement the proportion of SNPs as being benefitial taken from the prior;
- [X] Fix the time for first sample - Time of equilibrium teq = t1;
- [X] Work on the SLiM code to allow the simulation of Ne where Ne < ns;
- [X] Sample theta and Ne from prior: sample theta and use it to define the population size for M-D equilibrium - Ne0, then change the population size to Ne1;
- [X] Re-implement output path in slim code;
- [X] Re-think the utility of "convertToSubstitution" of neutral mutations: 
  - [X] model_0
  - [X] model_1