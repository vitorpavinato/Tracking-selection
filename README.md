## **Tracking Selection using Temporal Population Genomics Data**

**Authors:**

Vitor A. C. Pavinato, Stéphane de Mita &  Miguel Navascués

This repository contains the pipeline to run the simulations for the ABC-RF inference of demography and selection.

**INSTALLATION**

Dependencies:
- R packages [`moments`](https://cran.r-project.org/web/packages/moments/index.html) and [`ROCR`](https://ipa-tys.github.io/ROCR/);
- For parallel computing you also should install: `foreach`, `doParallel` and `parallel`;
- [SLiM 3.x](https://messerlab.org/slim/)
- [`bgzip 0.2.5 or higher`](http://www.htslib.org/download/)
- [`tabix 0.2.5 or higher`](http://www.htslib.org/download/)
- [`bcftools 1.6 or higher`](http://samtools.github.io/bcftools/)   
- `Python 2.7`
- [`egglib`](https://egglib.org)

**USAGE**
```
Rscript main.R
```

**CONFIGURATION**

In the main script you should speficy the path of the above mentioned dependencies. 
