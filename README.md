## **Tracking Selection using Temporal Population Genomics Data**

**Authors:**

Vitor A. C. Pavinato, Stéphane de Mita &  Miguel Navascués

This repository contains the implementation of a pipeline to run the simulations and to produce a reference table for the ABC-RF inference of demography and selection.

**CITATION**
```
@article{PAvinato:2021,
author = {Pavinato, Vitor A.C. and de Mita, S{\'e}phane and Marin, Jean-Michel and tNavascu{\'e}s, Miguel},
title = {{Joint inference of adaptive and demographic history from temporal population genomic data}},
journal = {bioRxiv},
year = {2021},
volume = {},
pages = {},
month = jan
}
```

**INSTALLATION**

Dependencies:
- R packages [`moments`](https://cran.r-project.org/web/packages/moments/index.html) and [`ROCR`](https://ipa-tys.github.io/ROCR/);
- For parallel computing you also should install: `foreach`, `doParallel` and `parallel`;
- [SLiM 3.x](https://messerlab.org/slim/)
- [`bgzip 0.2.5 or higher`](http://www.htslib.org/download/)
- [`tabix 0.2.5 or higher`](http://www.htslib.org/download/)
- [`bcftools 1.6 or higher`](http://samtools.github.io/bcftools/)   
- `Python 2.7`
- [`EggLib`](https://egglib.org)

**USAGE**
```
Rscript src/proof/main.R
```

**CONFIGURATION**

In the main script you should speficy the path of the above mentioned dependencies and other configurations that might be necessary (prior range, genome size, `\tau`, etc). 
