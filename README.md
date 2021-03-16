## **Tracking Selection using Temporal Population Genomics Data**

[![DOI](https://zenodo.org/badge/113848735.svg)](https://zenodo.org/badge/latestdoi/113848735)


**Authors:**

Vitor A. C. Pavinato, Stéphane de Mita &  Miguel Navascués

This repository contains the implementation of a pipeline to run the simulations and to produce a reference table for the ABC-RF inference of demography and selection.

**CITATION**
```
@article {Pavinato2021.03.12.435133,
	author = {Pavinato, Vitor A. C. and De Mita, St{\'e}phane and Marin, Jean-Michel and de Navascu{\'e}s, Miguel},
	title = {Joint inference of adaptive and demographic history from temporal population genomic data},
	elocation-id = {2021.03.12.435133},
	year = {2021},
	doi = {10.1101/2021.03.12.435133},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/03/12/2021.03.12.435133},
	eprint = {https://www.biorxiv.org/content/early/2021/03/12/2021.03.12.435133.full.pdf},
	journal = {bioRxiv}
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

In the main script you should speficy the path of the above mentioned dependencies and other configurations that might be necessary (prior range, genome size, sample size, etc). 
