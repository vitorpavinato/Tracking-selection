# TODO

## Model
- [ ] *Defined default SLiM version for the simulations* **High** The models and pipelines were extensively tested with SLiM v.3.1 to v3.2.1, so they should be used at least for this version of the program;
- [ ] *Refactor `model proof` to upload the .tree as a `late()` instead of `early()` event in period 2* **High** Keep it as it is now for the additional simulation with different recombination rates for the `model proof`; but change it if new RF training simulations with recombination would be required (**NOTE** the model for the analyses of bee data already has this changes, see line 50); This should avoid the warning message raised by SLiM when uploading the .tree file;
- [ ] *Add one more generation after the upload of the .tree file* in `model proof` **High** This is required to completely refactor the uploading of .tree file as a `late()` event (**SEE** the script for the analyses of bee data);
- [ ] *Add `initializeSLiMOptions(keepPedigrees = T)` at the header of `model proof` burning file* **Low** This is required in order to keep track of parents if run the model with SLiM version 3.3 or higher (**NOTE** This was already implemented for the analyses of bee data);
- [ ] *Change calculation of genetic load (see calculation used for bees data)* **Medium** It keeps cached fitness different from zero only, and implements one calculation for the rest of simulation instead of different instantiations sprinkled on the code;
- [ ] *Check if using a different equation for the Ne (latent variable) would change the results* **Low** The equation used seems to lack `-2` term in the numerator; the correct one should be $Ne =\frac{4N - 2}{var(k) + 2}$ (raised on 07/Nov/2019);
- [ ] *Fix the conditional for the Ne calculation and check how it would impact the simulations* **Low** Originally it checks for the size of the vector of `gametes` higher than 1, but it should check for values different from zero (raised on 07/Nov/2019);

## Pipeline
- [ ] *Remove from the pipelines (proof and bee) the command line calls for the recombination rate ("rr") and for the genome size ("genomeS")* **Low** It does not affect the call for SLiM within R, but it can make the code confusing;
- [ ] *Incorporate the bug fixes implemented in the pipeline for the analyses of bee data on the `model proof` pipeline* **Low** To make both pipelines consistent in how they handle raised errors.
- [ ] *Fix how the pipeline check if .vcf files were merged* **Low** This error was raised when running the pipeline for the analyses of bee data; apparently, even if the bcftools did not detect a merged file, it creates one with zero lines (raised on 07/Nov/2019);
- [ ] *Fix error handling by RADseq sampling* **Low** Apparently, the error raised above was a consequence of mishandling of the RADseq sampling check; In order to fix it, add a different exception/error message for the RADseq sampling (raised on 07/Nov/2019);
- [ ] *Fix how the number of fragments under selection is rounded* **Low** It was already fixed on the pipeline for the analyses of bee data; should be incorporated in the `model proof` pipeline.