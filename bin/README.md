[//]: # (Hello,)

[//]: # (I'm done with a complete version of the summstats.py script. As we discussed, it takes an input file in the format of the "example.txt" file attached,)
[//]: # (and generates, based on option values, an output file in the format of the "test.txt" file attached. I attached the current version of the source file) 
[//]: # (also - the program is complete with this single file.)

[//]: # (But it will be more interesting to test it with a more serious input file. Where are you with SLiM?)

[//]: # (To use it, you need EggLib version 3.0.0b18 because I've especially added a feature to make it easier to compute the SFS. It will be online soon.)

Below is the program manual:

**1. Command line arguments**

The program is a single Python script. It is invoked as this:

```bash
python summstats.py <OPTIONS>
```

Options are all in the form <KEY>=<VALUE> with keys as listed below:

    - input-file        name of the input file

    - output-file       name of the output file

    - LSS               comma-separated list of summary statistics computed on selected loci, taken from the list:

                            * He          heterozygosity

                            * Dj          Jost's D

                            * WCst        Weir and Cockerham's theta

    - WSS               list of statistics computed on a window around selected loci, taken LSS plus that list:

                            * S           number of polymorphic sites

                            * thetaW      Watterson's 4Nu estimator

                            * D           Tajima's D

                            * Da          net distance between populations

                            * ZZ          Rozas et al.'s ZZ

    - GSS               list of global summary statistics, taken from the WSS list, plus SFS

    - wspan             span of windows on each side of the focal locus (in bp, such that the window size is actually 2 * wspan + 1 bp)

    - SFS-bins          number of bins for SFS

    - select            method to select focal loci (to compute LSS and WSS), from the list:

                            * all         use all loci

                            * rand        use each locus randomly

                            * freq        use each nth locus

                            * list        use the "selection" column of the input file

    - select-num        number of loci to draw randomly (if select=rand)

    - select-freq       period of selected loci (if select=freq): if

                        freq=1, all loci are selected; if freq=2, each

                        second locus is selected; and so on


**2. Input file**

The input file is a text file with space/tab separated values. There isa header line 
giving column names and, then, one line per site (a site may be a SNP, and indel or 
whatever). No term may contain a space. The order of columns is unimportant, except that 
the individual columns must be after the columns with fixed names.

List of columns:

    "chrom"       a string identifying a chromosome, contig or whatever.
    
    "position"    site position, as an integer.

    "status"      site category, chosen from "IG" (intergenic), "S"

                 (synonymous), and "NS" (non-synonymous).

    "alleles"     a comma-separated list of allelic values. Each allele

                  is represented by a string or undefined length. There

                  must be at least two alleles.

    "selection"   a column of Y or N to indicate whether LSS statistics

                  should be reported for each site [OPTIONAL].

    "name@label"  one additional column per individual, where "name" is

                  the sample name and "label" is a population label

                  (both are strings).

    ...           and so on for all individuals.


The data (values for individual columns) are two-digit strings with integer values 
(diploid data only).

    - 0 for missing data.

    - 1 for the first allele.

    - 2 for the second allele, and so on.

**3. Output file**

The output file is a tab-separated table, with a header line and, then, one line per 
locus (matching the loci found in the input file. The first column is titled "ID" and 
gives the identifier of a loci, in the form "chrom:position", as taken from the input 
file. All subsequent columns provide the values for a statistic. The column name is 
the statistic code (with a "LSS", "WSS", or "GSS" prefix). For GSS, the value of the 
statistics is identical for all loci. Missing data are denoted by "NA".
