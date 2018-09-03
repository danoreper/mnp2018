# MNP
Code to generate results for 'Reciprocal F1 hybrids of two inbred mouse strains
reveal parent-of-origin and perinatal diet effects on behavior and expression' (Oreper, Schoenrock, 2018)


# Requirements
* OSX or linux
* jags
* R>=3.4.1
* Python>=2.7.12
* R CRAN packages: Cairo, car, coda, corrplot, data.table, evir, ggplot2, ggrepel, grid, gdata,  gtools, flextable, igraph, lattice, lme4, lmerTest, MASS,  mcmcplots,  mvtnorm, multcomp, nlme, officer, parallel, pracma, reshape, reshape2, rjags, stringr,tools, utils, yaml.
* R Bioconductor packages: affxparser, BiocGenerics, biomaRt, Biostrings, GenomicRanges, IRanges, limma, oligo, org.Mm.eg.db, Rsamtools, rtracklayer.
* Python packages: PyYAML


# Download
Git clone from the following url: https://github.com/danoreper/mnp2018.git

# Install
Enter the following at the command line, using working directory mnp2018\_LOCATION (the cloned, local mnp2018 repository):  
```mnp2018_LOCATION\$ bash install.sh``` 

For now, this install script only downloads data files (that are too large for github), and places them in their expected location relative to source code. jags, R, Python, and Python/R libraries need to be self-installed by the user.


# Run
1. To run locally:  
```mnp2018_LOCATION/src$ R CMD BATCH ./mnp/runAll.R```  
Note that running locally will not generate permutation testing threshholds, as it is too slow to be feasible.

2. Running locally will require up to 2 days, even without the permutation testing. If you have access to an LSF based cluster, computation will also finish in a few days, but for 400 permutations. To use an LSF based system, enter the following:  
```mnp2018_LOCATION/src$ bsub -M 20 -q week R CMD BATCH '--args ../config/defaultCluster.yaml' ./mnp/runAll.R```


# Directory structure

* mnp2018\_LOCATION/external: external tools that are packaged with this code; specifically the particular version of Affymetrix Power Tools that we used. 
* mnp2018\_LOCATION/src: All source code
* mnp2018\_LOCATION/src/genomerep: code to do general computations about the genome, such as where variants are
* mnp2018\_LOCATION/src/lm: code to do various sorts of linear modelling, including transformation, mixed models, multivariate models, etc. Also code to parse there various models.
* mnp2018\_LOCATION/src/parallel: code to parallelize jobs across killdevil/longleaf/multicore
* mnp2018\_LOCATION/src/util: code to load up system parameters via yaml, as well as some generally useful utilities
* mnp2018\_LOCATION/src/mnp: code that is very specific to the matnut pilot analysis
* mnp2018\_LOCATION/src/mnp/behavior: code for modeling and reporting mnp behavior
* mnp2018\_LOCATION/src/mnp/micro: code for modeling and reporting mnp microarray data
* mnp2018\_LOCATION/src/mnp/mediation: code for mnp mediation analysis
* mnp2018\_LOCATION/src/mnp/qpcr: code for modeling and reporting qPCR data
* mnp2018\_LOCATION/config: configuration files for running locally, on killdevil, or longleaf.
* mnp2018\_LOCATION/data: Data files for running analysis. Populated by install.sh which pulls from zenodo

* mnp2018\_LOCATION/output/: all output files from analyses
* mnp2018\_LOCATION/output/mnp/behavior: behavior analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/micro: microrray analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/qpcr: qpcr analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/mediation: mediation analysis, generated by runAll.R


