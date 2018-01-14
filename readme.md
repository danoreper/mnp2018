# MNP
Code to generate results for 'F1 reciprocal crosses of common inbred lines ..., Oreper, 2017 '

# Requirements
* OSX or linux
* jags
* R>=3.3.1
* Python>=2.7.12
* R CRAN packages: Cairo, car, coda, corrplot, data.table, evir, ggplot2, ggrepel, grid, gdata, igraph, reshape, reshape2, lattice, lme4, lmerTest, MASS,  mcmcplots,  mvtnorm, multcomp, nlme,  parallel, pracma, rjags, stringr,tools, utils, yaml.
* R Bioconductor packages: affxparser, BiocGenerics, biomaRt, Biostrings, GenomicRanges, IRanges, limma, oligo, org.Mm.eg.db, Rsamtools, rtracklayer.
* Python packages: PyYAML


# Download
Git clone from the following url: https://github.com/danoreper/mnp2018.git

# Install
Enter the following at the command line, using working directory mnp2018\_LOCATION (the cloned, local mnp2018 repository):
mnp2018\_LOCATION\$ bash install.sh 

For now, this install script only downloads data files (that are too large for github), and places them in their expected location relative to source code. jags, R, Python, and Python/R libraries need to be self-installed by the user.


# Run
1. To run locally: mnp2018\_LOCATION/src$ R CMD BATCH ./mnp/runAll.R
Note that running locally will not generate permutation testing threshholds, as it is too slow to be feasible.

2. Running locally will require up to 2 days, even without the permutation testing. If you have access to an LSF based cluster, computation will also finish in a few days, but for 400 permutations. To use an LSF based system, enter the following:

mnp2018\_LOCATION/src$ bsub -M 20 -q week R CMD BATCH '--args ../config/defaultCluster.yaml' ./mnp/runAll.R


# Directory structure
* mnp2018\_LOCATION/src: The source code that generates all analysis results
* mnp2018\_LOCATION/config: configuration files for running locally or on killdevil.
* mnp2018\_LOCATION/data: Data files for running analysis. Populated by install.sh

* mnp2018\_LOCATION/output/mnp/behavior: behavior analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/micro: microrray analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/qpcr: qpcr analyses, generated by runAll.R
* mnp2018\_LOCATION/output/mnp/mediation: mediation analysis, generated by runAll.R


