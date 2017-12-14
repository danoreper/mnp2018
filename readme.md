

Files within /data

b6_reference: data pertaining to the b6 reference genome

b6_reference/mm10.chrom.sizes:   the mm10 sizes of every chromosome, derived from http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&chromInfoPage=,   on  April 14, 2015

b6_reference/mus_musculus.GRCm38.75: folder for 38.75 reference in particular

b6_reference/mus_musculus.GRCm38.75/genome.fa: The fasta sequence of the B6 genome including all chromosome and MT

b6_reference/mus_musculus.GRCm38.75/Mus_musculus.GRCm38.75.gtf: The gtf which crucially contains exon locations

b6_reference/mus_musculus.GRCm38.75/snord: the Snord ensembl ids in csv files, with a separate csv per Snord family. Locations left off.


imprinted_genes_list/crowley_2015_brain_imprinted_genes.csv: The genes identified as imprinted in the crowley paper. The locations aren't used, but the names of the genes are.

imprinted_genes_list/Literature_Imprinted_Genes.csv: the mousebook imprinted genes, positions are unused.

imprinted_genes_list/feb2014.archive.ensembl.org_allImprintedGenes.csv: The collated set of all imprinted genes, (mousebook and crowley) complete with their locations as pulled from the ensembl 38.75 database (Feb2014). This is the KEY file describing which genes are imprinted.

mnp/2017-02-13_breederLog.csv: The true (non-behavior) covariates for all animals bred for this experiment (e.g., diet, strain, etc).The covariates from this breeding file are considered correct, and take precedence over any identical (in name or meaning) covariates in any of the phenotype files.

mnp/allcovariates_2011_12_08.csv: The covariates file used to decide which animals should be sequenced; there are some mistakes with respect to the breederLog, but we maintain this file for reproducibility.

mnp/behaviorModels.csv: The models employed per behavior for the most straightforward behaviors; a handful of other behaviors, in particular the PPI related ones, are entirely specified programmatically rather than in a document.

mnp/cel_files: The cel files from the Affy 1.1ST exon array. Each is named according to the animal sample id that was measured.

expression_choice/expression_choice_20120411.csv: The output of the algorithm to choose which animals would be sequenced

mnp/log2RMA-GENE-DEFAULT-Group1.txt:
APT applied with defaults and no masking; used primarily for its probeset annotations

mnp/MatNut_validation_samplelist_072516.csv: list of all samples and their status (pulverized or not, etc) as of 072516.
Primarily used for qPCR planning

mnp/microarray_lib_files: The affymetrix library files for the MoGene1.0 and MoGene1.1 exon ST arrays. Used to
generate expression matrix across all samples with proper normalization using affymetrix power tools,
as well as to identify probes to mask

mnp/MoGene-1_0-st-v1_probeinfo_b.txt: The MoGene1.0 exon st probe binding locations, as pulled from ensembl38.75 (see getEnsemblProbeAligns.sh)

mnp/nod.b6.variants.txt: the snp/indel variant positions (and the variant sequences) that vary between NOD and B6. Extracted from
ISVdb. Used for masking variants. 

mnp/phenotypes: the behavioral phenotype data. Not that all the covariates specified in the breederLog file are specified here too, but they may ben incorrect; thus, many of the covariates in these files are unused, and are pulled from the breederLog file instead.
mnp/BodyWeight_8.12.14.csv: the bodyweight data
mnp/phenotypes/Cocaine/cokedata_final_1-14-16.csv: the cocaine response data
mnp/phenotypes/ForcedSwim/newFSTdata_12-8-15.csv: the forced swim data
mnp/phenotypes/LightDark/dietstudy_light_dark.csv: the light dark test data
mnp/phenotypes/OpenField/dietstudy_openfield.csv: the open field test data
mnp/phenotypes/Restraint Stress/cortdata_1-14-16.csv: the restraint stress test data
mnp/phenotypes/SIH/dietstudy_allsih_1-14-16.csv: the SIH test data
mnp/phenotypes/sociability/SocialHabSoc.csv: the sociability data.
mnp/phenotypes/Startle_PPI/PPI_Rawdatafiles_combined.csv: the PPI data.
mnp/phenotypes/TailSuspension: dietstudy_tst.csv: the tail suspension test data.














