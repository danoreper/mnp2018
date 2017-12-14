datadir=./data/

rm -rf $datadir

mkdir $datadir
mkdir $datadir/b6_reference/
mkdir $datadir/b6_reference/mus_musculus.GRCm38.75/
mkdir $datadir/imprinted_genes_list/
mkdir $datadir/mnp
mkdir $datadir/mnp/microarray_lib_files
mkdir $datadir/mnp/qpcr
mkdir $datadir/mnp/phenotypes
mkdir $datadir/mnp/phenotypes/Cocaine
mkdir $datadir/mnp/phenotypes/ForcedSwim
mkdir $datadir/mnp/phenotypes/LightDark
mkdir $datadir/mnp/phenotypes/OpenField
mkdir $datadir/mnp/phenotypes/RestraintStress
mkdir $datadir/mnp/phenotypes/SIH
mkdir $datadir/mnp/phenotypes/Sociability
mkdir $datadir/mnp/phenotypes/Startle_PPI
mkdir $datadir/mnp/phenotypes/TailSuspension

#replace with wget or curl to zenodo
cp -r ~/Desktop/zenodo/ $datadir/zenodo
cd $datadir/zenodo
datadir=../

mv File_S2_mm10.chrom.sizes.txt       $datadir/b6_reference/mm10.chrom.sizes

unzip File_S3_Mus_musculus.GRCm38.75.gtf.zip; mv File_S3_Mus_musculus.GRCm38.75.gtf $datadir/b6_reference/mus_musculus.GRCm38.75/Mus_musculus.GRCm38.75.gtf

unzip File_S4_Snord.zip; cp -r File_S4_Snord $datadir/b6_reference/mus_musculus.GRCm38.75/snord

mv File_S5_crowley_2015_brain_imprinted_genes.csv $datadir/imprinted_genes_list/crowley_2015_brain_imprinted_genes.csv

mv File_S6_Literature_Imprinted_Genes.csv $datadir/imprinted_genes_list/Literature_Imprinted_Genes.csv

mv File_S7_feb2014.archive.ensembl.org_allImprintedGenes.csv $datadir/imprinted_genes_list/feb2014.archive.ensembl.org_allImprintedGenes.csv

unzip File_S8_nod.b6.variants.csv.zip; cp -r File_S8_nod.b6.variants.csv $datadir/mnp/nod.b6.variants.csv

mv File_S9_2017-02-13_breederLog.csv $datadir/mnp/2017-02-13_breederLog.csv

unzip File_S10_MoGene-1_1_libs.zip; mv File_S10_MoGene-1_1_libs $datadir/mnp/microarray_lib_files/MoGene-1_1-st-v1.r4.analysis-lib-files 

unzip File_S11_MoGene-1_0_libs.zip; mv File_S11_MoGene-1_0_libs $datadir/mnp/microarray_lib_files/MoGene-1_0-st-v1.r3.analysis-lib-files

unzip File_S12_MoGene-1_0_st-v1_probeinfo_b.txt.zip; mv File_S12_MoGene-1_0_st-v1_probeinfo_b.txt $datadir/mnp/MoGene-1_0-st-v1_probeinfo_b.txt 

unzip File_S13_cel_files.zip; mv File_S13_cel_files $datadir/mnp/cel_files

mv File_S14_log2RMA-GENE-DEFAULT-Group1.txt $datadir/mnp/log2RMA-GENE-DEFAULT-Group1.txt 

mv File_S15_MatNut_validation_samplelist_072516.csv $datadir/mnp/qpcr/MatNut_validation_samplelist_072516.csv:

mv File_S16_qpcrData.csv $datadir/mnp/qpcr/Matnut_pilot_alltaqmanplates_datacompletev3_100716.csv

mv File_S17_behaviorModels.csv $datadir/mnp/behaviorModels.csv

mv File_S18_bodyweight.csv $datadir/mnp/phenotypes/BodyWeight_8.12.14.csv

mv File_S19_cocaine.csv $datadir/mnp/phenotypes/Cocaine/cokedata_final_1-14-16.csv

mv File_S20_forced_swim.csv $datadir/mnp/phenotypes/ForcedSwim/newFSTdata_12-8-15.csv

mv File_S21_light_dark.csv $datadir/mnp/phenotypes/LightDark/dietstudy_light_dark.csv

mv File_S22_openfield.csv $datadir/mnp/phenotypes/OpenField/dietstudy_openfield.csv

mv File_S23_cort.csv $datadir/mnp/phenotypes/RestraintStress/cortdata_1-14-16.csv

mv File_S24_SIH.csv $datadir/mnp/phenotypes/SIH/dietstudy_allsih_1-14-16.csv

mv File_S25_sociability.csv $datadir/mnp/phenotypes/Sociability/SocialHabSoc.csv

mv File_S26_startle_PPI.csv $datadir/mnp/phenotypes/Startle_PPI/PPI_Rawdatafiles_combined.csv

mv File_S27_tail.csv $datadir/mnp/phenotypes/TailSuspension/dietstudy_tst.csv
