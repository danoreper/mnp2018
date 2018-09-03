datadir=./data2/

##rm -rf $datadir

mkdir $datadir
mkdir $datadir/b6_reference/
mkdir $datadir/b6_reference/mus_musculus.GRCm38.75/
mkdir $datadir/imprinted_genes_list/
mkdir $datadir/mpd
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


cd $datadir

##wget -m -nH --cut-dirs=1 -np -A 'File_*.zip, File_*.csv, File_*.txt' -r https://zenodo.org/record/1343994/

cd ./1343994/files
datadir=../../

cp File_S2_mm10.chrom.sizes      $datadir/b6_reference/mm10.chrom.sizes

unzip -o File_S3_Mus_musculus.GRCm38.75.gtf.zip; mv File_S3_Mus_musculus.GRCm38.75.gtf $datadir/b6_reference/mus_musculus.GRCm38.75/Mus_musculus.GRCm38.75.gtf

unzip -o File_S4_Snord.zip; mv File_S4_Snord $datadir/b6_reference/mus_musculus.GRCm38.75/snord

unzip -o File_S5_imprinted_genes.zip;

mv File_S5_imprinted_genes/File_S5a_crowley_2015_brain_imprinted_genes.csv $datadir/imprinted_genes_list/crowley_2015_brain_imprinted_genes.csv

mv File_S5_imprinted_genes/File_S5b_Literature_Imprinted_Genes.csv $datadir/imprinted_genes_list/Literature_Imprinted_Genes.csv

mv File_S5_imprinted_genes/File_S5c_feb2014.archive.ensembl.org_allImprintedGenes.csv $datadir/imprinted_genes_list/feb2014.archive.ensembl.org_allImprintedGenes.csv

unzip -o File_S6_nod.b6.variants.csv.zip;

mv File_S8_nod.b6.variants.csv $datadir/mnp/nod.b6.variants.csv

unzip -o File_S7_breeding.zip
mv offspring_breedingLog.csv $datadir/mnp/2017-02-13_breederLog.csv
mv damdata.csv $datadir/mnp/damdata.csv

unzip -o File_S8_MoGene_library_files

mv MoGene-1_1_libs $datadir/mnp/microarray_lib_files/MoGene-1_1-st-v1.r4.analysis-lib-files 

mv MoGene-1_0_libs $datadir/mnp/microarray_lib_files/MoGene-1_0-st-v1.r3.analysis-lib-files

unzip -o File_S9_MoGene-1_0_st-v1_probeinfo_b.txt.zip; mv MoGene-1_0_st-v1_probeinfo_b.txt $datadir/mnp/MoGene-1_0-st-v1_probeinfo_b.txt 

unzip -o File_S10_expression.zip;

mv File_S10_expression/cel_files $datadir/mnp/cel_files
mv File_S10_expression/log2RMA-GENE-DEFAULT-Group1.txt $datadir/mnp/log2RMA-GENE-DEFAULT-Group1.txt 

unzip -o File_S11_reciprocal_behavior_weight.zip

mv File_S11_reciprocal_behavior_weight/weight.txt $datadir/mnp/phenotypes/BodyWeight_8.12.14.csv

mv File_S11_reciprocal_behavior_weight/cocaine.txt $datadir/mnp/phenotypes/Cocaine/cokedata_final_1-14-16.csv

mv File_S11_reciprocal_behavior_weight/swim.txt $datadir/mnp/phenotypes/ForcedSwim/newFSTdata_12-8-15.csv

mv File_S11_reciprocal_behavior_weight/lightdark.txt $datadir/mnp/phenotypes/LightDark/dietstudy_light_dark.csv

mv File_S11_reciprocal_behavior_weight/openfield.txt $datadir/mnp/phenotypes/OpenField/dietstudy_openfield.csv

mv File_S11_reciprocal_behavior_weight/cort.txt $datadir/mnp/phenotypes/RestraintStress/cortdata_1-25-16.csv

mv File_S11_reciprocal_behavior_weight/SIH.txt $datadir/mnp/phenotypes/SIH/dietstudy_allsih_1-14-16.csv

mv File_S11_reciprocal_behavior_weight/sociability.txt $datadir/mnp/phenotypes/Sociability/SocialHabSoc.csv

mv File_S11_reciprocal_behavior_weight/PPI_Rawdatafiles_combined.csv $datadir/mnp/phenotypes/Startle_PPI

mv File_S11_reciprocal_behavior_weight/tail.txt $datadir/mnp/phenotypes/TailSuspension/dietstudy_tst.csv

cp File_S12_behaviorModels.csv $datadir/mnp/behaviorModels.csv

cp File_S13_MatNut_validation_samplelist_072516.csv $datadir/mnp/qpcr/MatNut_validation_samplelist_072516.csv:

unzip -o File_S14_qpcrData.zip
mv File_S14_qpcrData/qpcrData.csv $datadir/mnp/qpcr/Matnut_pilot_alltaqmanplates_datacompletev3_100716.csv


unzip -o File_S15_mouse_phenome_database.zip
mv File_S15_mouse_phenome_database/Chesler4* $datadir/mpd/
mv $datadir/mpd/Chesler4_original.csv $datadir/mpd/Chesler4.csv

mv File_S15_mouse_phenome_database/Wiltshire1* $datadir/mpd/
mv $datadir/mpd/Wiltshire1_original.csv $datadir/mpd/Wiltshire1.csv

