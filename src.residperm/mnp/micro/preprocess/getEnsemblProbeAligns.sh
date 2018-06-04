#
# loads all the probe mapping information for the affy 1.0 st r3 probes from the ensembl database. The probe ids can be remapped to affy 1.1 st r4 probe ids. 
# 
mysql --host='ensembldb.ensembl.org' --port=3306 -u anonymous -e "
#show DATABASES LIKE '%mus_musculus_func%'; 
USE mus_musculus_funcgen_75_38; 
#show TABLES; 
#select * FROM array_chip;

SELECT probe.name as "probename", seq_region.name as "seqname", 
       probe_feature.seq_region_start, probe_feature.seq_region_end, probe_feature.seq_region_strand, 
       probe_feature.analysis_id, probe_feature.mismatches, probe_feature.cigar_line

FROM probe, array_chip, probe_feature, seq_region 
WHERE array_chip.design_id='MoGene-1_0-st-v1' 
AND   probe.array_chip_id = array_chip.array_chip_id
AND   probe_feature.probe_id = probe.probe_id
AND   seq_region.seq_region_id = probe_feature.seq_region_id" > "./MoGene-1_0-st-v1_probeinfo_b.txt"