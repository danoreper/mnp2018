mysql -u root -p -e '

select 
v1.variant_id as variant_id, 
v1.chrom as chrom, 
v1.pos as pos,
strain.strain_name as strain_name, 
allele_index_1, 
allele_index_2, 
a1.allele as allele1, 
a2.allele as allele2, 
ref.allele as refallele, 
prob, 
is_max, 
c1.consequence_name as consequence_1, 
c2.consequence_name as consequence_2, 
g.gene_name as gene_name, 
t.transcript_name as transcript_name 

from variant v1

inner join genotype ila on v1.variant_id=ila.variant_id 
inner join allele a1 on a1.variant_id = ila.variant_id and a1.allele_index=ila.allele_index_1 
inner join allele a2 on a2.variant_id = ila.variant_id and a2.allele_index=ila.allele_index_2 
inner join strain on ila.strain_id=strain.strain_id 
inner join allele ref on ref.variant_id=v1.variant_id and ref.allele_index=0 

inner join consequence_info vinfo1 on vinfo1.variant_id = v1.variant_id and vinfo1.allele_index=ila.allele_index_1
inner join consequence c1 on c1.consequence_id = vinfo1.consequence_id
inner join gene g on vinfo1.gene_id = g.gene_id
inner join transcript t on vinfo1.transcript_id = t.transcript_id

inner join consequence_info vinfo2 on vinfo2.variant_id = v1.variant_id and vinfo2.allele_index=ila.allele_index_2 and vinfo1.transcript_id = vinfo2.transcript_id
inner join consequence c2 on c2.consequence_id = vinfo2.consequence_id

limit 10000000000;' paul_1410 > out.txt

