create.F1.table <- function(con, strain1, strain2)
{
    subq = function(i1, i2,name)
    {
        subquery =
            paste0("select * from (select ", 
                   "v1.variant_id as variant_id, ",
                   "v1.chrom as chrom, ",
                   "v1.pos as pos, ",
                   "s1.strain_name as strain_name_1, ", 
##                 "fa1.allele_index_",i1," as allele_index_1, ", 
##                 "fa2.allele_index_",i2," as allele_index_2, ", 
                   "a1.allele as allele1, ",
                   "a2.allele as allele2, ", 
                   "fa1.prob*fa2.prob/4 as prob, ", 
                   "c1.consequence_name as consequence_1, ",
                   "c2.consequence_name as consequence_2, ",
                   "g.gene_name as gene_name, ",
                   "t.transcript_name as transcript_name ",

                   "from variant v1 ",

                   "inner join founder_allele fa1 on v1.variant_id=fa1.variant_id ",
                   "inner join strain s1 on   fa1.strain_id=s1.strain_id ",
                   "inner join allele a1 on   a1.variant_id = fa1.variant_id and a1.allele_index=fa1.allele_index_",i1," ",

                   "inner join founder_allele fa2 on v1.variant_id=fa2.variant_id ",
                   "inner join strain s2 on   fa2.strain_id=s2.strain_id ",
                   "inner join allele a2 on   a2.variant_id = fa2.variant_id and a2.allele_index=fa2.allele_index_",i2," ", 

                   "inner join varinfo vinfo1 on vinfo1.variant_id = v1.variant_id and vinfo1.allele_index=fa1.allele_index_",i1, " ",
                   "inner join consequence c1 on c1.consequence_id = vinfo1.consequence_id ",
                   "inner join gene g on vinfo1.gene_id = g.gene_id ",
                   "inner join transcript t on vinfo1.transcript_id = t.transcript_id ",

                   "inner join varinfo vinfo2 on vinfo2.variant_id = v1.variant_id and vinfo2.allele_index=fa2.allele_index_",i2, " and vinfo1.transcript_id = vinfo2.transcript_id ",
                   "inner join consequence c2 on c2.consequence_id = vinfo2.consequence_id ",

                   "where s1.strain_name = '", strain1, "' and s2.strain_name = '", strain2, "' and a1.allele_index!=a2.allele_index) as ", name)
        return(subquery)
    }
    
     q = paste0(subq(1,1, "t1"),# " UNION ALL ",
              # subq(1,2, "t2"), " UNION ALL ",
              # subq(2,1, "t3"), " UNION ALL ",
               # subq(2,2, "t4"),
                ";")
    print(q)
    getout = dbGetQuery(con, q)
    return(getout)
}


strain1 = "NOD_ShiLtJ"
strain2 = "C57BL6J"
user     = "root"
password = "FrozenMixedVegetables"
db   = "full_1410"
host     = "127.0.0.1"
##an unused port. Someday, write code to check which ports are available (i.e. rotation student)
port     = 3308

##either get system call to work without hanging, or do this manually on the CL to open an SSH tunnel
##from port 3308 on this machine to port 3306 on the valdardb machine. Then when rmysql connects to port 3308 locally, it sees valdardb
##system("ssh -nNT -L 3308:localhost:3306 doreper@valdardb.its.unc.edu")

con = dbConnect(MySQL(), user=user, port = port, password =password, dbname=db, host=host)
dbSendQuery(con, 'set @@max_heap_table_size=4294967296;')
out = create.F1.table(con, strain1, strain2)
write.table(out, file.path(prop$mnp$data, prop$mnp$NOD.B6.variantsFile), sep="\t", row.names=F)
