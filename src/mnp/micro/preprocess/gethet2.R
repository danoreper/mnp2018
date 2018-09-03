source("./genomerep/variantdb2/dump_parser.R")

dbdir = "../../data/isvdb/full_1410"
db    = db_builder$get.db.lite(dbdir)

nodhet = list()
for(chr in db$genotype$getChrs())
{
    nodhet[[chr]] = db$read(strain1 = "NOD_ShiLtJ", strain2 = "C57BL6J", type = "genotype", phased = T, chr = chr)
}

nodhet = do.call(rbind, nodhet)
nodhet = nodhet[allele1!="None"]
nodhet = nodhet[allele1!=allele2]
fwrite(nodhet, "../data/mnp/d_b6_variants.csv")
rm("nodhet")
gc()


