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


hets = list()
for(s1 in c("CAST_EiJ", "PWK_PhJ", "WSB_EiJ"))
{
    het = list()
    for(chr in db$genotype$getChrs())
    {
        het[[chr]] = db$read(strain1 = s1, strain2 = "C57BL6J", type = "genotype", phased = T, chr = chr)
    }
    het = do.call(rbind, het)
    het = het[allele1!="None"]
    het = het[allele1!=allele2]
    hets[[s1]] = het
    rm("het")
    gc()
}


hets = do.call(rbind, hets)
hets = hets[!duplicated(variant_id)]
fwrite(hets, "../data/crowley15/fgh_b6_variants.csv")
rm("hets")
gc()
