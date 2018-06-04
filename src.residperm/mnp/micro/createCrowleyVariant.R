source("./loadParams.R")
source("./genomerep/variantdb2/dump_parser.R")

db = db_builder$get.db.lite(dat("isvdb/exon_1410"))

parseFunc <- function(strain, chr, df)
{
    df = df[as.character(consequence_1)!="reference"|
            as.character(consequence_2)!="reference"]
    
    df$chr = chr
    return(df)
}
    

out = db$iterate("genotype",
                 strain1s = c("CAST_EiJ", "PWK_PhJ", "WSB_EiJ"),
                 strain2s =rep("C57BL6J",3),
                 parseFunc = parseFunc)


## out= db$genotype$iterate(strains = c("CAST_EiJ", "PWK_PhJ", "WSB_EiJ"),
##                          parseFunc = parseFunc)


