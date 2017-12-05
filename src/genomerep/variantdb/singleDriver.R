source("./genomerep/buildGenomeData.R")
#msk = "/data/doreper/nd/output/yunjung2.bed"
#dbname = "yunjung_1410"

msk = "../data/vcf/rel1410/paul.bed"
dbname = "paul_1410"
buildGenomeData$buildSingleVariantDb(dbname, msk, rebuildVCF=F)

## msk = "/data/doreper/nd/vcf/rel1410/rachel_x.bed"
## buildGenomeData$buildSingleVariantDb("rachel_x", msk)
