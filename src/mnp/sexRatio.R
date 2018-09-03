##Sex ratio analysis.
library(data.table)
source("./loadParams.R")

## df = fread(dat("./mnp/male_female_breeding.csv"))
##  df = df[,c("Pipeline", "Breeding Batch", "Diet", "Strain", "Dam ID", "Litter", "Litter DOB", "Litter Size", "Pups Survived","# M pups wean", "# F pups wean")]
##  fwrite(df, dat("./mnp/damdata.csv"))

df = fread(dat("./mnp/damdata.csv"))
df[Strain=="B6J", Strain:="B6"]
df = df[Litter=="Yes"]
##df[,DamID:=paste0(Strain, DamID)]

setnames(df,
         c("Dam ID", "# M pups wean", "# F pups wean", "Breeding Batch", "Strain"),
         c("DamID", "pups.m", "pups.f", "BreedingBatch", "POE"))


afit = glm(df,
           formula = cbind(pups.m, pups.f)~BreedingBatch + Diet + POE + Diet:POE,
           family = binomial)
    
an = (anova(afit, test = "Chisq"))
rnames = rownames(an)
an = data.table(an)
an = cbind(variable=rnames, an)
fwrite(file = "../output/mnp/sexratio.txt", an)
print(an)
