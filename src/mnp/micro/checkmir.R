library(data.table)
curwd = getwd()
mirhubPath = "../../../resources/miRhub"
setwd(mirhubPath)


genelist = fp(curwd, "../output.resub.bak/mnp/micro/effect.table/p_0.05_Diet:Strain_fwer.csv")
df = fread(genelist)
genelist.file = "./dsgenes.txt"
    ##file.path(dirname(genelist),"ds.genes.list") 
fwrite(data.frame(df$gene_name), col.names=F, file = genelist.file) 

command = paste0("python3 miRhub.py ",  "-c0 ", genelist.file,"")
system(command)
setwd(curwd)
