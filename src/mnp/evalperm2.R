library(data.table)
source("./mnp/loadAllData.R")

load("../output/mnp.bak.2/intermediate/permsmall_Diet:Strain")
## inp  = loadAllData$createAllInputs()


nms = "Strain"
lrrc16a = "10408280"
Mir341  = "10398350"
Meg3    = "10398326"
Cnot2   = inp$probesetInfo[gene_name=="Cnot2"]$Probe.Set.ID
Airn    = inp$probesetInfo[gene_name=="Airn"]$Probe.Set.ID

    
dfsimple = (permOut$ident.full$results$dfsimple)
perm.lik   = permOut$perms$lik.rats
##perm.lik   = perm.lik[2:nrow(perm.lik),]

perm.anova.p = permOut$perms$perms
perm.anova.p = perm.anova.p[2:10001,]
##perm.anova.p = perm.anova.p[2:nrow(perm.anova.p),]

## perm.lik.p = 1 - pchisq(2*perm.lik, df=3)
## perm.lik.p = perm.lik.p[2:nrow(perm.lik.p)]

p.anova.emp = rep(NA, ncol(perm.anova.p))
names(p.anova.emp) = colnames(perm.anova.p)

p.lik.emp = rep(NA, ncol(perm.anova.p))
names(p.lik.emp) = colnames(perm.anova.p)

p.lik.p.emp = rep(NA, ncol(perm.anova.p))
names(p.lik.p.emp) = colnames(perm.anova.p)

for(i in 1:nrow(dfsimple))
{
    nm  = unlist(dfsimple[i,"Probe.Set.ID"])
    ident.anova.p = unlist(dfsimple[i, "p.value"])
    #ident.lik     = unlist(dfsimple[i, "lik.rat"])
    #ident.lik.p   = 1-pchisq(2*ident.lik, 3)

    p.anova.emp[nm] = (1+sum(perm.anova.p[,nm]<=ident.anova.p, na.rm =T))/nrow(perm.anova.p)
     #p.lik.emp[nm]   = sum(perm.lik[,nm]>=ident.lik, na.rm = T)/nrow(perm.lik)
 #   p.lik.p.emp[i] = 1 - sum(perm.lik.p[,nm]>=ident.lik.p, na.rm=T)/nrow(perm.lik.p)
}


fwer = apply(perm.anova.p[2:nrow(perm.anova.p),-4079], 1, min, na.rm=T)
##names(fwer) = names(p.anova.emp)

fdr =  (p.adjust(p.anova.emp, method = "fdr"))

print("fdr")
print(dfsimple[Probe.Set.ID %in% names(fdr)[fdr<=.05]])
genez = inp$probesetInfo[Probe.Set.ID %in% names(fdr)[fdr<=.05]]$gene_name
print(genez)

print("fwer")
genez = (dfsimple[p.value < quantile(fwer, .05)])
print(genez)
genez = inp$probesetInfo[Probe.Set.ID %in% genez$Probe.Set.ID]$gene_name
print(genez)

