library(data.table)
    
load("../output/mnp/intermediate/permsmall_Strain")

##nms = c("Strain", "Diet", "Strain:Diet")
##nms = "Diet"
nms = "Strain"
lrrc16a = "10408280"
Mir341  = "10398350"
Meg3    = "10398326"
Cnot2   = inp$probesetInfo[gene_name=="Cnot2"]$Probe.Set.ID
Airn    = inp$probesetInfo[gene_name=="Airn"]$Probe.Set.ID

    
dfsimple = (permOut$ident.full$results$dfsimple)
perm.lik   = permOut$perms$lik.rats

perm.anova.p = permOut$perms$perms
#perm.lik.p = 1 - pchisq(2*perm.lik, df=3)


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
    ident.lik     = unlist(dfsimple[i, "lik.rat"])
    ident.lik.p   = 1-pchisq(2*ident.lik, 3)

    p.anova.emp[i] = 1 - sum(perm.anova.p[,nm]>=ident.anova.p, na.rm =T)/nrow(perm.anova.p)
    p.lik.emp[i]   = sum(perm.lik[,nm]>=ident.lik, na.rm = T)/nrow(perm.lik)
 #   p.lik.p.emp[i] = 1 - sum(perm.lik.p[,nm]>=ident.lik.p, na.rm=T)/nrow(perm.lik.p)
}

