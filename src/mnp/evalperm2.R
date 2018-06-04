library(data.table)
source("./mnp/loadAllData.R")
source("./multipleTesting.R")


##load("../output/mnp.bak.2/intermediate/permsmall_Diet:Strain")

## load("~/Desktop/mnp.nullt/permsmall_Diet:Strain")
## permOut.ds = permOut

## load("~/Desktop/mnp.nullt/permsmall_Diet")
## permOut.d  = permOut

## load("~/Desktop/mnp.nullt/permsmall_Strain")
## permOut.s  = permOut

## inp  = loadAllData$createAllInputs()

plotgene <- function(psid, varabletype, inp)
{
    library(ggplot2)
    frm = inp$cov.data.full
    frm$ID = as.character(frm$ID)
##    frm$expression = NA
    frm = frm[ID %in% rownames(inp$exp.mat)]
    
    expfrm= data.table(ID = rownames(inp$exp.mat), expression = inp$exp.mat[,psid])
    ##browser()
    frm = frm[expfrm, on = "ID"]
    alldata$Strain.x = (as.integer(alldata$Strain)-(length(unique(alldata$Strain))+1)/2)
    aplot = ggplot(frm, aes(x = Diet, color=Strain, y=expression, size = 1.5))
    aplot = aplot + geom_jitter(width=.1)
    print(aplot)
}


ident = 1
##vr = "Diet:Strain"
limit = 10001
##vr = "Diet"
vr = "Strain"

if(vr=="Strain")
{
    permOut = permOut.s
} 
if(vr =="Diet")
{
    permOut = permOut.d
}
if(vr == "Diet:Strain")
{
    permOut = permOut.ds
}

    
lrrc16a = "10408280"
Mir341  = "10398350"
Meg3    = "10398326"
Cnot2   = inp$probesetInfo[gene_name=="Cnot2"]$Probe.Set.ID
Airn    = inp$probesetInfo[gene_name=="Airn"]$Probe.Set.ID

    
dfsimple = (permOut$ident.full$results$dfsimple)
perm.lik   = permOut$perms$lik.rats
##perm.lik   = perm.lik[2:nrow(perm.lik),]

perm.anova.p = permOut$perms$perms
perm.anova.p = perm.anova.p[2:limit,]
##perm.anova.p = perm.anova.p[2:nrow(perm.anova.p),]

## perm.lik.p = 1 - pchisq(2*perm.lik, df=3)
## perm.lik.p = perm.lik.p[2:nrow(perm.lik.p)]

p.anova.emp = rep(NA, ncol(perm.anova.p))
names(p.anova.emp) = colnames(perm.anova.p)

p.lik.emp = rep(NA, ncol(perm.anova.p))
names(p.lik.emp) = colnames(perm.anova.p)

p.lik.p.emp = rep(NA, ncol(perm.anova.p))
names(p.lik.p.emp) = colnames(perm.anova.p)

## for(i in 1:nrow(perm.anova.p))
## {
##     print(i)
##     perm.anova.p[i,] = multipleTesting$median.correct.2(perm.anova.p[i,])
## }
## dfsimple$p.value = multipleTesting$median.correct.2(dfsimple$p.value)


for(i in 1:nrow(dfsimple))
{
    if(i%%5000==0)
    {
        print(i)
    }
    nm  = unlist(dfsimple[i,"Probe.Set.ID"])
    ident.anova.p = unlist(dfsimple[i, "p.value"])
    #ident.lik     = unlist(dfsimple[i, "lik.rat"])
    #ident.lik.p   = 1-pchisq(2*ident.lik, 3)
    
    p.anova.emp[nm] = (ident+sum(perm.anova.p[,nm]<=ident.anova.p, na.rm =T))/(ident+nrow(perm.anova.p))
     #p.lik.emp[nm]   = sum(perm.lik[,nm]>=ident.lik, na.rm = T)/nrow(perm.lik)
 #   p.lik.p.emp[i] = 1 - sum(perm.lik.p[,nm]>=ident.lik.p, na.rm=T)/nrow(perm.lik.p)
}


setorder(dfsimple, p.value)

fwer = apply(perm.anova.p, 1, min)
##names(fwer) = names(p.anova.emp)

fdr =  (p.adjust(p.anova.emp, method = "fdr"))

print("fdr")
print(dfsimple[Probe.Set.ID %in% names(fdr)[fdr<=.05]])
genez.fdr = inp$probesetInfo[Probe.Set.ID %in% names(fdr)[fdr<=.05]]$gene_name
print(genez.fdr)

try(dev.off())
plot.new()
par(mfrow = c(1,2))
hist(p.anova.emp, main=paste("perm.pvals (",vr, ")"))

print("fwer")
print(paste("thresh = ", quantile(fwer, .05, na.rm=T)))
genez.fwer = (dfsimple[p.value < quantile(fwer, .05, na.rm=T)])
print(genez.fwer)
genez.fwer = inp$probesetInfo[Probe.Set.ID %in% genez.fwer$Probe.Set.ID]$gene_name
print(genez.fwer)
hist(dfsimple$p.value, main = paste("anova.pvals(",vr, ")"))




## fwer.mins = unlist(apply(perm.anova.p, 1, which.min))
## bads = tail(sort(table(names(p.anova.emp)[fwer.mins]), F),10)
## print(bads)
## print(sum(bads))

## bdf = dfsimple[Probe.Set.ID %in% names(bads)]

## inp$probesetInfo[Probe.Set.ID %in% names(bads)]



