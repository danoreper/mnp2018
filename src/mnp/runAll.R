print("starting")
## 
# Author: doreper
###############################################################################


try(rm("prop"))
sources <- function()
{
    if(!exists("prop"))
    {
        source("./loadParams.R")
    }
    
    source("./mnp/loadAllData.R")
##    source("./mnp/loadAllDataCrowley.R")
    source("./mnp/behavior/analysis.R")
    source("./mnp/micro/analysis.R")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")
}

sources()
dir.create(prop$mnp$output, showWarnings = F, recursive = T)
inp  = loadAllData$createAllInputs()
##browser()
## ps = inp$probesetInfo
## for(thresh in c(100, 500, 1000, 5000, 10000))
## {
##     ps1 = ps[minDistToImprinted<=thresh]
##     setorder(ps1, minDistToImprinted)
##     fwrite(ps1, file=outm("impThresh", paste0("bp_",thresh,".txt")), sep = "\t")
##     print(nrow(ps1))
## }
## inp$annot.data           = inp$annot.data[1:2000,]
## inp$annot.data.control   = inp$annot.data.control[1:2000,]
## inp$exp.mat              = inp$exp.mat[,1:2000]
## inp$exp.mat.control.full = inp$exp.mat.control.full[,1:2000]

alphas = c(.05)

fromFile = F
if(!fromFile)
{
    if(prop$mnp$SSVA.numperm==0)
    {
        print("running no perm")
        out = micro.analysis$run.noperm(inp)
        print("done running no perm")
    }  else {
        if(prop$mnp$mode=="ALL")
        {
            out.noperm = micro.analysis$run.noperm(inp)
            print("done running no perm")
            util$save(list =ls(), file = outm("noperm.RData"))
            
            print("running all perms")
            out = micro.analysis$runallPerms(inp)
            print("done all perms")
            util$save(file = outm("perm.RData"), list = ls())
        } else if (prop$mnp$mode == "MAIN") {
            print("running all perms MAIN")
            out = micro.analysis$runallPerms(inp)
            print("done all perms")
            util$save(file = outm("perm.MAIN.RData"), list = ls())
        } else if(prop$mnp$mode  == "INTERACTION") {
            print("running all perms INTERACTION")
            ##TODO consider running freedman?
            ##            out = micro.analysis$runallPermsFreed(inp)
            out = micro.analysis$runallPerms(inp)
            print("done all perms")
            util$save(file = outm("perm.INTERACTION.RData"), list = ls())
        }
    }
} else {
    if(prop$mnp$mode=="MAIN")
    {
        load(outm("perm.RData"))
        sources()
    }
}
    

##location of saved file
if(prop$mnp$mode == "COLLATE")
{
    load("../outputfinal/mnp/perm.RData")
    out1 = out
    load("../outputsynch/mnp/perm.INTERACTION.RData")
    out2 = out
    rm(out)
    sources()
    browser()
    ##TODO merge
    print(54)
}
if(prop$mnp$mode %in% c("ALL","MAIN", "COLLATE"))
{
    print("made it to collation")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")
    sources()
    if(prop$mnp$freed)
    {
        out$ident.full = try(out$results)
    }
    micro.report$reportAnalysis(inp$exp.mat,
                                inp$cov.data.full,
                                out$ident.full$results,
                                ##                            out$results$results,
                                inp$annot.data,
                                inp$probesetInfo,
                                out$threshholds,
                                inp$karyo)
}

qpcr.analysis$run(inp)




check.tf <- function()
{
    imprintedGenes = fread(dat(prop$genome$imprintedGenesFile))
    imp.genelist = imprintedGenes$ensembl_gene_id

    tf = fread(dat(prop$mnp$tf.file))
    tf = data.frame(tf)
    tf = data.table(tf)
    tf$Official.Full.Name
    tf$Transcription.factor

    ginfo = buildGenomeData$getGeneInfo()
    gID = ginfo$gene_id
    names(gID) = ginfo$gene_name
    tf.genelist = gID[tf$Transcription.factor]

    hits = intersect(imp.genelist, tf.genelist)
    ##chr13	24171246	24171366

    full    = fread("../output/sva/report/fullmerged.csv")
    lrrc16a = as.character(full[gene_name=="Lrrc16a"]$Probe.Set.ID[1])
    wt1     = as.character(full[gene_id=="ENSMUSG00000016458"]$Probe.Set.ID[1])

    x = out$ident.full$sv.info$exp.mat[,wt1]
    y = out$ident.full$sv.info$exp.mat[,lrrc16a]

}


beh.analysis$run(inp)

