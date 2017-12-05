print("starting")
####################3 
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
alphas = c(.05)

fromFile = F
if(!fromFile)
{
    print("running no perm")
    out = micro.analysis$run.noperm(inp)
    print("done running no perm")
    util$save(list =ls(), file = outm("noperm.RData"))
    
    if(prop$mnp$SSVA.numperm>0)
    {
        print("running all perms")
        out = micro.analysis$runallPerms(inp)
        print("done all perms")
        util$save(file = outm("perm.RData"), list = ls())
    }
} else {
    load(outm("perm.RData"))
    sources()
}

micro.report$reportAnalysis(inp$exp.mat,
                            inp$cov.data.full,
                            out$ident.full$results,
                            ##                            out$results$results,
                            inp$annot.data,
                            inp$probesetInfo,
                            out$threshholds,
                            inp$karyo)

qpcr.analysis$run(inp)

beh.analysis$run(inp)



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

