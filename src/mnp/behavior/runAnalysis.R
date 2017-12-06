source("./mnp/phen/analysis.R")
source("./mnp/getExpressionData2.R")



##requires gene expression analysis to have already been run; if it hasnt been run, no psids
psids = NA ##if nothing is loaded, no expression will be considered
load.flag = try(load(outb("/gold/all2.RData")))
if(class(load.flag)!="try-error")
{
    tophits = beh$.getTopHits(results, "Strain")
    ##perform analysis with no gene expression, with Lrrc16a, and with the best snord hit.
    psids = NA
    psids = c(psids, tophits[gene_name =="Lrrc16a"]$Probe.Set.ID[1])
    psids = c(psids, tophits[gene_name =="s115"]$Probe.Set.ID[1])
}


##read in all the phenotypes, covariates, and associated frames
phen = readPhens()

##iterate over all the probeset ids, rerunning the analysis.
for(i in 1:length(psids))
{
    psid = psids[i]
    gene_name = "nogene"
    geneExp = NULL
    if(!is.na(psid))
    {
        gene_name = tophits[Probe.Set.ID==psid]$gene_name
        geneExp = list(geneExp = exp.mat[,psid], name = gene_name)
    }
    print(gene_name)
    df  = beh.analysis$run(phen, geneExp = geneExp)
    df1 = beh.analysis$.adjust.pvals(df, phen)
    ##plotVarianceExplained(df1)

    ##    fname = paste(gene_name, psid, "phenEval2Old.csv", sep=".")
    fname = paste(gene_name, psid, "phenNew.csv", sep=".")

    df2= beh.analysis$.dispSignificantPhen(df1, outb(fname))
}

merged = beh$.mergeIntoPipelines(phen=phen)
##PCA analysis of the two pipelines
beh.analysis$.evalPCAphen(merged$pipel1, phen$breedLog, "pipeline_behavior_1")
beh.analysis$.evalPCAphen(merged$pipel2, phen$breedLog, "pipeline_behavior_2")
beh.analysis$.plotIntraPipelineCorrelations(merged)
