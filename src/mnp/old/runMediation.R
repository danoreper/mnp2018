##imp.genelist = fread(dat(prop$genome$imprintedGenesFile))$ensembl_gene_id

source("./loadParams.R")
source("./mnp/behavior/analysis.R")
source("./mnp/micro/analysis.R")
source("./mnp/micro/report.R")
source("./mnp/qpcr/analysis.R")
source("./mnp/loadAllData.R")
source("./mnp/mediationAnalysis.R")

source("./lm/mediation.R")



 inp  = loadAllData$createAllInputs()
## ##mnp surrogate variable corrected gene expression and covariates
 sv.info = mediation$get.sv.corrected.genes(inp, T)
 exp.mat  = sv.info$exp.mat
 cov.data = sv.info$cov.data


##mnp results for strain variable

load(datm("/perm.RData"))
strain.results = micro.report$postProcessResults(out$ident.full$results, inp$annot.data, cov.data, inp$probesetInfo)$full[variable == "Strain"]

strain.results = micro.report$postProcessResults(out$ident.full$results, inp$annot.data, inp$cov.data.full, inp$probesetInfo)$full[variable == "Strain"]                                                                                                         

setnames(strain.results, old = "anova.p.value", new = "strain.p.value")

strain.diet.results = micro.report$postProcessResults(out$ident.full$results, inp$annot.data, inp$cov.data.full, inp$probesetInfo)$full[variable == "Diet:Strain"]

strain.results$strain.diet.p.value = strain.diet.results$anova.p.value[match(strain.results$Probe.Set.ID, strain.diet.results$Probe.Set.ID)]                                                                                                              

fwrite(strain.results, datm("strainResults.csv"))


##outcome id
outcome.psid = strain.results[grepl(gene_name, pattern="Lrrc16a")]$Probe.Set.ID



##Mapping from gene name to geneid
ginfo = buildGenomeData$getGeneInfo()
gID = ginfo$gene_id
names(gID) = ginfo$gene_name
gName = ginfo$gene_name
names(gName) = ginfo$gene_id

outcome = exp.mat[, outcome.psid]
outcome = outcome[cov.data$ID]



##the subset of genes to focus on...
candidates = strain.results[minDistToImprinted<5000|gene_name=="Lrrc16a"]
#candidates = candidates[!grepl(gene_name, pattern="Lrrc16a")]
candidates = candidates[anova.q.value<.05]
setkey(candidates, "anova.p.value")
browser()


nullcandidates = strain.results[sample(1:nrow(strain.results),size = 1000)]


## ##base string-- hardcode for crowley
## mediatorModelString = sv.info$covariateModelString
## mediatedFactor      = "Strain"



source("./loadParams.R")
source("./mnp/behavior/analysis.R")
source("./mnp/micro/analysis.R")
source("./mnp/micro/report.R")
source("./mnp/qpcr/analysis.R")
source("./mnp/loadAllData.R")
source("./mnp/mediationAnalysis.R")
source("./lm/mediation.R")


dfs = list()
df.nulls = list()


df = mediation$checkMediationForCandidates(candidates=candidates,
                                           exp.mat=exp.mat,
                                           outcome=outcome,
                                           cov.data=cov.data,
                                           mediatorModelString=mediatorModelString,
                                           mediatedFactor=mediatedFactor)
df$cross = "BD_DB"
dfs = util$appendToList(dfs, df) 

## df.null =  checkMediationForCandidates(candidates=nullcandidates,
##                                        exp.mat=exp.mat,
##                                        outcome=outcome,
##                                        cov.data=cov.data,
##                                        mediatorModelString=mediatorModelString,
##                                        mediatedFactor=mediatedFactor)
## df.null$cross = df$cross
## df.nulls = util$appendToList(df.nulls, df.null)

##dfs[variable =="partial.mediation"]



##crowley inputs... TODO build a separate pipeline for this.
inp.crowley = mediation$getCrowleyMicroData()


##the subset of genes to focus on...
## candidates = strain.results[minDistToImprinted<5000]
## candidates = candidates[!grepl(gene_name, pattern="Lrrc16a")]
## ##candidates = candidates[anova.q.value<.05]
## setkey(candidates, "anova.p.value")



 mediatorModelString = "~ sex + tissue + cross + cross:sex + cross:tissue + cross:direction + cross:direction:tissue + (1 | mouse )"
 mediatedFactor = "cross:direction"

## mediatorModelString = "~ sex + tissue + direction + direction:sex + direction:tissue  + (1 | mouse )"
## mediatedFactor = "direction"


valid.crosses = list(c("FG","GF"), c("GH", "HG"), c("FH","HF"))

## for(valid.cross in valid.crosses)
##{
    exp.mat = inp.crowley$exp.mat
    cov.data = inp.crowley$cov.data ##[cross %in% valid.cross, drop=T]
    cov.data[,mouse:=factor(mouse)]

outcome = exp.mat[, outcome.psid]
names(outcome) = rownames(exp.mat)

outcome = outcome[cov.data$ID]

    df = mediation$checkMediationForCandidates(candidates=candidates,
                                               exp.mat=exp.mat,
                                               outcome=outcome,
                                               cov.data=cov.data,
                                               mediatorModelString=mediatorModelString,
                                               mediatedFactor=mediatedFactor)

##    df$cross = paste(valid.cross, collapse = "_")
    dfs = util$appendToList(dfs, df)


    ## df.null = checkMediationForCandidates(candidates=nullcandidates,
    ##                                       exp.mat=exp.mat,
    ##                                       outcome=outcome,
    ##                                       cov.data=cov.data,
    ##                                       mediatorModelString=mediatorModelString,
    ##                                       mediatedFactor=mediatedFactor)

    ## df.null$cross = paste(valid.cross, collapse = "_")
    ## df.nulls = util$appendToList(df.nulls, df.null) 
    
##}
df = rbindlist(dfs)

## df.null = rbindlist(df.nulls)
## save(list = ls(), file = "./mediationNull.RData")


regulatedProbes =df[variable == "factor.on.mediator" & anova.p.value<.05]$mediator_probeset
df.reg = df[df$mediator_probeset %in% regulatedProbes]

browser()




##comparison of chipbase hits to imprinted gene list.
tf.genelist = gID[tf$Transcription.factor]

hits = intersect(imp.genelist, tf.genelist)
##chr13	24171246	24171366

full    = fread("../output/sva/report/fullmerged.csv")
lrrc16a = as.character(full[gene_name=="Lrrc16a"]$Probe.Set.ID[1])
wt1     = as.character(full[gene_id=="ENSMUSG00000016458"]$Probe.Set.ID[1])

