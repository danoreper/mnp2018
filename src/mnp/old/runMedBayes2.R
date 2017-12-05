# TODO: Add comment
# 
# Author: doreper
###############################################################################

library(ggplot2)
library(data.table)

library(lme4)
library(utils)

##library(runjags)
source("./loadParams.R")
source("./mnp/loadAllData.R")
source("./utils.R")
source("./lm/fitBoxCoxModels.R")
source("./multipleTesting.R")
source("./bayes/processSamples.R")
source("./mnp/mediationBayes2.R")


strain.results = fread(datm("2017-05_strainResults.csv"))
strain.results$Probe.Set.ID = as.character(strain.results$Probe.Set.ID)

cands.lrrc = as.character(strain.results[gene_name=="Lrrc16a"]$Probe.Set.ID)
cands.top =  as.character(strain.results[anova.q.value<.05 & minDistToImprinted<5000]$Probe.Set.ID)
cands.all  = as.character(strain.results$Probe.Set.ID)
cands.imp  = as.character(strain.results[minDistToImprinted<5000]$Probe.Set.ID)
cands.airn = "10441787"
raw.data.BD = loadAllData$createAllInputs()
phenrepo = loadBehavior$getPhenotypeRepository()

#############  bd Phen   #############
print("running phen analysis")
pracma::tic()

allkeeps = list()
##Run phenotype analysis both without missing expression data samples.

mediator.ids.all = list();
##mediator.ids.all[[1]] = cands.all; mediator.ids.all[[2]] = cands.lrrc
##M.measures.all =  c("micro", "qpcr")

mediator.ids.all[[1]] = cands.lrrc
M.measures.all =  c("qpcr")
Y.measures = "behavior"
for(i in 1:length(M.measures.all))
{
    M.measures = M.measures.all[i]
    mediator.ids = mediator.ids.all[[i]]
    
for (discardMissingGeneExpression in c(T,F))
{
    ##the model specifications for every phenotype
    allModelSpecs = mnp.med$getAllMediationPhenModelSpecs(phenrepo)
    
    afunc = function(amodel, M.measures, Y.measures, mediator.ids, discardMissingGeneExpression, strain.results, raw.data.BD)
    {
        outcome.id = paste0(amodel$experiment, "_", amodel$phen)
        print(outcome.id)
        BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = M.measures, Y.measures =Y.measures, discardMissingGeneExpression = T, phenotypeSpec = amodel, raw.data = raw.data.BD)
        output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = mediator.ids, original.strain.results = strain.results)
        return(output)
    }
    force(afunc)
    
    accum = bsub$get.mc.accumulator(mc.cores=1)
##    accum = mnp.med$getAccumulator(1)
    
    accum$init(func= afunc,  otherGlobals = list(M.measures=M.measures, Y.measures=Y.measures, mediator.ids = mediator.ids, discardMissingGeneExpression = discardMissingGeneExpression, strain.results = strain.results, raw.data.BD = raw.data.BD))
    
    ##The model and the dataset changes dpending on the phenotype
    for(i in 1:length(allModelSpecs))
    {
        if( i %% 1000 ==0)
        {
            print(i)
        }
        accum$addCall(funcArgs = list(amodel = allModelSpecs[[i]]))
    }
    outs = accum$runAll()

    if(accum$outputs.files)
    {
        outs = bsub$getAllOutputs(outs, accum)
    }
    cleaned = list()
    for(out in outs)
    {
        validItem = suppressWarnings(length(out)>1)
        if(validItem)
        {
            cleaned = util$appendToList(cleaned, out)
        }
    }
    outputs = rbindlist(cleaned)
    outputs = outputs[moderators == "Diet=Ave"]

    froot = outm("mediation3", "behavior")
    postfix = paste0(discardMissingGeneExpression, "_", paste(M.measures, collapse=","), "_", paste(Y.measures, collapse = ","))
    tokeep = mnp.med$saveOutputs(outputs, froot = paste0(froot, "_", postfix))
    tokeep$postfix = postfix
    allkeeps = util$appendToList(allkeeps, tokeep)
}
}
allkeeps = rbindlist(allkeeps)
pracma::toc()
print("all done with BD Phen")

############################################
print("running BD analysis")


M.measures = "micro"
Y.measures.all = c("micro", "qpcr", "qpcr")
discards.all   = c(T, T, F)  
for (i in 1:length(Y.measures.all))
{
    discard = discards.all[i]
    pracma::tic()
    Y.measures = Y.measures.all[i]
    BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = M.measures, Y.measures = Y.measures, discardMissingGeneExpression = discard, phenotypeSpec = list(phen = cands.lrrc))
    
    output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = cands.all, original.strain.results = strain.results)
    output = output[moderators == "Diet=Ave"]

    froot = outm( "mediation2", "expression")
    postfix = paste0(discard, "_", paste(M.measures, collapse=","), "_", paste(Y.measures, collapse = ","))
    tokeep = mnp.med$saveOutputs(outputs, froot = paste0(froot, "_", postfix))
    
    print("all done")
}
###########################################


## #############  BD expression  #############
## print("running BD analysis")
## pracma::tic()
## BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = "micro", Y.measures =c("qpcr"), discardMissingGeneExpression = F,
##                                                phenotypeSpec = list(phen = cands.lrrc), merge.qpcr.plate = T)
## output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = cands.top, original.strain.results = strain.results)
## output = output[moderators == "Diet=Ave"]
## tokeep = mnp.med$saveOutputs(output, froot = BD.inputbuilder(cands.lrrc)$froot)
## browser()

stop()


#############   FGH  ################
print("running FGH analysis")
FGH.input = mnp.med$get.FGH.inputs()
outcome.id = FGH.input$allcandidates[gene_name=="Lrrc16a"]$Probe.Set.ID

output1 = mnp.med$run(input = FGH.input, outcome.id = outcome.id, mediator.ids = cands.all, original.strain.results = strain.results)
tokeep = mnp.med$saveOutputs(output1, froot = FGH.input$froot)
print("all done with FGH expression")





