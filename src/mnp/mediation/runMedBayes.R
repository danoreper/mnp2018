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
source("./mnp/mediationBayes.R")


strain.results = fread(datm("2017-05_strainResults.csv"))
strain.results$Probe.Set.ID = as.character(strain.results$Probe.Set.ID)

cands.lrrc = as.character(strain.results[gene_name=="Lrrc16a"]$Probe.Set.ID)
cands.top =  as.character(strain.results[anova.q.value<.05 & minDistToImprinted<5000]$Probe.Set.ID)
cands.all  = as.character(strain.results$Probe.Set.ID)
cands.imp  = as.character(strain.results[minDistToImprinted<5000]$Probe.Set.ID)
cands.airn = "10441787"
raw.data.BD = loadAllData$createAllInputs()
phenrepo = loadBehavior$getPhenotypeRepository()



#############   FGH  ################
print("running FGH analysis")
FGH.input = mnp.med$get.FGH.inputs()
outcome.id = FGH.input$allcandidates[gene_name=="Lrrc16a"]$Probe.Set.ID

output1 = mnp.med$run(input = FGH.input, outcome.id = outcome.id, mediator.ids = cands.all, original.strain.results = strain.results)
tokeep = mnp.med$saveOutputs(output1, froot = FGH.input$froot)
print("all done with FGH expression")



#############  BD expression  #############
print("running BD analysis")
pracma::tic()
BD.input  = mnp.med$get.BD.inputs()
output2 = mnp.med$run(input = BD.input, outcome.id = cands.lrrc, mediator.ids = cands.top, original.strain.results = strain.results)
output2 = output2$mediation[moderators == "Diet=Ave"]
pracma::toc()
tokeep = mnp.med$saveOutputs(output2, froot = BD.input$froot)
browser()
print(tokeep)
print("all done with BD expression") 
    
#############  BD   #############
    print("running BD analysis")
    pracma::tic()
    

    BD.input  = mnp.med$get.BD.inputs()
    outcome.id = BD.input$allcandidates[gene_name=="Lrrc16a"]$Probe.Set.ID
    cands.top =  as.character(BD.input$allcandidates[anova.q.value<.05 & minDistToImprinted<5000]$Probe.Set.ID)
    cands.all  = as.character(BD.input$allcandidates$Probe.Set.ID)
    cands.imp  = as.character(BD.input$allcandidates[minDistToImprinted<5000]$Probe.Set.ID)

    browser()
    output2 = mnp.med$run(input = BD.input, outcome.id = outcome.id, mediator.ids = cands.top, numShuffles = 0)
    
    save(file =fp(BD.input$froot, "out.big.RData"), list = ls())
    browser()
    
    output2$mediation[,q.value := p.adjust(p.value, method = "BH"), by = "moderators"]
    output2$mediation$outcome.id = NULL
    output2 = output2$mediation[p.value<.05]
    setcolorder(output2, c("mediator_name", "mediator.id", "moderators", "coef", "p.strain.on.mediator", "p.straindiet.on.mediator", "minDistToImprinted", "zero.mediation.lik", "FWER.thresh", "p.value", "q.value"))
    fwrite(file = "./matnutMediation.txt", output2, sep ="\t")

#############  BD Phen   #############
print("running phen analysis")
pracma::tic()
missingStrategies = rbind(c(F, T),
                          c(T, T))
colnames(missingStrategies) = c("discardMissingExpression", "discardMissingBehavior")
missingStrategies = data.frame(missingStrategies)

allkeeps = list()
##Run phenotype analysis both without missing expression data samples.
for (ro in 1:nrow(missingStrategies))
{
    print(missingStrategies[ro,])
    
    phenInputBuilder = mnp.med$getPhenInputBuilder(discardMissingExpression = missingStrategies$discardMissingExpression[ro],
                                                   discardMissingBehavior   = missingStrategies$discardMissingBehavior[ro])

    ##the model specifications for every phenotype
    allModelSpecs = mnp.med$getAllMediationPhenModelSpecs(phenInputBuilder$phenz)
    
    outputs  = list()

    ##The model and the dataset changes dpending on the phenotype
    for(i in c(1:length(allModelSpecs)))
    {
        amodel = allModelSpecs[[i]]
        outcome.id = paste0(amodel$experiment, "_", amodel$phen)
        print(outcome.id)
        
        phen.input = phenInputBuilder$get.phen.inputs(amodel)
        output = mnp.med$run(input = phen.input,
                              outcome.id = outcome.id,
                              mediator.ids = c(cands.all),
                             original.strain.results = strain.results)

        outputs = util$appendToList(outputs, output$mediation)

    }
    outputs = rbindlist(outputs)
    outputs = outputs[moderators == "Diet=Ave"]

    discardString = paste(unlist(unname(unlist(as.list(missingStrategies[ro,])))), collapse = "_")
    tokeep = mnp.med$saveOutputs(outputs, froot = paste0(phen.input$froot, "_",
                                                         discardString))
                                                            
    
    tokeep$discard = discardString
    allkeeps = util$appendToList(allkeeps, tokeep)
}
pracma::toc()
print("all done with BD Phen")












