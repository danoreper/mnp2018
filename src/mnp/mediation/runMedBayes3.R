# TODO: Add comment
# 
# Author: doreper
###############################################################################

library(ggplot2)
library(data.table)

library(lme4)
library(utils)

source("./util/loadParams.R")
source("./mnp/loadAllData.R")
source("./util/utils.R")
source("./lm/fitBoxCoxModels.R")
source("./multipleTesting.R")
source("./mnp/mediation/mediationBayes3.R")
source("./genomerep/buildGenomeData2.R")


##exons = data.table(data.frame(buildGenomeData$getExons()))

strain.results = fread(outm(fp("micro", "effect.table", "p_all_Strain.csv")))
strain.results$Probe.Set.ID = as.character(strain.results$Probe.Set.ID)

cands.all  = as.character(strain.results$Probe.Set.ID)
cands.imp  = as.character(strain.results[imprinted=="Y"]$Probe.Set.ID)

cands.lrrc = as.character(strain.results[gene_name=="Lrrc16a"]$Probe.Set.ID)
cands.airn = "10441787"
cands.mir341 = as.character(strain.results[gene_name=="Mir341"]$Probe.Set.ID)

tops = strain.results[["-log10.qval"]]> -log10(.05) & strain.results[["imprinted"]]=="Y"
cands.top =  as.character(strain.results[tops]$Probe.Set.ID)


    
mirtargs = fread(datm("TargetScan7.1__miR-341-3p.predicted_targets.txt"))
mirtargs$probeset = NA
for(i in 1:nrow(mirtargs))
{
    
    gn =mirtargs[i][["Target gene"]]
    print(gn)
    ps = strain.results[gene_name==gn|
                        grepl(gene_name, pattern = paste0(gn,","))]$Probe.Set.ID

    
    ## if(length(ps)>1)
    ## {
    ##     browser()
    ## }
    mirtargs$probeset[i]= ps[1]
}


outcomeProbesets = mirtargs$probeset
outcomeProbesets = setdiff(mirtargs$probeset, NA)

raw.data.BD = loadAllData$createAllInputs()
phenrepo    = raw.data.BD$phens 

print("Expression Mediation")


print("running BD mir341 analysis")

sharedVariables = list()
sharedVariables$Y.measures = "micro"
sharedVariables$M.measures = "micro"
sharedVariables$discard = T
sharedVariables$qpcrmerge = T
sharedVariables$raw.data.BD = raw.data.BD
sharedVariables$cands.mir341 = cands.mir341
sharedVariables$strain.results = strain.results

outputs = list()
f = function(outcomeProbeset, Y.measures, M.measures, discard, qpcrmerge, raw.data.BD, cands.mir341, strain.results)
{
    BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = M.measures, Y.measures = Y.measures, discardMissingGeneExpression = discard, phenotypeSpec = list(phen = outcomeProbeset), merge.qpcr.plate = qpcrmerge, raw.data = raw.data.BD)
    
    output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = cands.mir341, original.strain.results = strain.results)$interaction
    output$outcome = strain.results[Probe.Set.ID == outcomeProbeset]$gene_name
    return(output)
}


if(length(cands.all)/50<3)
{
    batchSize = min(3, length(cands.all))
} else {
    batchSize = 50
}

accum = parallel$get.mc.accum(func = f, mc.cores= prop$mnp$mc.cores, sharedVariables = sharedVariables)
if(prop$onCluster & length(cands.all)>1)
{
    accum = parallel$get.cluster.accum(system.type = prop$system.type,
                                       func = f,
                                       sharedVariables = sharedVariables,
                                       
                                       filesToSource = "./mnp/mediation/mediationBayes3.R",
                                       batchSize = batchSize,
                                       timeLimit.hours = ceiling(3 +(20*batchSize)/60),
                                       cpuMemLimit.GB = 2,
                                       outdir         = prop$tmpdir,
                                       saveProp       = T)
}

for(i in 1:2)##length(cands.all)) ##length(outcomeProbesets))
##for(i in 1100:1120) ##length(outcomeProbesets))
{
    print(i)
    outcomeProbeset = cands.all[i]
    accum$addCall(funcArgs = list(outcomeProbeset=outcomeProbeset))
              ##outcomeProbesets[i]
    ##outcomeProbeset = outcomeProbesets[i]
    ##outputs = util$appendToList(outputs, output)
##    output = output[moderators == "Diet=Ave"]
##    froot = outm("mediation")
##    postfix
##    postfix = paste0(discard, "_", qpcrmerge, "_",  paste(M.measures, collapse=","), "_", paste(Y.measures, collapse = ","))
  ##  tokeep = mnp.med$saveOutputs(output, froot = fp(froot,  postfix))
   ## print("all done")
}
out=accum$runAll()
out=accum$getAllOutputs(out)
browser()
out = rbindlist(out)

fwrite(out, outm("mir341_med.txt"))

stop()
############################################
print("running BD Lrrc expression analysis")


M.measures = "micro"
##Y.measures.all = c("micro", "qpcr")
##discards.all   = c(T, T)
Y.measures.all = c("micro")
discards.all = c(T)


for (i in 1:length(Y.measures.all))
{
    discard = discards.all[i]
    pracma::tic()
    Y.measures = Y.measures.all[i]
    qpcrmerges = F
    if("qpcr" %in% Y.measures)
    {
        qpcrmerges = c(T, F)
    }
    for(qpcrmerge in qpcrmerges)
    {
        BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = M.measures, Y.measures = Y.measures, discardMissingGeneExpression = discard, phenotypeSpec = list(phen = cands.lrrc), merge.qpcr.plate = qpcrmerge, raw.data = raw.data.BD)
        
        output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = cands.all, original.strain.results = strain.results)$mediation
        output = output[moderators == "Diet=Ave"]

        froot = outm("mediation")
        postfix = paste0(discard, "_", qpcrmerge, "_",  paste(M.measures, collapse=","), "_", paste(Y.measures, collapse = ","))
        tokeep = mnp.med$saveOutputs(output, froot = fp(froot,  postfix))
        print("all done")
    }
}


#############  bd Phen   #############
print("running BD phen analysis")
pracma::tic()

allkeeps = list()

afunc = function(amodel, M.measures, Y.measures, mediator.ids, discardMissingGeneExpression, strain.results, merge.qpcr.plate)
{
    outcome.id = paste0(amodel$experiment, "_", amodel$phen)
    print(outcome.id)
    BD.inputbuilder  = mnp.med$get.BD.inputBuilder(M.measures = M.measures,
                                                   Y.measures =Y.measures,
                                                   discardMissingGeneExpression = discardMissingGeneExpression,
                                                   phenotypeSpec = amodel,
                                                   merge.qpcr.plate = merge.qpcr.plate)
    output = mnp.med$run(inputBuilder = BD.inputbuilder, mediator.ids = mediator.ids, original.strain.results = strain.results)$mediation
    output = output[moderators == "Diet=Ave"]
    return(output)
}

##Run phenotype analysis both without missing expression data samples.
##missingGeneExpressions = c(F,T)
missingGeneExpressions = c(T)
mediator.ids.all = list();
mediator.ids.all[[1]] = cands.all; mediator.ids.all[[2]] = cands.lrrc
M.measures.all =  c("micro")##, "qpcr")
Y.measures = "behavior"

## mediator.ids.all[[1]] = cands.lrrc
## M.measures.all =  c("qpcr")
print("Behavior mediation")
for (discardMissingGeneExpression in missingGeneExpressions)
{

    for(i in 1:length(M.measures.all))
    {
        M.measures = M.measures.all[i]
        if("qpcr" %in% c(M.measures) &!discardMissingGeneExpression)
        {
            ##this doesn't currently work because missing qpcr data isnt handled correctly yet.
            next
        }

        mediator.ids = mediator.ids.all[[i]]
    
        qpcrmerges = F
        ##don't bother iterating over qpcrmerge = T if there arent any qpcr values to begin with.
        if("qpcr" %in% c(Y.measures, M.measures))
        {
            qpcrmerges = c(T, F)
        } 
        for(qpcrmerge in qpcrmerges)
        {            
            ##the model specifications for every phenotype
            allModelSpecs = mnp.med$getAllMediationPhenModelSpecs(phenrepo)
         
            sharedVariables = list(M.measures=M.measures, Y.measures=Y.measures, mediator.ids = mediator.ids, discardMissingGeneExpression = discardMissingGeneExpression, strain.results = strain.results, merge.qpcr.plate = qpcrmerge)

            accum = parallel$get.mc.accum(func = afunc, mc.cores=1, sharedVariables = sharedVariables)
            if(M.measures=="qpcr" && prop$onCluster)
            {
                accum = parallel$get.cluster.accum (system.type       = prop$system.type,
                                                    func = afunc,
                                                    sharedVariables   = sharedVariables,
                                                    
                                                    filesToSource     = c("./mnp/mediationBayes3.R"),
                                                    batchSize         = 1,
                                                    timeLimit.hours   = .5,
                                                    cpuMemLimit.GB   = 2,
                                                    coresPerJob       = 1,
                                                    saveProp          = T,
                                                    outdir = prop$tmpdir)

            } 

            print("******************************")
            print(paste0("discard missing expression: ", discardMissingGeneExpression))
            print(paste0("qpcr merge: ", qpcrmerge))
            print(paste0("M.measures:",M.measures))
            print(paste0("Y.measures:",Y.measures))
            
            ##The model and the dataset changes dpending on the phenotype
            for(j in 1:length(allModelSpecs))
            {
                print(paste0("Outcome:", allModelSpecs[[j]]$phen))

                accum$addCall(funcArgs = list(amodel = allModelSpecs[[j]]))
            }
   
            outs = accum$runAll()
            outs = accum$getAllOutputs(outs, removeFailing = T)
            outputs = rbindlist(outs)
            
            froot = outm("mediation")
            postfix = paste0(discardMissingGeneExpression, "_", qpcrmerge, "_", paste(M.measures, collapse=","), "_", paste(Y.measures, collapse = ","))
            tokeep = mnp.med$saveOutputs(outputs, froot = fp(froot,  postfix))
            tokeep$postfix = postfix

            allkeeps = util$appendToList(allkeeps, tokeep)
        }
    }
}

allkeeps = rbindlist(allkeeps)
pracma::toc()


stop()


#############   FGH  ################
print("running FGH analysis")
FGH.input = mnp.med$get.FGH.inputs()
outcome.id = FGH.input$allcandidates[gene_name=="Lrrc16a"]$Probe.Set.ID

output1 = mnp.med$run(input = FGH.input, outcome.id = outcome.id, mediator.ids = cands.all, original.strain.results = strain.results)
tokeep = mnp.med$saveOutputs(output1, froot = FGH.input$froot)
print("all done with FGH expression")





