source("./loadParams.R")
source("./bayes/processSamples.R")

source("./bsub.R")
source("./lm/fitBoxCoxModels.R")

source("./mnp/loadAllData.R")
source("./mnp/behavior/analysis.R")
source("./mnp/micro/analysis.R")
source("./mnp/micro/report.R")
source("./mnp/qpcr/analysis.R")

source("./lm/mediation.R")
source("./mnp/mediationjags.R")


mnp.med = new.env(hash=T)

##Stored in mnp.med
source("./mnp/BDmodel2.R")
##source("./mnp/FGHmodel.R")
##source("./mnp/PhenModel.R")

##TODO store outcome.id in inputBuilder?
mnp.med$run <- function(inputBuilder,  mediator.ids,  original.strain.results)
{
    mediation = mnp.med$runAllCandidateMediators(mediator.ids    = mediator.ids,
                                                 inputBuilder    = inputBuilder)


    matchind = match(mediation$mediator.id, original.strain.results$Probe.Set.ID)
    
    mediation$mediator_name            = original.strain.results$gene_name[matchind]
    mediation$p.strain.on.mediator     = original.strain.results$anova.p.value[matchind]
    mediation$imprinted                = original.strain.results$minDistToImprinted[matchind]<5000
    setkey(mediation, "p.value.ab")

    ##TODO: call saveoutputs in here?
    return(mediation)
}



##TODO: pass in named vectors/matrixes rather than IDs?
mnp.med$runAllCandidateMediators <- function(inputBuilder,
                                             mediator.ids)
{
    sampleInput = inputBuilder(mediator.ids[1])

    accum = mnp.med$getAccumulator(sampleInput$gibbsBatchSize)
    accum$init(mnp.med$get.gibbs.samples, otherGlobals = list(inputBuilder))
    for(i in 1:length(mediator.ids))
    {
        mediator.id = mediator.ids[i]
        accum$addCall(list(mediator.id = mediator.id))
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
    mediation = rbindlist(cleaned)
    
    return(mediation)
}        

mnp.med$saveOutputs <- function(froot, outputs)
{
    dir.create(froot, recursive = T, showWarnings = F)

##    outputs = outputs$mediation
    fwrite(outputs, file = fp(froot,"mediationAll.csv"))
    outputs$suppressor = outputs$coef.ab*outputs$coef.c < 0
    tokeep = c("outcome.id", "mediator_name",  "imprinted", "moderators", "suppressor", "coef.c",  "coef.ab", "coef.a", "coef.b", "p.value.c", "p.value.ab")
    tokeep = outputs[,tokeep, with = F]
    fwrite(tokeep, file = fp(froot,"mediationHighlights.csv"))
    return(tokeep)
}

mnp.med$get.gibbs.samples <- function(inputBuilder,
                                      mediator.id)
{
    input = inputBuilder(mediator.id)
    nchains = 1 #4
    thin = 1 ##prop$mnp.med$thin
    burninFrac = .2

    iterAdapt  = burninFrac*input$iter
    iterSample = (1-burninFrac)*input$iter

    print("about to call jags.model")

    mytic()
   
    jags <- suppressWarnings(jags.model(file = textConnection(input$jagsmodel), data=input$dataForJags, n.chains = nchains, n.adapt = iterAdapt,quiet=!showTiming))
    mytoc("total time for jags burnin")

    
    mytic()
    vals = coda.samples(model=jags, variable.names=input$observeInfo, n.iter = iterSample, thin=thin, quiet=!showTiming,progress.bar="none")[[1]]
    mytoc("Total time for jags recorded sampling")


    mediation             = input$mediationFunc(vals, input$dataForJags)
    mediation$mediator.id = mediator.id
    mediation$outcome.id  = input$outcome.id
    
##    print(mediation)
    return(mediation)
}


mnp.med$getAccumulator <- function(batchSize)
{
##   accum = bsub$get.mc.accumulator(mc.cores= 1)
    
    if(prop$onCluster)
    {
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = 3, queue="hour")
        accum = bsub$get.bsub.accumulator("./mnp/mediationBayes2.R", bsubCommand, batchSize=batchSize)
    } else {
        
        batchSize = prop$mnp$mc.cores*100
        accum = bsub$get.mc.accumulator(mc.cores= prop$mnp$mc.cores)
    }
    return(accum)
}

mnp.med$samplesToMediation <- function(a, b, c.prime)
{
    ##just to make code terser; c is cprime.
    c = c.prime
    ab = a*b
    hpd.ab = HPDinterval(ab)
    ## if(class(hpd.ab)=="try-error")
    ## {

    ## }
    
    out = data.frame(
        p.value.ab  = mnp.med$pval(ab),
        p.value.c   = mnp.med$pval(c),
        p.value.a   = mnp.med$pval(a),
        p.value.b   = mnp.med$pval(b),


        lb.ab       = hpd.ab[1,1],
        ub.ab       = hpd.ab[1,2],
        coef.ab     = round(median(ab)*100)/100,
        coef.c      = round(median(c)*100)/100,
        coef.a      = round(median(a)*100)/100,
        coef.b      = round(median(b)*100)/100)

##        zero.mediation.lik = density(ab, from =0, to=0, n=1)$y)
   
    return(out)
}

##TODO: put into utils or multple test correction or bayes?
mnp.med$pval <- function(samplez)
{
    p.value  = (min(sum(samplez>0), sum(samplez<0)))/length(samplez)
    return(p.value)
}

mnp.med$get.sv.corrected.genes <- function(inp, fromFile = T)
{

    residualizeOutCovariates = c()
    if(!fromFile)
    {
        pracma::tic()
        svinfo = surrogatcalc$generate.svinfo(svFunc = micro.analysis$get.SV.func(),
                                              exp.mat = inp$exp.mat,
                                              exp.mat.control = inp$exp.mat.control,
                                              cov.data = inp$cov.data.full,
                                              covariateModelString = micro.analysis$formCovariateFullModelString(),
                                              residualizeOutCovariates = residualizeOutCovariates,
                                              residualizeOutSV = T,
                                              accum = bsub$get.mc.accumulator(mc.cores = prop$mnp$mc.cores))
                
        print(paste0("got surrogate variables:",pracma::toc()))
        gc()
        
        print("done getting resids")
        save(file= outm("sv.c.corrected.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outm("sv.c.corrected.RData"))
    }

    return(svinfo)
}
