library(rjags)
source("./loadParams.R")

source("./parallel/accumulator.R")
source("./lm/fitBoxCoxModels.R")
source("./lm/mediation.R")

source("./mnp/loadAllData.R")
source("./mnp/behavior/analysis.R")
source("./mnp/micro/analysis.R")
source("./mnp/micro/report.R")
source("./mnp/qpcr/analysis.R")



mnp.med = new.env(hash=T)

##Stored in mnp.med
source("./mnp/mediation/BDmodel2.R")
##source("./mnp/FGHmodel.R")

##TODO store outcome.id in inputBuilder?
mnp.med$run <- function(inputBuilder,  mediator.ids,  original.strain.results)
{
    print("calling run!")
    mediation.batch = mnp.med$runAllCandidateMediators(mediator.ids    = mediator.ids,
                                                       inputBuilder    = inputBuilder)


    mediation = mediation.batch$mediation
    interacts = mediation.batch$interaction

    matchind = match(mediation$mediator.id, original.strain.results$Probe.Set.ID)
    mediation$mediator_name            = original.strain.results$gene_name[matchind]
    mediation$p.strain.on.mediator     = original.strain.results$anova.p.value[matchind]
    mediation$imprinted                = original.strain.results$imprinted[matchind]
    setkey(mediation, "p.value.ab")

    matchind = match(interacts$mediator.id, original.strain.results$Probe.Set.ID)
    interacts$gene_name                = original.strain.results$gene_name[matchind]
    interacts$imprinted                = original.strain.results$imprinted[matchind]
    
    

    ##TODO: call saveoutputs in here?
    return(list(mediation=mediation, interaction = interacts))
}



##TODO: pass in named vectors/matrixes rather than IDs?
mnp.med$runAllCandidateMediators <- function(inputBuilder,
                                             mediator.ids)
{
    sampleInput = inputBuilder(mediator.ids[1])
    if(length(mediator.ids)/sampleInput$gibbsBatchSize<3)
    {
        batchSize = min(3, length(mediator.ids))
    } else {
        batchSize = sampleInput$gibbsBatchSize
    }
    
    sharedVariables = list(inputBuilder = inputBuilder)
    accum = parallel$get.mc.accum(func = mnp.med$get.gibbs.samples,mc.cores= prop$mnp$mc.cores, sharedVariables = sharedVariables)
    if(prop$onCluster & length(mediator.ids)>1)
    {
        accum = parallel$get.cluster.accum(system.type = prop$system.type,
                                           func = mnp.med$get.gibbs.samples,
                                           sharedVariables = sharedVariables,

                                           filesToSource = "./mnp/mediation/mediationBayes3.R",
                                           batchSize = batchSize,
                                           timeLimit.hours = ceiling(3 +(20*batchSize)/60),
                                           cpuMemLimit.GB = 2,
                                           outdir         = prop$tmpdir,
                                           saveProp       = T)
    }

    for(i in 1:length(mediator.ids))
    {
        mediator.id = mediator.ids[i]
        accum$addCall(list(mediator.id = mediator.id))
    }
        
    outs = accum$runAll()
    outs = accum$getAllOutputs(outs, removeFailing = T)
    
    df = rbindlist(lapply(outs, "[[", "df") )

    F = data.frame(
        mediator.id = mediator.ids,
        outcome.id  = df$outcome.id[1],
        F1  = unlist(lapply(outs, "[[", "F1")),
        F2  = unlist(lapply(outs, "[[", "F2")))
    
##    mediation = rbindlist(outs)
    return(list(mediation = df, interaction = F))
}        

mnp.med$saveOutputs <- function(froot, outputs)
{
    dir.create(froot, recursive = T, showWarnings = F)

##    outputs = outputs$mediation
    fwrite(outputs, file = fp(froot,"mediationAll.csv"))
    outputs$suppressor = outputs$coef.ab*outputs$coef.c < 0
    tokeep = c("outcome.id", "mediator_name",  "imprinted", "moderators", "suppressor", "coef.c",  "coef.ab", "coef.a", "coef.b", "p.value.c", "p.value.ab")
    tokeep = outputs[,tokeep, with = F]
    setorder(tokeep, p.value.ab)
    fwrite(tokeep, file = fp(froot,"mediationHighlights.csv"))
    return(tokeep)
}

mnp.med$get.gibbs.samples <- function(inputBuilder,
                                      mediator.id)
{

    input = inputBuilder(mediator.id)
    print(input$outcome.id)
    nchains = 1 #4
    thin = 1 ##prop$mnp.med$thin
    burninFrac = .2

    iterAdapt  = burninFrac*input$iter
    iterSample = (1-burninFrac)*input$iter

##    print("about to call jags.model")

    mytic()
##    browser()
    jags <- (jags.model(file = textConnection(input$jagsmodel), data=input$dataForJags, n.chains = nchains, n.adapt = iterAdapt)) 
    ## jags <- suppressWarnings(jags.model(file = textConnection(input$jagsmodel), data=input$dataForJags, n.chains = nchains, n.adapt = iterAdapt,quiet=!showTiming))
    mytoc("total time for jags burnin")

    
    mytic()
    
    vals = coda.samples(model=jags, variable.names=input$observeInfo, n.iter = iterSample, thin=thin, quiet=!showTiming,progress.bar="none")[[1]]
    mytoc("Total time for jags recorded sampling")


    mediation             = input$mediationFunc(vals, input$dataForJags)


    
    mediation$df$mediator.id = mediator.id
    mediation$df$outcome.id  = input$outcome.id
    
##    print(mediation)
    return(mediation)
}


mnp.med$samplesToMediation <- function(a, b, c.prime)
{
    ##just to make code terser; c is cprime.
    c = c.prime
    ab = a*b
    hpd.ab = HPDinterval(ab)
    ## if(class(hpd.ab)=="try-error")
    ## {
    ##     browser()
    ## }
    
    out = data.frame(
        p.value.ab  = mnp.med$pval(ab),
        p.value.c   = mnp.med$pval(c),
        p.value.a   = mnp.med$pval(a),
        p.value.b   = mnp.med$pval(b),


        lb.ab       = hpd.ab[1,1],
        ub.ab       = hpd.ab[1,2],
        coef.ab     = median(ab),
        coef.c      = median(c),
        coef.a      = median(a),
        coef.b      = median(b))
##        zero.mediation.lik = density(ab, from =0, to=0, n=1)$y)
    return(out)
}

##TODO: put into utils or multple test correction or bayes?
mnp.med$pval <- function(samplez)
{
    p.value  = (min(sum(samplez>0), sum(samplez<0)))/length(samplez)
    return(p.value)
}

mnp.med$getLocalParArgs <- function()
{
    parallelArgs = list(mc.cores = prop$mnp$mc.cores,
                        mclBatch = prop$mnp$mc.cores*100)
    
}


mnp.med$get.sv.corrected.genes <- function(inp, fromFile = T)
{
    if(!fromFile)
    {
        parallelArgs = mnp.med$getLocalParArgs()
        
        if(prop$onCluster)
        {
            parallelArgs = list(system.type     = prop$system.type,
                                batchSize       = 100,
                                timeLimit.hours = 2,
                                cpuMemLimit.GB  = 4,
                                outdir          = prop$tmpdir,
                                saveProp        = T,
                                coresPerJob     = prop$mnp$mc.cores)
        }

        residualizeOutCovariates = c()
        nullModelString = NULL
        strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = T) 
        transformParams          = fit.model.bc$getDefaultTransformParams()

        print("getting surrogate variable info")
        pracma::tic()
        svinfo = surrogatcalc$generate.svinfo(svFunc = micro.analysis$get.SV.func(),
                                              exp.mat = inp$exp.mat,
                                              exp.mat.control = inp$exp.mat.control,
                                              cov.data = inp$cov.data.full,
                                              covariateModelString = micro.analysis$formCovariateFullModelString(),
                                              residualizeOutCovariates = residualizeOutCovariates,
                                              residualizeOutSV = T,
                                              transformParams = transformParams,
                                              strategy = strategy,
                                              parallelArgs = parallelArgs)
                
        print(paste0("got surrogate variables:",pracma::toc()))
        gc()
        
        print("done getting resids")
        save(file=outm("sv.c.corrected.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outm("sv.c.corrected.RData"))
        source("./mnp/mediation/mediationBayes3.R")
        source("./mnp/mediation/BDmodel2.R")
    }

    return(svinfo)
}



##TODO move this elsewhere-- mediationBayes?
mnp.med$buildTrainingData <- function(trainingData, colsToIndex)
{
    dataForJags = list()
    trainingData = mnp.med$formIntegerRepresentationOfFactors(trainingData, colsToIndex)
    dataForJags = c(dataForJags, trainingData)
    dataForJags = c(dataForJags, mnp.med$countLevels(trainingData, colsToIndex))    
    return(dataForJags)
}


mnp.med$countLevels <- function(trainingData, colsToIndex)
{
    dataForJags = list()
    for(cname in colsToIndex)
    {
        ##jags needs to know the count of levels of each factor
        dataForJags[[paste0(cname, ".n")]]       = length(levels(trainingData[[cname]]))
    }
    return(dataForJags)
}

mnp.med$formIntegerRepresentationOfFactors <- function(trainingData, colsToIndex)
{
    for(cname in colsToIndex)
    {
        if(cname %in% colnames(trainingData))
        {
            ##jags needs an integer representation of any levels of a factor- store these explictly in the data frame
            trainingData[[cname]] = as.factor(trainingData[[cname]])
            valz = trainingData[[cname]]		
            trainingData[[paste0(cname, ".k")]]  = as.integer(valz)
        }
    }
    
    return(trainingData)
}
