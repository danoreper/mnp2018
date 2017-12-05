source("./loadParams.R")
source("./bayes/processSamples.R")
source("./bayes/fitjags.R")
source("./bsub.R")
source("./lm/fitBoxCoxModels.R")


##source("./mnp/mediationAnalysis.R")
source("./mnp/loadAllData.R")
source("./mnp/behavior/analysis.R")
source("./mnp/micro/analysis.R")
source("./mnp/micro/report.R")
source("./mnp/qpcr/analysis.R")

source("./lm/mediation.R")



mnp.med = new.env(hash=T)

##Stored in mnp.med
source("./mnp/BDmodel.R")
source("./mnp/FGHmodel.R")
source("./mnp/PhenModel.R")

mnp.med$run <- function(input, outcome.id,  mediator.ids,  original.strain.results, numShuffles = 0)
{
    ## Only need to run the prior model once, as it is completely independent of the actual mediator data
    ## gibbs.samples.prior = mnp.med$runAllCandidateMediators(outcome.id      = outcome.id,
    ##                                                        mediator.ids    = mediator.ids[1],
    ##                                                        outcome.mat     = input$exp.mat,
    ##                                                        cov.data        = input$cov.data,
    ##                                                        mediationFunc   = input$mediationFunc,
    ##                                                        lmer.modelString= input$lmer.modelString,
    ##                                                        jags.modelname  = input$priorModel,
    ##                                                        colsToIndex     = input$colsToIndex,
    ##                                                        observeInfo     = input$observeInfo,
    ##                                                        gibbsBatchSize  = input$gibbsBatchSize,
    ##                                                        additionalData  = input$additionalData,
    ##                                                        numShuffles     = 0)

    gibbs.samples = mnp.med$runAllCandidateMediators(outcome.id      = outcome.id,
                                                     mediator.ids    = mediator.ids,
                                                     cov.data        = input$cov.data,
                                                     outcome.mat     = input$outcome.mat,
                                                     
                                                     mediationFunc   = input$mediationFunc,
                                                     lmer.modelString= input$lmer.modelString,
                                                     jags.modelname  = input$jags.modelname,
                                                     colsToIndex     = input$colsToIndex,
                                                     observeInfo     = input$observeInfo,
                                                     gibbsBatchSize  = input$gibbsBatchSize,
                                                     additionalData  = input$additionalData,
                                                     numShuffles     = numShuffles)
    
    matchind = match(gibbs.samples$mediation$mediator.id, original.strain.results$Probe.Set.ID)
    
    gibbs.samples$mediation$mediator_name            = original.strain.results$gene_name[matchind]
    gibbs.samples$mediation$p.strain.on.mediator     = original.strain.results$anova.p.value[matchind]
    gibbs.samples$mediation$imprinted                = original.strain.results$minDistToImprinted[matchind]<5000
    setkey(gibbs.samples$mediation, "p.value.ab")
    
    return(gibbs.samples)
}



##TODO: pass in named vectors/matrixes rather than IDs?
mnp.med$runAllCandidateMediators <- function(outcome.id,
                                             mediator.ids,
                                             cov.data,
                                             outcome.mat,

                                             mediationFunc,
                                             lmer.modelString,
                                             jags.modelname,
                                             colsToIndex,
                                             observeInfo,
                                             gibbsBatchSize,
                                             additionalData,
                                             numShuffles = 0)
{
    
    force(mediationFunc)

    outcome = outcome.mat[,outcome.id]

    
    shuffles = util$generateShuffles(1:nrow(outcome.mat), numShuffles = numShuffles, identityFirst =T)


    cleaned = list()
    for(s in 1:ncol(shuffles))
    {
        if(ncol(shuffles)>1)
        {
            print(paste0("on shuffle: ", s))
        }
        ashuffle = shuffles[,s]
        shuffledOutcomes = outcome.mat[ashuffle,]

        accum = mnp.med$getAccumulator(gibbsBatchSize)
        accum$init(mnp.med$get.gibbs.samples, otherGlobals = list(cov.data = cov.data,
                                                                  mediationFunc = mediationFunc,
                                                                  lmer.modelString = lmer.modelString,
                                                                  jags.modelname   = jags.modelname,
                                                                  colsToIndex      = colsToIndex,
                                                                  observeInfo      = observeInfo,
                                                                  additionalData   = additionalData))

        for(i in 1:length(mediator.ids))
        {
            mediator.id = mediator.ids[i]
            
            mediator = shuffledOutcomes[, mediator.id]
            mediator = mediator[as.character(cov.data$ID)]
            accum$addCall(list(outcome = outcome,
                               outcome.id = outcome.id,
                               mediator = mediator,
                               mediator.id = mediator.id,
                               s.id = s))
        }
        outs = accum$runAll()
        
        if(accum$outputs.files)
        {
            outs = bsub$getAllOutputs(outs, accum)
        }
        
        for(out in outs)
        {
            validItem = suppressWarnings(length(out))[1]
            if(validItem)
            {
                cleaned = util$appendToList(cleaned, out)
            }
        }
    }
    mediation = rbindlist(lapply(cleaned, "[[", "mediation"))
    
    if(ncol(shuffles)>1)
    {
        df = mediation[s.id %in% 2:ncol(shuffles), list(minp =min(p.value)), by = c("s.id", "moderators")]
        t1 = df[,list(FWER.thresh=multipleTesting$minp.gev(minp, alpha = .05)), by = "moderators"]
        mediation = t1[mediation[s.id==1], on="moderators"]
    } 
    mediation$s.id = NULL
    
    samples   = NULL
    ##    samples = lapply(cleaned, "[[", "samples")
    
    return(list(mediation = mediation, samples = samples))
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

mnp.med$get.gibbs.samples <- function(cov.data,
                                      outcome.id,
                                      outcome,
                                      mediator.id,
                                      mediator,
                                      s.id = NULL,

                                      lmer.modelString,
                                      mediationFunc,

                                      jags.modelname,
                                      colsToIndex,
                                      observeInfo,
                                      additionalData)
{

    cov.data$mediator = mediator
    

    ##mnp.med$transformPhenotype(outcome, cov.data, lmer.modelString)                               
    cov.data$y = outcome

    iter= 8000   ##prop$mnp.med$iter                                                               
    nchains = 1 #4
    thin = 1 ##prop$mnp.med$thin
    burninFrac = .2

    samples = suppressWarnings(fitjags$runJags(trainingData= cov.data,
                                               modelname   = textConnection(jags.modelname),
                                               iter        = iter,
                                               burninFrac  = burninFrac,
                                               nchains     = nchains,
                                               thin        = thin,
                                               phen        = "y",
                                               observeInfo = observeInfo,
                                               colsToIndex = colsToIndex,
                                               additionalData = additionalData))


    mediation = mediationFunc(samples, cov.data)
    mediation$mediator.id = mediator.id
    mediation$outcome.id  = outcome.id
    mediation$s.id        = s.id

##    print(mediation)
    return(list(samples=NULL, mediation=mediation))
}


mnp.med$getAccumulator <- function(batchSize)
{
##   accum = bsub$get.mc.accumulator(mc.cores= 1)
    browser()
    if(prop$onCluster)
    {
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = 3, queue="hour")
        accum = bsub$get.bsub.accumulator("./mnp/mediationBayes.R", bsubCommand, batchSize=batchSize)
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
    ## if(class(hpd)=="try-error")
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
        save(file=outf( "sv.c.corrected.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outf( "sv.c.corrected.RData"))
    }

    return(svinfo)
}


##TODO: put all functions below here into a more general mediation module?
mnp.med$getModelString <- function(originalModel, mediatorLikelihood, outcomeLikelihood)
{
    ##if this method is being called, it expects the keyword stubs are in there. check to make sure
    if(!grepl(pattern = "MEDIATOR_LIKELIHOOD", originalModel))
    {
        stop()
    }
    if(!grepl(pattern = "OUTCOME_LIKELIHOOD", originalModel))
    {
        stop()
    }
    
    newModel = gsub(originalModel, pattern = "MEDIATOR_LIKELIHOOD",
                    replacement = mediatorLikelihood)

    newModel = gsub(newModel, pattern = "OUTCOME_LIKELIHOOD",
                    replacement = outcomeLikelihood)

##    cat(newModel)
    return(newModel)
}

mnp.med$getEffectDistribution <- function(simpleMainEffect)
{
    p = paste0
    m = simpleMainEffect

    astring = p(        "for(i in 1:2)\n")
    astring = p(astring,"{\n")
    astring = p(astring,"   precision.",m,"[i] <- pow(bound*sd.multiple, -2)\n")
    astring = p(astring,"   f.",m,"[i,1] <- 0\n")
    astring = p(astring,"   for(j in 2:n.",m,")\n")
    astring = p(astring,"   {\n")
    astring = p(astring,"        f.",m,"[i,j] ~ dnorm(0, precision.",m,"[i])\n")
    astring = p(astring,"   }\n")
    astring = p(astring,"}\n")
    astring = p(astring,"\n\n")
    
    return(astring)
}


mnp.med$getRandomEffectDistribution <- function(randomEffect)
{
    p = paste0
    m = randomEffect

    astring = p(        "for(i in 1:2)\n")
    astring = p(astring,"{\n")
    astring = p(astring,"   std.",m,"[i] ~ dunif(0, bound*sd.multiple)\n")
    astring = p(astring,"   precision.",m,"[i] <- pow(std.",m,"[i],-2)\n")
    astring = p(astring,"   for(j in 1:n.",m,")\n")
    astring = p(astring,"   {\n")
    astring = p(astring,"        ranef.",m,"[i,j] ~ dnorm(0, precision.",m,"[i])\n")
    astring = p(astring,"   }\n")
    astring = p(astring,"}\n")
    astring = p(astring,"\n\n")
    
    return(astring)
}


mnp.med$getLikContribution <- function(simpleMainEffect, type)
{
    m = simpleMainEffect
    p = paste0
    astring = p("      ",type,".",simpleMainEffect,"[i,      ",m,".j[s]]")
    return(astring)
}



##TODO: delete?

## mnp.med$caterplot <- function(gibbs.samples, effect, exact = F, thetitle=phen, limz = NULL)
## {
##     if(!exact)
##     {
##         colz = processSamples$getColsForEffect(gibbs.samples, effect)
##     } else {
##         colz = effect
##     }
##     toPlot = gibbs.samples[,colz]
##     if(!exact)
##     {
##         colnames(toPlot) = processSamples$getShortName(effect, colnames(toPlot))
##     }
    
##     if(!is.null(limz))
##     {
##         caterplot(toPlot, col="black", val.lim = limz)
##     } else
##     {
##         caterplot(toPlot, col="black")
##     }
## ##    caterplot(toPlot,style = "plain", col="black")
##     title(thetitle)
## }

## mnp.med$plotEffects <- function(contr.strain,  effect, exact, atitle, lambdasPerPhen)
## {
##     froot = fp(prop$mnp$output, "ovx", atitle)
##     dir.create(froot, recursive = T, showWarnings=F)
##     gibbs.samples = contr.strain$gibbs.samples
    
##     for(phen in names(gibbs.samples))
##     {
##         samplez = mnp.med$getSamples.for.effect.type(contr.strain, phen, effect, exact = exact, lambdasPerPhen)
##         afile = fp(froot, paste0("cat_", phen, "_",effect,".pdf"))
##         print(afile)
##         pdf(afile)
##         mnp.med$caterplot(samplez, colnames(samplez), paste0(atitle,":", phen), exact = T)
##         dev.off()
##     }
## }


