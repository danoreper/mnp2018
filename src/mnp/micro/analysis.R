library(multcomp)
library(nlme)
library(MASS)
library(data.table)
library(lme4)
library(lattice)
library(reshape)
library(ggplot2)
library(parallel)
library("biomaRt")
##library(lmerTest) #for lmerTest, caTools, bitOps, and possibly others need to be intentionally installed, automatic dependency checking doesnt work right. really annoying.
options(warn=1)

## if(!exists("prop"))
## {
##     stop("master script must load default params")
## }

source("./utils.R") 
source("./genomerep/buildGenomeData2.R")
source("./multipleTesting.R")
source("./lm/fitBoxCoxModels.R")

source("./lm/multivariate.lm.R")
source("./enrichmentTesting.R")
source("./surrogatcalc.R")

source("./mnp/plotting.R")
source("./mnp/micro/preprocess/evalprobes2.R")
source("./mnp/micro/preprocess/extractFromProbes.R")
source("./mnp/loadAllData.R")

micro.analysis <- new.env(hash=T)

micro.analysis$get.SV.func <- function(local = F)
{
    surrogat        = prop$mnp$surrogat
    scaleOnForPCA   = T ##prop$mnp$scaleOnForPCA
    transformParams = fit.model.bc$getDefaultTransformParams()

    
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = F)

    if(surrogat == "SSVA") # uses all control probes
    {
        parallelArgs = micro.analysis$getBestParArgs(100, 3)
        if(local)
        {
            parallelArgs = micro.analysis$getLocalParArgs()
        }
        if(!is.null(parallelArgs$system.type) && parallelArgs$system.type %in% c("killdevil", "longleaf"))
        {
            parallelArgs$timeLimit.hours = .95
        }
        print("using SSVA")
        num.sv = prop$mnp$num.sv
        f = surrogatcalc$get.mm.SSVA.func(num.sv = num.sv,
                                          transformParams = transformParams,
                                          strategy,
                                          scaleOnForPCA = scaleOnForPCA,
                                          residualize   = prop$mnp$ssva.resid,
                                          parallelArgs  = parallelArgs)
        
    } else if(surrogat == "SSVA_PERM_GENERATE") {
        
        parallelArgs = micro.analysis$getBestParArgs(100, 3)

        numPerm       = prop$mnp$SSVA.numperm
        pvalThresh    = prop$mnp$pvalThresh
        f = surrogatcalc$get.mm.SSVA.perm.func(transformParams = transformParams,
                                               pvalThresh = pvalThresh,
                                               numPerm = numPerm,
                                               scaleOnForPCA = scaleOnForPCA,
                                               parallelArgs = parallelArgs)
    } else if(surrogat == "SVA") {
        parallelArgs = micro.analysis$getBestParArgs(100, 3)
        num.sv = prop$mnp$num.sv
        topGenes = prop$mnp$SVA.topgenes
        f  = surrogatcalc$get.mm.SVA.func(num.sv = num.sv,
                                          transformParams,
                                          topGenes = topGenes,
                                          scaleOnForPCA = scaleOnForPCA)
    } 
   return(f)
}

micro.analysis$formCovariateFullModelString <- function(includeStrainByDiet = T, dietBeforeStrain =T)
{
    pipelineString = ""

    modelPipeline   = T
##    modelPipeline = length(prop$mnp$pipelines)>1
    if(modelPipeline)
    {
        pipelineString   = "Pipeline +"
    }
    
    pipelinestrainString = ""
    if(modelPipeline & prop$mnp$modelStrainPipeline)
    {
        pipelinestrainString ="Pipeline:Strain + "  
    }
    
    strainDietString = ""
    if(includeStrainByDiet)
    {
        strainDietString = "Diet:Strain + " 
    }
    
    dietStrainAdditive = "Diet + Strain + "
    if(!dietBeforeStrain)
    {
        dietStrainAdditive = "Strain + Diet + "
    }
    
    covariates = paste0( " ~ 1 + Batch + ", pipelineString, dietStrainAdditive, pipelinestrainString, strainDietString, " (1|Dam.ID)")

    return(covariates)
}

micro.analysis$getResidualizeOutCovariates <- function()
{
    residualizeOutCovariates = c()
    if(prop$mnp$residualize.expression.pre.straindiet)
    {
        residualizeOutCovariates = c("Batch", "Pipeline")
    }
    return(residualizeOutCovariates)
}

micro.analysis$getMainEffectContrasts <- function()
{
    return(c("Diet"))
}

##TODO use this in surrogatcalc
micro.analysis$getBestParArgs <- function(batchSize=1, mem.gb = 22, timeLimit.hours = 23)
{
    if(prop$onCluster)
    {
        parallelArgs = list(system.type     = prop$system.type,
                            filesToSource   = "./mnp/micro/analysis.R",
                            batchSize       = batchSize,
                            timeLimit.hours = timeLimit.hours,
                            cpuMemLimit.GB  = mem.gb,
                            outdir          = prop$tmpdir,
                            saveProp        = T,
                            coresPerJob     = prop$mnp$mc.cores)

    } else {
        parallelArgs = micro.analysis$getLocalParArgs()
    }
    
    return(parallelArgs)
}

micro.analysis$getLocalParArgs <- function()
{
    parallelArgs = list(mc.cores = prop$mnp$mc.cores,
                        mclMult = 10)
    
    return(parallelArgs)
}

micro.analysis$run.noperm <- function(inp)
{
    svFunc                   = micro.analysis$get.SV.func()
    covariateModelString     = micro.analysis$formCovariateFullModelString()
    nullModelString          = NULL
    residualizeOutCovariates = micro.analysis$getResidualizeOutCovariates()
    modelParser              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                               mainEffectContrasts  = micro.analysis$getMainEffectContrasts(),
                                                               medianAdjust.p.value = prop$mnp$medianAdjust.p.values)
        
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = T) 
    
    transformParams          = fit.model.bc$getDefaultTransformParams()
    ##TODO remove
    if(prop$mnp$tiny)
    {
        print("eliminating lambdas")
        transformParams$lambdasToTry = 1
    }



    ident.full = surrogatcalc$runAnalysis(svFunc                   = svFunc,
                                          exp.mat                  = inp$exp.mat,
                                          exp.mat.control          = inp$exp.mat.control.full,
                                          cov.data                 = inp$cov.data.full,
                                          covariateModelString     = covariateModelString,
                                          residualizeOutCovariates = residualizeOutCovariates,
                                          modelParser              = modelParser,
                                          strategy                 = strategy,
                                          transformParams          = transformParams,
                                          parallelArgs             = micro.analysis$getBestParArgs(100, 3))
                                          ##parallelArgs                    = micro.analysis$getLocalParArgs())


    out = list(threshholds = NULL, ident.full = ident.full)
    return(out)
}

micro.analysis$runallPerms <- function(inp, local = F)
{
    
    svFunc                   = micro.analysis$get.SV.func()
    covariateModelString     = micro.analysis$formCovariateFullModelString()
    nullModelString          = NULL
    residualizeOutCovariates = micro.analysis$getResidualizeOutCovariates()

    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = T) 
    transformParams          = fit.model.bc$getDefaultTransformParams()

    ##TODO remove
    if(prop$mnp$tiny)
    {
        print("eliminating lambdas")
        transformParams$lambdasToTry = 1
    }
    parg = micro.analysis$getBestParArgs(100, 3)
    if(local)
    {
        parg = micro.analysis$getLocalParArgs()
    }
    
    print("about to submit sv generation jobs")
    sv.info = surrogatcalc$generate.svinfo(svFunc = svFunc,
                                           exp.mat = inp$exp.mat,
                                           exp.mat.control = inp$exp.mat.control.full,
                                           cov.data = inp$cov.data.full,
                                           covariateModelString     = covariateModelString,
                                           nullModelString          = nullModelString,
                                           residualizeOutCovariates = residualizeOutCovariates,
                                           residualizeOutSV = (length(residualizeOutCovariates)>0),
                                           transformParams  = transformParams,
                                           strategy         = strategy,
                                           parallelArgs     = parg)
##             

    print("got surrogateVariables")
    modelParser              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                               mainEffectContrasts  = micro.analysis$getMainEffectContrasts(),
                                                               medianAdjust.p.value = prop$mnp$medianAdjust.p.values)
    
    
    numPerm = prop$mnp$SSVA.numperm
    
    cov.data             = sv.info$cov.data
    covariateModelString = sv.info$covariateModelString
    ##    exp.mat              = sv.info$exp.mat
##################
    
    cov.data$coded.interaction = factor(paste0(cov.data$Diet, ".", cov.data$Strain))
    ##Special string just for permutation testing of diet by strain effect.
    interactionModelString = unlist(gsub(covariateModelString, pattern = "Strain\\:Diet", replacement = "coded.interaction"))
    
    shuffleVariable = "Dam.ID"
    ##Strain and Diet:Strain permutation
    ##only shuffle diets WITHIN the same strain 
    permTypes = list(
        list(variableOfInterest = "Strain",             covariateModelString = covariateModelString, constraint = as.integer(as.factor(cov.data$Diet))),
        list(variableOfInterest = "coded.interaction",  covariateModelString = interactionModelString, identVar = "Diet:Strain", identString = covariateModelString),
         list(variableOfInterest = "Diet",              covariateModelString = covariateModelString,  constraint=as.integer(as.factor(cov.data$Strain))))
    
    allthresh = list()
    
    for(p in c(2,3,1))
    {
        
        permStatistics = list()
        permType = permTypes[[p]]
        print(paste0("working on permutation for ", permType$variableOfInterest))
        varsToSwap           = permType$variableOfInterest
        if("varsToSwap" %in% names(permType))
        {
            varsToSwap = permType$varsToSwap
        }
        variableOfInterest   = permType$variableOfInterest
        constraint           = permType$constraint
        covariateModelString = permType$covariateModelString
        identVar             = permType$identVar
        identString          = permType$identString
        
        permOut = micro.analysis$generatePerms(exp.mat = sv.info$exp.mat,
                                               cov.data = cov.data,
                                               covariateModelString = covariateModelString,
                                               shuffleVariable = shuffleVariable,
                                               varsToSwap = varsToSwap,
                                               variableOfInterest = variableOfInterest,
                                               constraint = constraint,
                                               alphas = alphas,
                                               modelParser = modelParser,
                                               transformParams = transformParams,
                                               strategy = strategy,
                                               identVar = identVar,
                                               identString = identString)
##        browser()
        save(list = ls(), file = outm("intermediate", paste0("perm_", variableOfInterest)))
        threshholds = permOut$threshholds
        allthresh   = util$appendToList(allthresh, threshholds)
    }

    ##ident.full is the same for every permutation test... a bit of redundancy above, but its not worth optimizing.
    ident.full      = permOut$ident.full
    threshholds     = do.call(rbind, allthresh)

    threshholds$variable[threshholds$variable == "coded.interaction"] =  "Diet:Strain"
    out = list(ident.full = list(results = ident.full, sv.info = sv.info), threshholds = threshholds)
    
  
    return(out)
}
                          

micro.analysis$generatePerms <- function(exp.mat,
                                         cov.data,
                                         covariateModelString,
                                         shuffleVariable,
                                         varsToSwap,
                                         variableOfInterest,
                                         constraint,
                                         alphas,
                                         modelParser,
                                         transformParams,
                                         strategy,
                                         identVar        = NULL,
                                         identString     = NULL)
{
    parallelArgs = micro.analysis$getBestParArgs(1,10, 40)
    parallelArgs$func = fit.model.bc$fit

    parallelArgs$sharedVariables = list(y.mat         = exp.mat,
                                        modelParser     = modelParser,
                                        transformParams = transformParams,
                                        checkAnova = T,
                                        strategy   = strategy)

    accum = parallel$getAccum(parallelArgs)

    ##TODO pass in another accumulator
    
    trueLabels = as.integer(as.factor(cov.data[[shuffleVariable]]))

    print("checking for bug")
    print(shuffleVariable)
    print(dim(trueLabels))
    shuffles = util$generateShuffles(trueLabels = trueLabels, numShuffles = prop$mnp$SSVA.numperm, subjectTo = constraint, identityFirst = T, seed = 3) 
    numPerm  = ncol(shuffles)

    ##The parallel args used within each permutation.
    withinPermParallelArgs = micro.analysis$getLocalParArgs ()
    for(shuffleIndex in 1:numPerm)
    {
        print(paste0("adding shuffle:", shuffleIndex))
        
        covString = covariateModelString
        if(!is.null(identString)&&shuffleIndex==1)
        {
            covString = identString
        }

        cov.shuffled       = cov.data
        shuffle            = shuffles[,shuffleIndex]
        for(varToSwap in varsToSwap)
        {
            cov.shuffled[[varToSwap]] = cov.data[[varToSwap]][shuffle]
        }

        addargz = list(funcArgs = list(cov.data = cov.shuffled, covariateModelString = covString, parallelArgs = withinPermParallelArgs))
        if(prop$onCluster)
        {
            propObj = as.environment(as.list(prop))
            propObj$onCluster = F
            addargz$propObj = propObj
        }
        do.call(accum$addCall, addargz)

            
    }
    
    outs = accum$runAll()
    
    perms = micro.analysis$collatePerms(accum, outs, numPerm, variableOfInterest, identVar)
        
    return(perms)
}


##TODO: move to permutation testing?
micro.analysis$collatePerms <- function(accum, outs, numPerm, variableOfInterest, identVar)
{
    iter = accum$getOutputIterator(outs)
    badind = c()
    i= 0
    while(iter$hasNext())
    {
        i = i+1
        out = iter$nextItem()
        varOfInterest = variableOfInterest
        if(!is.null(identVar) && i==1)
        {
            varOfInterest = identVar
        }

        if (i %% 1000 == 0) { cat(sep="","[",i,"]") }

        if(class(out)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }

        pvalues          = try(out$per.variable[variable == varOfInterest]$anova.p.value)
        if(class(pvalues)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }
            
        if(i==1)
        {
            perm.pvalues = matrix(NA, nrow = numPerm, ncol = length(pvalues))
            ident.full = out
        }
        
        perm.pvalues[i,] = pvalues
    }

    if(length(badind)>0)
    {
        print(paste("bads:", badind))
        perm.pvalues = perm.pvalues[,-badind]
    }
    perms = perm.pvalues[1:nrow(perm.pvalues),]
    ident = perm.pvalues[1,]
    alphas = .05
    threshholds = multipleTesting$evaluate.perms(perms,
                                                 ident,
                                                 alpha = alphas,
                                                 direction = 1)

    save(file = outm("rawperms_block", variableOfInterest), list = ls())
    
    threshholds = data.table(variable       = variableOfInterest,
                             alpha          = alphas,
                             threshhold     = threshholds$empirical.global.threshhold,
                             threshhold.gev = threshholds$empirical.global.gev.threshhold)

    return(list(threshholds= threshholds,
                ident.full = ident.full,
                perm.pvalues = perm.pvalues))
}



