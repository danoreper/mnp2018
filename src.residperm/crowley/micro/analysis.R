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
options(warn=1)

source("./utils.R") 
source("./genomerep/buildGenomeData2.R")
source("./multipleTesting.R")
source("./lm/fitBoxCoxModels.R")

source("./lm/multivariate.lm.R")
source("./enrichmentTesting.R")
source("./surrogatcalc.R")

source("./mnp/micro/preprocess/evalprobes2.R")
source("./mnp/micro/preprocess/extractFromProbes.R")
source("./mnp/loadAllData.R")

crowley.analysis <- new.env(hash=T)

crowley.analysis$get.SV.func <- function()
{
    surrogat        = prop$crowley$surrogat
    scaleOnForPCA   = T 
    transformParams = fit.model.bc$getDefaultTransformParams()

    parallelArgs = crowley.analysis$getBestParArgs(50, 3)
    ##parallelArgs = crowley.analysis$getNoClusterAccum()

    
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = F)

    if(surrogat == "SSVA") # uses all control probes
    {
        print("using SSVA")
        num.sv = prop$crowley$num.sv
        f = surrogatcalc$get.mm.SSVA.func(num.sv = num.sv,
                                          transformParams = transformParams,
                                          strategy,
                                          scaleOnForPCA = scaleOnForPCA,
                                          parallelArgs = parallelArgs)

    } else if(surrogat == "SSVA_PERM_GENERATE") {
        parallelArgs = crowley.analysis$getBestParArgs(100, 3)
        numPerm       = prop$mnp$SSVA.numperm
        pvalThresh    = prop$mnp$pvalThresh
        f = surrogatcalc$get.mm.SSVA.perm.func(transformParams = transformParams,
                                               pvalThresh = pvalThresh,
                                               numPerm = numPerm,
                                               scaleOnForPCA = scaleOnForPCA,
                                               parallelArgs = parallelArgs)

    } else if(surrogat == "SVA") {
        parallelArgs = crowley.analysis$getBestParArgs(100, 3)
        num.sv = prop$mnp$num.sv
        topGenes = prop$mnp$SVA.topgenes
        f  = surrogatcalc$get.mm.SVA.func(num.sv = num.sv,
                                          transformParams,
                                          topGenes = topGenes,
                                          scaleOnForPCA = scaleOnForPCA)
        
    } 
   return(f)
}

crowley.analysis$formCovariateFullModelString <- function()
{
     modelString = paste0("~ 1 + sex + cross + tissue + ",
                          "sex:cross + sex:tissue + cross:tissue  + cross:direction + ",
                          "sex:cross:tissue + sex:cross:direction + tissue:cross:direction +",
                          "sex:cross:tissue:direction +", 
                          " (1|mouse)")

    return(modelString)
}


crowley.analysis$getMainEffectContrasts <- function()
{
    return(c("Diet"))
}

##TODO use this in surrogatcalc
crowley.analysis$getBestParArgs <- function(batchSize=1, mem.gb = 22)
{
    if(prop$onCluster)
    {
        parallelArgs = list(system.type = prop$system.type,
                            filesToSource = "./mnp/micro/analysis.R",
                            batchSize = batchSize,
                            timeLimit.hours = 24,
                            cpuMemLimit.GB = mem.gb,
                            coresPerJob   = prop$mnp$mc.cores)

    } else {
        parallelArgs = crowley.analysis$getLocalParArgs()
    }

    return(parallelArgs)
}

crowley.analysis$getLocalParArgs <- function()
{
    parallelArgs = list(mc.cores = prop$mnp$mc.cores,
                        mclBatch = prop$mnp$mc.cores*100)
    return(parallelArgs)
}


crowley.analysis$run.noperm <- function(inp)
{
    svFunc                   = crowley.analysis$get.SV.func()
    covariateModelString     = crowley.analysis$formCovariateFullModelString()
    nullModelString          = NULL
    residualizeOutCovariates = c() #crowley.analysis$getResidualizeOutCovariates()
    modelParser              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                               mainEffectContrasts  = crowley.analysis$getMainEffectContrasts(),
                                                               medianAdjust.p.value = prop$mnp$medianAdjust.p.values)
        
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = F) ##For some reason, lme fails a lot with this particular model 
    
    ident.full = surrogatcalc$runAnalysis(svFunc                   = svFunc,
                                          exp.mat                  = inp$exp.mat,
                                          exp.mat.control          = inp$exp.mat.control.full,
                                          cov.data                 = inp$cov.data.full,
                                          covariateModelString     = covariateModelString,
                                          residualizeOutCovariates = residualizeOutCovariates,
                                          modelParser              = modelParser,
                                          strategy                 = strategy,
                                          transformParams          = fit.model.bc$getDefaultTransformParams(),
                                          parallelArgs             = crowley.analysis$getBestParArgs(100, 3))


    out = list(threshholds = NULL, ident.full = ident.full)
    return(out)
}

crowley.analysis$runallPerms <- function(inp)
{
    
    svFunc                   = crowley.analysis$get.SV.func()
    covariateModelString     = crowley.analysis$formCovariateFullModelString()
    nullModelString          = NULL
    residualizeOutCovariates = c()##crowley.analysis$getResidualizeOutCovariates()

    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = T) 
    transformParams          = fit.model.bc$getDefaultTransformParams()

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
                                           parallelArgs             = crowley.analysis$getBestParArgs(100, 3))
##                                           accum                    = crowley.analysis$getNoClusterAccum())

    print("got surrogateVariables")
    modelParser              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                               mainEffectContrasts  = crowley.analysis$getMainEffectContrasts(),
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
        list(variableOfInterest = "Strain",            covariateModelString = covariateModelString),
        list(variableOfInterest = "coded.interaction", covariateModelString = interactionModelString, identVar = "Diet:Strain", identString = covariateModelString),
         list(variableOfInterest = "Diet",              covariateModelString = covariateModelString,  constraint=as.integer(as.factor(cov.data$Strain))))
    
    allthresh = list()
    for(p in 1:length(permTypes))
    {
        
        permStatistics = list()
        permType = permTypes[[p]]
        print(paste0("working on permutation for ", permType))
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
        
        permOut = crowley.analysis$generatePerms(exp.mat = sv.info$exp.mat,
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

        save(list = ls(), file = outm("intermediate", paste0("perm_", permType)))
        threshholds = permOut$threshholds
        allthresh   = util$appendToList(allthresh, threshholds)
    }
    
    ##ident.full is the same for every permutation test... a bit of redundancy above, but its not worth optimizing.
    ident.full      = permOut$ident.full
    threshholds     = do.call(rbind, allthresh)
    threshholds$variable[threshholds$variable == "coded.interaction"] = identVar
    out = list(ident.full = list(results = ident.full, sv.info = sv.info), threshholds = threshholds)
    return(out)
}
                          

crowley.analysis$generatePerms <- function(exp.mat,
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
    parallelArgs = crowley.analysis$getBestParArgs(1,10)
    parallelArgs$func = fit.model.bc$fit
    if(prop$onCluster)
    {
        parallelArgs$sharedVariables = list(y.mat         = exp.mat,
                                            modelParser     = modelParser,
                                            transformParams = transformParams,
                                            checkAnova = T,
                                            strategy   = strategy)
    }
    accum = parallel$getAccum(parallelArgs)

    
    trueLabels = as.integer(as.factor(cov.data[[shuffleVariable]]))

    print("checking for bug")
    print(shuffleVariable)
    print(dim(trueLabels))
    shuffles = util$generateShuffles(trueLabels = trueLabels, numShuffles = prop$mnp$SSVA.numperm, subjectTo = constraint, identityFirst = T) 
    numPerm  = ncol(shuffles)

    withinPermParallelArgs = crowley.analysis$getLocalParArgs ()
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

        addargz = list(cov.data = cov.shuffled, covariateModelString = covString, parallelArgs = withinPermParallelArgs)
        if(prop$onCluster)
        {
            propObj = as.environment(as.list(prop))
            propObj$onCluster = F
            addargz$propObj = propObj
        } 
        do.call(accum$addCall, addargz)

    }
    
    outs = accum$runAll()
    perms = crowley.analysis$collatePerms(accum, outs, variableOfInterest, identVar)
        
    return(perms)
}


##TODO: move to permutation testing?
crowley.analysis$collatePerms <- function(accum, outs, variableOfInterest, identVar)
{
    badind = c()
    for(i in 1:length(outs)) ##the first shuffle is the identity permutation
    {
        print(i)
        varOfInterest = variableOfInterest
        if(!is.null(identVar) && i==1)
        {
            varOfInterest = identVar
        }

##        if (i %% 1000 == 0) { cat(sep="","[",i,"]") }
        out             = try(accum$getOut(outs[[i]]))
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
            perm.pvalues = matrix(NA, nrow = length(outs), ncol = length(pvalues))
            ident.full = out
        }
        perm.pvalues[i,] = pvalues
    }

    if(length(badind)>0)
    {
        print(paste("bads:", badind))
        perm.pvalues = perm.pvalues[,-badind]
    }
    perms = perm.pvalues[2:nrow(perm.pvalues),]
    ident = perm.pvalues[1,]
    alphas = .05
    threshholds = multipleTesting$evaluate.perms(perms,
                                                 ident,
                                                 alpha = alphas,
                                                 direction = 1)
    
    threshholds = data.table(variable       = variableOfInterest,
                             alpha          = alphas,
                             threshhold     = threshholds$empirical.global.threshhold,
                             threshhold.gev = threshholds$empirical.global.gev.threshhold)

    return(list(threshholds= threshholds,
                ident.full = ident.full))
}



