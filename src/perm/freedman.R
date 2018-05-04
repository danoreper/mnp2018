source("./lm/formulaWrapper.R")

micro.analysis$runallPermsFreed <- function(inp)
{
    permTypes   = c("Strain", "Diet", "Diet:Strain")
##    permTypes = c("Strain")

    
    covariateModels = list(Strain = micro.analysis$formCovariateFullModelString(includeStrainByDiet = F, dietBeforeStrain = T))
#,
#                           Diet   = micro.analysis$formCovariateFullModelString(includeStrainByDiet = F, dietBeforeStrain = F))
#    covariateModels[["Diet:Strain"]] = micro.analysis$formCovariateFullModelString(includeStrainByDiet = T, dietBeforeStrain = T)

    allthresh = list()
    allperms = list()
    for(permType in names(covariateModels))
    {
        print(paste0("working on permutation for ", permType))
        print(covariateModels[[permType]])
        
        permOut = micro.analysis$generatePermsFreed(exp.mat              = inp$exp.mat,
                                                    exp.mat.control      = inp$exp.mat.control.full,
                                                    cov.data             = inp$cov.data.full,
                                                    covariateModelString = covariateModels[[permType]],
                                                    variableOfInterest   = permType,
                                                    alphas               = .05)

        print("generated all")
        ##save(list = ls(), file = outm("intermediate", paste0("perm_", permType)))
        allthresh   = util$appendToList(allthresh, permOut$perms$threshholds)
        allperms    = util$appendToList(allperms, permOut$perms$perms)

    }
    threshholds     = do.call(rbind, allthresh)

    out = list(results = permOut$ident.full, threshholds = threshholds, perms = allperms)
    return(out)
}

micro.analysis$getPermResidParser <- function()
{
    parser = new.env(hash = T)
    
    parser$parse <- function(fit)
    {
        if(is.null(fit$fit))
        {

            y.mu    = rep(NA, length(fit$y.transformed))
            epsilon = rep(NA, length(fit$y.transformed))
            inversionParams = NA
            rands = NA
        } else {
            
            epsilon = resid(fit$fit)
            rands = ranef(fit$fit)
                        
            y.mu = fit$y.transformed - epsilon - rands[fit$fit$data$Dam.ID, "(Intercept)"]

            inversionParams = list(lambda1 = fit$lambda,
                                   lambda2 = fit$lambda2,
                                   center1 = fit$center1,
                                   scale1  = fit$scale1,
                                   center2 = fit$center2,
                                   scale2  = fit$scale2)

        }
        return(list(y.mu = y.mu, epsilon = epsilon, rands = rands, inversionParams=inversionParams))
    }

    parser$collate <- function(outs, accum)
    {
        out = list()
        out$epsilon      = do.call(cbind, lapply(outs, "[[", "epsilon"))
        out$y.mu         = do.call(cbind, lapply(outs, "[[", "y.mu"))

        out$rands        = do.call(cbind, lapply(outs, "[[", "rands"))
        colnames(out$rands) = colnames(out$epsilon)
        out$invertParams = lapply(outs, "[[", "inversionParams")

        return(out)
    }
    return(parser)

}

micro.analysis$generatePermsFreed <- function(exp.mat,
                                              exp.mat.control,
                                              cov.data,
                                              covariateModelString,
                                              variableOfInterest,
                                              alphas)
{
    svFunc                   = micro.analysis$get.SV.func()
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = T) 
    modelParser              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                               mainEffectContrasts  = micro.analysis$getMainEffectContrasts(),
                                                               medianAdjust.p.value = prop$mnp$medianAdjust.p.values)

    ##    reducedModelString = formulaWrapper$removeEffect(variableOfInterest, covariateModelString)
    reducedModelString = formulaWrapper$removeEffectAndInteractions(variableOfInterest, covariateModelString)$modified.string

##    reducedModelString = covariateModelString
    print(paste0("original model: ", covariateModelString))
    print(paste0("reduced model: ", reducedModelString))
    print("getting resids from reduced")
    
    transformParams          = fit.model.bc$getDefaultTransformParams()
    ##TODO remove
    if(prop$mnp$tiny)
    {
        transformParams$lambdasToTry = 1
    }

    prs = micro.analysis$getBestParArgs(100, 3)
##    prs = micro.analysis$getLocalParArgs()
    parser = micro.analysis$getPermResidParser()
##    debug(parser$parse)
    toperm = surrogatcalc$runAnalysis(svFunc                   = svFunc,
                                      exp.mat                  = exp.mat,
                                      exp.mat.control          = exp.mat.control,
                                      cov.data                 = cov.data,
                                      covariateModelString     = reducedModelString,
                                      modelParser              = parser,
                                      strategy                 = strategy,
                                      residualizeOutSV         = F,
                                      transformParams          = transformParams,
                                      parallelArgs             = prs)$results

    
##    browser()
    print("running full model")
    ident.full = surrogatcalc$runAnalysis(svFunc                   = svFunc,
                                          exp.mat                  = exp.mat,
                                          exp.mat.control          = exp.mat.control,
                                          cov.data                 = cov.data,
                                          covariateModelString     = covariateModelString,
                                          modelParser              = modelParser,
                                          strategy                 = strategy,
                                          transformParams          = transformParams,
                                          parallelArgs             = micro.analysis$getBestParArgs(100, 3))
  ##  browser()
    ## getPermAccum <- function()
    ## {
    ##     parallelArgs      = micro.analysis$getBestParArgs(1, 10)
    ##     parallelArgs$func = surrogatcalc$runAnalysisHelper
    ##     parallelArgs$sharedVariables = list(sv.info = ident.full$sv.info,
    ##                                         nullModelString = NULL,
    ##                                         modelParser     = modelParser,
    ##                                         transformParams = NULL, ## apply no transform at all to the y_pi generated by fitting the reduced model
    ##                                         strategy        = strategy,
    ##                                         parallelArgs    = micro.analysis$getLocalParArgs()) ##dont parallelize within permutation, or too many jobs
        
    ##     accum = parallel$getAccum(parallelArgs)
    ##     return(accum)
    ## }

    #GARBAGE
    getPermAccum <- function()
    {
        parallelArgs      = micro.analysis$getBestParArgs(1, 10)
        parallelArgs$func = surrogatcalc$runAnalysis
        parallelArgs$sharedVariables = list(
            svFunc          = svFunc,
            exp.mat.control = exp.mat.control,
            cov.data        = cov.data,
            covariateModelString = covariateModelString,
            modelParser     = modelParser,
            strategy        = strategy,
            transformParams = NULL,
            parallelArgs    = micro.analysis$getLocalParArgs())
        
        accum = parallel$getAccum(parallelArgs)

        return(accum)
    }

    
    accum = getPermAccum()
    shuffles     = util$generateShuffles(trueLabels = 1:nrow(exp.mat), numShuffles = prop$mnp$SSVA.numperm, identityFirst = T)
    shuffles.dam = util$generateShuffles(trueLabels = rownames(toperm$rands), numShuffles = prop$mnp$SSVA.numperm, identityFirst = T)

    numPerm  = ncol(shuffles)

    
    for(shuffleIndex in 1:numPerm)
    {
        print(paste0("adding shuffle:", shuffleIndex))
        shuffle                = shuffles[,shuffleIndex]

##        browser()
        
        shuffle.dam            = shuffles.dam[,shuffleIndex]
        names(shuffle.dam)     = rownames(toperm$rands)
        
        y.shuffled             = toperm$y.mu +
                                 toperm$epsilon[shuffle,] +
                                 toperm$rands[shuffle.dam[cov.data$Dam.ID],]

        rownames(y.shuffled)   = rownames(cov.data)
##        accum$addCall(funcArgs = list(y.mat = y.shuffled))
        accum$addCall(funcArgs = list(exp.mat = y.shuffled))
    }

    outs = accum$runAll()
    perms = micro.analysis$collatePermsFreed(accum, outs, numPerm, variableOfInterest)
        
    return(list(perms=perms,ident.full = ident.full))
}


##TODO: move to permutation testing?
micro.analysis$collatePermsFreed <- function(accum, outs, numPerm, variableOfInterest)
{
##    browser()
    print("collating perms")
    iter = accum$getOutputIterator(outs)
    badind = c()
    i= 0
    perm.pvalues = NULL
    while(iter$hasNext())
    {
        i = i+1
        out = iter$nextItem()
        ##GARBAGE
        out = try(out$results)

        if(class(out)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }

        varOfInterest = variableOfInterest

        if (i %% 1000 == 0) { cat(sep="","[",i,"]") }

        if(class(out)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }

        pvalues          = try(out$per.variable[variable == varOfInterest]$anova.p.value)
        namez            = try(out$per.variable[variable == varOfInterest]$Probe.Set.ID)
        if(class(pvalues)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }

        if(is.null(perm.pvalues))
        {
            perm.pvalues = matrix(NA, nrow = numPerm, ncol = length(pvalues))
            colnames(perm.pvalues) = namez
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


    save(file = outm("rawperms", variableOfInterest), list = ls())
    
    threshholds = multipleTesting$evaluate.perms(perms,
                                                 ident,
                                                 alpha = alphas,
                                                 direction = 1)
    
    threshholds = data.table(
        
        variable       = variableOfInterest,
        alpha          = alphas,
        threshhold     = threshholds$empirical.global.threshhold,
        threshhold.gev = threshholds$empirical.global.gev.threshhold)

    return(list(threshholds = threshholds, perms = perms))
}


