source("./lm/formulaWrapper.R")

micro.analysis$runallPermsFreed <- function(inp)
{
##    permTypes   = c("Strain", "Diet", "Diet:Strain")
##    permTypes = c("Strain")

##    covariateModels = list()
    covariateModels = list(Strain = micro.analysis$formCovariateFullModelString(includeStrainByDiet = F, dietBeforeStrain = T),
                            Diet   = micro.analysis$formCovariateFullModelString(includeStrainByDiet = F, dietBeforeStrain = F))
    covariateModels[["Diet:Strain"]] = micro.analysis$formCovariateFullModelString(includeStrainByDiet = T, dietBeforeStrain = T)

    allthresh = list()
    allperms  = list()
    alllik    = list()

    dfsimples = list() 
    for(permType in names(covariateModels))
    {
        print(paste0("working on permutation for ", permType))
        print(covariateModels[[permType]])
        
        permOut = micro.analysis$generatePermsFreed(exp.mat              = inp$exp.mat,
                                                    exp.mat.control      = inp$exp.mat.control.full,
                                                    cov.data             = inp$cov.data.full,
                                                    covariateModelString = covariateModels[[permType]],
                                                    variableOfInterest   = permType,
                                                    alpha               = .05)

        
##        sto,p()
        print("generated all")
        ##save(list = ls(), file = outm("intermediate", paste0("perm_", permType)))
        allthresh   = util$appendToList(allthresh, permOut$perms$threshholds)
        allperms[[permType]] = permOut$perms$perms
        alllik  [[permType]] = permOut$perms$lik.rats


        result = permOut$ident.full$results
        dfsimple = result$dfsimple 
##        browser()
        dfsimple$variable      = permType
        
        dfsimples = util$appendToList(dfsimples, dfsimple) 
        
        ##dfalt    = result$dfalt      
##        per.variable = dfalt$per.variable[variable == permType]
        
        ## browser()
        
        if(permType=="Diet:Strain")
        {
            dfalt = result$dfalt
        }
    }
   
    threshholds     = rbindlist(allthresh)
    dfsimples       = rbindlist(dfsimples)
    per.variable    = dfalt$per.variable
    per.variable    = per.variable[variable %in% names(covariateModels)]
    per.variable$anova.p.value.old = per.variable$anova.p.value
    per.variable$anova.p.value = NULL
    per.variable[dfsimples, on = c("Probe.Set.ID", "variable")]

    per.variable = dfsimples[per.variable, on = c("Probe.Set.ID", "variable")]
    setnames(per.variable, "p.value", "anova.p.value")
    dfalt$per.variable = per.variable
    
    out = list(results = dfalt, threshholds = threshholds, perms = allperms, liks = alllik)
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

            ##garbage

##            browser()

##            browser()
            ##garbage
            if(class(fit$fit)=="lme")
            {
                frm = fit$fit$data
            } else if(class(fit$fit)=="lmerMod")
            {
                frm = attr(fit$fit, "frame")
                rands = rands[[1]]
            }


            y.mu = fit$y.transformed - epsilon - rands[frm$Dam.ID, "(Intercept)"]

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
    modelParser.ident              = lm.parsing$getHeavySingleParser(include.resid = F,
                                                                     mainEffectContrasts  = micro.analysis$getMainEffectContrasts(),
                                                                     medianAdjust.p.value = prop$mnp$medianAdjust.p.values)

    ##garbage
    modelParser = lm.parsing$getHeavyDoubleParser(include.resid = F,
                                                  mainEffectContrasts = c())
                                                  
    
    ##    reducedModelString = formulaWrapper$removeEffect(variableOfInterest, covariateModelString)
    reducedModelString = formulaWrapper$removeEffectAndInteractions(variableOfInterest, covariateModelString)$modified.string

##    reducedModelString = covariateModelString
    print(paste0("original model: ", covariateModelString))
    print(paste0("reduced model: ", reducedModelString))
    print("getting resids from reduced")
    
    transformParams          = fit.model.bc$getDefaultTransformParams()
    ##TODO remove
    ##transformParams$lambdasToTry = -3
    if(prop$mnp$tiny)
    {
        transformParams$lambdasToTry = 1
    }

    prs = micro.analysis$getBestParArgs(100, 3)
##    prs = micro.analysis$getLocalParArgs()
    parser      = micro.analysis$getPermResidParser()
    identParser =
##    browser()
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

    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = T, prefer.lme = T)
    modelParser = lm.parsing$getHeavyDoubleParser(include.resid = F, mainEffectContrasts = c("Diet"),includeSingle = T)
    
    ident.full = (surrogatcalc$runAnalysis(svFunc                   = svFunc,
                                          exp.mat                  = exp.mat,
                                          exp.mat.control          = exp.mat.control,
                                          cov.data                 = cov.data,
                                          covariateModelString     = covariateModelString,
                                          nullModelString          = formulaWrapper$removeEffect(variableOfInterest, covariateModelString),
                                          modelParser              = modelParser,
                                          strategy                 = strategy,
                                          transformParams          = transformParams,
                                          parallelArgs             = micro.analysis$getBestParArgs(100, 3)))

    ## strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = T)
    ## modelParser = lm.parsing$getHeavySingleParser(include.resid = F, mainEffectContrasts = c())

    ## Z = system.time(surrogatcalc$runAnalysis(svFunc                   = svFunc,
    ##                                       exp.mat                  = exp.mat,
    ##                                       exp.mat.control          = exp.mat.control,
    ##                                       cov.data                 = cov.data,
    ##                                       covariateModelString     = covariateModelString,
    ##                                       modelParser              = modelParser,
    ##                                       strategy                 = strategy,
    ##                                       transformParams          = transformParams,
    ##                                       parallelArgs             = micro.analysis$getBestParArgs(100, 3)))

    
##    browser()
##    summaryRprof("prof.out")
    perv = ident.full$results$dfalt$per.variable
##    browser()    
    
    ## pdf("dietplot.pdf")
    ## hist(perv[variable=="Diet"]$anova.p.value)
    ## dev.off()

    ## pdf("strainplot.pdf")
    ## hist(perv[variable=="Strain"]$anova.p.value)
    ## dev.off()


    
    ## pdf("dietplot2.pdf")
    ## hist(multipleTesting$median.correct.2(perv[variable=="Diet"]$anova.p.value))
    ## dev.off()

    ## pdf("strainplot2.pdf")
    ## hist(multipleTesting$median.correct.2(perv[variable=="Strain"]$anova.p.value))
    ## dev.off()


    
    ## browser()
    save(list=ls(), file = outm("identperm"))
    ##added
    modelParser = lm.parsing$getHeavyDoubleParser(include.resid = F, mainEffectContrasts = c())
##    debug(modelParser$collate)
##    debug(modelParser$parse)
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = T, prefer.lme = T) 
    getPermAccum <- function()
    {
        parallelArgs      = micro.analysis$getBestParArgs(13, mem.gb = 2, timeLimit.hours = 23)
        parallelArgs$func = surrogatcalc$runAnalysisHelper
        parallelArgs$sharedVariables = list(sv.info = ident.full$sv.info,
##                                            nullModelString = NULL,
                       ##                     nullModelString =  formulaWrapper$removeEffect(variableOfInterest, sv.info$covariateModelString),
                                            modelParser     = modelParser,
                                            transformParams = NULL, ## apply no transform at all to the y_pi generated by fitting the reduced model
                                            strategy        = strategy,
                                            parallelArgs    = micro.analysis$getLocalParArgs()) ##dont parallelize within permutation, or too many jobs
        				       
        accum = parallel$getAccum(parallelArgs)
        return(accum)
    }


##    browser()
    accum = getPermAccum()
    shuffles     = util$generateShuffles(trueLabels = 1:nrow(exp.mat), numShuffles = prop$mnp$SSVA.numperm, identityFirst = T)
    shuffles.dam = util$generateShuffles(trueLabels = rownames(toperm$rands), numShuffles = prop$mnp$SSVA.numperm, identityFirst = T)
    numPerm  = ncol(shuffles)
    
    for(shuffleIndex in 1:numPerm)
    {
        print(paste0("adding shuffle:", shuffleIndex))
        shuffle                = shuffles[,shuffleIndex]

##      arser  browser()
        
        shuffle.dam            = shuffles.dam[,shuffleIndex]
        names(shuffle.dam)     = rownames(toperm$rands)
        
        y.shuffled             = toperm$y.mu +
                                 toperm$epsilon[shuffle,] +
                                 toperm$rands[shuffle.dam[cov.data$Dam.ID],]

        rownames(y.shuffled)   = rownames(cov.data)
        accum$addCall(funcArgs = list(y.mat = y.shuffled))
       ## accum$addCall(funcArgs = list(exp.mat = y.shuffled))
    }

    ##    browser()
##    debug(accum$runAll)
    outs = accum$runAll()
##    browser()
    perms = micro.analysis$collatePermsFreed(accum, outs, numPerm, variableOfInterest)
   
    return(list(perms=perms,ident.full = ident.full))

}


##TODO: move to permutation testing?
micro.analysis$collatePermsFreed <- function(accum, outs, numPerm, variableOfInterest)
{
##    browser()
##    browser()
    print("collating perms")
    iter = accum$getOutputIterator(outs)
    badind = c()
    i= 0
    perm.pvalues = NULL

##    browser()
##    debug(iter$nextItem)
    
    while(iter$hasNext())
    {
        i = i+1
        out = iter$nextItem()

         ##GARBAGE
        ##out = try(out$results)

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

        ##garbage
        ## pvalues          = try(out$per.variable[variable == varOfInterest]$anova.p.value)
        ## namez            = try(out$per.variable[variable == varOfInterest]$Probe.Set.ID)


##        browser()
##        dfs = try(out$dfsimple)
        ## out$dfsimple$pvalues          = try(unlist(out$dfsimple$p.value))
        ## out$dfsimple$namez            = try(unlist(out$dfsimple$Probe.Set.ID))
        ## out$dfsimple$lik.rat          = try(unlist(out$dfsimple$lik.rat))
        
        if(class(out$dfsimple$pvalues)=="try-error")
        {
            badind = c(badind, i)
            print(paste0("permutation failed!!!", i))
            next
        }

        if(is.null(perm.pvalues))
        {
            lik.rats     = matrix(NA, nrow = numPerm, ncol = nrow(out$dfsimple))
            perm.pvalues = matrix(NA, nrow = numPerm, ncol = nrow(out$dfsimple))
            colnames(perm.pvalues) = out$dfsimple$namez
            colnames(lik.rats) = out$dfsimple$namez
        }
        perm.pvalues[i,] = out$dfsimple$p.value
        lik.rats[i,]     = out$dfsimple$lik.rat
        
    }

    if(length(badind)>0)
    {
        print(paste("bads:", badind))
        perm.pvalues = perm.pvalues[,-badind]
    }
##    perms = perm.pvalues[1:nrow(perm.pvalues),]
##    ident = perm.pvalues[1,]
    alphas = .05


    save(file = outm("rawperms", variableOfInterest), list = ls())
    
    threshholds = multipleTesting$evaluate.perms(perm.pvalues,
##                                                 ident,
                                                 alpha = alphas,
                                                 direction = 1)
    
    threshholds = data.table(
        
        variable       = variableOfInterest,
        alpha          = alphas,
        threshhold     = threshholds$empirical.global.threshhold,
        threshhold.gev = threshholds$empirical.global.gev.threshhold)

    return(list(threshholds = threshholds, perms = perm.pvalues, lik.rats = lik.rats))
}


