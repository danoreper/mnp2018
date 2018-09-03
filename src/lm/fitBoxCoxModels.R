source("./util/loadParams.R")
source("./lm/formulaWrapper.R")
source("./lm/lm.parsing.R")
source("./utils.R")
source("./lm/fit.models.R")
source("./parallel/accumulator.R")
DEBUG = F

fit.model.bc = new.env(hash=T)


fit.model.bc$getIdentityTransformParams <- function()
{
    out = list()
    out$lambdasToTry    = NULL
    out$lambdaPerY      = NULL
    out$normalizeBeforeTransform = F
    out$normalizeAfterTransform  = F
    out$extremelb      = -3
    out$extremeub      = 3
    return(out)
}

fit.model.bc$getDefaultTransformParams <- function()
{
    out = list()
    out$lambdasToTry = c(-3, -2, -1, -.5, 0, .5, 1,2,3)
    out$lambdaPerY      = NULL
    out$normalizeBeforeTransform = T
    out$normalizeAfterTransform  = T
    out$extremelb      = -3
    out$extremeub      = 3
    return(out)
}

fit.model.bc$getDefaultParser <- function(covariateModelString,
                                          nullModelString = NULL)
{
    if(is.null(nullModelString))
    {
        parser = lm.parsing$getFullSingleParser(T)
    } else {
        parser = lm.parsing$getFullDoubleParser(T)
    }
    return(parser)
}


fit.model.bc$fit <- function(y.mat,
                             cov.data,
                             covariateModelString,
                             nullModelString  = NULL,
                             modelParser      = fit.model.bc$getDefaultParser(covariateModelString, nullModelString),
                             transformParams  = fit.model.bc$getDefaultTransformParams(),
                             
                             checkAnova       = T,
                             strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString),
                                                                           prefer.lme = T),
                             parallelArgs = parallel$getDefaultLocalArgs())
{
    
    if(is.vector(y.mat))
    {
        y.mat = as.matrix(ncol = 1, y.mat)
    }
    


    notlocal = parallel$useCluster(parallelArgs)
    parallelArgs$func = fit.model.bc$fit.single
    if(notlocal)
    {
        parallelArgs$saveProp = T
        parallelArgs$filesToSource = c(parallelArgs$filesToSource, "./lm/fitBoxCoxModels.R")
    }
    parallelArgs$sharedVariables = list(cov.data = cov.data,
                                        covariateModelString = covariateModelString,
                                        nullModelString      = nullModelString,
                                        checkAnova     = checkAnova,
                                        strategy       = strategy,
                                        modelParser    = modelParser)

    accum = parallel$getAccum(parallelArgs)

    for(j in 1:ncol(y.mat))
    {
        y = y.mat[,j]
        if(!is.null(transformParams))
        {
            if(!is.null(transformParams$lambdaPerY))
            {
                lambdasToTry = transformParams$lambdaPerY[j]
            } else {
                lambdasToTry = transformParams$lambdasToTry
            }
            transformParams$lambdasToTry = lambdasToTry
        }
        accum$addCall(list(y=y, transformParams = transformParams))
    }
    
    outs  = accum$runAll()
    outs  = accum$getAllOutputs(outs, removeFailing = F)

    ##name every output... this seems like it should go elsewhere perhaps
    cnamez = colnames(y.mat)
    if(is.null(cnamez))
    {
        cnamez = paste0("phen_",1:ncol(y.mat))
    }
    names(outs) = cnamez
    outs = modelParser$collate(outs, accum)

    return(outs)
}

fit.model.bc$fit.single <- function(y,
                                    cov.data,
                                    covariateModelString,
                                    nullModelString         = NULL,
                                    modelParser             = fit.model.bc$getDefaultParser(covariateModelString, nullModelString),
                                    transformParams          = fit.model.bc$getDefaultTransformParams(),
                                    checkAnova               = T,
                                    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString),
                                                                                  prefer.lme = T))
{
    loopOverLambdas <- function(checkAnovaInLoop, lambdasToTry)
    {
         bestFit     = NULL
         best.pval   = 0
         best.lambda = 1
         bestY       = NULL
         bestAnova   = NULL
    
         if(DEBUG) {print("*********************")}
         if(DEBUG) {print(paste0("Looping over lambdas, using model:", covariateModelString))}
         if(DEBUG) { print(paste0("checkAnova=", checkAnovaInLoop))}
         for(lambda in lambdasToTry)
         {
             y.transformed = transform(y=y,lambda=lambda)

             ## if(DEBUG){ print(paste("y.transformed:", paste0(sprintf("%.3f", head(y.transformed, 5)), collapse = ",")))}
             afit = callfit(y.transformed$y, covariateModelString, checkAnova = checkAnovaInLoop)
             if(!is.null(afit) && class(afit)!="try-error")
             {
                 pval = afit$normality.pval
                 if(DEBUG) { print(paste0("fit lambdaToTry=", lambda, ", shapiro.p = ", pval))}
                 if(pval>best.pval)
                 {
                     bestFit     = afit
                     best.pval   = pval
                     best.lambda = lambda
                     bestY       = y.transformed
                     bestAnova   = afit$anovaWrapper 
                 }
             } else { if(DEBUG){print(paste0("fit failed, skipping lambdaToTry=", lambda))}}
         }
         
         out = (list(fit = bestFit$fit,
                    best.pval = best.pval,
                    lambda = best.lambda,
                    lambda2       = bestY$lambda2,
                    y.transformed = bestY$y,
                    center1       = bestY$center1,
                    scale1        = bestY$scale1,
                    center2       = bestY$center2,
                    scale2        = bestY$scale2,
                    anovaWrapper = bestAnova))
         return(out)
    }

    
    transform <- function(y, lambda)
    {
        y.transformed = fit.model.bc$transform.by.lambda(y = y,
                                                         lambda = lambda,
                                                         normalizeBeforeTransform = normalizeBeforeTransform,
                                                         normalizeAfterTransform = normalizeAfterTransform)
        return(y.transformed)
    }
    
    callfit <- function(y.transformed, modelString, checkAnova, checkNormality = T)
    {
        afit = try(fit.modelg$fit.single(y = y.transformed,
                                        cov.data             = cov.data,
                                        covariateModelString = modelString,
                                        checkAnova           = checkAnova,
                                        checkNormality       = checkNormality,
                                        strategy             = strategy))

        return(afit)
    }

    ##NA transform params means don't transform the data at all, just fit the model as is. If it fails, don't parse, return failure
    if(is.null(transformParams))
    {
        if(DEBUG){print("boxcox: no transform, so running on raw data, without looping")}
           
        bestFit = callfit(y, covariateModelString, checkAnova = checkAnova, checkNormality = F)

        if(class(bestFit)!="try-error")
        {
            out = list(fit = bestFit$fit, best.pval = NA, lambda = NA, y.transformed = y, anovaWrapper = bestFit$anovaWrapper)
            if(!is.null(nullModelString))
            {
                out$fit.null = callfit(out$y.transformed, nullModelString, checkAnova = F)$fit
            }
            parsed   = modelParser$parse(out)
            return(parsed)
        } else {
            return(bestFit)
        }
    }

    tp = transformParams
    normalizeBeforeTransform = tp$normalizeBeforeTransform
    normalizeAfterTransform  = tp$normalizeAfterTransform
    extremelb                = tp$extremelb
    extremeub                = tp$extremeub                                   
    lambdasToTry             = tp$lambdasToTry

    out = loopOverLambdas(checkAnovaInLoop=F, lambdasToTry = lambdasToTry)
    orig.lambda = out$lambda
    if(DEBUG) { print(paste0("best lambda without checking anova: ", orig.lambda))}
    
    if(out$lambda <= extremelb || out$lambda >= extremeub || is.null(out$fit))
    {
        if(DEBUG) { print(paste0(orig.lambda, " failed or is beyond the bounds", "(",extremelb,",", extremeub, "), use inverse normal instead"))}
        out = loopOverLambdas(checkAnovaInLoop = checkAnova, lambdasToTry = "inverse_normal")
    } else {
        if(checkAnova)
        {
            if(DEBUG) {print("now checking if selected lambda allows anova to be taken")}
            out = loopOverLambdas(checkAnovaInLoop=T, lambdasToTry = orig.lambda)
            if(is.null(out$fit)||is.null(out$anovaWrapper))
            {
                if(DEBUG)
                {
                    print(paste0("anova failed on selected lambda=", orig.lambda,
                                 ", so loop over lambdas again (with the exception of the failing lambda, and extrema),", 
                                 "but this time check anova along the way"))
                }                
                lambdasToTry = setdiff(lambdasToTry, c(orig.lambda, extremelb, extremeub))
                out = loopOverLambdas(checkAnovaInLoop=T, lambdasToTry)
            }
        }
    }
    if(DEBUG) {print(paste0("Best lambda, in the end: ", out$lambda))}
    if(!is.null(nullModelString))
    {
        if(DEBUG) {print("after transforming Y by the ALT derived model box cox parameter, do anova againt the null model")}
        out$fit.null = callfit(out$y.transformed$y, nullModelString, checkAnova = F)$fit
    }

    
    parsed   = modelParser$parse(out)
    
    return(parsed)
}
    
fit.model.bc$computeResids <-  function(y.mat,
                                        cov.data,
                                        covariateModelString,
                                        transformParams = fit.model.bc$getDefaultTransformParams(),
                                        strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F,
                                                                                      prefer.lme = T),
                                        parallelArgs = parallel$getDefaultLocalArgs())  
{
    parser = new.env(hash=T) 

    parser$parse <- function(fit)
    {
        if(is.null(fit$fit))
        {
            return(rep(NA, nrow(y.mat)))
        } else {
            return(resid(fit$fit))
        }
    }

    parser$collate <- function(outs, accum)
    {
        out     = do.call(cbind, outs)
        rownames(out) = rownames(y.mat)
        colnames(out) = colnames(y.mat)
        return(out)
    }
    
    out = fit.model.bc$fit(y.mat = y.mat,
                           cov.data = cov.data,
                           covariateModelString = covariateModelString,
                           modelParser          = parser,
                           transformParams = transformParams,
                           checkAnova     = F,
                           strategy = strategy,
                           parallelArgs = parallelArgs)
    
    return(out)
}


fit.model.bc$residOutCovariates <- function(y.mat,
                                            cov.data,
                                            covariateModelString,
                                            nullModelString = NULL,
                                            residualizeOutCovariates,
                                            transformParams = fit.model.bc$getDefaultTransformParams(),
                                            strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme=T), 
                                            parallelArgs = parallel$getDefaultLocalArgs())
{
    nuis.covariateModelString =
        formulaWrapper$onlyKeepEffectAndInteractions(residualizeOutCovariates,
                                                     covariateModelString)$modified.string
    
    exp.mat = fit.model.bc$computeResids(y.mat          = y.mat,
                                         cov.data       = cov.data, 
                                         covariateModelString = nuis.covariateModelString,
                                         transformParams= transformParams,
                                         strategy = strategy,
                                         parallelArgs = parallelArgs)

    ##After we're done residulizing out terms, remove them from the covariate model string too.
    covariateModelString  = formulaWrapper$removeEffectAndInteractions(residualizeOutCovariates,
                                                                       covariateModelString)$modified.string
    
    if(!is.null(nullModelString))
    {
        nullModelString =
            formulaWrapper$removeEffectAndInteractions(residualizeOutCovariates,
                                                       nullModelString)$modified.string
    }
    return(list(residualizedData = exp.mat,
                residualizedModelString = covariateModelString,
                residualizedNullString = nullModelString))
}


##TODO move to separate file
fit.model.bc$transform.by.lambda <- function(y,
                                             lambda,
                                             normalizeBeforeTransform,
                                             normalizeAfterTransform)
{
    output =list()
    if(normalizeBeforeTransform)
    {
        y = scale(y,center = T)
        output$center1 = attr(y, "scaled:center")
        output$scale1  = attr(y, "scaled:scale")
    }

    if(lambda == "inverse_normal")
    {
        output$lambda2 = NA
        y = fit.model.bc$inverseNormal(y)
        
    } else {
        y = fit.model.bc$boxcox.for.lambda(y, lambda)
        output$lambda2 = y$lambda2
        y = y$y
    }

    if(normalizeAfterTransform)
    {
        y = scale(y, center = T)
        output$center2 = attr(y, "scaled:center")
        output$scale2  = attr(y, "scaled:scale")
    }

    y = y[1:length(y)]
    output$y = y
    if(lambda<0)
    {
        output$y = -output$y
    }
    output$lambda1 = lambda
    return(output)
}

fit.model.bc$inverseNormal <- function(y)
{
    allRanks = rank(y, na.last="keep")
    maxRank = max(allRanks, na.rm=TRUE);
    return(qnorm(allRanks/(maxRank+1)));
}

fit.model.bc$boxcox.for.lambda <- function(y, lambda1)
{
    offset = .1  #we cant have lambda2==-y_i, so we will make sure lambda2>=-y_i+offset
    lambda2 = max(0, max(-y, na.rm=T)+offset) #dont bother having a lambda2 unless there is at least one zero value

    if(lambda1 == 0)
    {
        out = log(y + lambda2)
    }
    else
    {
        out = (((y + lambda2)^lambda1)-1)/lambda1
    }

    return(list(y = out, lambda2 = lambda2))
}

fit.model.bc$invert <- function(y.transformed, perm=T)
{
    if(y.transformed$lambda1 == "inverse_normal")
    {
        stop("cant invert inverse_normal")
    }
    
    y = y.transformed$y
    if(!is.null(y.transformed$scale2))
    {
        y = y*y.transformed$scale2 + y.transformed$center2
    }

    if(y.transformed$lambda1!="inverse_normal")
    {
        lambda1 = y.transformed$lambda1
        lambda2 = y.transformed$lambda2

        
        y = y*lambda1+1
        if(perm && lambda1>1 && min(y)<=0)
        {
            print(min(y))
            print(mean(y))
            y = y - min(y)+.1
        }
        
        y = y^(1/lambda1)-lambda2
    }

    if(!is.null(y.transformed$scale1))
    {
        y = y*y.transformed$scale1 + y.transformed$center1
    }

    if(any(is.na(y)))
    {
        browser()
    }
    return(y)
}
