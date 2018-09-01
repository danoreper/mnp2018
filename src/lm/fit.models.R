source("./lm/formulaWrapper.R")
source("./lm/lm.parsing.R")
source("./utils.R")

##TODO: rename environment
fit.modelg = new.env(hash=T)

fit.modelg$fit.single <- function(y,
                                  cov.data,
                                  covariateModelString,
                                  checkAnova     = T,
                                  
                                  checkNormality = F,
                                  gurka          = F,
                                  strategy       = NULL)
{
    if(is.null(strategy))
    {
        strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = T)
    }
    
    cov.data[["y"]]=y


    ##undebug(fit.modelg$.fit)
    fit = fit.modelg$.fit(pheno = "y",
                          covariateModelString = covariateModelString,
                          cov.data = cov.data,
                          checkAnova = checkAnova,
                          strategy = strategy)

##    browser()
    
    if(checkNormality)
    {
        noiseVector = lm.parsing$getNoiseVector(fit$fit, gurka=gurka)
        normality.pval = try(shapiro.test(noiseVector)$p.value)
    }

    pvals = NULL
    anovaWrapper = NULL

    fit$normality.pval = normality.pval
    ## output = list(fit = fit$fit,
    ##               anovaWrapper = fit$anovaWrapper,
    ##               ##       pvals = pvals,
    ##               normality.pval = normality.pval)
    return(fit)
}

fit.modelg$getDefaultModelStrategy <- function(anovaComparison,
                                               prefer.lme)
{
    force(anovaComparison)
    force(prefer.lme)
    
    strategyFunc = function(covariateModelString, checkAnova)
    {
        covInfo = formulaWrapper$parseCovariateString(covariateModelString)
        if(length(covInfo$ranef)==0)
        {
            funcName = "lm"
        } else {
            if(!checkAnova)
            {
                funcName = "lme4::lmer"  # "lme" #
            } else {
                if(length(covInfo$ranef)>1 |!prefer.lme)
                {
                    funcName = "lmerTest::lmer"
                } else {
                    funcName = "lme"
                }
            }
        }
        func = eval(parse(text=funcName))
        
        opts = list()
        opts$checkAnova = checkAnova
        ## if(!is.null(nullModelString))
        ## {
        ##     opts$checkAnova = F
        ## }

        ##TODO: refactor to have lm specific opts
        opts$method = "REML"
        opts$REML   = T
        opts$na.action = na.omit
        if(anovaComparison)## && funcName == "lme")
        {
            opts$method = "ML"
            opts$REML = F
        }
        return(list(func = func, funcName = funcName, opts=opts, anovaComparison = anovaComparison))
    }
    return(strategyFunc)
}



##This may seem like pointless indirection, but its actually critical because of the way lme works and its interaction with environment;
##reassigns variable name to dataset.
fit.modelg$.fit <- function(pheno,
                            cov.data,
                            covariateModelString,
                            checkAnova,
                            strategy)
    
{
    strategyParams = strategy(covariateModelString = covariateModelString, checkAnova = checkAnova)
    ##don't delete this-- its actually critical!!!
    dataSet = cov.data
    func                 = strategyParams$func
    funcName             = strategyParams$funcName
    opts                 = strategyParams$opts
    ##don't delete this!!!

    out = try(fit.modelg$.fithelper(pheno                = pheno,
                                dataSet              = dataSet,
                                covariateModelString = covariateModelString,
                                checkAnova           = checkAnova,
                                func                 = func,
                                funcName             = funcName,
                                opts                 = opts))

    ## TODO consider integrating recovery
    ## if(class(out)=="try-error" & funcName == "lme")
    ## {
    ##     print("recovering from failed lme by trying lme4 instead")
    ##     strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = strategyParams$anovaComparison, prefer.lme=F)

    ##     strategyParams = strategy(covariateModelString = covariateModelString, checkAnova = checkAnova)
    ##     ##don't delete this-- its actually critical!!!
    ##     dataSet = cov.data
    ##     func                 = strategyParams$func
    ##     funcName             = strategyParams$funcName
    ##     opts                 = strategyParams$opts
    ##     ##don't delete this!!!
    ##     out = try(fit.modelg$.fithelper(pheno                = pheno,
    ##                                     dataSet              = dataSet,
    ##                                     covariateModelString = covariateModelString,
    ##                                     checkAnova           = checkAnova,
    ##                                     func                 = func,
    ##                                     funcName             = funcName,
    ##                                     opts                 = opts))
    ## }
    return(out)
}
fit.modelg$.fithelper <- function(pheno,
                                  dataSet,
                                  covariateModelString,
                                  checkAnova = checkAnova,
                                  func,
                                  funcName,
                                  opts)

{
    ##lmer spams console with partial matches.
    ##To avoid this, temporarily turning it off, before restoring state
    ##just before fithelper function concludes. Not thread safe, so silencing won't
    ##work properly if running on multiple cores on the same machine; not sure how to fix.
    pmatch = options("warnPartialMatchDollar")[[1]]
    options(warnPartialMatchDollar = F) 
    
    if(checkAnova)
    {
##        browser()
    }
    covInfo = formulaWrapper$parseCovariateString(covariateModelString)
    
    if(funcName == "lm")
    {
        
        if(length(covInfo$ranef)>0)
        {
            stop("lm assigned, but at least 1 random effect")
        }

        astring = paste0(pheno, covInfo$modified.string)
        ##        fit = func(as.formula(paste0(pheno, paste0("~ ",paste(covInfo$fixef, collapse=" + ")))), data=dataSet)
        fit = func(as.formula(astring), data=dataSet)
    } else if (funcName == "lmerTest::lmer" || funcName == "lme4::lmer")
    {
        formla = as.formula(paste0(pheno, covariateModelString))

        fit = func(formla, data=dataSet, REML =opts$REML)
    } else if (funcName == "lme") {
        ## we are doing this rigamarole with substitute so that the "fixed" and "random"
        ## call fields of the resulting lme object make sense in the environment of dataSet,
        ## and not just in this functions scope.
        ## NB anova comparing 2 lme models needs the fixed and random fields.
        ## ; Put another way, LME does something wrong with environmental variables such that we need to do
        ## a lot of work to call it properly within a function.

        ##TODO correct this such that it uses the new formula 
        fixef  = covInfo$fixef
        ranef  = covInfo$ranef

        fixedd = paste(pheno, "~", formulaWrapper$mergeTerms(covInfo$fixef))
        ranedd = list()
        for(aranef in covInfo$ranef)
        {
            ranterms.left  = formulaWrapper$mergeTerms(aranef$components)
            ranterms.right = formulaWrapper$mergeTerms(aranef$group)
            ranterms = paste(ranterms.left, "|", ranterms.right)
            ranedd = util$appendToList(ranedd, ranterms)
        }
        ranedd = paste("~", do.call(paste, c(ranedd, sep=" + ")))
        
        fit = (eval(substitute(
            func (data   = dataSet,
                  fixed  = as.formula(fixedd),
                  random = as.formula(ranedd),
                  method = opts$method,
                  na.action = opts$na.action))))
        
    } else {
        stop("unimplemented")
    }
    
    
    ##TODO: why on earth is this necessary?? data is supposed to be stored
    ##fit$data = dataSet

    ##check anova MUST happen in the scope of this function rather than higher up because of some sort of internal bug in
    ##lmerTest::anova internal environment

    anovaWrapper = NULL
    if(checkAnova)
    {
        ##TODO move into parsing
        anovaWrapper = try(lm.parsing$getAnovaWrapper(fit))
        
        pvals = try(anovaWrapper$an[,  anovaWrapper$pvalueCol])
       ## print(pvals)
        if(class(pvals)=="try-error"|| !anovaWrapper$pvalueCol %in% colnames(anovaWrapper$an))
        {
            print("anova failure")
            return(NULL)
        }
    }
    out = list(fit = fit, anovaWrapper = anovaWrapper)

    try(options(warnPartialMatchDollar = pmatch))
    return(out)
}

