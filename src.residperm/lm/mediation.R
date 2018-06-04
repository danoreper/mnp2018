source("./lm/fitBoxCoxModels.R")
source("./utils.R")
lm.mediation = new.env(hash=T)

lm.mediation$simpleTest = function(outcome.vec,
                                   mediator.vec,
                                   cov.data,
                                   mediatorModelString,
                                   mediatedFactor,
                                   transformParams = fit.model.bc$getDefaultTransformParams(),
                                   strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = T,
                                                                                 prefer.lme = T))
{

##gsub(m.alt, pattern = paste0(mediatedFactor, replacement
    
    cov.data$mediator = mediator.vec


    m.null.mediator   = formulaWrapper$removeEffectAndInteractions(mediatedFactor, mediatorModelString)$modified.string
    
    m.alt             = formulaWrapper$appendEffect("mediator", mediatorModelString)$modified.string
    m.null.complete   = formulaWrapper$removeEffectAndInteractions(mediatedFactor, m.alt)$modified.string
    m.null.partial    = formulaWrapper$removeEffectAndInteractions("mediator",     m.alt)$modified.string

    ##These don't change regardless of the model we are fitting
    baseargs = list(cov.data = cov.data, transformParams = transformParams, strategy = strategy)

    
    argz = list()

    argz[["complete.mediation"]] = list(y.mat = outcome.vec,  covariateModelString = m.alt,               nullModelString = m.null.complete)
    argz[["factor.on.mediator"]] = list(y.mat = mediator.vec, covariateModelString = mediatorModelString, nullModelString = m.null.mediator)
    argz[["partial.mediation"]]  = list(y.mat = outcome.vec,  covariateModelString = m.alt,               nullModelString = m.null.partial)
    

    dfs = list()
    for(i in 1:length(argz))
    {

        comparison.type = names(argz)[[i]]
        call.args = argz[[comparison.type]]

        afit.mediator = do.call(fit.model.bc$fit, c(baseargs, call.args))
        ## if(comparison.type == "complete.mediation")
        ## {
        ##     browser()
        ## }
        anova.p.value = afit.mediator$phen_1$aw$p.value
        df = data.frame(comparison.type = comparison.type, anova.p.value = anova.p.value)
        dfs = util$appendToList(dfs, df)
    }

    df = rbindlist(dfs)
    df = dcast(df, . ~ comparison.type, value.var = "anova.p.value")
    df[["."]] = NULL
    
    return(df)
}
