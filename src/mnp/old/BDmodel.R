mnp.med$get.BD.inputs <- function()
{
    raw.data = loadAllData$createAllInputs()
    
    ##only for debugging
    ## observeInfo = data.frame(variable = c("precision.Strain", "precision.epsilon", "intercept", "beta", 
    ##                                       "f.Strain", "f.Diet", "f.Strain.Diet", "f.Batch", "f.Pipeline" ),
    ## factorforindexconversion = c(NA, NA, NA, NA,
    ##                              "Strain", "Diet", NA, "Batch", "Pipeline"),
    ## keep = T)

    observeInfo = data.frame(variable = c("f.Strain[1,2]",
                                          "f.Strain[2,2]",
                                          "f.Strain.Diet",
                                          "beta"),
                             factorforindexconversion = c("f.Strain", NA, NA, NA),
                             keep = T)
        
    observeInfo$variable = as.character(observeInfo$variable)
##    observeInfo$factorForIndexConversion = as.character(observeInfo$factorForIndexConversion)
    
    froot = outm("mediation", "BD")
    dir.create(froot, showWarnings = F, recursive = T)
    

    sv.info = mnp.med$get.sv.corrected.genes(raw.data, T)
    cov.data = sv.info$cov.data

    ##TODO extract this to be its own integer index building module
    n.Strain = length(levels(cov.data$Strain))
    n.Diet   = length(levels(cov.data$Diet))
    nonzero = as.integer(cov.data$Strain)!=1 & as.integer(cov.data$Diet)!=1
    cov.data$Strain.Diet = 0
    cov.data$Strain.Diet[nonzero] = (as.integer(cov.data$Diet[nonzero])-2)*(n.Strain -1) + (as.integer(cov.data$Strain[nonzero]) - 1)
    cov.data$Strain.Diet = cov.data$Strain.Diet + 1

    lmer.modelString = formulaWrapper$appendEffect("mediator", sv.info$covariateModelString)$modified.string
##   lmer.modelString = sv.info$covariateModelString 


    s.info = formulaWrapper$parseCovariateString(lmer.modelString)
    simpleFixedEffects = unlist(lapply(FUN=paste, X= s.info$fixef, collapse = "."))
    simpleFixedEffects = setdiff(simpleFixedEffects, "1")
    simpleFixedEffects = setdiff(simpleFixedEffects, "mediator")

     randomEffects = unlist(lapply(FUN = "[[", X = s.info$ranef, "group"))

    additionalData = list()
##    additionalData = list(n.Strain.Diet = 4)

    return(list(froot            = froot,
                lmer.modelString = lmer.modelString,
                cov.data         = cov.data,
                outcome.mat      = sv.info$exp.mat,
                mediationFunc    = mnp.med$BD.mediationFunc,
                observeInfo      = observeInfo,
                ##                jags.modelname = "./mnp/BD.bug",
                ##                jags.modelname = textConnection(modelstring),
                jags.modelname   = mnp.med$get.BD.models(setdiff(simpleFixedEffects, "mediator"), randomEffects)$postModel,
                priorModel       = mnp.med$get.BD.models(simpleFixedEffects, randomEffects)$priorModel,
                colsToIndex      = setdiff(c(simpleFixedEffects, randomEffects), "mediator"),
                ##gibbsBatchSize = 45,
                gibbsBatchSize   = 45,
                additionalData   = additionalData))
}

mnp.med$BD.mediationFunc <- function(mcobj, dataForJags)
{

    dfs = list()
    as = list()
    bs = list()
    cs = list()
    b =  mcobj[,"beta"]
    for(i in 1:length(levels(dataForJags$Diet)))
    {
        interactionOffset = paste0("f.Strain.Diet[2,",i,"]")
        a = (mcobj[,"f.Strain[2,2]"] + mcobj[,interactionOffset])

        interactionOffset = paste0("f.Strain.Diet[1,",i,"]")
        c.prime = (mcobj[,"f.Strain[1,2]"] + mcobj[,interactionOffset])

        df = mnp.med$samplesToMediation(a, b, c.prime)
        df$moderators = paste0("Diet=", levels(dataForJags$Diet)[i])

        
        dfs = util$appendToList(dfs, df)
        as[[i]] = a
        bs[[i]] = b
        cs[[i]] = c.prime
    }
    
    df = mnp.med$samplesToMediation(mcmc(do.call(c, as)),
                                    mcmc(do.call(c, bs)),
                                    mcmc(do.call(c, cs)))

    df$moderators = paste0("Diet=Ave")

    dfs = util$appendToList(dfs, df)
    dfs = rbindlist(dfs)

    return(dfs)
}


mnp.med$get.BD.basemodel <- function(simpleMainEffects, randomEffects)
{
    p = paste0
    modelString ="model {

    ##TODO pass this in
    bound <- 5

    sd.multiple <-1
    ##################################
    ##  Priors for intercept
    ##############################
    intercept[1] ~ dnorm(0, pow(bound*sd.multiple, -2))
    intercept[2] ~ dnorm(0, pow(bound*sd.multiple, -2))

"
    for(simpleMainEffect in simpleMainEffects)
    {
        modelString = p(modelString, mnp.med$getEffectDistribution(simpleMainEffect))
    }

    for(randomEffect in randomEffects)
    {
        modelString = p(modelString, mnp.med$getRandomEffectDistribution(randomEffect))
    }

   
    modelString = p(modelString,
                    
    "
    ############################################################
    ## Beta- the effect of the mediator m on y.
    ## Assumed to be independent of all other covariates.
    ############################################################
    beta ~ dnorm(0, pow(bound*sd.multiple, -2))
    ##beta <- 0

    ############################
    ## Sigmas for y and mediator.
    ############################
    for(i in 1:2)
    {
        std.epsilon[i] ~ dunif(0,25)
        precision.epsilon[i] <- pow(std.epsilon[i],-2)
    }

    ###############################################
    # Likelihood functions, indexed by s, the sample
    ###############################################
    for(s in 1:n)
    {
       ##compute the means first
       for(i in 1:2)
       {

       u[i,s] <- intercept[i]                    + \n")

    me = lapply(FUN = mnp.med$getLikContribution, simpleMainEffects, type = "f")
    me = do.call(paste, c(me, sep = "        +\n"))
    modelString = paste(modelString, me, "         +\n")
    
    re = lapply(FUN = mnp.med$getLikContribution, randomEffects, type = "ranef")
    re = do.call(paste, c(re, sep = "        +\n"))
    modelString = paste(modelString, re, "

       }","

       ##The mediator likelihood.
       MEDIATOR_LIKELIHOOD

       ##Data likelihood
       OUTCOME_LIKELIHOOD
    }
}")

    return(modelString)
}

mnp.med$get.BD.models <- function(fixedMainEffects, randomEffects)
{
    modelString =mnp.med$get.BD.basemodel(fixedMainEffects, randomEffects)
    
    postModel = mnp.med$getModelString(modelString,
                                       "mediator[s] ~ dnorm(u[2,s], precision.epsilon[2])",
                                       "y[s] ~ dnorm(u[1,s] + beta*mediator[s], precision.epsilon[1])")

    priorModel = mnp.med$getModelString(modelString, "", "")

    out = list(postModel = postModel, priorModel = priorModel)
    return(out)
}



