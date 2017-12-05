##TODO encode brain and female

source("./crowley/loadAllDataCrowley.R")
source("./crowley/micro/analysis.R")


mnp.med$get.crowley.sv.corrected.genes <- function(inp=loadDataCrowley$createAllInputs(),
                                                   fromFile = T)
{
    

    if(!fromFile)
    {
        residualizeOutCovariates = c()
        strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString), prefer.lme = F) ##For some reason, lme fails a lot with this particular model 

        print("getting surrogate variable info")
        pracma::tic()
        svinfo = surrogatcalc$generate.svinfo(svFunc = crowley.analysis$get.SV.func(),
                                              exp.mat = inp$exp.mat,
                                              exp.mat.control = inp$exp.mat.control,
                                              cov.data = inp$cov.data.full,
                                              covariateModelString = crowley.analysis$formCovariateFullModelString(),
                                              residualizeOutCovariates = residualizeOutCovariates,
                                              residualizeOutSV = T,
                                              strategy = strategy,
                                              accum  = crowley.analysis$getBestAccumulator(100, 3))))
                
        print(paste0("got surrogate variables:",pracma::toc()))
        gc()
        
        print("done getting resids")
        save(file=outm("sv.crowley.corrected.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outm("sv.crowley.corrected.RData"))
    }

    return(svinfo)
}


mnp.med$get.FGH.inputs <- function()
{

    inp.crowley = loadDataCrowley$createAllInputs()
    micro.data  = mnp.med$get.crowley.sv.corrected.genes(inp.crowley, F)$exp.mat
    
    observeInfo = data.frame(variable = c("intercept", "f.sex", "f.cross", "f.tissue",
                                          "f.sex.cross", "f.sex.tissue", "f.cross.tissue", "f.cross.direction",
                                          "f.sex.cross.tissue", "f.sex.cross.direction", "f.cross.tissue.direction",
                                          "f.sex.cross.tissue.direction", "beta" ),
                             factorForIndexConversion = NA,
                             keep = T)
    observeInfo$variable = as.character(observeInfo$variable)
    observeInfo$factorForIndexConversion = as.character(observeInfo$factorForIndexConversion)
    
    froot = outm("mediation", "BD")
    dir.create(froot, showWarnings = F, recursive = T)

    strain.results = fread(datm("strainResults.csv"))
    strain.results$Probe.Set.ID = as.character(strain.results$Probe.Set.ID)

    ##raw.data = loadAllData$createAllInputs()
##    sv.info = mediation$get.sv.corrected.genes(raw.data, T)

    FGH.lmer.modelString = paste0("~ 1 + sex + cross + tissue + ",
                                  "sex:cross + sex:tissue + cross:tissue + cross:direction + ",
                                  "sex:cross:tissue + sex:cross:direction +tissue:cross:direction +",
                                  "sex:cross:tissue:direction +", 
                                  "mediator + (1|mouse)")

    additionalData = list(n.cross.tissue = 1 + (length(levels(inp.crowley$cov.data$cross))-1) * (length(levels(inp.crowley$cov.data$tissue))-1))

    ##TODO put this into load all data crowley?
    cov.data     = inp.crowley$cov.data
    cov.data$sex = as.integer(cov.data$sex) - 1
    
    cov.data$cross.tissue = 0
    
    
    n.cross  = length(levels(cov.data$cross))
    n.tissue = length(levels(cov.data$tissue))
    nonzero = as.integer(cov.data$cross)!=1 & as.integer(cov.data$tissue)!=1
    cov.data$cross.tissue[nonzero] = (as.integer(cov.data$tissue[nonzero])-2)*(n.cross -1) + (as.integer(cov.data$cross[nonzero]) - 1)
    cov.data$cross.tissue = cov.data$cross.tissue + 1

    
    return(list(froot            = froot,
                allcandidates    = strain.results,
                lmer.modelString = FGH.lmer.modelString,
                cov.data         = cov.data,
                exp.mat          = micro.data,
                mediationFunc    = mnp.med$FGH.mediationFunc,
                observeInfo = observeInfo,
                ##                jags.modelname = "./mnp/FGH.bug",
                jags.modelname   = mnp.med$get.FGH.models()$postModel,
                priorModel       = mnp.med$get.FGH.models()$priorModel,
                colsToIndex      = c("cross","sex", "tissue",  "mouse"),
                gibbsBatchSize   = 30,
                additionalData = additionalData))
}


mnp.med$FGH.mediationFunc <- function(mcobj, dataForJags)
{
    samplez = list()
    dfs = list()
    
    for(i in 1:length(levels(dataForJags$cross)))
    {
        across = levels(dataForJags$cross)[i]
        if(across %in% c("FF","GG","HH"))
        {
            next
        }
        interactionOffset = paste0("f.cross.direction[2,",i,"]")
        samplez[[i]]  = (mcobj[,"beta"]*(mcobj[,interactionOffset]))
        
        df = mnp.med$samplesToMediation(samplez[[i]])
        df$moderators = paste0("Cross=", across)
        dfs = util$appendToList(dfs, df)
    }

    dfs = rbindlist(dfs)
    return(dfs)

}


mnp.med$get.FGH.models <- function()
{
    modelString = 
    "model {


    ###################################
    ## intercept, sex,
    ##################################
    ## i=1 for direct effects, i=2 for effects on mediator
    for(i in 1:2)
    {
	intercept[i] ~ dnorm(0, .0000001)
	f.sex[i]     ~ dnorm(0, .0000001)
    }


    ###################################
    ## f.cross 
    ###################################
    precision.cross <- .0000001
    for(i in 1:2)
    {
        f.cross[i,1]  <- 0
    	for(j in 2:n.cross)
    	{
             f.cross[i,j] ~ dnorm(0, precision.cross)
        }
    }


    ###########################################
    ## f.tissue  (consider making random effect)
    ###########################################
    precision.tissue <- .0000001
    for(i in 1:2)
    {
        f.tissue[i,1]  <- 0
    	for(j in 2:n.tissue)
    	{
             f.tissue[i,j] ~ dnorm(0, precision.tissue)
        }
    }


    ###################################
    ## f.sex-by-f.cross-effect
    ##################################
    precision.sex.cross <- .0000001
    for(i in 1:2)
    {
        f.sex.cross[i,1]  <- 0
    	for(j in 2:n.cross)
    	{
             f.sex.cross[i,j] ~ dnorm(0, precision.sex.cross)
        }
    }


    ###################################
    ## f.sex-by-f.tissue-effect (consider making random effect)
    ##################################
    precision.sex.tissue <- .0000001
    for(i in 1:2)
    {
        f.sex.tissue[i,1]  <- 0
    	for(j in 2:n.tissue)
    	{
             f.sex.tissue[i,j] ~ dnorm(0, precision.sex.tissue)
        }
    }


    ###################################
    ## f.cross-by-f.tissue-effect (consider making random effect)
    ##################################
    precision.cross.tissue <- .0000001
    for(i in 1:2)
    {
        f.cross.tissue[i,1]  <- 0
    	for(j in 2:n.cross.tissue)
    	{
             f.cross.tissue[i,j] ~ dnorm(0, precision.cross.tissue)
        }
    }

    ###################################
    ## f.cross-by-f.direction-effect 
    ##################################
    precision.cross.direction <- .0000001
    for(i in 1:2)
    {
        f.cross.direction[i,1]  <- 0	
    	for(j in 2:n.cross)
    	{
             f.cross.direction[i,j] ~ dnorm(0, precision.cross.direction)
        }
    }

    ###################################
    ## f.sex-by-f.cross-by-f.tissue-effect (consider making random effect)
    ##################################
    precision.sex.cross.tissue <- .0000001
    for(i in 1:2)
    {
        f.sex.cross.tissue[i,1]  <- 0
    	for(j in 2:n.cross.tissue)
    	{
             f.sex.cross.tissue[i,j] ~ dnorm(0, precision.sex.cross.tissue)
        }
    }

    
    ###################################
    ## f.sex-by-f.cross-by-f.direction-effect (consider making random effect)
    ##################################
    precision.sex.cross.direction <- .0000001
    for(i in 1:2)
    {
        f.sex.cross.direction[i,1]  <- 0
    	for(j in 2:n.cross)
    	{
             f.sex.cross.direction[i,j] ~ dnorm(0, precision.sex.cross.direction)
        }
    }


    ###################################
    ## f.cross-by-f.tissue-by-f.direction-effect (consider making random effect)
    ##################################
    precision.cross.tissue.direction <- .0000001
    for(i in 1:2)
    {
        f.cross.tissue.direction[i,1]  <- 0
    	for(j in 2:n.cross.tissue)
    	{
             f.cross.tissue.direction[i,j] ~ dnorm(0, precision.cross.tissue.direction)
        }
    }

    ###################################
    ## f.sex-by-f.cross-by-f.tissue-by-f.direction-effect (consider making random effect)
    ##################################
    precision.sex.cross.tissue.direction <- .0000001
    for(i in 1:2)
    {
        f.sex.cross.tissue.direction[i,1]  <- 0
    	for(j in 2:n.cross.tissue)
    	{
             f.sex.cross.tissue.direction[i,j] ~ dnorm(0, precision.sex.cross.tissue.direction)
        }
    }



    ###################################
    ## mouse effect
    ##################################
    for(i in 1:2)
    {
        std.mouse[i] ~ dunif(0,1000)
        precision.mouse[i] <- pow(std.mouse[i],-2)
        for(j in 1:n.mouse)
    	{
    	     ranef.mouse[i,j] ~ dnorm(0, precision.mouse[i])
    	}
    }


    ############################################################
    ## Beta- the effect of the mediator m on y.
    ## Assumed to be independent of all other covariates.
    ############################################################
    beta ~ dnorm(0, .0000001)
    ##beta <- 0

    ############################
    ## Sigmas for y and mediator.
    ############################
    for(i in 1:2)
    {
        std.epsilon[i] ~ dunif(0,1000)
        precision.epsilon[i] <- pow(std.epsilon[i],-2)
    }


    #########################################
    # Likelihood functions, indexed by s, the sample
    ###############################################
    for(s in 1:n)
    {
       ##compute the means first
       for(i in 1:2)
       {
       u[i,s] <- intercept[i]                                       +
                 f.sex[i]*sex[s]                                    +
		 f.cross[i,                     cross.j[s]]         +	
		 f.tissue[i,                   tissue.j[s]]         +

		 f.sex.cross[i,                  cross.j[s]]*sex[s]  +
		 f.sex.tissue[i,                tissue.j[s]]*sex[s]  +
		 f.cross.tissue[i,               cross.tissue[s]] +
		 f.cross.direction[i,            cross.j[s]]*direction[s] +
		 
		 f.sex.cross.tissue[i,           cross.tissue[s]]*sex[s] +
		 f.sex.cross.direction[i,        cross.j[s]]*sex[s]*direction[s] +
		 f.cross.tissue.direction[i,     cross.tissue[s]]*direction[s] +
		 f.sex.cross.tissue.direction[i, cross.tissue[s]]*sex[s]*direction[s] +

		 ranef.mouse[i,     mouse.j[s]]
       }

       ##The mediator likelihood.
       MEDIATOR_LIKELIHOOD
       ##Data likelihood
       OUTCOME_LIKELIHOOD

##       epsilon[1,s] <- y[s] - u[1,s]
    }
}"

    postModel = mnp.med$getModelString(modelString,
                                       "mediator[s] ~ dnorm(u[2,s], precision.epsilon[2])",
                                       "y[s] ~ dnorm(u[1,s] + beta*mediator[s], precision.epsilon[1])")

    priorModel = mnp.med$getModelString(modelString, "", "")
    
    out = list(postModel = postModel, priorModel = priorModel)
    return(out)

}
