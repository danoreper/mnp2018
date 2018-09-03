source("./util/loadParams.R")
source("./util/loadParams.R")
source("./mnp/loadAllData.R")
source("./lm/formulaWrapper.R")
source("./mnp/behavior/processPhenData.R")
source("./lm/fitBoxCoxModels.R")


mnp.med$get.BD.inputBuilder <- function(M.measures, Y.measures, discardMissingGeneExpression,phenotypeSpec = NULL, merge.qpcr.plate = F, raw.data = NULL)
{
    measures = c(M.measures, Y.measures)
    modelBehavior = "behavior" %in% union(M.measures, Y.measures)
    outcomeID = phenotypeSpec$phen
    if(modelBehavior)
    {
        outcomeID = paste0(phenotypeSpec$experiment, "_", phenotypeSpec$phen)
    }
    
    froot = outm("mediation")
    dir.create(froot, showWarnings = F, recursive = T)
    observeInfo = c("Strain.f[1,2]", "Strain.f[2,2]", "Strain.Diet.f", "beta")

    if(is.null(raw.data))
    {
        raw.data = loadAllData$createAllInputs()
    }

    ##The data building blocks building mediation models
    cov.data    = raw.data$phens$breedLog
    qpcr.data   = getRawTaqData(merge.qpcr.plate)

    fromFile = F
    micro.data  = mnp.med$get.sv.corrected.genes(raw.data, fromFile)$exp.mat
    phen.data   = NULL
    if(modelBehavior)
    {
        phen.data = getPhenData(phenotypeSpec, raw.data$phens)
        cov.data = phen.data[!duplicated(phen.data$ID)]
    }

    
    ##get the samples to keep after accounting for missing data
    tokeep =  getSamplesToKeep(cov.data, qpcr.data, micro.data, phen.data, c(measures), discardMissingGeneExpression)

    
    ##Filter out the bad samples, add any additional covariates
    cov.data   = cov.data[ID %in% tokeep]
    micro.data = data.table(ID = rownames(micro.data), micro.data, key = "ID")
    micro.data = micro.data[J(tokeep)]
    setkey(qpcr.data, "qpcr.ID")
    qpcr.data   = qpcr.data[J(tokeep)]
    qpcr.data$qpcr.s   = match(qpcr.data$qpcr.ID,  cov.data$ID)
    micro.data$micro.s = match(micro.data$ID, cov.data$ID)

    if(modelBehavior)
    {
        phen.data = phen.data[ID %in% tokeep]
        ##TODO may need to rename id column here
        phen.data$behavior.s   = match(phen.data$ID,  cov.data$ID)
    }
    
    ##explicitly encode Strain.Diet interaction as 4 factor levels; the first is the 0 reference
    cov.data = updateInteractionVariables(cov.data)

    outcome.modelString = "~ 1 + Pipeline + Batch + Diet + Strain + Strain.Diet + (1|Dam.ID)"
##    outcome.modelString = formulaWrapper$appendEffect("mediator", outcome.modelString)$modified.string
    colsToIndex         = c("Pipeline", "Batch", "Strain", "Diet", "Strain.Diet", "Dam.ID")

    if("behavior" %in% c(measures))
    {
        outcome.modelString = phenotypeSpec$lmerformula
        fixed.random = getSimpleFixedAndRandomEffects(outcome.modelString)
        colsToIndex  = union("Pipeline", c(fixed.random$fixedEffects, fixed.random$randomEffects))
    }
    
    dataForJags         =  mnp.med$buildTrainingData(cov.data, colsToIndex)
    dataForJags$y.n = nrow(cov.data)

    jagsmodel = buildModelString(outcome.modelString, M.measures, Y.measures, merge.qpcr.plate)

    
    
    input.template = list(
        iter             = 16000,
        froot            = froot,
        outcome.id       = outcomeID,
        dataForJags      = dataForJags,
        mediationFunc    = mnp.med$BD.mediationFunc,
        observeInfo      = observeInfo,
        jagsmodel        = jagsmodel,
        gibbsBatchSize   = 50)

    rm("raw.data")
    
    getInput <- function(mediatorID)
    {
        input = copy(input.template)
        measures = list()
        measures[[1]] = M.measures
        measures[[2]] = Y.measures
        ids = c(mediatorID, outcomeID)
        
        ## build jags data for behavioral measurements
        for(r in 1:length(measures))
        {
            if ("behavior" %in% measures[[r]])
            {
                beh = data.frame(behavior.obs = phen.data[[ids[r]]],
                                 behavior.s   = phen.data$behavior.s,
                                 behavior.j = r)
                
                newnames = gsub(colnames(beh), pattern = "behavior", replacement = paste0("behavior",r))
                setnames(beh, old = colnames(beh), new = newnames)
                tdata = mnp.med$buildTrainingData(beh, colsToIndex = c())
                input$dataForJags = c(input$dataForJags, tdata)
                input$dataForJags$behavior.n = nrow(beh)
            }
        }

        
        ##build jags data for qpcr measurements
        for(r in 1:length(measures))
        {
            if("qpcr" %in% measures[[r]])
            {
                qpcr.copy = copy(qpcr.data)
                newnames = gsub(colnames(qpcr.copy), pattern = "qpcr", replacement = paste0("qpcr",r))
                setnames(qpcr.copy, old = colnames(qpcr.copy), new = newnames)
                tdata = mnp.med$buildTrainingData(qpcr.copy, colsToIndex = c(paste0("qpcr",r,".Plate")))
                input$dataForJags = c(input$dataForJags, tdata)
                input$dataForJags$qpcr.n = nrow(qpcr.copy)
            }
        }
        input$dataForJags$qpcr.Plate.n = length(levels(qpcr.data$qpcr.Plate))

        ##build jags data for microarray measurements

        for(r in 1:length(measures))
        {
            if("micro" %in% measures[[r]])
            {
                micro = data.frame(micro.obs = micro.data[[ids[r]]],
                                    micro.s   = micro.data$micro.s,
                                    micro.j = r)

                newnames = gsub(colnames(micro), pattern = "micro", replacement = paste0("micro",r))
                setnames(micro, old = colnames(micro), new = newnames)
                tdata = mnp.med$buildTrainingData(micro, colsToIndex = c())
                input$dataForJags = c(input$dataForJags, tdata)
                input$dataForJags$micro.n = nrow(micro)
            }
        }

        return(input)

    }
                
    return(getInput)
}

getPhenData <- function(phenotypeSpec, phens)
{
    useBasic = suppressWarnings(is.na(phenotypeSpec$derivedSubset)|phenotypeSpec$derivedSubset=="")
    if(useBasic)
    {
        
        ##phen.data = phens$getExperiment(phenotypeSpec$experiment, all=T, useBreedLog = T)
        phen.data = phens$getExperiment(phenotypeSpec$experiment, all=F, useBreedLog = T)
    } else {
        phen.data = phenotypeSpec$derivedSubset
    }
     
    phen = phen.data[[phenotypeSpec$phen]]
    
    lmerformula = formulaWrapper$removeEffectAndInteractions("Pipeline", phenotypeSpec$lmerformula )$modified.string
    phen = fit.model.bc$fit(phen, phen.data, lmerformula)$phen_1
    lambda = phen$lambda
    print(lambda)
    phen = phen$y.transformed
    if(lambda<0)
    {
         phen = -1*phen
         print(lambda)
    }

    phen.data[[paste0(phenotypeSpec$experiment, "_", phenotypeSpec$phen)]] = phen
    return(phen.data)
}

getRawTaqData <- function(merge.qpcr.plate)
{
    qpcr.data = getTaqmanData("Lrrc16a")
    qpcr.data = qpcr.data[,c("ID", "Plate", "Delta.Ct"), with = F]
    
    qpcr.data$Plate = as.factor(as.character(qpcr.data$Plate))
    if(merge.qpcr.plate)
    {
        qpcr.data[ ,list(Delta.Ct = mean(Delta.Ct), Plate = 4),by="ID"] 
    }
    
    setnames(qpcr.data, old = c("Plate", "ID"), new = c("qpcr.Plate", "qpcr.ID"))

    qpcr.data$qpcr.obs = qpcr.data$Delta.Ct *-1
    qpcr.data$qpcr.obs = fit.model.bc$transform.by.lambda(y=qpcr.data$qpcr.obs, lambda = .5, T, T)$y
    qpcr.data$Delta.Ct = NULL

    
    return(qpcr.data)
}
        
 getSamplesToKeep <- function(cov.data, qpcr.data, micro.data, phen.data, measures, discardMissingGeneExpression)
{
    ##Filter out samples missing gene expression in the modalities that are being modelled
    tokeep = cov.data$ID
    ## always remove samples missing behavior if behavior is being modelled.
    if("behavior" %in% measures) { tokeep = intersect(tokeep, phen.data$ID)}
    
    if(discardMissingGeneExpression)
    {
        if("qpcr"  %in% measures) {tokeep = intersect(tokeep, qpcr.data$qpcr.ID) } 
        if("micro" %in% measures) {tokeep = intersect(tokeep, rownames(micro.data))}
    } else{
        if("qpcr"  %in% measures) {tokeep = union(tokeep, qpcr.data$qpcr.ID) }
        if("qpcr"  %in% measures) {tokeep = union(tokeep, rownames(micro.data))}  
    }
    return(tokeep)
}

##TODO extract this to be its own integer index building module=
updateInteractionVariables <- function(cov.data)
{
    n.Strain = length(levels(cov.data$Strain))
    n.Diet   = length(levels(cov.data$Diet))
    nonzero = as.integer(cov.data$Strain)!=1 & as.integer(cov.data$Diet)!=1
    cov.data$Strain.Diet = 0
    cov.data$Strain.Diet[nonzero] = (as.integer(cov.data$Diet[nonzero])-2)*(n.Strain -1) + (as.integer(cov.data$Strain[nonzero]) - 1)
    cov.data$Strain.Diet = cov.data$Strain.Diet + 1
    return(cov.data)
}


##TODO incorporate phenotype modelling
buildModelString <- function(outcome.modelString, M.measures, Y.measures, merge.qpcr.plate)
{
    measures = list()
    measures[[1]] = M.measures
    measures[[2]] = Y.measures
    
    qpcr.modelString     = "~ 1 + qpcr.Plate"
    micro.modelString    = "~ 1 "
    behavior.modelString = "~ 1 "

    ##The stuff that is in every mediation model, regardless of covariates
    jagsmodel ="model {\n\n"
    jagsmodel = paste0(jagsmodel, mnp.med$basicMedModel())    
    jagsmodel  = paste0(jagsmodel, addEffectsToModel(outcome.modelString, numMediators = 1))

    ## qpcr measures if they are needed
    if("qpcr" %in% c(measures))
    {jagsmodel  = paste0(jagsmodel, addEffectsToModel(qpcr.modelString, numMediators = 1))}

    ##micro measures if they are needed.
    if("micro" %in% union(M.measures, Y.measures))
    {jagsmodel  = paste0(jagsmodel, addEffectsToModel(micro.modelString,   numMediators = 1))}
        
    
    if("behavior" %in% union(M.measures, Y.measures))
    {jagsmodel  = paste0(jagsmodel, addEffectsToModel(behavior.modelString,   numMediators = 1))}

    outcomeContribution = mnp.med$getLatentContribution(outcome.modelString)
    outcomeLikelihood = gsub(pattern = "OUTCOME_CONTRIBUTION", replacement = outcomeContribution,
       "for(s in 1:y.n)
       {
           for(j in 1:2)
           {
               y.u[j,s] <- y.intercept[j] + 
               OUTCOME_CONTRIBUTION
           }
           y[1,s] ~ dnorm(y.u[1,s],               y.precision.epsilon[1])
           y[2,s] ~ dnorm(y.u[2,s] + beta*y[1,s], y.precision.epsilon[2])
       }\n\n")
    

    qpcr.likelihood= c("", "")
    for (r in c(1,2))
    {
        measure = measures[[r]]
        varianceModel = ifelse(merge.qpcr.plate & !"micro" %in% measure, "1000000000",paste0(" qpcr.precision.epsilon[",r,"]"))
        plateEffect   = ifelse(merge.qpcr.plate, "",  paste0(" + qpcr.Plate.f[",r,",qpcr",r,".Plate.k[o]]"))
        ##TODO: this falls apart in the presence of missing data-- qpcr.s is na for missing data, can't have NA as an index. Need to explicitly create Fake qpcr id data for the missing animals, complete with a plate, or simulate the plate in jags.
        qpcr.likelihood[r] = paste0(
    "   for(o in 1:qpcr.n)
        {
            qpcr",r,".u[o] <- y[",r,", qpcr",r,".s[o]] ", plateEffect, "\n    
            qpcr",r,".obs[o] ~  dnorm(qpcr",r,".u[o], ", varianceModel, ") 
        }\n\n")
    }

    
    micro.likelihood = c("", "")
    for(r in c(1,2))
    {
        measure = measures[[r]]
        
        errorModelling = "qpcr" %in% measure & "micro" %in% measure
        varianceModel = ifelse(errorModelling, paste0("micro.precision.epsilon[micro",r,".j[o]]"),"100000000")
        
        micro.likelihood[r] = paste0(
        "   for(o in 1:micro.n)
        {
            micro",r,".obs[o] ~ dnorm(y[micro",r,".j[o], micro",r,".s[o]], ",varianceModel, ")
        }\n\n")
    }

    behavior.likelihood = c("", "")
    for(r in c(1,2))
    {
        measure = measures[[r]]
        
        varianceModel = "100000000"
        
        behavior.likelihood[r] = paste0(
        "   for(o in 1:behavior.n)
        {
            behavior",r,".obs[o] ~ dnorm(y[behavior",r,".j[o], behavior",r,".s[o]], ",varianceModel, ")
        }\n\n")
    }   
        
    jagsmodel = paste0(jagsmodel, outcomeLikelihood)
    for(r in c(1,2))
    {
        measure = measures[[r]]
        if("micro" %in% measure)
        {
            jagsmodel = paste0(jagsmodel, micro.likelihood[r])
        }
        if("qpcr" %in% measure)
        {
            jagsmodel = paste0(jagsmodel, qpcr.likelihood[r])
        }
        if("behavior" %in% measure)
        {
            jagsmodel = paste0(jagsmodel, behavior.likelihood[r])
        }
        
    }


    jagsmodel = paste0(jagsmodel, "}")
    
    cat(jagsmodel)
    return(jagsmodel)
}


addEffectsToModel <- function(modelString, numMediators)
{
    p = paste
    fixed.random = getSimpleFixedAndRandomEffects(modelString)
    modelString = ""
    for(simpleMainEffect in fixed.random$fixedEffects)
    {
        modelString = p(modelString, mnp.med$getEffectDistribution(simpleMainEffect, numMediators))
    }
    
    for(randomEffect in fixed.random$randomEffects)
    {
        modelString = p(modelString, mnp.med$getRandomEffectDistribution(randomEffect, numMediators))
    }
    return(modelString)
}


getSimpleFixedAndRandomEffects <- function(modelString)
{
    s.info = formulaWrapper$parseCovariateString(modelString)
    simpleFixedEffects = unlist(lapply(FUN=paste, X= s.info$fixef, collapse = "."))
    simpleFixedEffects = setdiff(simpleFixedEffects, c("1", "mediator"))
    simpleFixedEffects = setdiff(simpleFixedEffects, "mediator")
    randomEffects = unlist(lapply(FUN = "[[", X = s.info$ranef, "group"))
    return(list(fixedEffects = simpleFixedEffects, randomEffects = randomEffects))
}



mnp.med$BD.mediationFunc <- function(mcobj, dataForJags)
{

    dfs = list()
    as = list()
    bs = list()
    cs = list()
    b =  mcobj[,"beta"]
    for(k in 1:length(levels(dataForJags$Diet)))
    {
        interactionOffset = paste0("Strain.Diet.f[1,",k,"]")
        a = (mcobj[,"Strain.f[1,2]"] + mcobj[,interactionOffset])

        interactionOffset = paste0("Strain.Diet.f[2,",k,"]")
        c.prime = (mcobj[,"Strain.f[2,2]"] + mcobj[,interactionOffset])

        df = mnp.med$samplesToMediation(a, b, c.prime)
        df$moderators = paste0("Diet=", levels(dataForJags$Diet)[k])

        
        dfs = util$appendToList(dfs, df)
        as[[k]] = a
        bs[[k]] = b
        cs[[k]] = c.prime
    }

    for(k in 1:4)
    {
        margsum = as[[k]]*bs[[k]]
        if(k==1)
        {
            tots = margsum
        } else {
            tots = tots + margsum
        }
    }
    tots = tots/4


    for(k in 1:4)
    {
        margsum = (tots - as[[k]]*bs[[k]])^2
        if(k==1)
        {
            newtot = margsum
        } else {
            newtot = newtot + margsum
        }
    }
    newtot = newtot/4

##    browser()
    F1 = median(newtot)/sd(newtot)
    
    Y = c(as[[1]]*bs[[1]],
          as[[2]]*bs[[2]],
          as[[3]]*bs[[3]],
          as[[4]]*bs[[4]])

    N = length(as[[1]])
    X = factor(c(rep(1, N), rep(2,N), rep(3,N), rep(4,N)))
    an = (anova(lm(Y~X)))
    F2  = an["X", "F value"] 
    
    
    df = mnp.med$samplesToMediation(mcmc(do.call(c, as)),
                                    mcmc(do.call(c, bs)),
                                    mcmc(do.call(c, cs)))

    
    df$moderators = paste0("Diet=Ave")

    dfs = util$appendToList(dfs, df)
    dfs = rbindlist(dfs)

    return(list(df = dfs, F1 = F1, F2 = F2))
}


mnp.med$getEffectDistribution <- function(simpleMainEffect, numMediators=1)
{
    p = paste0
    m = simpleMainEffect
    M = 1 + numMediators
    
    astring = p(        "for(j in 1:",M,")\n")
    astring = p(astring,"{\n")
    astring = p(astring,"  ",m,".precision[j] <- pow(bound*sd.multiple, -2)\n")
    astring = p(astring,"  ",m,".f[j,1] <- 0\n")
    astring = p(astring,"   for(k in 2:",m,".n)\n")
    astring = p(astring,"   {\n")
    astring = p(astring,"        ",m,".f[j,k] ~ dnorm(0, ",m,".precision[j])\n")
    astring = p(astring,"   }\n")
    astring = p(astring,"}\n")
    astring = p(astring,"\n\n")
    
    return(astring)
}


mnp.med$getRandomEffectDistribution <- function(randomEffect, numMediators=1)
{
    p = paste0
    m = randomEffect
    M = 1 + numMediators
    
    astring = p(       "for(j in 1:",M,")\n")
    astring = p(astring,"{\n")
    astring = p(astring,"  ",m,".std[j] ~ dunif(0, bound*sd.multiple)\n")
    astring = p(astring,"  ",m,".precision[j] <- pow(",m,".std[j],-2)\n")
    astring = p(astring,"   for(k in 1:",m,".n)\n")
    astring = p(astring,"   {\n")
    astring = p(astring,"       ",m,".ranef[j,k] ~ dnorm(0, ",m,".precision[j])\n")
    astring = p(astring,"   }\n")
    astring = p(astring,"}\n")
    astring = p(astring,"\n\n")
    
    return(astring)
}

mnp.med$basicMedModel <- function()
{
    jagsmodel = "
        ##TODO pass this in
        bound <- 5

        sd.multiple <-1
        ##################################
        ##  Priors for intercept
        ##############################
        for(j in 1:2)
        {
            y.intercept[j] ~ dnorm(0, pow(bound*sd.multiple, -2))
            qpcr.kappa[j]  ~ dnorm(0, pow(bound*sd.multiple, -2))
##          micro.kappa[j] ~ dnorm(0, pow(bound*sd.multiple, -2))
        }
        
        ############################################################
        ## Beta- the effect of the mediator m on y.
        ## Assumed to be independent of all other covariates.
        ############################################################
        beta ~ dnorm(0, pow(bound*sd.multiple, -2))

        ############################
        ## Sigmas for y and mediator.
        ############################
        for(j in 1:2)
        {
            y.std.epsilon[j]           ~ dunif(0,bound*bound)
            y.precision.epsilon[j]     <- pow(y.std.epsilon[j],-2)

            qpcr.std.epsilon[j]       ~ dunif(0, bound*bound)
            qpcr.precision.epsilon[j] <- pow(qpcr.std.epsilon[j],-2)

            ##TODO: necessary?
            micro.std.epsilon[j]       ~ dunif(0, bound*bound)
            micro.precision.epsilon[j] <- pow(micro.std.epsilon[j],-2)
        }


        ##beta <- 0
        \n\n"
    return(jagsmodel)
}

mnp.med$getLatentContribution <- function(modelString)
{
    jagsmodel = ""
    fixed.random = getSimpleFixedAndRandomEffects(modelString)
    me = lapply(FUN = mnp.med$getLikContribution, fixed.random$fixedEffects, type = "f")
    me = do.call(paste, c(me, sep = "        +\n"))
    jagsmodel = paste(jagsmodel, me, "         +\n")
    
    re = lapply(FUN = mnp.med$getLikContribution, fixed.random$randomEffects, type = "ranef")
    re = do.call(paste, c(re, sep = "        +\n"))
    jagsmodel = paste(jagsmodel, re, "\n")
    return(jagsmodel)
}

mnp.med$getLikContribution <- function(simpleMainEffect, type)
{
    m = simpleMainEffect
    p = paste0
    astring = p("      ",simpleMainEffect,".", type, "[j,      ",m,".k[s]]")
    return(astring)
}


##TODO: refactor to another phenotype input files? 
mnp.med$getAllMediationPhenModelSpecs <- function(phenz)
{
    ##The simple models that can be described just in a csv file
    allmodels  = fread(datm("behaviorModels.csv"))

    allmodels2 = list()
    for(i in 1:nrow(allmodels)) { allmodels2 = util$appendToList(allmodels2, allmodels[i,])}
    allmodels = allmodels2

    ##The Startle models that are more complicated to describe, so we'll create the necessary datasets here.
    startle = phenz$getExperiment("startle")
    per.pp  = processBehavior$buildStartleByGroupMean(startle, phenz$breedLog, all = F)
    for(pp in c("PP74", "PP78", "PP82", "PP86", "PP90"))
    {

        phen = paste0("mean_", startlechoice, "_",  pp)
        
        modelSpec = list(experiment = "startle", derivedSubset = per.pp, phen = phen,
                      lmerformula = " ~ 1 + Batch + Diet + Strain + Strain:Diet  + (1 | Dam.ID) ")
        allmodels = util$appendToList(allmodels, modelSpec)
        
    }
    
    modelSpec = list(experiment = "startle",
                  derivedSubset = processBehavior$getAS50(startle, phenz$breedLog, all =F),
                  phen = "as50.Average_normalized",
                  lmerformula = " ~ 1 + Batch + Diet + Strain + Strain:Diet  + (1 | Dam.ID) ")
    allmodels = util$appendToList(allmodels, modelSpec)

    ## modelSpec = list(experiment = "startle",
    ##                  derivedSubset = processBehavior$getAS50(startle, phenz$breedLog, all =F),
    ##                  phen = "as50.Latency_normalized",
    ##                  lmerformula = " ~ 1 + Batch + Diet + Strain + Strain:Diet  + (1 | Dam.ID) ")
##    allmodels = util$appendToList(allmodels, modelSpec)
    return(allmodels)
}
