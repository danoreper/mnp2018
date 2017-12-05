source("./mnp/BDmodel.R")
source("./mnp/behavior/processPhenData.R")
       
mnp.med$getPhenInputBuilder <- function(discardMissingExpression = T, discardMissingBehavior = T)
{
    inputBuilder = new.env(hash=T)
    ##only for debugging
    ## observeInfo = data.frame(variable = c("precision.Strain", "precision.epsilon", "intercept", "beta", 
    ##                                       "f.Strain", "f.Diet", "f.Strain.Diet", "f.Batch", "f.Pipeline" ),
    ## factorforindexconversion = c(NA, NA, NA, NA,
    ##                              "Strain"[, "Diet", NA, "Batch", "Pipeline"),
    ## keep = T)
    
    observeInfo = data.frame(variable = c("f.Strain[1,2]",
                                          "f.Strain[2,2]",
                                          "f.Strain.Diet",
                                          "beta"),
                             factorforindexconversion = c("f.Strain", NA, NA, NA),
                             keep = T)
    observeInfo$variable = as.character(observeInfo$variable)
    
    froot = outm( "mediation", "behavior")
    dir.create(froot, showWarnings = F, recursive = T)

    raw.data   = loadAllData$createAllInputs()
    exp.mat    = mnp.med$get.sv.corrected.genes(raw.data, T)$exp.mat
    
    inputBuilder = new.env(hash=T )
    inputBuilder$phenz = raw.data$phens
    
    inputBuilder$get.phen.inputs <- function(modelSpec) 
    {

        additionalData = list()

        useBasic = suppressWarnings(is.na(modelSpec$derivedSubset)|modelSpec$derivedSubset=="")
        if(useBasic)
        {
            cov.data = inputBuilder$phenz$getExperiment(modelSpec$experiment, all = T)
        } else {
            cov.data = modelSpec$derivedSubset
        }

        keep = !is.na(cov.data[[modelSpec$phen]])
        allnames = setdiff(union(cov.data[keep]$ID, rownames(exp.mat)),NA)
        setkey(cov.data, "ID")
        cov.data = cov.data[J(ID=allnames)]
        exp.mat.copy = data.table(ID = rownames(exp.mat), exp.mat, key = "ID")
        exp.mat.copy  = exp.mat.copy[J(ID=allnames)]

        if(discardMissingBehavior)
        {
            subcov = cov.data[is.na(cov.data[[modelSpec$phen]])]
            toRemove = unique(subcov$ID)
            cov.data = cov.data[!(ID %in% toRemove)]
            exp.mat.copy  = exp.mat.copy[!(ID %in% toRemove)]
        }
        if(discardMissingExpression)
        {
            subex  = exp.mat.copy[as.vector(is.na(exp.mat.copy[,2,with=F]))]
            toRemove = unique(subex$ID)
            cov.data = cov.data[!(ID %in% toRemove)]
            exp.mat.copy  = exp.mat.copy[!(ID %in% toRemove)]
        }


        ##TODO extract this to be its own integer index building module
        n.Strain = length(levels(cov.data$Strain))
        n.Diet   = length(levels(cov.data$Diet))
        nonzero = as.integer(cov.data$Strain)!=1 & as.integer(cov.data$Diet)!=1
        cov.data$Strain.Diet = 0
        cov.data$Strain.Diet[nonzero] = (as.integer(cov.data$Diet[nonzero])-2)*(n.Strain -1) + (as.integer(cov.data$Strain[nonzero]) - 1)
        cov.data$Strain.Diet = cov.data$Strain.Diet + 1

        exp.mat.copy = exp.mat.copy[J(ID = cov.data$ID), on="ID"]
        outcome.mat = cov.data[[modelSpec$phen]]
        
        ####transform the data
        ####Need to remove pipeline for lmer to work.
        lmer.behavior = formulaWrapper$removeEffectAndInteractions("Pipeline", modelSpec$lmerformula)$modified.string
        
        outcome.mat = fit.model.bc$fit(outcome.mat, cov.data, lmer.behavior)$phen_1

        lambda = outcome.mat$lambda
        print(lambda)
        outcome.mat = outcome.mat$y.transformed
        if(lambda<0)
        {
            outcome.mat = -1*outcome.mat
            print(lambda)
        }

        
        outcome.mat = as.matrix(outcome.mat)
        outcome.id = paste0(modelSpec$experiment, "_", modelSpec$phen)
        colnames(outcome.mat) = outcome.id
        exp.mat.copy = cbind(exp.mat.copy, outcome.mat)
        cov.data[[modelSpec$phen]] = NULL

        
        s.info = formulaWrapper$parseCovariateString(modelSpec$lmerformula)
        simpleFixedEffects = unlist(lapply(FUN=paste, X= s.info$fixef, collapse = "."))
        simpleFixedEffects = setdiff(simpleFixedEffects, "1")


        randomEffects = unlist(lapply(FUN = "[[", X = s.info$ranef, "group"))

        exp.mat.copy = as.matrix(exp.mat.copy[,2:ncol(exp.mat.copy)])
        rownames(exp.mat.copy) = cov.data$ID
        cov.data = data.frame(cov.data)


        return(list(froot            = froot,
                    lmer.modelString = lmer.behavior,
                    cov.data         = cov.data,
                    outcome.mat      = exp.mat.copy,
                    mediationFunc    = mnp.med$BD.mediationFunc,
                    observeInfo      = observeInfo,
                    ##                jags.modelname = "./mnp/BD.bug",
                    jags.modelname   = mnp.med$get.BD.models(simpleFixedEffects, randomEffects)$postModel,
                    priorModel       = mnp.med$get.BD.models(simpleFixedEffects, randomEffects)$priorModel,
                    colsToIndex      = union("Pipeline", c(simpleFixedEffects, randomEffects)),
                    gibbsBatchSize   = 40,
                    additionalData   = additionalData))
    }
    return(inputBuilder)
}


##TODO: refactor to another phenotype input files? 
mnp.med$getAllMediationPhenModelSpecs <- function(phenz)
{
    ##The simple models that can be described just in a csv file
    allmodels  = fread("../data/sva/behaviorModels.csv")

    allmodels2 = list()
    for(i in 1:nrow(allmodels)) { allmodels2 = util$appendToList(allmodels2, allmodels[i,])}
    allmodels = allmodels2

    ##The Startle models that are more complicated to describe, so we'll create the necessary datasets here.
    startle = phenInputBuilder$phenz$getExperiment("startle")
    per.pp  = processBehavior$buildStartleByGroupMean(startle, phenz$breedLog, all = T)
    for(pp in c("PP74", "PP78", "PP82", "PP86", "PP90"))
    {

        phen = paste0("mean_", startlechoice, "_",  pp)
        
        modelSpec = list(experiment = "startle", derivedSubset = per.pp, phen = phen,
                      lmerformula = " ~ 1 + Batch + Diet + Strain + Strain:Diet  + (1 | Dam.ID) ")
        allmodels = util$appendToList(allmodels, modelSpec)
        
    }
    
    modelSpec = list(experiment = "startle",
                  derivedSubset = processBehavior$getAS50(startle, phenz$breedLog, all = T),
                  phen = "as50.Average_normalized",
                  lmerformula = " ~ 1 + Batch + Diet + Strain + Strain:Diet  + (1 | Dam.ID) ")

    allmodels = util$appendToList(allmodels, modelSpec)
    return(allmodels)
}

