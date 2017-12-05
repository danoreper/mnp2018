library("MASS")
library("rjags") 
library("coda")
library("mvtnorm")
library("pracma")
##library("coefplot2")
library(ggplot2)
load.module("glm")
options(warnPartialMatchDollar=T)

showTiming = F
##library("gdata")
##install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source")
##library("caterplot")
##caterplot(mcobj[,rownames(sigterms)])
##coefplot2(mcobj[,rownames(sigterms)])

library(mcmcplots)

fitjags = new.env(hash=T)


fitjags$runJags <- function(trainingData,
                            modelname,
                            observeInfo,
                            iter=10000,
                            burninFrac = .4,
                            nchains= 1,
                            thin = 1,
                            phen = phen,
                            colsToIndex,
                            additionalData = list()) 
{

    if(showTiming)
    {
        print("about to run model")
    }
    burnin = round(burninFrac*iter) #number of burnin iterations
    dir.create("../output", showWarnings=F)
    
    allOut = fitjags$.runJags(trainingData = trainingData,
                             modelname = modelname,
                             observeInfo=observeInfo,
                             iterAdapt = burnin,
                             iterSample = (iter-burnin),
                             nchains,
                             phen = phen,
                             colsToIndex,
                             thin = thin,
                             additionalData = additionalData)
    if(showTiming)
    {
        print("done running jags")
    }
    
    gc()
    vals  = allOut$vals
    mcobj = mcmc(vals)
    return(mcobj)
}


fitjags$.runJags <- function(trainingData,
                            modelname,
                            observeInfo,
                            iterAdapt,
                            iterSample,
                            nchains=1,
                            thin=1,
                            phen,
                            colsToIndex,
                            additionalData = list()) 
{
    if(showTiming)
    {
        tic()
    }

    browser()
    jags = fitjags$.buildModel(trainingData = trainingData,
                               modelname = modelname,
                               nchains = nchains,
                               iterAdapt = iterAdapt,
                               phen=phen,
                               colsToIndex = colsToIndex,
                               additionalData = additionalData)
    if(showTiming)
    {
        print("total time for jags burnin:")
        toc()
    }
    startTime = proc.time()
    ##use the training data coming out of build model- it has been modified to include the appropriate factor levels
    vals= fitjags$.buildSamples(jags = jags$jags, observeInfo = observeInfo, iterSample = iterSample, trainingData = jags$trainingData, thin=thin)
    totalTime = proc.time() - startTime
    if(showTiming)
    {
        print("total time for jags recorded sampling:")
        print(totalTime)
    }
    
    return(list(vals=vals, trainingData = jags$trainingData, dataForJags = jags$dataForJags, jagsModel = jags$jags))
}


fitjags$.buildModel <- function(trainingData, modelname, nchains, iterAdapt, phen, colsToIndex, additionalData = list())
{
    ##get an integer per factor level for passing into jags. will variables of the form: varname.j
    trainingData = fitjags$.formIntegerRepresentationOfFactors(trainingData, colsToIndex = colsToIndex)

    
    ##Get num levels per factor
    countPerFactor = list()
  ##  covmatPerFactor = list()
    for(cname in colsToIndex)
    {
        ##jags needs to know the count of levels of each factor
	countPerFactor[[paste0("n.", cname)]]       = length(levels(trainingData[[cname]]))
##        covmatPerFactor[[paste0("covmat.", cname)]] = diag(length(levels(trainingData[[cname]])))
    }

    
    ##Jags fit
    dataForJags   = as.list(trainingData)
    dataForJags$n = nrow(trainingData)
    dataForJags   = c(dataForJags, countPerFactor)
    dataForJags$y = dataForJags[[phen]]
    dataForJags = c(dataForJags, additionalData)
    browser()
    
    cat(modelname)
##    print(paste0("running: ", modelname))
    jags <- jags.model(modelname, data=dataForJags, n.chains = nchains, n.adapt = iterAdapt,quiet=!showTiming)
    if(showTiming)
    {
        print("done fitting")
    }

    return(list(jags=jags, dataForJags = dataForJags, trainingData = trainingData))
}


fitjags$.buildSamples <- function(jags, observeInfo, iterSample, trainingData, thin=1)
{
    vals = coda.samples(model=jags, 
			variable.names=observeInfo$variable, 
			n.iter =iterSample,
			thin=thin, quiet=!showTiming,progress.bar="none")

    ## z = gelman.diag(vals, multivariate = F)
    ## if(any(z$psrf[,2]>1.01, na.rm =T))
    ## {
    ##     print(z$psrf)
    ##     browser()
    ## }
    vals = vals[[1]] #single chain for now
    gc()
    if(showTiming)
    {
        print("done sampling")
    }
    toSpecifyName = which(!is.na(observeInfo$factorForIndexConversion))
    for(i in toSpecifyName)
    {
        
        observedVariable = observeInfo[i, "variable"] 
        relevantFactor   = observeInfo[i,"factorForIndexConversion"]
        levelz = levels(trainingData[[relevantFactor]])
        oldcolnames = paste0(observedVariable,"[", 1:length(levelz),   "]")
        newcolname  = paste0(observedVariable, "[", as.character(levelz), "]")
        replaceInds = match(oldcolnames, colnames(vals))
        if(all(is.na(replaceInds)))
        {
            next
        }
        colnames(vals)[replaceInds]=newcolname
    }

    return(vals)
}

fitjags$.formIntegerRepresentationOfFactors <- function(trainingData, colsToIndex)
{
    for(cname in colsToIndex)
    {
        ##print(cname)
        trainingData[[cname]] = as.factor(trainingData[[cname]])
        ##jags needs an integer representation of any levels of a factor- store these explictly in the data frame
        valz = trainingData[[cname]]		
        ## trainingData[[paste0(cname, ".j")]]  = match(valz, levels(valz))
        trainingData[[paste0(cname, ".j")]]  = as.integer(valz)
    }
    
    return(trainingData)
}


###########
##Prediction-- consider putting this in another file


fitjags$fitDataAndCreatePredictive <- function(trainingData, modelname, iter, burninfrac, observeInfo, toPredict) 
{
    observeInfo    = rbind(observeInfo, toPredict$predictObserveInfo)
    trainingData   = rbind(trainingData, toPredict$predictInputFrame)
    reconstructed  = fitjags$runModel(trainingData = trainingData, 
                                      modelname = modelname, 
                                      observeInfo = observeInfo, 
                                      iter=iter,
                                      burninFrac = burninfrac, 
                                      nchains = 1)
    return(reconstructed)
}

fitjags$.getPredictVariablesToObserve <- function(n.training, predictCovariates, predictPhenotype) 
{
    numExist = n.training
    numPredict = nrow(predictCovariates)
    startInd = numExist+1
    variablesToObserve = paste0(predictPhenotype, "[", startInd:(startInd + numPredict - 1), "]")
    return(variablesToObserve)
}

##TODO: refactor to use jags instead
fitjags$getPredictionParams <- function(n.training, predictCovariates, predictPhenotypes) 
{
    counter = 1
    predictObserveInfo = list()
    for(predictPhenotype in predictPhenotypes)
    {
        print(predictPhenotype)
        variablesToObserve = fitjags$.getPredictVariablesToObserve(n.training, 
                                                                   predictCovariates = predictCovariates,
                                                                   predictPhenotype = predictPhenotype)
        predictObserveInfo[[counter]] = data.frame(variable = variablesToObserve, keep=T, factorForIndexConversion=NA)	
        counter = counter + 1
    }
    predictObserveInfo = do.call(rbind, predictObserveInfo)
        
    ##TODO:fill in programatically reather than hardcoding mom, pop, etc? probably not worth it 
    predictInputFrame = cbind(predictCovariates)
    
    for(phen in predictPhenotypes)
    {
        predictInputFrame[[phen]] = NA
    }
    predictInputFrame = readFertilityData$addDerived(predictInputFrame)
    
    return(list(predictInputFrame = predictInputFrame, predictObserveInfo = predictObserveInfo))
}

fitjags$extractPredictionFrame <- function(toPredict, reconstructed, trainingData) 
{
    phens = c("littersize", "wean","fwean")
    newdf = rep(1:nrow(toPredict$predictInputFrame), each=nrow(reconstructed))
    newdf = toPredict$predictInputFrame[newdf,]			
    
    for(phen in phens)
    {
        varToObserve = posthocModel$getPredictVariablesToObserve(nrow(trainingData),toPredict$predictInputFrame,phen)
        print(phen)
        phenvals = NA
        if(all(!varToObserve %in% colnames(reconstructed)))
        {
            next
        }
        phenvals = try(c(reconstructed[,varToObserve]))
        newdf[[phen]] = phenvals
    }
    return(newdf)
}
