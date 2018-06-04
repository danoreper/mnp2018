run <- function()
{
    ##source("./mnp/mediation/mediationPlots2.R")

    ## if(!exists("prop"))
    ## {
    ##     source("./loadParams.R")
    ## }

    ## source("./mnp/loadAllData.R")
    ## source("./mnp/behavior/analysis.R")
    ## source("./mnp/micro/analysis.R")
    ## source("./mnp/micro/report.R")
    ## source("./mnp/qpcr/analysis.R")

    ## ## out = micro.analysis$runallPerms(inp)
    ## micro.report$reportAnalysis(inp$exp.mat,
    ##                             inp$cov.data.full,
    ##                             out$ident.full$results,
    ##                             ##                            out$results$results,
    ##                             inp$annot.data,
    ##                             inp$probesetInfo,
    ##                             out$threshholds,
    ##                             inp$karyo)

##    source("./mnp/behavior/cleanData.R");
    ## source("mnp/behavior/loadData.R");
    ## source("mnp/behavior/analysis.R");
    ## phenz = loadBehavior$getPhenotypeRepository();
    ## beh.analysis$runAll(phenz)

    ## source("./perm/freedman.R")
    ## out = micro.analysis$runallPermsFreed(inp)

    source("./lm/fitBoxCoxModels.R")
    failed = c()
    clusterOut = list()
    for(i in 1:length(funcArgs))
    {
        argz = c(funcArgs[[i]], otherGlobals)
        clusterOut[[i]] = try(do.call(func, argz))
        if(length(clusterOut[[i]])==1 && class(clusterOut[[i]])=="try-error")
        {
            ##print(clusterOut[[i]])
            failed = c(failed, i)
        }
    }

}


