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
    first = F
    source("./mnp/behavior/loadData.R");
    source("./mnp/behavior/analysis.R");
    source("./mnp/loadAllData.R")
    source("./mnp/behavior/analysis.R")
    source("./mnp/micro/analysis.R")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")
    if(first)
    {
        load(outm(fp("micro","perm.RData")))
        inp  = loadAllData$createAllInputs()
    }
    source("./mnp/behavior/loadData.R");
    source("./mnp/behavior/analysis.R");
    source("./mnp/loadAllData.R")
    source("./mnp/behavior/analysis.R")
    source("./mnp/micro/analysis.R")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")

    ## qpcr.analysis$run(inp)

    micro.report$reportAnalysis(inp$exp.mat,
                            inp$cov.data.full,
                            out$ident.full$results,
                            ##                            out$results$results,
                            inp$annot.data,
                            inp$probesetInfo,
                            out$threshholds,
                            inp$karyo)

    
    ## phenz = loadBehavior$getPhenotypeRepository();
    ## beh.analysis$runAll(phenz)
}


