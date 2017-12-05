run <- function()
{
    ##source("./mnp/runAll.R")
    ##    source("./mnp/mediationPlots2.R")

    ## source("./mnp/plotting.R")
    ## source("./mnp/micro/report.R")
    ## source("./mnp/qpcr/analysis.R")
    ## micro.report$reportAnalysis(inp$exp.mat,
    ##                             inp$cov.data.full,
    ##                             out$ident.full$results,
    ##                             inp$annot.data,
    ##                             inp$probesetInfo,
    ##                             out$threshholds,
    ##                             inp$karyo)

    ## source("./loadParams.R")
    ## source("./genomerep/variantdb2/dump_parser.R")
    
    ## db = db_builder$get.db.lite("../../data/isvdb/exon_1410/")                                  
    ## df = db$read("CC001", "CC002", chr=1, type="diplotype", storeKnownFields= F, phased = T)

    ## df1 = db$read("CC001", "CC002", chr=1, type="genotype", storeKnownFields= F, phased = F)
    ## df2 = db$read("CC001", "CC002", chr=1, type="genotype", storeKnownFields= F, phased = T)
    ##    browser()

##    source("./mnp/micro/preprocess/gethet2.R")


    ## source("./crowley/loadAllDataCrowley.R")
    ## source("./crowley/micro/analysis.R")
    ## inp = loadDataCrowley$createAllInputs()
    ## crout = crowley.analysis$run.noperm(inp)

    ## source("./mnp/runAll.R")
    ## source("./mnp/micro/analysis.R")
    
    ## out = micro.analysis$runallPerms(inp)
    
    if(!exists("prop"))
    {
        source("./loadParams.R")
    }

   ##  source("./mnp/loadAllData.R")
   ##  source("./mnp/behavior/analysis.R")
   ##  source("./mnp/micro/analysis.R")
   ##  source("./mnp/micro/report.R")
   ##  source("./mnp/qpcr/analysis.R")
   ## load(outm("perm.RData"))
    source("./mnp/loadAllData.R")
    source("./mnp/behavior/analysis.R")
    source("./mnp/micro/analysis.R")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")

    outq = data.table(qpcr.analysis$run(inp))
    setnames(outq, old = "assay", new ="effect")
    
    outq[effect=="Lrrc16a"]$effect = "POE on Lrrc16a"
    outq[effect=="Meg3"]$effect    = "Diet-by-POE on Meg3" 
    outq$pvalue = format.pval(outq$pvalue, digits = 2)

    outq = outq[phen=="Delta.Ct"]
    outq$phen = NULL
  
##    setcolorder
    fwrite(file = outm("report", "qpcr.txt"), outq)
    
    micro.report$reportAnalysis(inp$exp.mat,
                                inp$cov.data.full,
                                out$ident.full$results,
                                ##                            out$results$results,
                                inp$annot.data,
                                inp$probesetInfo,
                                out$threshholds,
                                inp$karyo)
    


}
