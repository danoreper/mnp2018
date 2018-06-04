source("./loadParams.R")
source("./mnp/micro/analysis.R")
source("./mnp/loadAllData.R")
source("./mnp/micro/preprocess/evalprobes2.R")
library(ggplot2)

qpcr.analysis = new.env(hash=T)

## reportDir = outm("report")
## dir.create(reportDir, showWarnings = F)
## mainReport = fp(reportDir, "report.txt")
## sink(mainReport, append=T)

toReport = function(str)
{
     print(str)
     ##sink(str, con = mainReport)
}


qpcr.analysis$run <- function(inp)
{
    modes = c(2)
    exp.mat                 = inp$exp.mat
    probesetInfo            = inp$probesetInfo
    goi = c("Lrrc16a", "Meg3", "Rfng", "Mir341")
    ##get microarray expression info for just genes of interest
    ##the expression info is currently unused in this analysis, but it could be useful
    expression.mic = getMicroarrayData(probesetInfo, goi, exp.mat)

############ Covariate Info ####################
    cov.data = inp$cov.data.allSamples 

    ## the microarray and covariate data which doesnt change regardless of the taqman assay we are inspecting.
    tabulated = expression.mic[cov.data]
    setkey(tabulated, "ID")

    assay="Mir341"
    mode = 2
    toplot = tabulated
    toplot$y = tabulated[[paste0(assay, ".mic")]]
    atitle = bquote(italic(.(assay))*": microarray")
    alabel = bquote("log"[2] * "(expression)")
    aplot = plot.poe.data(alldata=toplot, ylabel = alabel, atitle = atitle, mode = mode)
    show.and.write(aplot, atitle, mode, fname = fp("micro", paste0(assay,".","micro")))
    

    setnames(tabulated, old="Lrrc16a.mic", new = "Carmil1.mic")
    resultsAll = list()
###########  Analysis and plotting  ####################
## repeat analysis for Lrrc16a and Meg3 data

    for(assay in c("Carmil1", "Meg3"))
    {
        expression.taq = getTaqmanData(assay)
        print(paste0("Assay=", assay))
        
        expression.taq$Strain = NULL
        alldata = merge(tabulated, expression.taq, all=T)
        alldata$batchplate = as.factor(paste(alldata$Batch, alldata$Plate, sep="_"))
        alldata = alldata[!is.na(alldata$Delta.Ct)]

        m.string = " ~ 1  + Pipeline + Diet + Strain + Diet:Strain + (1|Dam.ID) + (1|batchplate)"

        regresstypes = list(full = alldata,
                            old = alldata[!is.na(Carmil1.mic)],
                            new = alldata[is.na(Carmil1.mic)])

        results = list()
        for(regtype in names(regresstypes))
        {
            print(paste0("regtype=", regtype))
            regressData = regresstypes[[regtype]]
            ##            for(phen in c("goi.taq", "control.taq", "Delta.Ct"))
            for(phen in c("control.taq", "Delta.Ct"))
            {
                print(paste0("phen=", phen))
                toReport(paste0(assay, ":", phen, ":", regtype))
     
                m = suppressWarnings(qpcr.analysis$.callFit(trainingData = regressData,
                                                    phen=phen,
                                                    m.string=m.string))

                if(assay == "Carmil1")
                {
                    p.val = m$anovaWrapper$an["Strain", m$anovaWrapper$pvalueCol]
                } else if (assay=="Meg3"){
                    p.val = m$anovaWrapper$an["Diet:Strain", m$anovaWrapper$pvalueCol]
                } else {
                    stop("unimplemented")
                }
                
                ##                print(m$anovaWrapper$an)
                result = data.frame(assay    = assay,
                                    dataset  = regtype,
                                    n.pups   = length(unique(regressData$ID)),
                                    phen     = phen,
                                    pvalue   = p.val)

                
                results = util$appendToList(results, result)
                
            }
        }

        results = data.frame(do.call(rbind, results))

        results$pvalue = as.character(signif(results$pvalue, 2))        
        print(results)

##        results$logp = -log10(results$pvalue)
        
        write.table(file = outm(fp("qPCR", paste0(assay,"_qpcr.regressOnStrain.csv"))), results, row.names = F)

        
        
        
        ##TODO: move plotting out of here?
        toplot = alldata
        toplot$y = -toplot$Delta.Ct
        for(mode in modes)
        {
            atitle = bquote(italic(.(assay))*": qPCR")
            legend.position = NULL
            if(assay == "Carmil1" & mode == 2)
            {
                legend.position = c(.99, .3)
            }


            aplot = plot.poe.data(toplot,
                                  ylabel = expression(paste("-", Delta,"Ct")),
                                  atitle = atitle, mode = mode, legend.position = legend.position)
            show.and.write(aplot, atitle, mode, fname = fp("qPCR", paste0(assay,".","qPCR")))
        }
        
        
        toplot = alldata
        toplot$y = alldata[[paste0(assay, ".mic")]]
        for(mode in modes)
        {
            print(paste0("mode=",mode))
            atitle = bquote(italic(.(assay))*": microarray")
            alabel = bquote("log"[2] * "(expression)")
            aplot = plot.poe.data(toplot, ylabel = alabel, atitle = atitle, mode = mode)
            show.and.write(aplot, atitle, mode, fname = fp("micro", paste0(assay,".","micro")))
        }
        resultsAll = util$appendToList(resultsAll, results)
    }
    resultsAll = rbindlist(resultsAll)

    toflex = resultsAll[phen == "Delta.Ct"]
    toflex$Effect = ""
    toflex[assay =="Carmil1"]$Effect = "POE on Carmil1"
    toflex[assay =="Meg3"]$Effect    = "DietxPOE on Meg3"
    toflex = toflex[,c("Effect", "dataset", "n.pups", "pvalue")]
    
    browser()
    mytab = regulartable(data = toflex)
    mytab = bold(mytab, part = "header")
    mytab = align( mytab, align = "center", part = "all")

    rep = list()
    rep[["x"]] = mytab
    rep[["dataset"]]="Dataset"
    rep[["n.pups"]] = "# Pups"
    rep[["pvalue"]] = "p value"
    mytab = do.call(set_header_labels, rep)
    mytab = autofit(mytab)

    doc = read_docx(fp("./mnp/template.docx"))
    doc = body_add_flextable(doc, mytab)
    print(doc, target = fp(outdir, paste0("qPCR/qpcr.docx")))
    
    return(resultsAll)
}

qpcr.analysis$.callFit <- function(trainingData, phen, m.string)
{
    transformParams = fit.model.bc$getDefaultTransformParams()
    strategy        = fit.modelg$getDefaultModelStrategy(anovaComparison = F, prefer.lme = F)
    checkAnova      = T
    trainingData$y = trainingData[[phen]]
    mA = fit.model.bc$fit(trainingData[[phen]],
                          cov.data = trainingData,
                          covariateModelString = m.string,
                          transformParams = transformParams,
                          checkAnova = checkAnova,
                          strategy = strategy)[[1]]
    print(mA$lambda)

    return(mA)
}


