source("./mnp/loadAllDataCrowley.R")
mediation = new.env(hash=T)




##TODO: refactor to use only Probe.Set.ID... conisder a version that uses gene_id 
mediation$checkMediationForCandidates <- function(candidates, exp.mat, outcome.vec, mediator.vec, cov.data, mediatorModelString, mediatedFactor)
{
    strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = T, prefer.lme = F)
    dfs = list()
    for(i in 1:nrow(candidates))
    {
        thegene = candidates[i]$gene_name
        ## if(thegene=="Lrrc16a")
        ## {
        ##     browser()
        ## }
        print(thegene)
        
        ps.id = candidates[i]$Probe.Set.ID
        if(!(ps.id %in% colnames(exp.mat))) { next }
        
        mediator.vec = exp.mat[,ps.id]
        mediator.vec = mediator.vec[cov.data$ID]
        df = try(lm.mediation$simpleTest(outcome.vec         = outcome.vec,
                                         mediator.vec        = mediator.vec,
                                         cov.data            = cov.data,
                                         mediatorModelString = mediatorModelString,
                                         mediatedFactor      = mediatedFactor,
                                         strategy            = strategy))

        if(class(df)=="try-error")
        {
            next
        }
        df$mediator = thegene
        df$mediator_probeset = ps.id
        dfs = util$appendToList(dfs, df)
        
    }
    df = rbindlist(dfs)
    setcolorder(df, c("mediator_probeset", "mediator", "factor.on.mediator", "complete.mediation", "partial.mediation"))
    return(df)
}



mediation$get.sv.corrected.genes <- function(inp, fromFile = T)
{

    residualizeOutCovariates = c()
    if(!fromFile)
    {
        pracma::tic()
        svinfo = surrogatcalc$generate.svinfo(svFunc = micro.analysis$get.SV.func(),
                                              exp.mat = inp$exp.mat,
                                              exp.mat.control = inp$exp.mat.control,
                                              cov.data = inp$cov.data.full,
                                              covariateModelString = micro.analysis$formCovariateFullModelString(),
                                              residualizeOutCovariates = residualizeOutCovariates,
                                              residualizeOutSV = T,
                                              accum = bsub$get.mc.accumulator(mc.cores = 8))
                
        print(paste0("got surrogate variables:",pracma::toc()))
        gc()
        
        print("done getting resids")
        save(file=outm( "sv.corrected.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outm("sv.corrected.RData"))
    }

    return(svinfo)
    ## exp.mat = svinfo$exp.mat
    ## lrrc16a  = getProbesetId(inp$probesetInfo, "Lrrc16a")
    ## geneExp = exp.mat[,lrrc16a]
    ## geneExp = list(geneExp = geneExp, name = "Lrrc16a")

}
                                             




##TODO refactor
mediation$.model.SEM.behaviorWithMicroarray <- function (inp)
{
    residualizeOutCovariates = c("Batch", "Strain", "Diet")
    
    if(!fromFile)
    {
        pracma::tic()
##        f = function(exp.mat, exp.mat.control, cov.data, covariateModelString)
        
        svinfo = surrogatcalc$generate.svinfo(svFunc = micro.analysis$get.SV.func(),
                                              exp.mat = inp$exp.mat,
                                              exp.mat.control = inp$exp.mat.control,
                                              cov.data = inp$cov.data.full,
                                              covariateModelString = micro.analysis$formCovariateFullModelString(),
                                              lambdasToTry = prop$mnp$lambdasToTry,
                                              residualizeOutCovariates = residualizeOutCovariates,
                                              residualizeOutSV = T)
                
        print(paste0("got surrogate variables:",pracma::toc()))
        gc()
        
        print("done getting resids")
        save(file=outm("behresids2.RData"), list=c("svinfo"))
    }
    if(fromFile)
    {
        load(file = outm("behresids2.RData"))
    }
    
    exp.mat = svinfo$exp.mat
    lrrc16a  = getProbesetId(inp$probesetInfo, "Lrrc16a")
    geneExp = exp.mat[,lrrc16a]
    geneExp = list(geneExp = geneExp, name = "Lrrc16a")
    out = mediation$.parSemFunc(phen=inp$phens, geneExp)
    
    
    
    
    numSamp = 46*10
    samplez = sample(1:ncol(resids), numSamp)
    samplez = util$generateShuffles(1:nrow(resids), numShuffles = numSamp)
    
    accum = getAccumulator()
    accum$init(mediation$.parSemFunc, otherGlobals =list(phen=inp$phens))
    outs = list()
    for(i in 1:numSamp)
    {
        ## gm = samplez[i]
        ## geneExp = resids[,gm]
        ## geneExp = list(geneExp = geneExp, name = gm)

        gm = i
        names(geneExp) = names(geneExp)[samplez[,i]]
        geneExp = list(geneExp = geneExp, name = gm)
        accum$addCall(list(geneExp = geneExp))
        
##        outs = util$appendToList(outs, out)
    }
    
    outs=accum$runAll()

    if(accum$outputs.files)
    {
        collated = collateFiles(outs)
    } else {

        collated = outs;

        for(i in 1:length(collated))
        {
            if(class(collated[[i]])=="try-error")
            {
                collated[[i]] = NA
            }
        }
        collated = do.call(rbind, collated)
    }
    write.table(outs, file ="outs.txt", sep=",", row.names = F)
}


mediation$.parSemFunc <- function(phen, geneExp)
{
    df = run(phen, geneExp = geneExp)
    df1 = beh$.beh$.adjust.pvals(df, phen)
    out = df1
    ## out = data.table(combined =enrichmentTesting$fisherCombined(df$gene.pval)$chisq.analytic.p,
    ##                  combined.lowp = enrichmentTesting$fisherCombined(df$gene.pval[df.orig$strain.pval<.05])$chisq.analytic.p,
    ##                  combined.lowq = enrichmentTesting$fisherCombined(df$gene.pval[df1.orig$strain.pval.qval.fdr<.05])$chisq.analytic.p)
    ## out$geneName = geneExp$name

    ##out = enrichmentTesting$fisherCombined(df$gene.pval)$chisq.analytic.p
    return(out)
}

## covariates = " ~ 1 +  as.factor(Batch)  + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID)"
mediation$regressResiduals.exp.vs.beh <- function(geneResid.mat,
                                                  probesetInfo,
                                                  allPhen = readPhens(),
                                                  covariates = " ~ 1 +  as.factor(Batch)  + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID)",
                                                  geneName   = "Lrrc16a",
                                                  experiment = "light.dark",
                                                  pheno      = "Total.Distance",
                                                  lambdas    = prop$mnp$lambdasToTry)
{
     psid  = getProbesetId(probesetInfo, geneName)

     dataSet = allPhen$frame[[experiment]]
     
     fit.lambda = getBestLambda(lambdas = lambdas,
                                        pheno = pheno,
                                        covariateModelString = covariates,
                                        dataSet = dataSet,
                                        uselme = F,
                                        normalizeBeforeTransform = T,
                                        normalizeAfterTransform = T,
                                        checkAnova = F)
     beh = resid(fit.lambda$fit)
     names(beh) = paste0("Mouse.", dataSet$ID)
     beh = beh[names(beh) %in% rownames(geneResid.mat)]
     geneExp = exp.mat[names(beh), psid]
     fit  = lm(geneExp, beh)
     wrap = lm.parsing$getTsWrapper(fit)
     p.value = wrap$ts["beh",wrap$pvalueCol]
     return(p.value)
}



  ## ##TODO bring back
  ##   if(!is.null(behaviorExpCors))
  ##   {
  ##       ##merge in annotations info TODO get rid of this if we merge probesetinfo and annot.data
  ##       annot.descrip   = getAnnotDescriptions(annot.data)
  ##       annot.descrip   = data.table(annot.descrip, key ="Probe.Set.ID")
  ##       behaviorExpCors = behaviorExpCors[annot.descrip]
  ##       write.table(behaviorExpCors, file = outm("behaviorExp.csv"), row.names=F)
  ##       plotBehaviorExpCors(behaviorExpCors, prop$mnp$output)
  ##       print("done correlating with behavior!")
  ##   }

mediation$.runBehaviorCorAnalysis <- function(exp.mat.full, probesetInfo, annot.data)
{
    pracma::tic()
    print("starting correlation with behavior")
    ##read in phenotypes
    phens = readPhens()
    ##get the behavior correlations data table
    behaviorExpCors = mediation$.Calc.Cors.wrapper(phens, exp.mat.full)
    behaviorExpCors = data.table(behaviorExpCors, key="Probe.Set.ID")
    
    ##merge in probeset info
    behaviorExpCors = probesetInfo[behaviorExpCors]
    pracma::toc()

    ##merge in annotations info TODO get rid of this if we merge probesetinfo and annot.data
    return(behaviorExpCors)
}

mediation$.Calc.Cors <- function(phen, phen.data, exp.mat, annot.data) 
{
    x <- phen.data[[phen]]
    ok <- !is.na(x)
    x <- x[ok]
    y.all = exp.mat[ok,]
    cat(sep="", "\n", phen, ":\n")

    lfunc <- function(indexes, phen, phen.data, exp.mat, annot.data)
    {

        d <- data.frame(Probe.Set.ID = annot.data$Probe.Set.ID[indexes])
        d$Est.Cor       <- NA
        d$LogP.Cor      <- NA
        d$Est.SpearCor  <- NA
        d$LogP.SpearCor <- NA
        d$Est.KendCor   <- NA
        d$LogP.KendCor  <- NA
        d$phen          <- phen
        
        for (i in 1:length(indexes))
        {
            
            Probe.Set.ID = as.character(d$Probe.Set.ID[i])
            index = indexes[i]
            if (index %% 1000 == 0) { cat(sep="","[",index,"]") } 
            y <- y.all[, index]
            if(any(y!=y.all[,Probe.Set.ID]))
            {
                
            }
            
            ct <- cor.test(x,y, method=c("pearson"))
            d$Est.Cor[i] <- ct$estimate
            d$LogP.Cor[i] <- -log10(ct$p.value)
            
            ct <- suppressWarnings(cor.test(x,y, method=c("spearman"), exact=TRUE))
            d$Est.SpearCor[i] <- ct$estimate
            d$LogP.SpearCor[i] <- -log10(ct$p.value)
            
            ct <- suppressWarnings(cor.test(x,y, method=c("kendall"), exact=TRUE))
            d$Est.KendCor[i] <- ct$estimate
            d$LogP.KendCor[i] <- -log10(ct$p.value)
        }
        return(data.table(d))
    }

    alist = mnp.lapply.wrapper(len = ncol(exp.mat), FUN = lfunc, phen = phen, phen.data = phen.data, exp.mat = exp.mat, annot.data = annot.data)
    df = rbindlist(alist)
    return(df)
}

mediation$.Calc.Cors.wrapper <- function(phens, exp.mat.full)
{
    merged = mergeIntoPipelines(phens)
    cor.frames = list()
    for(pipelind in 1:2)
    {
        pipel            = paste0("pipel",pipelind)
        phenz            = merged[[pipel]]
        phenz            = phenz[,lapply(.SD,mean),by="ID"]
        phen.cols        = setdiff(colnames(phenz), "ID")
        exp.mat          = exp.mat.full[rownames(exp.mat.full) %in% paste0("Mouse.",phenz$ID),]
        phenz            = phenz[J(ID = unlist(lapply(strsplit(rownames(exp.mat),"\\."),"[", 2)))]
        
        ##cor.list <- mclapply(phen.cols, Calc.Cors, phen.data = phenz, exp.mat=exp.mat, annot.data=annot.data,mc.cores=mc.cores)
        cor.list <-   lapply(phen.cols, Calc.Cors, phen.data = phenz, exp.mat=exp.mat, annot.data=annot.data)
        cor.frame = do.call(rbind, cor.list)
        cor.frame$pipeline.LogQ.Cor      = -log10(p.adjust(10^-(cor.frame$LogP.Cor), method="fdr"))
        cor.frame$pipeline.LogQ.SpearCor = -log10(p.adjust(10^-(cor.frame$LogP.SpearCor), method="fdr"))
        cor.frame$pipel = pipelind
        cor.frames[[pipelind]] = cor.frame
    }

    cor.frame = data.table(do.call(rbind, cor.frames))
    cor.frame$Probe.Set.ID = as.character(cor.frame$Probe.Set.ID)

    return(cor.frame)
}



##TODO refactor
qpcr.resid <- function()
{
    followup = read.table(datm("followup/Matnut_pilot_alltaqmanplates_datacompletev3_100716.csv"), sep=",", header=T)
    followup = data.table(followup)
    toExamine = "Lrrc16a"
    ##toExamine = "Meg3"
    followup = followup[followup$Assay == toExamine]
    followup = followup[!is.na(FAM.Ct)]

    setnames(followup, old = c("FAM.Ct", "VIC.Ct"), new = c(paste0("taqman.", toExamine), "taqman.Rfng"))

    followup$ID = paste0("Mouse.",followup$ID)
    setkey(followup, "ID")

    cov.data.2 = data.table(cov.data, key = "ID")

    followup = cov.data.2[followup]
    followup = followup[Strain %in% c("NOD.B6", "B6.NOD")]
    followup$batchplate = as.factor(paste(followup$Batch, followup$Plate, sep="_"))


    phen.control = "Delta.Ct"
    m.string = " ~ 1  + Pipeline + Diet + Strain + Strain:Diet + (1|DamID) + (1|batchplate)"
    lambdasToTry = c(-3, -2, -1, -.5, 0, .5, 1, 2, 3)
    m.phen         = callFit(lambdasToTry, followup, phen.control, m.string)
    resids         = data.table(res = resid(m.phen$fit), ID  = followup$ID)
    resids = resids[j=list(res = mean(res)), by = "ID"]
    rownames(resids) = resids$ID
    geneExp = resids$res
    names(geneExp) = resids$ID
    return(geneExp)
}

mediation$runAnalysis <- function(inp)
{
    mediation$.model.SEM.behaviorWithMicroarray(inp)
    mediation$.regressResiduals.exp.vs.beh()
    ## TURN This off in testing as it takes a while, and generally don't find anything.
    behaviorExpCors = NULL
    if(prop$mnp$computeBehaviorCor)
    {
        behaviorExpCors = mediation$.runBehaviorCorAnalysis(exp.mat.full, probesetInfo, annot.data)
    }
}
