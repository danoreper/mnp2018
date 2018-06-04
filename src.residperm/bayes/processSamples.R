processSamples = new.env(hash=T)

processSamples$get.p.value <- function(y)
{
    pval = sum(y<0)/length(y)
    if(pval<1-pval)
    {
        direction = 1
    } else {
        pval = min(1-pval, pval)
        direction = -1
    }

    return(data.frame(pval = pval, direction = direction))
}

processSamples$getRegexName <- function(effectName)
{
    regexEffectName = gsub(effectName, pattern = "\\.", replacement = "\\\\.")
    return(regexEffectName)
}

processSamples$getShortName <- function(effectName, originalName)
{
    regexEffectName = processSamples$getRegexName(effectName)
    out = gsub(originalName, pattern = paste0(regexEffectName,"|\\[|\\]"), replacement = "")
    return(out)
}

processSamples$getColsForEffect <- function(gibbs.samples, effectName)
{
    regexEffectName = gsub(effectName, pattern = "\\.", replacement = "\\\\.")
    colstring = paste0(regexEffectName,"\\[.*\\]")
    colsForEffect = grepl(pattern = colstring, perl=T, x = colnames(gibbs.samples))
    colsForEffect = colnames(gibbs.samples)[colsForEffect]
    return(colsForEffect)
}

processSamples$getComparisonFrame <- function(gibbs.samples, effectName, shortname = F)
{
    compareCols = processSamples$getColsForEffect(gibbs.samples, effectName)
    dfs = list()

    c1 = paste0(effectName, ".1")
    c2 = paste0(effectName, ".2")

    for(i in 1:length(compareCols)) ##1:(length(compareCols)-1))
    {
        col.i = compareCols[i]
        for(j in 1:length(compareCols))##(i+1):length(compareCols))
        {
            col.j = compareCols[j]

            browser()
            colvec = unname(gibbs.samples[,col.j] - gibbs.samples[,col.i])
            ##hpd = HPDinterval(mcmc(colvec))
            if(i==j)
            {
                pval = data.frame(pval=1, direction=0)
            } else {
                pval = processSamples$get.p.value(colvec)
            }
            df = pval
            df[[c1]] = col.i
            df[[c2]] = col.j
            
            dfs = util$appendToList(dfs, df)
        }
    }
    
    dfs = do.call(rbind, dfs)
    if(shortname)
    {
        dfs[[c1]] = processSamples$getShortName(effectName, dfs[[c1]])
        dfs[[c2]] = processSamples$getShortName(effectName, dfs[[c2]])
    }
    return(dfs)
}

processSamples$getSingleFrame <- function(gibbs.samples, effectName, shortname=F )
{

    cols = processSamples$getColsForEffect(gibbs.samples, effectName)
    dfs = list()
    for(col.i in cols)
    {
        pval = processSamples$get.p.value(gibbs.samples[,col.i])
        df = pval
        df[[effectName]] = col.i
        dfs = util$appendToList(dfs, df)
    }
    dfs = do.call(rbind, dfs)
    if(shortname)
    {
        dfs[[effectName]] = processSamples$getShortName(effectName, dfs[[effectName]])
    }
    return(dfs)
}


processSamples$computeCriticalIntervalInfo <- function(gibbs.samples)
{
    asum = summary(gibbs.samples)
    aqt   = (asum$quantiles)
    ast   = (asum$statistics)

    dfs = list()
    for(i in 1:nrow(aqt))
    {
        df = data.frame(effect.name = rownames(aqt)[i],
                        mean = ast[i,"Mean"],
                        l2.5 = aqt[i, "2.5%"],
                        u97.5= aqt[i, "97.5%"])
        
        dfs = util$appendToList(dfs, df)
    }
    dfs = do.call(rbind, dfs)
    return(dfs)

}

processSamples$getSignificantVars <- function(mcobj)
{
    hpd = coda::HPDinterval(mcobj)
    sigterms = (hpd[(hpd[,1]*hpd[,2])>0,])
    sigall = mcobj[,rownames(sigterms)]
    return(sigall)
}
