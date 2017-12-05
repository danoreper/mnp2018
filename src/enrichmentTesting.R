enrichmentTesting = new.env(hash=T)

enrichmentTesting$fisherCombined <- function(pvals)
{
    chisq.stat.part    = -2*log(pvals)
    chisq.stat         = sum(chisq.stat.part)
    chisq.analytic.p   = 1 - pchisq(chisq.stat, df=2*length(pvals))
    df = list(
      #chisq.stat.part = chisq.stat.part,
                    chisq.stat = chisq.stat,
                    chisq.analytic.p = chisq.analytic.p)
    
    return(df)
}


## trueSubset-- vector of booleans indicating which element of pvals together form a subset
##permSubset-- matrix of booleans indicating a premutation chosen set of boolean denoting subset
##membership. As many rows as there are permutations, as many cols as there are pvals
enrichmentTesting$fisherCombinedWrapper <- function(pvals, trueSubset, permSubset)
{
    true.stat= fisherCombined(pvals[trueSubset])
    perm.stat = list()
    for(i in 1:ncol(permSubset))
    {
        permedPvals    = pvals[permSubset[,i]]
        perm.stat[[i]] = fisherCombined(permedPvals)
        ##qqplot(x = qunif(ppoints(length(permedPvals))), y = permedPvals)
    }

    numGreater = sum(true.stat$chisq.stat < unlist(lapply(perm.stat , "[[", "chisq.stat")))
    df = data.frame(analytic.p = true.stat$chisq.analytic.p,
                    perm.p     = numGreater/length(perm.stat))
    #hist(chisq.stat.perm)
    return(df)
}
