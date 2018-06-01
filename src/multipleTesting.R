library(evir)
multipleTesting = new.env(hash=T)

## multipleTesting$median.correct <- function(pvals)
## {
##     zscores        = qnorm(pvals)
##     med            = median(zscores)
##     zscores        = zscores - med
##     pvals.adjusted = pnorm(zscores) 
##     return(pvals.adjusted)
## }

multipleTesting$median.correct.2 <- function(pvals)
{
##    browser()
    zscores        = qchisq(pvals, df=1)
    med            = median(zscores)/.4549
    zscores        = zscores/med
    pvals.adjusted = pchisq(zscores,df=1) 
    return(pvals.adjusted)
}

multipleTesting$minp.gev <- function(pvals, alphas)
{

    pvals = na.omit(pvals)
    evd.pars <- gev(-log10(pvals))$par.ests
    gev.thresh <- qgev(p=(1-alphas), xi=evd.pars["xi"], sigma=evd.pars["sigma"], mu=evd.pars["mu"])
    gev.thresh = 10^(-gev.thresh)
    return(gev.thresh)
}
    

multipleTesting$get.empirical.p.value.for.q <- function(pvals, qvals, alphalevel)
{
    tokeep = qvals<alphalevel
    if(sum(tokeep)>0)
    {
        pval = pvals[tokeep]
        qval = qvals[tokeep]
        out = pval[which.max(qval)]
    } else {
        return(NA)
    }
}

multipleTesting$get.empirical.one.tail.p.value <- function(nulldist, stat, direction)
{
    nulldist = nulldist*direction
    out = sum(stat>nulldist)
    return(out)
}

multipleTesting$get.empirical.one.tail.thresh.foralpha <- function(nulldist, alpha=.05, direction)
{
    nulldist = nulldist*direction
    
    out = quantile( probs = 1-alpha)
    return(out)
}


## TODO: direction here is backwards... redo. also rename to avoid implication that these have to be p.values, it can be any arbitrary statistic.
##perm.p.value is a matrix represnting multiple tests under permutation of labels, where perm.p.value[i,j] corresponds to the p-value computed for the ith permutation of labels, and for the jth phenotype.
##noperm.p.value is a vector of length equal to the number of phenotype which corresponds to the pvalues computed with
## the original labels and without permutation.
multipleTesting$evaluate.perms <- function(perm.p.value, alpha, direction)
{
    perm.p.value   = perm.p.value * direction
    ##noperm.p.value = noperm.p.value * direction

    ## filter out failing permutations
    badPerms = apply(perm.p.value, 1, function(x){all(is.na(x))})
    if(any(badPerms))
    {
        warning(paste0("some permutations failed: ", paste(which(badPerms), collapse = ",")))
        perm.p.value = perm.p.value[!badPerms,]
    }

    local.ecdf = apply(perm.p.value, 2, ecdf)
    
    out = list()
    
    ## out$empirical.local.p.value  = unlist(lapply(1:length(local.ecdf),
    ##                                              FUN = function(j, local.ecdf, noperm.p.value){local.ecdf[[j]](noperm.p.value[j])},
    ##                                              local.ecdf, noperm.p.value))
    
    minp = apply(perm.p.value, 1, min, na.rm = T)
    global.ecdf = ecdf(minp)
    ## out$empirical.global.p.value =
    ##     unlist(lapply(noperm.p.value,
    ##                   FUN = function(noperm.p.value, global.ecdf){global.ecdf(noperm.p.value)}, global.ecdf))
     
    out$empirical.global.threshhold = unname(quantile(minp, alpha, na.rm =T))
    out$empirical.global.gev.threshhold = try( multipleTesting$minp.gev(minp, alpha), silent =T)
    if(class(out$empirical.global.gev.threshhold)=="try-error")
    {
        out$empirical.global.gev.threshhold = NA
    }
    
    out$alpha = alpha

    ## out$perm.p.value = perm.p.value
    ## out$noperm.p.value = noperm.p.value 
    return(out)
}

##Y is a matrix of outcomes, each outcome is a column
multipleTesting$get.perm.corr.p.value <- function(Y, x, shuffles=NULL, numPerm=NULL, alpha)
{
    noperm.p.value = multivariate.lm$getPvalues(Y,x)
    perm.p.value = matrix(NA, ncol = ncol(Y), nrow =numPerm)
    for(shuffle in 1:numPerm)
    {
        xPerm = x[shuffles[,shuffle]]
        perm.p.value[shuffle,] = multivariate.lm$getPvalues(Y,xPerm)
    }
    out = multipleTesting$evaluate.perms(perm.p.value, noperm.p.value, alpha)
    return(out)
}

