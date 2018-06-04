library(lmerTest)

S = 50
J = 5
A = 2
a = rbind(rnorm(J), rnorm(J))
a[1,1] = 10#10*a[1,1]  

B = 4
b = rbind(rnorm(J), rnorm(J), rnorm(J), rnorm(J)) 

n = 100

lsize = 4
ndam = ceiling(n/lsize)


d.i = factor(1:ndam)

a.i = factor(sample(c(rep(1,ceiling(ndam/A)), rep(2,ceiling(ndam/A))))[1:ndam])

b.i = factor(sample(c(rep(1,ceiling(ndam/B)),
               rep(2,ceiling(ndam/B)),
               rep(3,ceiling(ndam/B)),
               rep(4,ceiling(ndam/B))))[1:ndam])

sigma.j = runif(J, 0,1)
#sigma.j[1] = .05
sigma.j[] = .5


df = data.frame(a.x = rep(a.i, each=lsize),
                b.x = rep(b.i, each=lsize),
                d.x = rep(d.i, each=lsize))

ys       = matrix(nrow = n, ncol = J)
true.res = matrix(nrow = n, ncol = J)

for(j in 1:J)
{
    i = 1
    for(d in 1:ndam)
    {
        for(pup in 1:lsize)
        {
            true.res[i,j] = rnorm(1, 0, sigma.j)
            ys[i,j] = a[a.i[d],j] + b[b.i[d],j] + true.res[i,j]
            i = i + 1
        }
    }
}








getpval <- function(df)
{
    
    ## fit = lm(y ~ b.x + a.x, df)
    ## an = data.frame(anova(fit))
    ## pval = an["a.x", "Pr..F."]
    ## return(list(pval=pval, fit=fit, resids = resid(fit)))

    fit = lmerTest::lmer(formula = y ~ b.x + a.x + (1|d.x), data = df)
    an = anova(fit, type = 1)
    pval = an["a.x", "Pr(>F)"]
    return(list(pval=pval, fit = fit))
}

getredpval <- function(df)
{
    ##fit = lm(y ~ b.x , df)
    fit = lmerTest::lmer(y ~ b.x + (1|d.x), df)
    return(list(fit=fit, resids = resid(fit)))
}


y.mu   = matrix(nrow = n, ncol = J)
resids = matrix(nrow = n, ncol = J)
dms    = matrix(nrow = ndam, ncol = J)

for(j in 1:J)
{
    ##browser()
    df$y    = ys[,j]
    afit    = getredpval(df)
    ress = afit$resids
    browser()
    rownames(dms) = rownames(ranef(afit$fit)$d.x)
    dm   = ranef(afit$fit)$d.x[["(Intercept)"]]
    dms[,j] = dm
    ##browser()
    
    y.mu[,j] = df$y - afit$resids - dm[as.integer(df$d.x)]
    resids[,j] = afit$resids
}

shuffles = matrix(nrow = S, ncol = n)
shuffles[1,] = 1:n

shuffles.dam = matrix(nrow = S, ncol = ndam)
shuffles.dam[1,] = 1:ndam
for(s in 2:S)
{
    shuffles[s,] = sample(1:n)
    shuffles.dam[s,] = sample(1:ndam)
}


pvals = matrix(nrow = S, ncol = J)
for(s in 1:S)
{
    print(s)
    for(j in 1:J)
    {
        ##browser()
        df$y = y.mu[,j]+resids[shuffles[s,],j] + rep(dms[shuffles.dam[s,],j], each = lsize)
        afit    = getpval(df)
        pvals[s,j] = afit$pval
    }
}


minps = (apply(pvals, 1, min)[-1])
thresh = quantile(minps, .05)

print(sum((p.adjust(pvals[1,], method = "fdr")<.05) & sum(pvals[1,]>thresh)))

## pvals.0 = rep(NA, J)
## for(j in 1:J)
## {
##     df$y    = y[,j]
##     afit    = getpval(df)
##     pvals.0[j] = afit$pval
## }
