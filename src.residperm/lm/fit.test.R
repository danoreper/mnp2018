source("./lm/fitBoxCoxModels.R")

nsamp = 200
x1 = 1:nsamp
x2 = sample(x1, nsamp, replace =T)
x3 = factor(sample(1:5,nsamp, replace=T))


tp = fit.model.bc$getDefaultTransformParams()
tp$lambdasToTry = c(-3,-2,-1,-.5, -.3333, 0, .3333, .5, 1,2,3) 
tp$normalizeBeforeTransform = F
tp$normalizeAfterTransform = F
##tp = NULL

df = data.frame(x1 = x1, x2 = x2, x3 = x3)
df$y.system = (1010 + 5*x1 + -4*x2 + as.integer(x3))
df$err = rnorm(length(x1))

lamb = 2
df$y = "^"(((df$y.system+df$err)*(1/lamb)),(1/lamb))

fit = fit.model.bc$fit(df$y, df, covariateModelString = "~ x1 + x2 + x3", transformParams = tp)

z1 = df$y


z2 = fit.model.bc$transform.by.lambda(z1, lambda = 2, normalizeBeforeTransform = T, normalizeAfterTransform=T)
z1.p = fit.model.bc$invert(z2)

modelParser = new.env(hash=T)

modelParser$parse=function(fit)
{
    pval = fit$anovaWrapper$an["x3", "Pr(>F)"]
    return(pval)
}

modelParser$collate=function(parsedouts)
{
    return(parsedouts)
}

pvalz = fit.model.bc$fit(df$y, df, covariateModelString = "~ x1 + x2 + x3", transformParams = tp, modelParser=modelParser)
