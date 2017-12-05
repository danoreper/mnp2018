multivariate.lm = new.env(hash=T)

multivariate.lm$beta <- function(X,Y)
{
    beta    = solve(t(X)     %*% X)      %*% (t(X)%*%Y)
    ## if(class(beta)=="try-error")
    ## {
    ##     browser()
    ## }
    return(beta)
}

multivariate.lm$residuals <- function(X, Y)
{
    R = try(Y - X %*% multivariate.lm$beta(X,Y))
    if(class(R)=="try-error")
    {

        R = matrix(nrow = nrow(Y), ncol = ncol(Y))
        qr = lm.fit(X, Y[,1])$qr
        for(i in 1:ncol(Y))
        {
            y = Y[,i]
            R[,i] = qr.resid(qr, y)
        }
    }
    return(R)
}
multivariate.lm$computeRSS <- function(X, Y)
{
    R = multivariate.lm$residuals(X, Y)
    RSS = colSums(R^2)
    return(RSS)
}

multivariate.lm$batchPvals <- function(Y, X, X.null)
{
    RSS2 = multivariate.lm$computeRSS(X,Y)
    RSS1 = multivariate.lm$computeRSS(X.null, Y)
    
    p2 = ncol(X)
    p1 = ncol(X.null)
    n = nrow(Y)
    df1 = p2 - p1
    df2 = n - p2 
    F = ((RSS1 - RSS2)/(df1))/(RSS2/(df2))
    ##print(paste0(F, ",", RSS2, ",", RSS1))
    pValues = 1 -pf(F,df1 = df1, df2 = df2)
    return(pValues)
}

multivariate.lm$get.p.values <- function(Y,cov.data, fullmodel, nullmodel)
{
    cov.data$y =Y[,1]
    formla         = as.formula(paste0("y ", fullmodel))
    formla.null    = as.formula(paste0("y ", nullmodel))
    X              = model.matrix(formla,      data = cov.data)
    X.null         = model.matrix(formla.null, data = cov.data)
    parsedfit      = multivariate.lm$batchPvals(Y=Y, X, X.null)
    return(parsedfit)
    
}

multivariate.lm$batch.shapiro <- function(Y,cov.data,model)
{
    cov.data$y =Y[,1]
    formla         = as.formula(paste0("y ", model))
    X              = model.matrix(formla, data =cov.data)
    R = multivariate.lm$residuals(X, Y)
    ##TODO-- is there a faster hand written way to do this?
    out = apply(R, 2, function(r) { shapiro.test(r)$p.value})
    return(out)
}

multivariate.lm$getPvalues <- function(Y,x)
{
    y = Y[,1]
    afit          = lm(y~ 1+x)
    afit.null     = lm(y~ 1)
    X      = model.matrix(afit)
    X.null = model.matrix(afit.null)
    
    ##Comment out when not debugging
    ##lmfit = lm(y~x)
    ##an = anova(lm (y~x))
        ##print(an["x","Pr(>F)"])
    
    pValues = multivariate.lm$batchPvals(Y, X, X.null)
    return(pValues)
}


