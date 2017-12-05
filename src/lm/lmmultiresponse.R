SS<-function(x)
{
  cov(x,x)*(length(x)-1)
}

lm.multiresponse <- function(formula, response.matrix, data,
        null.formula  = NULL,
        null.fit      = NULL,
        rsquared      = FALSE,
        pvalue        = FALSE,
        logP          = FALSE,
		    LOD			      = FALSE, 
		    total.ss      = FALSE,
		    weights       = rep(1, nrow(data)),
        model.args    = list(),
        verbose.at.every = 0)
# fits a linear model with a constant design matrix to a matrix of
# responses. Return F-tests and pvalues for comparison with a null model
# specified either as a formula or a fitted lm object.
{
  if (any(1!=weights) & (rsquared | total.ss | any(0==weights)))
  {
    # note: the present function gets tss and all derived statistics wrong
    # and will probably choke if any(weights==0)
    # However, it is several times faster than wlm and does get rss correct, which means that it also
    # gets logP and LOD correct.
    # Nonetheless, wlm is preferred for predictability until lm.multi has its weighted tss fixed.
    return (
        wlm.multiresponse(formula, response.matrix, data,
                null.formula  = null.formula,
                rsquared      = rsquared,
                pvalue        = pvalue,
                logP          = logP,
        		    LOD			      = LOD, 
        		    weights       = weights,
                model.args    = model.args,
                verbose.at.every = verbose.at.every)
        )
  }
  
  formula         <- as.formula(formula)
  response.matrix <- as.matrix(response.matrix)

  # check response is univariate
  terms.object <- terms(formula)
  if (0==attr(terms.object, "response"))
  {
    stop("Must specify response in formula\n")
  }
  if (1!=attr(terms.object, "response")
        | length(all.vars(terms.object))!=nrow(attr(terms.object, "factors")))
  {
    stop("Multivariate response not allowed\n")
  }

  # transform response if requested
  response.name <- all.vars(terms.object)[1]
  response.expr <- rownames(attr(terms.object, "factors"))[1]
  if (response.name != response.expr)
  {
    FUN2 <- eval(parse(text=paste(
        "FUN <- function(", response.name,"){",
                response.expr,
        "}",
        sep="")))
    response.matrix <- apply(response.matrix, 2, FUN2)
  }
  if (!all(is.finite(response.matrix)))
  {
    stop("Response must be finite\n")
  }

  # fit model to all responses
  data[,response.name] <- response.matrix[,1]
  fit <- do.call("lm", args=c(
      model.args,
      list(
          formula=as.formula(formula),
          data=quote(data),
          weights=weights
          )
      ))
  qr  <- fit$qr
  tss <- numeric(ncol(response.matrix))
  rss <- numeric(ncol(response.matrix))

  rootw <- sqrt(weights)
  for (i in 1:ncol(response.matrix))
  {
    # show progress
    if (0!=verbose.at.every)
    {
        if (0 == i %% verbose.at.every)
        {
            cat(i,"/",ncol(response.matrix),"\n")
        }
    }
    # rate limiting step
    y <- response.matrix[,i]
    tss[i] <- sum((y-mean(y))^2*weights)  # incorrect when any(1!=weights), but close
#    rss[i] <- sum( qr.resid(qr, y) ^ 2 ) # unweighted version for speed
    rss[i] <- sum( qr.resid(qr, rootw*y) ^ 2 ) # always correct
  }

  retval <- list(
  		n           = length(resid(fit)),
      rss         = rss,
      rank        = fit$rank,
      df.residual = fit$df.residual)

  if (total.ss | TRUE) retval$total.ss <- tss
  if (rsquared) retval$rsquared <- (tss - rss) / tss
  

  #------------------------------
  # Fit null model for comparison
  #------------------------------
  if (!is.null(null.fit) | !is.null(null.formula))
  {
    if (is.null(null.formula))
    {
      formula.as.string(null.fit$terms)
    }
    if (split.formula(null.formula)$response!=split.formula(null.formula)$response)
    {
      stop("Response expression in formula and null formula differ: ",
          split.formula(null.formula)$response, " vs ",
          split.formula(formula)$response,
          "\n")
    }
    if (is.null(null.fit))
    {
      null.fit <- do.call("lm", args=c(
          model.args,
          list(
          formula=as.formula(null.formula),
          data=quote(data),
          weights=weights
          )
      ))
    }

    # get statistics from null fit
    qr0  <- null.fit$qr
    rss0 <- numeric(ncol(response.matrix))
    for (i in 1:ncol(response.matrix))
    {
      y <- response.matrix[,i]
      rss0[i] <- sum( qr.resid(qr0, rootw*y) ^ 2 ) # correct
    }
    retval$null.rss  <- rss0
    retval$null.rank <- null.fit$rank

    if (pvalue | logP)
    {
     # calculate F-test for comparison of models
      dfr <- fit$df.residual
      delta.dfp  <- fit$rank - null.fit$rank

      fss <- rss0 - rss
      f   <- fss / rss * dfr / delta.dfp

      if (pvalue | logP)
      {
        pval <- pf(f, delta.dfp, dfr, lower.tail=F)
        if (pvalue) retval$pvalue <- pval
        if (logP)   retval$logP   <- -log10(pval)
      }
    }
  	if (LOD)
  	{
  		retval$LOD <- (retval$n/2) *(log10(rss0) - log10(rss))
  	}
  }
  retval
}

wlm.multiresponse <- function(formula, response.matrix, data,
        null.formula  = NULL,
        rsquared      = FALSE,
        pvalue        = FALSE,
        logP          = FALSE,
		    LOD			      = FALSE, 
		    weights       = rep(1, nrow(data)),
        model.args    = list(),
        verbose.at.every = 0)
# fits a linear model with a constant design matrix to a matrix of
# responses. Return F-tests and pvalues for comparison with a null model
# specified either as a formula or a fitted lm object.
{
  formula         <- as.formula(formula)
  response.matrix <- as.matrix(response.matrix)

  # check response is univariate
  terms.object <- terms(formula)
  if (0==attr(terms.object, "response"))
  {
    stop("Must specify response in formula\n")
  }
  if (1!=attr(terms.object, "response")
        | length(all.vars(terms.object))!=nrow(attr(terms.object, "factors")))
  {
    stop("Multivariate response not allowed\n")
  }

  # transform response if requested
  response.name <- all.vars(terms.object)[1]
  response.expr <- rownames(attr(terms.object, "factors"))[1]
  if (response.name != response.expr)
  {
    FUN2 <- eval(parse(text=paste(
        "FUN <- function(", response.name,"){",
                response.expr,
        "}",
        sep="")))
    response.matrix <- apply(response.matrix, 2, FUN2)
  }
  if (!all(is.finite(response.matrix)))
  {
    stop("Response must be finite\n")
  }

  mlm.formula = paste("response.matrix ~", split.formula(formula)$predictor.string)
  mlm.fit = do.call("lm", args=c(
      model.args,
      list(
      formula=as.formula(mlm.formula),
      data=quote(data),
      weights=weights
      )))
  mlm.sum =summary(mlm.fit)

  # get statistics from fit
  rss = colSums(sapply(mlm.sum, function(x){x$residuals})^2)
  r2  = sapply(mlm.sum, function(x){x$r.squared})
  tss = rss/(1-r2)

  result <- list(
  		n           = NROW(resid(mlm.fit)),
      tss         = tss,
      rss         = rss,
      rank        = mlm.fit$rank,
      df.residual = mlm.fit$df.residual)

  if (rsquared) result$rsquared <- r2

  #------------------------------
  # Fit null model for comparison
  #------------------------------
  if (!is.null(null.formula))
  {
    if (split.formula(formula)$response!=split.formula(null.formula)$response)
    {
      stop("Response expression in formula and null formula differ: ",
          split.formula(null.formula)$response, " vs ",
          split.formula(formula)$response,
          "\n")
    }
    mlm.formula0 = paste("response.matrix ~", split.formula(null.formula)$predictor.string)
    mlm.fit0 <- do.call("lm", args=c(
        model.args,
        list(
        formula=as.formula(mlm.formula0),
        data=quote(data),
        weights=weights
        )
    ))
    mlm.sum0=summary(mlm.fit0)
    
    # get statistics from null fit
    result$null.rss = colSums(sapply(mlm.sum0, function(x){x$residuals})^2)
    result$null.rank = mlm.fit0$rank
    if (pvalue | logP)
    {
     # calculate F-test for comparison of models
      dfr <- result$df.residual
      delta.dfp  <- result$rank - result$null.rank

      fss <- result$null.rss - result$rss
      f   <- fss / rss * dfr / delta.dfp

      if (pvalue | logP)
      {
        pval <- pf(f, delta.dfp, dfr, lower.tail=F)
        if (pvalue) result$pvalue <- pval
        if (logP)   result$logP   <- -log10(pval)
      }
    }
  	if (LOD)
  	{
  		result$LOD <- (result$n/2) *(log10(result$null.rss) - log10(rss))
  	}
  }
  result
}
