##Function for computing surrogate variables, as in Leek and Storey, 2007

surrogatcalc = new.env(hash=T)


surrogatcalc$refine.ssva <- function(exp.ssva.mat, num.sv, scaleOnForPCA)
{
#    print("refining!")
    ##Step 1.2
    if(num.sv == 0)
    {
        return(matrix(0,0,0))
    }

    exp.resid.pca <- prcomp(exp.ssva.mat, scale.=scaleOnForPCA)
    ##Modified form of step 1.3-1.10. N.B. we are fitting a mixed model, so permutation doesnt make
    ##sense as a way to pick the number of principal components. Instead, we just eyeball the 
    ##eigenvalues and see where the inflection point is.

    ##TODO bring this back
    ##plot(exp.resid.pca$sdev^2, pch=as.character(c(1:9,0))) # suggests 7 pcs
    num.mice <- nrow(exp.ssva.mat)
    sv.mat <- matrix(nrow=num.mice, ncol=num.sv) # holds refined sv's
    for (ii in 1:num.sv) #Step 2.3 is the result of 1.0-1.10 
    {
        ##step 2.2 (step 2.1 is cached from 1.1)
        rough.sv      <- exp.resid.pca$x[,ii]
        sv.mat[ , ii] <- rough.sv#refined.pca$x[ , indx]
    }
    rownames(sv.mat) = rownames(exp.ssva.mat)
    return(sv.mat)
}


surrogatcalc$refine.sv <- function(exp.mat, exp.resid.mat, num.sv, topgenes, scaleOnForPCA)
{

    if(num.sv == 0)
    {
        return(matrix(0,0,0))
    }  
    ##Step 1.2
    exp.resid.pca <- prcomp(exp.resid.mat, scale.=scaleOnForPCA)
    ##Modified form of step 1.3-1.10. N.B. we are fitting a mixed model, so permutation doesnt make
    ##sense as a way to pick the number of principal components. Instead, we just eyeball the 
    ##eigenvalues and see where the inflection point is.
    plot(exp.resid.pca$sdev^2, pch=as.character(c(1:9,0))) # suggests 7 pcs
    num.mice <- nrow(exp.resid.mat)
    
    sv.mat <- matrix(nrow=num.mice, ncol=num.sv) # holds refined sv's
    for (ii in 1:num.sv) #Step 2.3 is the result of 1.0-1.10 
    {
        ##step 2.2 (step 2.1 is cached from 1.1)
        rough.sv      <- exp.resid.pca$x[,ii]
        ##getting the highest correlated genes to the eigengene seems to be equivalent,
        ##but just in case, following the exact text and computing a p value for every linear model.
        ##		gene.cors     <- apply(exp.mat, 2, function(y){ abs(cor(rough.sv, y)) })
        ##step 2.4
        pval = rep(0, ncol(exp.mat))
        for(m in 1:ncol(exp.mat))
        {
            mthGeneExpresion =  exp.mat[,m]
            fit = lm(rough.sv ~ mthGeneExpresion )
            pval[m] = coefficients(summary(fit))["mthGeneExpresion","Pr(>|t|)"]
        }
		
        ##modified form of step 2.5 TODO look into this
        top.thou.jj   <- order(pval, decreasing=F)[1:topgenes]
        ##		top.thou.jj   <- order(gene.cors, decreasing=TRUE)[1:topgenes]
        ##		
        ##step 2.6
        refined.pca   <- prcomp(exp.mat[ , top.thou.jj], scale.=scaleOnForPCA)
        
        ##we want to find the eigengene of the reduced gene set, that is most correlated with
        ## the original eigengene step 2.7 Note that we use the absolute value of correlation,
        ## since flipping a principal component positive or nega5tive should account for the
        ## same amount of variance, it will just change the regression coefficient for that
        ## component. i.e. corr = -.95 is a better pc than corr = .02
        pca.cors      <- apply(refined.pca$x, 2, function(y){ abs(cor(rough.sv, y)) })
        indx          = which.max(pca.cors)
        print(indx)
        ##		sv.mat[ , ii] <- refined.pca$x[ , 1]
        sv.mat[ , ii] <- refined.pca$x[ , indx]
    }
    rownames(sv.mat) = rownames(exp.resid.mat)
    return(sv.mat)
}


surrogatcalc$refine.ssva_perm <- function(exp.mat, exp.resid.mat, cov.data, numPerm, pvalThresh, scaleOnForPCA)
{
    num.mice = nrow(exp.mat)
    ##Step 1.2
    exp.resid.pca <- prcomp(exp.resid.mat, scale.=scaleOnForPCA)
    pdf(file.path(output, "./ssva_scree.pdf"))
    plot(exp.resid.pca$sdev^2, pch=as.character(c(1:9,0))) # suggests 7 pcs
    dev.off()
    
    T_k  = exp.resid.pca$sdev^2
    T_bk = matrix(-1,numPerm, length(T_k)) 
    
    for(b in 1:numPerm)
    {
        if(b%%100==0)
        {
            print(paste0("permutation ",b))
        }
        toPerm = exp.resid.mat;
        for (gene_i in 1:ncol(exp.resid.mat)) #gene_i is the gene index, we are going to permute the the sample labels, per gene
        {
            toPerm[,gene_i] = sample(toPerm[,gene_i], length(toPerm[,gene_i]),replace=F)
        }
		
        ##TODO comment out
        exp.resid.null = toPerm 
        T_bk[b,] = (prcomp(exp.resid.null, scale.=scaleOnForPCA)$sdev)^2
    }
    print("done permuting")
    
    pvals = rep(NA, num.mice)
    for(ii in 1:num.mice)
    {
        pvals[ii] = sum(T_k[ii]<T_bk[,ii])/numPerm
        if(ii>1)
        {
            pvals[ii] = max(pvals[ii-1], pvals[ii])
        }
    }
    num.sv = sum(pvals<pvalThresh)
    print(pvals)
    print(paste0("num.sv is ", num.sv))
                                           
    num.mice <- nrow(exp.resid.mat)
    sv.mat <- matrix(nrow=num.mice, ncol=num.sv) # holds refined sv's
    for (ii in 1:num.sv) #Step 2.3 is the result of 1.0-1.10 
    {
        ##step 2.2 (step 2.1 is cached from 1.1)
        rough.sv      <- exp.resid.pca$x[,ii]
        sv.mat[ , ii] <- rough.sv ##refined.pca$x[ , indx]
    }
    rownames(sv.mat) = rownames(exp.resid.mat)
    return(sv.mat)
}


## Convenient functions for generating functions that fullfill the run.sva interface, which use mixed model (mm) 
## The surrogatecalc functions that are called are agnostic as to the model form, and only depend on the residual matrix,
## whereas these functions below all assume a very specific form of model; i.e., f(y) = Xb + Zu + eps

##returns a function satisfying the interface f(traininedExpMat, cov.data, covariateModelString) for SSVA
surrogatcalc$get.mm.SSVA.func <- function(num.sv,
                                          transformParams = fit.model.bc$getDefaultTransformParams(),
                                          strategy = NULL,
                                          scaleOnForPCA=T,
                                          residualize = F,
                                          parallelArgs = parallel$getDefaultLocalArgs())
{
    
    f = function(exp.mat, exp.mat.control, cov.data, covariateModelString)
    {
##        print("running SSVA func")
        exp.ssva.mat = exp.mat.control
        
        if(residualize)
        {
            exp.ssva.mat = surrogatcalc$.computeResidsWrapper(exp.mat              = exp.mat.control,
                                                           cov.data             = cov.data,
                                                           covariateModelString = covariateModelString,
                                                           transformParams      = transformParams,
                                                           strategy             = strategy,
                                                           parallelArgs         = parallelArgs)
        }
##        print("got SSVA residuals")
        sv.mat   = surrogatcalc$refine.ssva(exp.ssva.mat   = exp.ssva.mat,
                                            num.sv          = num.sv,
                                            scaleOnForPCA   = scaleOnForPCA)
        return(sv.mat)
    }
    return(f)
}


surrogatcalc$get.mm.SVA.func <- function(num.sv,
                                         transformParams = fit.model.bc$getDefaultTransformParams,
                                         strategy = NULL,                                         
                                         topGenes = 1000,
                                         scaleOnForPCA = T,
                                         parallelArgs = parallel$getDefaultLocalArgs)
{
    f = function(exp.mat, exp.mat.control, cov.data, covariateModelString)
    {
        exp.resid.mat = surrogatcalc$.computeResidsWrapper(exp.mat              = exp.mat,
                                                           cov.data             = cov.data,
                                                           covariateModelString = covariateModelString,
                                                           transformParams      = transformParams,
                                                           strategy             = strategy,
                                                           parallelArgs         = parallelArgs)
        
        sv.mat = surrogatcalc$refine.sv(exp.mat       = exp.mat,
                                        exp.resid.mat = exp.resid.mat,
                                        num.sv        = num.sv,
                                        topgenes      = topGenes,
                                        scaleOnForPCA = scaleOnForPCA)
        return(sv.mat)
    }
    return(f)
}


surrogatcalc$get.mm.SSVA.perm.func <- function(transformParams = fit.model.bc$getDefaultTransformParams(),
                                               strategy = NULL,
                                               pvalThresh = .05,
                                               numPerm = 200,
                                               scaleOnForPCA = T,
                                               parallelArgs = parallel$getDefaultLocalArgs())
{
    f = function(exp.mat, exp.mat.control, cov.data, covariateModelString)
    {
        exp.resid.mat = surrogatcalc$.computeResidsWrapper(exp.mat              = exp.mat.control,
                                                           cov.data             = cov.data,
                                                           covariateModelString = covariateModelString,
                                                           transformParams      = transformParams,
                                                           strategy             = strategy,
                                                           parallelArgs         = parallelArgs)
        
        sv.mat   = surrogatcalc$refine.ssva_perm(exp.mat.control = exp.mat.control,
                                                 exp.resid.mat = exp.resid.mat,
                                                 cov.data = cov.data,
                                                 numPerm = numPerm,
                                                 pvalThresh = pvalThresh,
                                                 scaleOnForPCA = scaleOnForPCA)
        return(sv.mat)
    }
    return(f)
}


##compute and prep resids for surrogate variable construction
surrogatcalc$.computeResidsWrapper <- function(exp.mat,
                                               cov.data,
                                               covariateModelString,
                                               transformParams,
                                               strategy,
                                               rescale=F,
                                               parallelArgs = parallel$getDefaultLocalArgs())
{
    resids = fit.model.bc$computeResids(y.mat = exp.mat,
                                        cov.data = cov.data,
                                        covariateModelString = covariateModelString,
                                        transformParams = transformParams,
                                        strategy = strategy,
                                        parallelArgs = parallelArgs)


    na.cols <- which(apply(resids, 2, function(x){all(is.na(x))}))
    if(length(na.cols)>0) { resids <- resids[ , -na.cols] }
    ##RESIDUALS ARE SCALED!!! IS THIS GOOD? Should it be done in here?
    if(rescale) { resids <- scale(resids) }
    return(resids)
}

surrogatcalc$generate.svinfo <- function(svFunc,
                                         exp.mat,
                                         exp.mat.control,
                                         cov.data,
                                         covariateModelString,
                                         nullModelString = NULL,
                                         residualizeOutCovariates = c(),
                                         residualizeOutSV         = (length(residualizeOutCovariates)>0),
                                         transformParams          = fit.model.bc$getDefaultTransformParams(),
                                         strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(covariateModelString),
                                                                                       prefer.lme=T),
                                         parallelArgs = parallel$getDefaultLocalArgs())
{

    if(!residualizeOutSV && length(residualizeOutCovariates)>0)
    {
        warning("not residualizing out surrogate variables, but residualizing out other covariates... this is allowed but seems unlikely to be the user's intent.")
    }

    pracma::tic()

    if(any(rownames(cov.data)!= rownames(exp.mat)))
    {
        stop()
        browser()
    }
       
    sv.mat = svFunc(exp.mat, exp.mat.control, cov.data, covariateModelString)
    print(paste0("got surrogate variables:",pracma::toc()))
    gc()
    
    
    ##if(prop$mnp$saveIntermediates){save(list=ls(), file = fp(prop$mnp$output, "sv.RData"))}

###### merge covariates and surrogate variables, remove the samples that aren't being modelled in exp.mat, 
###### and reorder cov.data to follow same order as exp.mat
    colnames(sv.mat) = paste0("sv.",1:ncol(sv.mat))
    sv.mat = data.table(sv.mat)
    rownames(sv.mat) = rownames(exp.mat.control)
    sv.mat$ID = rownames(sv.mat)
    cov.data = sv.mat[cov.data, on = "ID"]
   
    rownames(cov.data) = cov.data$ID
    setkey(cov.data,  "ID")

    cov.data = cov.data[rownames(exp.mat)]
    rownames(cov.data) = cov.data$ID
    sv.mat$ID = NULL
    svString = paste(colnames(sv.mat), collapse = " + ")

    inserted = formulaWrapper$insertEffect(svString, 2, covariateModelString)
    covariateModelString = inserted$modified.string

    if(!is.null(nullModelString))
    {
        nullModelString = formulaWrapper$insertEffect(svString, 2, covariateModelString)$modified.string
    }
    
    ## Residualize out nuisance variables first if necessary.
    ##TODO: is this correct? should we be residulaizing with a box cox transform?
    ##TODO: also, should we pass that lambda on to the next step? Or just choose a new lambda as we have done?    
    
    if(residualizeOutSV)
    {
        residualizeOutCovariates  = c(colnames(sv.mat), residualizeOutCovariates)
    }

    
    if(length(residualizeOutCovariates)>0)
    {
        if(is.null(residualizeOutCovariates))
        {
            browser()
        }

        if(any(rownames(cov.data)!= rownames(exp.mat)))
        {
            stop()
            browser()
        }

        pracma::tic()
        
        residualized = fit.model.bc$residOutCovariates(y.mat    = exp.mat,
                                                       cov.data = cov.data,
                                                       covariateModelString = covariateModelString,
                                                       residualizeOutCovariates = residualizeOutCovariates,
                                                       transformParams = transformParams,
                                                       strategy = strategy,
                                                       parallelArgs    = parallelArgs)
        
        print(paste("residualized out nuisance covariates: ", pracma::toc()))

        exp.mat              = residualized$residualizedData;
        covariateModelString = residualized$residualizedModelString
        nullModelString      = residualized$residualizedNullString
    }
    
    out = list(covariateModelString = covariateModelString, nullModelString = nullModelString, cov.data = cov.data, sv.mat = sv.mat, exp.mat = exp.mat)
    return(out)
}


surrogatcalc$runAnalysis <- function(svFunc,
                                     exp.mat,
                                     exp.mat.control,
                                     cov.data,
                                     covariateModelString,
                                     nullModelString = NULL,
                                     modelParser, ##TODO give this a better default arg.
                                     strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString),
                                                                                   prefer.lme = T),
                                     residualizeOutCovariates = c(),
                                     residualizeOutSV         = (length(residualizeOutCovariates)>0),
                                     
                                     transformParams = fit.model.bc$getDefaultTransformParams(),
                                     parallelArgs    = parallel$getDefaultLocalArgs())
{
    ##create covariates frame that includes surrogate variables,
    ## updates the covariate string to include surrogate variables, and
    ## updates covariate string and expresssion mat to have residualized out all covariates (if any residualized covariates are
    ##specified.

    sv.info = surrogatcalc$generate.svinfo(svFunc                   = svFunc,
                                           exp.mat                  = exp.mat,
                                           exp.mat.control          = exp.mat.control,
                                           cov.data                 = cov.data,
                                           covariateModelString     = covariateModelString,
                                           nullModelString          = nullModelString,
                                           residualizeOutCovariates = residualizeOutCovariates,
                                           residualizeOutSV         = residualizeOutSV,
                                           transformParams          = transformParams,
                                           strategy                 = strategy,
                                           parallelArgs             = parallelArgs)

    ## print("got svinfo")

    results = fit.model.bc$fit(y.mat                = sv.info$exp.mat,
                               cov.data             = sv.info$cov.data,
                               covariateModelString = sv.info$covariateModelString,
                               nullModelString      = nullModelString,
                               modelParser          = modelParser,
                               transformParams      = transformParams,
                               checkAnova           = T,
                               strategy             = strategy,
                               parallelArgs         = parallelArgs)
    
    ## print(paste0("got fits:", pracma::toc()))
    ## if(T){save(list=ls(), file = fp(prop$mnp$output, "results.RData"))}

############### Generate permutation thresholds
    pracma::tic()
    out = (list(results = results,
                sv.info         = sv.info))

    return(out)
}

##wrapper around compute variance components that additionally computes the variance explaned by surrogate variables
surrogatcalc$computeVarExplained <- function(fit)
{
    varexp = lm.parsing$varexp(fit)
    response = varexp$response
    components = varexp$components
    ## for reasons unknown, variable names are not consisent between anova summary and coeficient fit summary, and we'd like them to be so we can easily join tables
    names(components)[names(components)=="intercept"] = "(Intercept)"
    
    num.sv = sum(grepl(pattern = "sv\\.[0-9]+", names(components)))
    if(num.sv>0)
    {
        svcols = paste0("sv.", 1:num.sv)
        components$sv.all = rowSums(components[,svcols, drop=F])
      ##  components = components[,-which(colnames(components) %in% svcols)]
    }
    varexp = (diag(var(components)/var(as.vector(response))))
    return(varexp)
}

