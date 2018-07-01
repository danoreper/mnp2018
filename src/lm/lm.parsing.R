library(nlme)
library(lmerTest)
library(data.table)

source("./utils.R")


lm.parsing = new.env(hash=T)

lm.parsing$getHeavySingleParser <- function(include.resid = F,
                                            mainEffectContrasts = c(),
                                            interactionContrasts = list(),
                                            varianceComputerList = lm.parsing$getDefaultVarianceComputer(),
                                            medianAdjust.p.value = T,
                                            computeFDR = T)
{
    force(include.resid)
    force(mainEffectContrasts)
    force(varianceComputerList)
    force(medianAdjust.p.value)
    force(computeFDR)
    
    parser = new.env(hash=T)
    
    parser$parse <- function(fit)
    {
        out = try(lm.parsing$parseModelFit(fit,
                                       include.resid = include.resid,
                                       mainEffectContrasts = mainEffectContrasts,
                                       interactionContrasts = interactionContrasts,
                                       varianceComputerList = varianceComputerList))
        if(class(out)=="try-error")
        {
            debug(lm.parsing$parseModelFit)
            lm.parsing$parseModelFit(fit,
                                     include.resid = include.resid,
                                     mainEffectContrasts = mainEffectContrasts,
                                     interactionContrasts = interactionContrasts,
                                     varianceComputerList = varianceComputerList)
            
        }
        return(out)
    }

    parser$collate <- function(outs, accum)
    {
        res = outs
        res.filt = list()
        cnamez = names(outs)
        failing = c()
        for(i in 1:length(res))
        {
            if(class(outs[[i]])!="try-error")
            {
                res.filt[[i]] = res[[i]]
                res.filt[[i]]$per.probe$Probe.Set.ID    = cnamez[i]
                res.filt[[i]]$per.variable$Probe.Set.ID = cnamez[i]
                res.filt[[i]]$per.level$Probe.Set.ID    = cnamez[i]
            } else {
                failing = c(failing,cnamez[[i]])
            }
        }

        results = list()
        results$failing      = failing
        results$per.probe    = rbindlist(lapply(res.filt, "[[", "per.probe"))
        results$per.variable = rbindlist(lapply(res.filt, "[[", "per.variable"))
        results$per.level    = rbindlist(lapply(res.filt, "[[", "per.level"))

        if(medianAdjust.p.value)
        {
            ## results$old.p.value = results$p.value
            results$per.variable[,anova.p.value := multipleTesting$median.correct.2(anova.p.value), by = variable]
        }
        if(computeFDR)
        {
            results$per.variable[,anova.q.value := p.adjust(anova.p.value , method = "fdr"), by=variable]
        }
        
        return(results)
    }
##    debug(parser$collate)
    return(parser)
}


lm.parsing$getHeavyDoubleParser  <- function(include.resid = F,
                                             mainEffectContrasts = c(),
                                             varianceComputerList = lm.parsing$getDefaultVarianceComputer())
{
    force(include.resid)
    force(mainEffectContrasts)
    force(varianceComputerList)
    parser = new.env(hash=T)

    parser$parse <- function(fit)
    {
        alt.vs.null = lm.parsing$getAnovaWrapper(fit$fit, fit$nullFit, type=1)
        p.value = try(alt.vs.null$an[[alt.vs.null$pvalueCol]][2])
        if(class(p.value)=="try-error")
        {
            p.value = NA
        }
        p.value = data.frame(p.value = p.value, lambda = fit$selectedLambda)
        return(p.value)
    }

    parser$collate <- function(outs, accum)
    {
        res = outs
        res              = do.call(rbind, res)
        cnamez           = names(outs)
        res$Probe.Set.ID = names(outs)
        return(res)    
    }

    return(parser)        
}


lm.parsing$getFullSingleParser <- function(includeData=F)
{
    force(includeData)
    parser = new.env(hash=T)

    parser$parse <- function(fit)
    {
        if(!includeData)
        {
            theclass = class(fit)
            if(theclass=="lme")
            {
                fit$data = NULL
            } else if(theclass %in% c("lmerMod", "merModLmerTest")){
                model.frame(fit) = NULL
            } else if(theclass == "lm")
            {
                fit$model = NULL
            } else {
                print(theclass)
                stop("unimplemented")
            }
        }

        return(fit)
    }

    parser$collate <- function(outs, accum)
    {
        outs1 = outs
        names(outs1) = names(outs)
        return(outs1)
    }
    return(parser)
}


lm.parsing$getFullDoubleParser <- function(includeData =F)
{
    force(includeData)
    parser = new.env(hash=T)

    parser$parse <- function(lm.out)
    {
        if(!includeData)
        {
            lm.out$fit$data = NULL
            lm.out$fit.null$data = NULL
        }
        aw = lm.parsing$getAnovaWrapper(lm.out$fit, lm.out$fit.null)
        lm.out$aw = aw
        return(lm.out)
    }

    parser$collate <- function(outs, accum)
    {
        outs1 = outs
        names(outs1) = names(outs)
        return(outs1)
    }
    return(parser)
}


## if(any(is.na(dt.pervariable[["anova.p-value"]])))
        ## {
        ##     print("calc signif failed at some point!!!")
        ## }
  
##TODO: figure out all main effect factors with more than 2 levels automatically, and avoid having to pass anything in at all. And or, pass in the contrasts per variable name.
lm.parsing$parseModelFit <- function(fullfit,
                                     include.resid = F,
                                     mainEffectContrasts=c(),
                                     interactionContrasts = list(),
                                     varianceComputerList=lm.parsing$getDefaultVarianceComputer())
{
    per.probe    = data.frame(lambda = fullfit$lambda)

    if(is.null(fullfit$fit))
    {
        return(NULL)
    }

    if(include.resid)
    {
        per.probe$resid  = paste(resid(fullfit$fit), collapse = ",")
##        per.probe$sample = paste(names(resid(fullfit$fit)), collapse = ",")
    }


    aw  = fullfit$anovaWrapper
    if(is.null(aw))
    {
        aw = lm.parsing$getAnovaWrapper(fullfit$fit)
        if(is.null(aw))
        {
            return(NULL)
        }
    }
    
    per.variable = data.table(aw$an)
    ##TODO why is this needed suddenly??
    setnames(per.variable, aw$pvalueCol, "p.value")

    
    setnames(per.variable, old = colnames(per.variable), new = paste0("anova.", colnames(per.variable)))
    per.variable$variable = rownames(aw$an)

        
    w = lm.parsing$getTsWrapper(fullfit$fit)
    per.level  = data.table(w$ts)
    setnames(per.level, old = colnames(per.level), new = paste0("coef.", colnames(per.level)))
    per.level$variable.level = rownames(w$ts)

###########################################
##compute variance explained
    for(varianceType in names(varianceComputerList))
    {
        varFunc = varianceComputerList[[varianceType]]
        varexp = varFunc(fullfit$fit)
        per.variable$varexp = varexp[match(per.variable$variable, names(varexp))]
        names(varexp)= paste0(varianceType, names(varexp))
        per.probe = cbind(per.probe, data.frame(as.list(varexp)))
    }

    ##refactor into a method
    ##contrasts
    for(contrastvar in mainEffectContrasts)
    {
        amcp = mcp(avar="Tukey")
        names(amcp) = contrastvar
        comps         = suppressWarnings(glht(fullfit$fit, linfct = amcp))
        tukey.p       = summary(comps)$test$pvalues
        nmz           = names((summary(comps))$test$tstat)
        nmz = gsub(pattern = " - ", replacement = ".vs.", nmz)
        nmz = paste0(nmz, ".p.value")
        names(tukey.p) = nmz
        tukey.p = data.frame(as.list(tukey.p))
        per.probe = cbind(per.probe, tukey.p)
    }

    ##TODO bring this back
    ##Refactor into a method.
    ## TODO reimplement strain-by-diet contrasts
    for(i in 1:length(interactionContrasts))
    {
        var1 = interactionContrasts[[i]][1]
        var2 = interactionContrasts[[i]][2]

        browser()
        fit.with.interaction = fullfit$fit
        contrast.mat = lm.parsing$form.interaction.contrast.mat(fit.with.interaction, var1, var2)
        glht.out = glht(fit.with.interaction, linfct = contrast.mat)
        tukeyDietByStrains =  (summary(glht.out)$test$pvalues)
        names(tukeyDietByStrains) = names(((summary(glht.out))$test)$coefficients)
        
        names(tukeyDietByStrains) = gsub(names(tukeyDietByStrains), pattern = var2, replacement = "")
        names(tukeyDietByStrains) = gsub(names(tukeyDietByStrains), pattern = var1, replacement = "")
        tukeyDietByStrains = data.frame(as.list(tukeyDietByStrains), check.names = T)
        ## varExp = lm.parsing$varexp(fit.with.interaction)
        ## varExp = var(varExp$components)/var(as.vector(varExp$response))
        per.probe = cbind(per.probe, tukeyDietByStrains)
    }

    ##clean up column names
    for(fram in list(per.variable, per.level))
    {
        if(ncol(fram)>0)
        {
            setnames(x = fram,
                     old = colnames(fram),
                     new = sub(x = colnames(fram), pattern = "-", replacement ="."))
        }
    }
    ret = list(per.probe = per.probe, per.variable = per.variable, per.level = per.level)
    return(ret)
}

lm.parsing$getDefaultVarianceComputer <- function()
{
    f <- function(fit)
    {
        varexp = lm.parsing$varexp(fit)
        response = varexp$response
        components = varexp$components
        ## for reasons unknown, variable names are not consisent between anova summary and coeficient fit summary, and we'd like them to be so we can easily join tables
        names(components)[names(components)=="intercept"] = "(Intercept)"
        varexp = (diag(var(components)/var(as.vector(response))))
        return(varexp)
    }
    return(list(var=f))
}


lm.parsing$get.Z.U <- function(afit, frame = lm.parsing$getCovFrame(afit), sum=F)
{
    if(class(afit)=="lme")
    {
        randomformulas = formula(afit$modelStruct$reStruct)
        Z = list()
        U = list()
        for(i in 1:length(randomformulas))
        {
            formName = names(randomformulas)[i]
            
            randomformula = randomformulas[[i]]
            randomformula = as.character(randomformula)
            if(length(randomformula)>2)
            {
                stop("unimplemented")
            }
            randomformula = randomformula[2]
            termz = strsplit(randomformula, split = "\\+")[[1]]
            termz = str_trim(termz)
            for(j in 1:length(termz))
            {

                term = termz[j]

                if(term =="1")
                {
                    randomformula = formName
                } else {
                    randomformula = paste0(formName,":",term)
                }
                
                randomformula = paste("~", randomformula)
                contrs = list()
                contrs[[formName]] = contrasts(eval(parse(text=formName), envir = frame), contrasts=F)
                if(term!="1")
                {
                    contrs[[term]]     = contrasts(eval(parse(text=term),     envir = frame), contrasts=F)
                }
                ##modelmat = model.matrix(formula(randomformula),  frame, contrasts.arg = contrs)
                modelmat = model.matrix(formula(randomformula),  frame, contrasts.arg = contrs)
                modelmat = modelmat[,-which(colnames(modelmat) %in% "(Intercept)")]

                ##remove unnecessary levels from model matrix
                if(term!="1")
                {
                    modelmatorig = modelmat
                    modelmat = list()

                    for(k in 2:length(colnames(ranef(afit))))
                    {
                        colm = paste0(formName, rownames(ranef(afit)), ":", colnames(ranef(afit))[k]) 
                        modelmat = util$appendToList(modelmat, modelmatorig[,colm])
                    }
                    modelmat = do.call(cbind, modelmat)
                }

                fullname = paste0(formName, ".", names(ranef(afit))[j])
                U[[fullname]] = ranef(afit)[,j]
                Z[[fullname]] = modelmat
            }
        }

        if(length(names(randomformulas))>1)
        {
            stop("unimplemented")
        }
    }
    else {
        Z     = lme4::getME(afit, "Ztlist")
        Z     = lapply(FUN=t, Z)
        ## Znew  = list()
        ## for(zname in names(Z))
        ## {
        ##     nesting = strsplit(zname, "\\.")
        ##     nest.1 = nesting[[1]]
        ##     nest.2 = nesting[[2]]
        ##     if(!(nest.1 %in% names(Znew)))
        ##     {
        ##         Znew[[nest.1]] = list()
        ##     }
        ##     Znew[[nest.1]][[nest.2]] = Z[[zname]]
        ## }
        ## Z = Znew

        U     = ranef(afit)
        Unew = list()
        for(u.name in names(U))
        {
            usub = U[[u.name]]
            
            for(u.name.sub in colnames(usub))
            {
                Unew[[paste0(u.name, ".", u.name.sub)]] = usub[[u.name.sub]]
            }
        }
        U = Unew
    }


    
    out = (list(Z=Z, U=U))
    if(sum)
    {
        firstname = names(U)[1]
        total = 0*(Z[[firstname]] %*% U[[firstname]]) 
        for(u.name in names(U))
        {
            total = total + Z[[u.name]] %*% U[[u.name]]
        }
        out$total = as.vector(total)
    }


    return(out)
}

lm.parsing$getLHS <- function(afit, fixedterms = lm.parsing$getFixedTerms(afit))
{
    out = stringr::str_trim(as.character(fixedterms)[2])    
    return(out)
}

lm.parsing$getFixedTerms <- function(afit)
{
    theclass = class(afit)
    
    if(theclass=="lme")
    {
        termz = afit$terms
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        termz = terms(afit, fixed.only=T)
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(termz)
}


lm.parsing$getCovFrame <- function(afit)
{
    theclass = class(afit)
    if(theclass=="lme")
    {
        frame = afit$data
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        frame = model.frame(afit)
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(frame)

}

lm.parsing$getX <- function(afit, frame = lm.parsing$getCovFrame(afit))
{
    theclass = class(afit)
    
    if(theclass=="lme")
    {
        termz = afit$terms
        X = model.matrix(as.formula(termz), frame)
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        X     = lme4::getME(afit, "X")
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(X)
}


##doreper version: works on lme or lmer object
lm.parsing$varexp <- function(afit)
{
    frame = lm.parsing$getCovFrame(afit)
    Z.U   = lm.parsing$get.Z.U(afit)
    Z = Z.U[["Z"]]
    U = Z.U[["U"]]

    termz = lm.parsing$getFixedTerms(afit)
    X     = lm.parsing$getX(afit)
    lhs   = lm.parsing$getLHS(afit, termz)
    Beta  = fixef(afit)
  
    resids = resid(afit)
    fixedTerms = unlist(attr(termz, "term.labels"))

    effectsPerSample = X*matrix(rep(as.vector(Beta),nrow(X)), byrow=T,nrow=nrow(X))


    newdf = data.frame(ID = rownames(effectsPerSample), intercept=Beta[["(Intercept)"]])
    for (fixt in fixedTerms)
    {
  ##      print(fixt)
        newTokens = list()
        interactTokens = unlist(strsplit(fixt, ":"))
        for(interactToken in interactTokens)
        {
            acol = try(frame[[interactToken]])
            if(class(acol)=="try-error"|is.null(acol))
            {
                acol = eval(parse(text=interactToken), envir=frame)
            }
            if(is.factor(acol))
            {
                levs = levels(acol)
                newTokens = util$appendToList(newTokens, paste0(interactToken, levs))
            } else if(is.logical(acol))
            {
                newTokens = util$appendToList(newTokens, paste0(interactToken, "TRUE"))
            } else if(is.character(acol))
            {
                stop("unimplemented")
            }
            else
            {
                newTokens = util$appendToList(newTokens, interactToken)
            }
        
        }

        interactions = levels(do.call(interaction, c(newTokens, sep=":"))) 
        interactions = intersect(interactions, colnames(effectsPerSample))
        groupedCol = effectsPerSample[,interactions]
        if(length(interactions)>1)
        {
            groupedCol = rowSums(groupedCol)
        }
        groupedCol = matrix(groupedCol,ncol=1)
        colnames(groupedCol) = fixt

        newdf = cbind(newdf, groupedCol)
    }
    rownames(newdf) = newdf$ID
    newdf$ID = NULL

    newraneff = data.frame(ID = rownames(newdf))
    for(aranef in names(U))
    { 
        ranefmat = U[[aranef]]
        zmat     = Z[[aranef]]

        persample = try(matrix(zmat%*%ranefmat,ncol=1))
        if(class(persample)=="try-error")
        {
            browser()
        }
            
        colnames(persample) = aranef
        newraneff = cbind(newraneff, persample)
    }
    rownames(newraneff)=newraneff$ID
    newraneff$ID = NULL

    newdf = cbind(newdf, newraneff)
    newdf$resid = resids

    ##    totalVar = var(frame[[lhs]])
    response = eval(parse(text=lhs), envir=frame)

    vars  = cov(newdf)    
    vars  = vars/var(as.vector(response))
    if(any(diag(vars)<0))
    {
        browser()
    }
    return(list(response = response, components = newdf))
}


##WV version.. doesn't seem to work right
lm.parsing$lmer.vartable <- function(fit.lmer, pctvar=TRUE)
{

     y      <- fitted(fit.lmer)+resid(fit.lmer)
     n      <- length(y)
     totvar <- var(y)
     tss    <- totvar*(n-1)
                                         # variance components
     vc    <- VarCorr(fit.lmer)
     out.df <- NULL
     for (comp.name in names(vc))
     {
         comp      <- vc[[comp.name]]
         comp.sd   <- attr(comp, "stddev") #comp@factors$correlation@sd
         comp.var  <- comp.sd^2
         comp.type <- colnames(comp) #@factors$correlation)
        
         df <- data.frame(
             Name     = rep(comp.name,length(comp.var)),
             Type     = comp.type,
             Variance = comp.var)
         out.df <- rbind(out.df, df)
     }
     scale <- attr(vc, "sc")
     df <- data.frame(
         Name     = "Residual",
         Type     = "(Intercept)",
         Variance = scale^2)
     out.df <- rbind(out.df, df)
     out.df$PctVar <- 100 * out.df$Variance / totvar
    
    ## fixed components

    an <- anova(fit.lmer, type=1)
    print(an)
     if (0!=nrow(an))
 	{
             fix.df <- data.frame(
                 Name     = rownames(an),
                 Type     = rep("Fixed", nrow(an)),
                 Variance = an$"Sum Sq"/n,
                 PctVar   = 100 * an$"Sum Sq"/tss)
             out.df <- rbind(out.df, fix.df)
 	}
     rownames(out.df) <- 1:nrow(out.df)
    
     return (out.df)
} 

lm.parsing$.getLevz <- function(fit.with.interaction, var1)
{
    if(class(fit.with.interaction)=="lme")
    {
        adat = fit.with.interaction$data   
    } else {
        adat = attr(fit.with.interaction, "frame")
    }
    
    var1.levelz = var1
    acol = adat[[var1]]
    
    if(is.factor(acol[1]))
    {
        var1.levelz = paste0(var1, levels(factor(acol)))
    } else if(is.logical(acol[1]))
    {
        var1.levelz = paste0(var1,c("FALSE", "TRUE"))
    } else if(is.character(acol[1]))
    {
        stop("unimplemented...")
    }
    return(var1.levelz)
}

##TODO, implment for case where var2 is also a factor, and var1 is not. Make things more general
##TODO, make this work for lme
lm.parsing$form.interaction.contrast.mat <- function( fit.with.interaction, var1, var2)
{

    ##browser()
    var1.levelz = lm.parsing$.getLevz(fit.with.interaction, var1)
    var2.levelz = lm.parsing$.getLevz(fit.with.interaction, var2)
##    var2.levelz = var2.levelz[-length(var2.levelz)]

                              
    alleffectnames = names(fixef(fit.with.interaction))
    pairz = expand.grid(var1.levelz, var2.levelz)
    pairz$sep = ":"
    pairz = do.call(paste, pairz)

    refpair = paste0(var1.levelz[1],":", var2.levelz[2])##setdiff(pairz, alleffectnames)
    pairz = intersect(pairz, alleffectnames)
    pairz.inds = match(pairz, alleffectnames)
    ## allDeltas = t(combn(contrast.inds,2))
    allDeltaStrings = t(combn(pairz,2))
    allDeltas = cbind(match(allDeltaStrings[,1], alleffectnames),
                      match(allDeltaStrings[,2], alleffectnames))
    nvar = length(names(fixef(fit.with.interaction)))
    ncontrast = length(pairz) + nrow(allDeltas)
    contrast.mat = matrix(0, nrow = ncontrast, ncol = nvar)
    rownames(contrast.mat) = rep("", nrow(contrast.mat))
    
    counter = 1
    for(c1 in pairz.inds)
    {
        contrast.mat[counter,c1]=1
        contrast.mat[counter,which(alleffectnames == var2.levelz[2])]=-1
        rownames(contrast.mat)[counter] = paste(pairz[counter], "-", refpair[1])
        counter = counter+1
    }
    for(i in 1:nrow(allDeltas))
    {
        contrast.mat[counter, allDeltas[i,2]] = 1
        contrast.mat[counter, allDeltas[i,1]] = -1
        rownames(contrast.mat)[counter] = paste(allDeltaStrings[i,2], "-",
                                                allDeltaStrings[i,1])

        counter = counter+1
    }

    return(contrast.mat)
}

lm.parsing$getFixefVec <- function(fit, varname)
{
    if(class(fit)=="lm")
    {
        return(coef(fit)[varname])
    } else {
        return(fixef(fit)[varname])
    }
}
               
lm.parsing$getAnovaWrapper <- function(fit, fit.null=NULL, type=1)
{
    if(is.null(fit.null))
    {
        if(class(fit) =="lmerModLmerTest")
        {
            ##browser()
            return(list(an = anova(fit, type=type),
                        pvalueCol = "Pr(>F)"))
            
        } else if(class(fit)=="lm")
        {
            return(list(an = stats::anova(fit),
                        pvalueCol = "Pr(>F)"))
        } else {
            ##Default anova does NOT allow specification of type... not sure what to do here
            return(list(an = stats::anova(fit),
                        pvalueCol = "p-value"))
        }
    ## } else {
    ##     if(class(fit) =="merModLmerTest")
    ##     {
    ##         return(list(an = lmerTest::anova(fit, fit.null, type=type),
    ##                     pvalueCol = "Pr(>F)"))
    ##     } else if(class(fit)=="lm")
    ##     {
    ##         return(list(an = stats::anova(fit, fit.null),
    ##                     pvalueCol = "Pr(>F)"))
    ##     } else {
    ##         ##Default anova does NOT allow specification of type... not sure what to do here
    ##         return(list(an = stats::anova(fit, fit.null),
    ##                     pvalueCol = "p-value"))
    ##     }
    ## }
    } else {
        if(class(fit) =="lmerModLmerTest")
        {
            an = lmerTest::anova(fit, fit.null, type=type)
            return(list(an = an,
                        p.value = an[2,"Pr(>F)"],
                        F.stat = an[2, "F"]))
         } else if(class(fit)=="lm")
         {
             an = stats::anova(fit, fit.null)
             return(list(an = an, 
                         p.value = an[2,"Pr(>F)"],
                         F.stat = an[2, "F"]))
         } else {
             ##Default anova does NOT allow specification of type... not sure what to do here
             an = stats::anova(fit, fit.null)
             return(list(an = an, 
                         p.value = an[2,"Pr(>F)"],
                         F.stat = an[2, "F"]))
             
         }
    }
}

lm.parsing$getTsWrapper <- function(fit)
{

    if(class(fit) =="lmerModLmerTest")
    {
        return(list(
            ts        = coefficients(summary(fit)),
            estCol    = "Estimate",
            tvalueCol = "t value",
            pvalueCol = "Pr(>|t|)",
            seCol     = "Std. Error",
            dfCol     = "df"))
    } else if (class(fit) =="lm") {
        return(list(
            ts        = coefficients(summary(fit)),
            estCol    = "Estimate",
            tvalueCol = "t value",
            pvalueCol = "Pr(>|t|)",
            seCol     = "Std. Error",
            dfCol     = "df"))
    } else {
        return(list(
            ts        = summary(fit)$tTable,
            estCol    = "Value",
            tvalueCol = "t-value",
            pvalueCol = "p-value",
            seCol    = "Std.Error",
            dfCol    = "DF"
        ))
    }
}



lm.parsing$compareModels <- function(modelfull, modelrestricted)
{
    an.compare = try(anova(modelfull, modelrestricted))
    out = try(list(p.value = an.compare["modelrestricted", "p-value"],
                   L.Ratio = an.compare["modelrestricted", "L.Ratio"],
                   num.df  = an.compare["modelfull",       "df"],
                   denom.df= an.compare["modelrestricted", "df"]))
                   
    return(out)
}


lm.parsing$getNoiseVector <- function(fit.with.interaction, gurka)
{
    noiseVector           = try(resid(fit.with.interaction))
    if(class(noiseVector)=="try-error")
    {
        browser()
    }
    if(gurka && class(fit.with.interaction) %in% c("lmerMod", "merModLmerTest", "lme"))
    {
        noiseVector = noiseVector + lm.parsing$get.Z.U(fit.with.interaction, sum = T)$total
    }
    return(noiseVector)
}

lm.parsing$get.lik.ratio.stat <- function(modelAlt, modelNull)
{
    out = as.numeric(2*(logLik(modelAlt)-logLik(modelNull)))
    return(out)
}


lm.parsing$getY <- function(fitWrapper)
{
    ##TODO implement
    return(0)
}



## ##TODO: excise, refactor heavy parser
## compareModels <- function(modelfull, modelrestricted)
## {
##     if(length(modelfull)!=length(modelrestricted))
##     {
##         stop("mismatched lengths")
##     }
##     comparisons = rep(NA, length(modelfull))
    
##     for(i in 1:length(modelfull))
##     {
##         pval = try(lm.parsing$compareModels(modelfull[[i]]$fit, modelrestricted[[i]]$fit))
        
##         if(class(pval)=="try-error")
##         {
##             pvals[i] = NA
##         } else {
##             pvals[i] = pval
##         }
##     }
##     return(pvals)
## }
