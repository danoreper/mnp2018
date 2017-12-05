library(WVmisc)
library(nlme)
library(MASS)
library(data.table)
library(lme4)
library(lattice)
library(reshape)
library(ggplot2)
library(parallel)
library("biomaRt")
##library(lmerTest) #for lmerTest, caTools, bitOps, and possibly others need to be intentionally installed, automatic dependency checking doesnt work right. really annoying.
options(warn=1)

if(!exists("prop"))
{
    stop("master script must load default params")
}
source("./genomerep/buildGenomeData.R")


##lme (data =cov.data, as.formula(paste0(lhs, " ~ 1 + ", svString, "as.factor(Batch) +", pipelineString, "Diet + Strain +  Strain:Diet")), random = ~ 1 | DamID)

##used in batching mclapply jobs- memory blows up otherwise, since garbage collection apparently never happens untill a given call to mclapply finishes
##fits a model, walking over possible box cox constants, using the best one. Allows the user to pass in a surrogate variable matrix optionally.
fit.model <- function(js = NULL, exp.mat, cov.data, sv.mat=matrix(0,0,0), lambdaPerGene=NULL, includeStrainByDiet = T, checkAnova = T, numBatch = 500, scanType="regular", method = "REML")
{
    runsinglefit <- function( jj, lambdaPerGene= NULL, checkAnova, exp.mat, lambdasToTry, pheno, covariates, method = "REML") 
    {
        ##        print(jj)
        
        if (jj %% 100 == 0) 
        { cat(sep="","[",jj,"]") } 
        
        y <- exp.mat[ , jj]
        cov.data$y=y
        
        if(!is.null(lambdaPerGene))
        {
            lambdasToTry = lambdaPerGene[jj]
        }
        fit.bestLambda = try( getBestLambda(lambdasToTry, pheno, covariates, cov.data,checkAnova = checkAnova, scanType=scanType, method=method))
        if (caught.error(fit.bestLambda))
        {
            fit.bestLambda = list(fit=NULL, selectedLambda=NULL, anovaWrapper=NULL)
        } 
        
        fit.bestLambda$Probe.Set.ID = colnames(exp.mat)[jj]
        return(fit.bestLambda)
    }
    
    if(is.null(js))
    {
        js = 1:ncol(exp.mat)
    }
    
    exp.resid.mat <- matrix(NA, nrow=nrow(exp.mat), ncol=ncol(exp.mat)) 
    
    svString = ""
    for(k in seq(1, ncol(sv.mat), length=ncol(sv.mat)))
    {
        colnam = paste0("sv.", k)
        cov.data[[colnam]] =sv.mat[,k]
        svString =  paste0(svString, colnam, "+")
    }

    pipelineString = ""

    modelPipeline = length(unique(cov.data$Pipeline))>1
    if(modelPipeline)
    {
        pipelineString ="Pipeline +"
    }
    
    pipelinestrainString = ""
    if(modelPipeline & prop$mnp$modelStrainPipeline)
    {
        pipelinestrainString ="Pipeline:Strain + "  
    }
    
    strainDietString = ""
    if(includeStrainByDiet)
    {
        strainDietString = "Strain:Diet + " 
    }

    if(scanType=="regular")
    {
        covariates = paste0( " ~ 1 + ", svString, "as.factor(Batch) + ", pipelineString, "Diet + Strain + ", pipelinestrainString, strainDietString, " (1|DamID)")
    }
    else if (scanType =="nostrain")
    {
        covariates = paste0( " ~ 1 + ", svString, "as.factor(Batch) + ", pipelineString, "Diet + ", " (1|DamID)")
    }
    else if (scanType  =="nodiet")
    {
        covariates = paste0( " ~ 1 + ", svString, "as.factor(Batch) + ", pipelineString, "Strain + ", " (1|DamID)")    
    }
    else if (scanType =="nostrainnodiet")
    {
        covariates = paste0( " ~ 1 + ", svString, "as.factor(Batch) + ", pipelineString, " (1|DamID)")
    }

    ##print(scanType)
    ##print(covariates)
    pheno ="y"
    
    
    inds = js
    
    batchCuts = rep(factor(T), length(inds))
    if(numBatch>1)
    {	
        numBatch = min(numBatch, length(inds))
        batchCuts = cut(inds, numBatch)
    }
    fitsAll = list()
    kk=1
    for(levl in levels(batchCuts))
    {
        indsGroup = inds[which(batchCuts==levl)]
        ##print(paste0("group is:",indsGroup))
        fits = lapply(as.list(indsGroup), runsinglefit, exp.mat = exp.mat, lambdaPerGene = lambdaPerGene, checkAnova = checkAnova, lambdasToTry = prop$mnp$lambdasToTry, pheno = pheno, covariates = covariates, method = method)
        fitsAll[[kk]] = fits
        kk= kk + 1
        ##gc()
    }
    fitsAll = do.call(c, fitsAll)
    date()
    return(fitsAll)
}

runModelsWithSv <- function(exp.mat, sv.mat, cov.data) 
{
    ##func to be lapplied
    lfunc <- function(indexes, exp.mat, sv.mat, cov.data)
    {
        ##print(paste0("indexes:", indexes))
        resids          = list()
        ## for(scanType in prop$mnp$scantypes) {resids[[scanType]] = list()}

        fits = fit.model(js=indexes, exp.mat = exp.mat, sv.mat = sv.mat, cov.data = cov.data, scanType = "regular", numBatch=1, checkAnova = T)
        
        lambdaPerGene = rep(NA, max(indexes))
        lambdaPerGene[indexes] = unlist(lapply(fits, "[[", "selectedLambda"))
        
        ##build a ML fit for both the full model, and the various restricted models
        mlfits = list()
        ##
        for(scanType in prop$mnp$scantypes)
        {
            ## print(scanType)
            ## yes, this does rebuild full model as well, but using ML method, which will be necessary later to compare models with anova.
            ## Also, we wont bother computing anova on individual models, as these will be used solely to compare models
            mlfits[[scanType]] = fit.model(js=indexes, exp.mat = exp.mat, sv.mat = sv.mat, cov.data = cov.data, scanType = scanType, numBatch=1, lambdaPerGene = lambdaPerGene, method ="ML",checkAnova=F)
        }
        results = parseModelFits(fits, mlfits)
        return(results)
    }
    

    ## func to accumulate results
    accumfunc <- function(res)
    {
        ##
        crossEffects = plyr::rbind.fill(lapply(res,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
        ## crossEffectsAll[[kk]] = data.frame(Probe.Set.ID = colnames(exp.mat)[indsGroup])
        crossEffects = cbind(Probe.Set.ID = crossEffects$Probe.Set.ID, crossEffects[,-which(colnames(crossEffects)=="Probe.Set.ID")])
        ##        for(scanType in prop$mnp$scantypes) {resids[[scanType]] = do.call(cbind, resids[[scanType]])} 
        ##      results         = list(results=crossEffectsAll, resids = residsAll)
        return(crossEffects)
    }
    
    ##Collate results across batches
    alist = mnp.lapply.wrapper(len = ncol(exp.mat), FUN = lfunc, exp.mat = exp.mat, sv.mat = sv.mat, cov.data = cov.data)

    return(accumfunc(alist))
}

mnp.lapply.wrapper <- function(len, FUN, ...)
{
    return(util$lapply.wrapper(len, FUN, batchSize = prop$mnp$batchSize, mc.cores = prop$mnp$mc.cores, ...)
}

formulaWrapper.boxCoxString <- function(phenotype, bestLambda1, dataSet)
{
    offset = .1  #we cant have lambda2==-y_i, so we will make sure lambda2>=-y_i+offset
    lambda2 = max(0, max(-dataSet[[phenotype]])+offset) #dont bother having a lambda2 unless there is at least one zero value
    
    if(bestLambda1 == 0)
    {
        return(paste0(    "(log(", phenotype, " + ", lambda2, "))"   
                      ))
    }
    else
    {
        return (paste0(    "((((", phenotype, " + ", lambda2, ")", "^", bestLambda1, ") -1 )/", bestLambda1 ,")"   ));
    }
}

getFormla <- function(lambda, pheno, covariates, dataSet)
{
    formla = (paste0(formulaWrapper.boxCoxString(pheno, lambda, dataSet), covariates))
    return(formla)
}

inverseNormal <- function(phenCol)
{
    allRanks = rank(phenCol, na.last="keep")
    maxRank = max(allRanks, na.rm=TRUE);
    return(qnorm(allRanks/(maxRank+1)));
}

getBestLambda <- function(lambdas, pheno, covariates, dataSet, checkAnova=T, method="REML", scanType="regular") 
{
    fitMixedModel <- function(pheno, lambda, dataSet, covariates, method = "REML") 
    {
        ##print(lambda)
        if(is.na(lambda))
        {
            
        }
        if(lambda=="inverse_normal")
        {
            lhs = paste0("inverseNormal(",pheno,")")
        } else {
            lhs = formulaWrapper.boxCoxString(pheno, lambda, dataSet)
        }
        if(prop$mnp$normalizeExpressionAfterTransform)
        {
            lhs = paste0("scale(",lhs,")")
        }
        if(!prop$mnp$uselme)
        {
            formla = (paste0(lhs, covariates))
            ##fit.with.interaction = try(lme4::lmer(formla, data=dataSet))
            ##TODO: this is probably wrong now and needs to make use of substitute and eval.
            fit.with.interaction = try(lmerTest::lmer(formla, data=dataSet))
            print(formla)
        } else
        {
            tokens = stringr::str_trim(unlist(strsplit(covariates,"\\+")))
            isRandom = grepl("\\|", tokens)
            fixedPart = paste(tokens[!isRandom], collapse="+")
            randomPart = paste(tokens[isRandom], collapse="+")
            randomPart = as.formula(paste0("~",gsub("\\(|\\)", "", randomPart)))
            fixedPart  = as.formula(paste0(lhs,fixedPart))

            ## we are doing this rigamarole with substitute so that the "fixed" and "random"
            ## call fields of 
            ##the resulting lme object make sense in the environment of dataSet,
            ##and not just in this functions scope. NB anova comparing 2 lme models requires use of the fixed and random fields.
            fit.with.interaction = try(eval(substitute(lme (data =dataSet, fixedPart, random = randomPart, method = method))))
        }

        return(fit.with.interaction)
    }
    
    pvals = rep(NA, length(lambdas))
    for(l in 1:length(lambdas))
    {
        fit.with.interaction = NULL
        lambda = lambdas [l]
        ##print(lambda)
        
        fit.with.interaction = fitMixedModel(pheno = pheno, lambda = lambda, dataSet = dataSet, covariates = covariates, method = method)
        
        if(class(fit.with.interaction)=="try-error")
        {
            next
        }
        resid.fit = NULL
        resid.fit            = try(resid(fit.with.interaction))
        if(class(resid.fit)=="try-error")
        {
            
        }
        
        anovaWrapper = NULL
        
        if(checkAnova)
        {
            anovaWrapper = try(getAnovaWrapper(fit.with.interaction))
            if(class(anovaWrapper)=="try-error"| !anovaWrapper$pvalueCol %in% colnames(anovaWrapper$an))
            {
                print("anova failure")
                ##			
                next
            }
        }	
        pval                 = NULL
        pval                 = try(shapiro.test(resid.fit)$p.value)
        if(class(pval)=="try-error")
        {
            
        }
        
        pvals[l] = pval
    }
    selectedLambda  = lambdas[which.max(pvals)]
    if(length(selectedLambda)==0|is.na(selectedLambda))
    {
        return(list(fit=NULL, selectedLambda = NULL, anovaWrapper = NULL))
    } 
    ## In principle we could cache the model and avoid some computation,
    ## but really, most of the time it turns out that we require
    ## recomputing with the inverse normal transform anyway. for simplicity, just recompute
    if((selectedLambda>2) || (selectedLambda<(-2)))
    {
        selectedLambda = "inverse_normal"
    }
    fit.with.interaction = fitMixedModel(pheno = pheno, lambda = selectedLambda, dataSet = dataSet, covariates = covariates, method = method)
    
    
    if(class(fit.with.interaction)=="try-error")
    {
        
    }
    if(checkAnova)
    {
	anovaWrapper = try(getAnovaWrapper(fit.with.interaction))
    }
    ##print(selectedLambda)
    return(list(fit=fit.with.interaction, selectedLambda = selectedLambda, anovaWrapper = anovaWrapper))
}

getLambdas <- function(fits)
{
    unlist(lapply(fits, "[[", "selectedLambda"))
}

compareModels <- function(modelfull, modelrestricted)
{
    an.compare = try(anova(modelfull, modelrestricted))
    log.pvalue.restricted = -log10(an.compare["modelrestricted","p-value"])
    return(log.pvalue.restricted)
}

## Test for effects of cross at any diet
##Step 2.8

parseModelFits <- function(fullfits, mlfits)
{

    lambdas  = getLambdas(fullfits)
    for (i in 1:length(fullfits))
    {
        fullfit      = fullfits[[i]]
        Probe.Set.ID = fullfit$Probe.Set.ID 
        aw  = fullfit$anov
        if(is.null(aw))
        {
            
        }
        
        includesPipe = "Pipeline" %in% rownames(aw$an)
        pipeP = NA 
        if(includesPipe)
        {
            pipeP = -log10(aw$an["Pipeline",         aw$pvalueCol])
        }
        out <- data.frame(
            LogP.Anova.SV          = -log10(aw$an["sv.mat",           aw$pvalueCol]),
            LogP.Anova.Batch       = -log10(aw$an["as.factor(Batch)", aw$pvalueCol]),
            LogP.Anova.Pipeline    =  pipeP,
            LogP.Anova.Diet        = -log10(aw$an["Diet",             aw$pvalueCol]),
            LogP.Anova.Cross       = -log10(aw$an["Strain",            aw$pvalueCol]),
            LogP.Anova.CrossByDiet = -log10(aw$an["Diet:Strain",       aw$pvalueCol]))
	
        if(prop$mnp$modelStrainPipeline)
        {
            out$LogP.Anova.CrossByPipeline = -log10(aw$an["Pipeline:Strain",       aw$pvalueCol])
        }
        
        if(is.null(out$LogP.Anova.Cross))
        {
            
        }
        if(is.na(out$LogP.Anova.Cross))
        {
            print("calc signif failed at some point!!!")
        }
        
        w = getTsWrapper(fullfit$fit)
        out[ , paste0("Est.CrossInStdCtrl") ]  <- w$ts[ "StrainNOD.B6", w$estCol ]
        out[ , paste0("Tval.CrossInStdCtrl") ] <- w$ts[ "StrainNOD.B6", w$tvalueCol ]
        out[ , paste0("LogP.CrossInStdCtrl") ] <- -log10(w$ts[ "StrainNOD.B6", w$pvalueCol ])
        out[ , paste0("SE.CrossInStdCtrl")]    = w$ts[ "StrainNOD.B6", w$seCol]
        out[ , paste0("DF.CrossInStdCtrl")]    = w$ts[ "StrainNOD.B6", w$dfCol]


        alternateDietz = grepl("Diet.*:", rownames(w$ts))
        alternateDietz = rownames(w$ts)[alternateDietz]
        alternateDietz = sub(pattern = ":.*", x = alternateDietz, replacement = "")
        alternateDietz = sub(pattern = "Diet",x = alternateDietz, replacement = "")

        
        for (diet in alternateDietz) 
        {
            interaction.string <- paste("Diet", diet, ":StrainNOD.B6",sep="")
            ##pred.string <- paste(diet, ":StrainNOD.B6",sep="")
            out[ , paste0("Est.CrossBy", diet, ".vs.StdCtrl") ]  <-  w$ts[ interaction.string, w$estCol ]
            out[ , paste0("Tval.CrossBy", diet, ".vs.StdCtrl") ] <-  w$ts[ interaction.string, w$tvalueCol ]
            out[ , paste0("LogP.CrossBy", diet, ".vs.StdCtrl") ] <- -log10(w$ts[ interaction.string, w$pvalueCol ])
        }

        
        out = cbind(out, computeVarExplained(fullfits[[i]]$fit,prop$mnp$num.sv))

        for(scanType in prop$mnp$scantypes)
        {
            if(scanType!="regular")
            {
                
                modelfull       = mlfits[["regular"]][[i]]$fit
                modelrestricted = mlfits[[scanType]][[i]]$fit
                cmpr = try(compareModels(modelfull, modelrestricted))
                if(class(cmpr)=="try-error")
                {
                    cmpr = NA
                }
                
                if(scanType =="nostrain")
                {
                    out$LogP.Anova.Cross.with.interaction = cmpr
                }
                else if(scanType=="nodiet")
                {
                    out$LogP.Anova.Diet.with.interaction = cmpr
                }
            }
        }
        out$selectedLambda = fullfits[[i]]$selectedLambda
        out$Probe.Set.ID = fullfits[[i]]$Probe.Set.ID
    }
    ##
    return(out)
}



getAnnotDescriptions <- function(annot.data)
{
    annotSub = data.table(copy(annot.data), "Probe.Set.ID")
    descrip.cols <- c(which(colnames(annot.data)=="Probe.Set.ID"), which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data))
    annotSub = annotSub[,descrip.cols, with=F]
    annotSub$Probe.Set.ID = as.character(annotSub$Probe.Set.ID)
    setkey(annotSub, "Probe.Set.ID")
    return(annotSub)
}
    
stackedPQ <- function(annot.data, results)
{
    annotSub = getAnnotDescriptions(annot.data)
    acop = data.table(copy(results))
    acop$Probe.Set.ID = as.character(acop$Probe.Set.ID)
    setkey(acop, "Probe.Set.ID")
    
    acop = annotSub[acop]
    
    k = melt.data.table(acop,
                        measure.vars = c("LogP.Anova.Cross",
                                         "LogP.Anova.Diet",
                                         "LogP.Anova.CrossByDiet",
                                         "LogP.Anova.Cross.with.interaction",
                                         "LogP.Anova.Diet.with.interaction"),
                        value.name="LogP.Anova",
                        variable.factor=F)
    
    k$variable[k$variable=="LogP.Anova.Cross"]                  = "strain"
    k$variable[k$variable=="LogP.Anova.Diet"]                   = "diet"
    k$variable[k$variable=="LogP.Anova.CrossByDiet"]            = "diet:strain"
    k$variable[k$variable=="LogP.Anova.Cross.with.interaction"] = "strain+diet:strain"
    k$variable[k$variable=="LogP.Anova.Diet.with.interaction"]  = "diet+diet:strain"
    
    k[,pval:=10^-k$LogP.Anova]
    k[,qval:=p.adjust(pval, method="fdr"), by = variable]
    k[,Log.qval:= -log10(qval), by= variable]
    
    setkey(k, "Probe.Set.ID", "variable")
    return(k)
}


adjust.pval.zscores <- function(stackedPQdata)
{
    for(vartype in unique(stackedPQdata$variable))
    {
        print(vartype)
        pvals = stackedPQdata[variable==vartype]$pval
        zscores = qnorm(pvals)
        med = median(zscores)
        zscores = zscores - med
        stackedPQdata$pval[variable==vartype] = pnorm(zscores)
    }
    return(stackedPQdata)
}

adjustSnordNames <- function(snords, stackedPQdata)
{
    g1 = unlist(lapply(strsplit(stackedPQdata$gene_name, split=","), "[",1))
    g2 = unlist(lapply(strsplit(stackedPQdata$gene_name, split=","), "[",2))
    snordId1 = snords$snord[match(g1, snords$name)]
    snordId2 = snords$snord[match(g2, snords$name)]
    
    gene_name_sn= g1
    gene_name_sn[!is.na(snordId1)] = paste0("s",snordId1[!is.na(snordId1)])

    gene_name_sn2 = g2
    gene_name_sn2[!is.na(snordId2)] = paste0("s", snordId2[!is.na(snordId2)])
    
    
    gene_name_sn[stackedPQdata$twogenes] = paste0(gene_name_sn[stackedPQdata$twogenes],
                                             ",",
                                              gene_name_sn2[stackedPQdata$twogenes])
    stackedPQdata$gene_name_orig = stackedPQdata$gene_name
    stackedPQdata$gene_name = gene_name_sn
    stackedPQdata$isSnord = !is.na(snordId1)|!is.na(snordId2)
    return(stackedPQdata)
}

getAnovaWrapper <- function(fit)
{
    fit.lmer           = (class(fit) =="merModLmerTest")
    if(!fit.lmer)
    {
        return(list(an = stats::anova(fit),
                    pvalueCol = "p-value"))
    }
    else
    {
        return(list(
            an = lmerTest::anova(fit, type=1),
            pvalueCol = "Pr(>F)"))
    }
}

getTsWrapper <- function(fit)
{
    fit.lmer = (class(fit) =="merModLmerTest")
    
    if(!fit.lmer)
    {
        return(list(
            ts        = summary(fit)$tTable,
            estCol    = "Value",
            tvalueCol = "t-value",
            pvalueCol = "p-value",
            seCol    = "Std.Error",
            dfCol    = "DF"
        ))
    }
    else
    {
        return(list(
            ts        = coefficients(summary(fit)),
            estCol    = "Estimate",
            tvalueCol = "t value",
            pvalueCol = "Pr(>|t|)",
            seCol     = "Std. Error",
            dfCol     = "df"))
    }
}


##					methylsuff = "MethylSuff",
##					vitddef    = "VitDDef",
##					stdctrl    = "StdCtrl",
##					lowpro     = "LowPro"))

Calc.Cors <- function(phen, phen.data, exp.mat, annot.data) 
{
    d <- data.frame(Probe.Set.ID = annot.data$Probe.Set.ID)
    d$Est.Cor       <- NA
    d$LogP.Cor      <- NA
    d$Est.SpearCor  <- NA
    d$LogP.SpearCor <- NA
    d$Est.KendCor   <- NA
    d$LogP.KendCor  <- NA
    d$phen          <- phen
    x <- phen.data[[phen]]
    ok <- !is.na(x)
    x <- x[ok]
    cat(sep="", "\n", phen, ":\n")
    
    for (i in 1:ncol(exp.mat))
    {
        
        if (i %% 1000 == 0) { cat(sep="","[",i,"]") } 
        y <- exp.mat[ok, i]
        
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
    return(d)
}

Calc.Cors.wrapper <- function(phens, exp.mat.full)
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

getSnordResults <- function(results, snordfile, annot.data) 
{
    snord.data       <- read.csv(snordfile, stringsAsFactors=FALSE)
                                        #snord.data$Location <- sub("\xca", "", snord.data$Location)
    snord.data$start <- as.integer(gsub(",", "", (sub("-.*", "", snord.data$Location))))
    snord.data$end   <- as.integer(gsub(",", "", sub(".*-", "", snord.data$Location)))
    
    snord.results <- dfapply(snord.data, snord.data$mRNA.Accession, results.add.FUN = rbind,
                             FUN = function(d, ...){
                                 ai <- which(as.character(d$mRNA.Accession)==annot.data$mRNA.Accession)
                                        #				print(ai)
                                 lp <- NA
                                 beta <- NA
                                 n <- NA
                                 se <- NA
                                 degf <- NA
                                 if (1<=length(ai)) {
                                     lp   <- results$LogP.CrossInStdCtrl[ai]
                                     beta <- results$Est.CrossInStdCtrl[ai]
                                     se   <- results$SE.CrossInStdCtrl[ai]
                                     degf <- results$DF.CrossInStdCtrl[ai] 
                                     n <- annot.data$Num.Pure.Probes[ai]      
                                 }
                                 d <- cbind(d, data.frame(LogP.Cross=lp, Est.Cross=beta, SE.Cross=se, DF.Cross=degf, Num.Probes=n))
                                 
                                        # d$Anova.LogP.Cross <- results$Anova.LogP.Cross[ai]
                                        # d$Est.MainOnly.Cross <- results$Est.MainOnly.Cross[ai]
                                 return (d)
                             })
    snord.results <- snord.results[!is.na(snord.results$LogP.Cross),]
    
    snord.results$Lower95 <- snord.results$Est.Cross - qt(0.975, df=snord.results$DF.Cross)*snord.results$SE.Cross
    snord.results$Upper95 <- snord.results$Est.Cross + qt(0.975, df=snord.results$DF.Cross)*snord.results$SE.Cross
                                        #the midpoint in megabases of the relevant mRNA
    snord.results$Mid.Mb <- 0.5*(snord.results$start + snord.results$end)/10^6
    return(snord.results)
}

getExpressSummary <- function(annot.data, cov.data)
{
    b6id = cov.data$ID[cov.data$Strain =="B6.NOD"]
    nodid = cov.data$ID[cov.data$Strain =="NOD.B6"]
    
    expressSummary = data.table(data.frame(Probe.Set.ID = as.character(annot.data$Probe.Set.ID),
                                           B6.NOD.expression = apply(FUN =sum, 1, X = annot.data[,b6id])/length(b6id),
                                           NOD.B6.expression = apply(FUN =sum, 1, X = annot.data[,nodid])/length(nodid)), key = "Probe.Set.ID")
    return(expressSummary)
}

computeVarExplained <- function(fit, num.sv)
{
    lmervar = lmer.var(fit)
    response = lmervar$response
    components = lmervar$components
    if(num.sv>0)
    {
        svcols = paste0("sv.", 1:num.sv)
        components$newsv = rowSums(components[,svcols, drop=F])
        components = components[,-which(colnames(components) %in% svcols)]
    }
    df = data.frame(as.list(diag(var(components)/var(as.vector(response)))))
    colnames(df) = paste0("var.", colnames(df))
    return(df)
}

appendVarExplained <- function(fits, results, num.sv)
{
    fvar = function(fit)
    {
        fit = fit$fit
        afit = computeVarExplained(fit,num.sv)
        return(afit)
    }
    
    vars = do.call(rbind, lapply(fits,fvar))
    colnames(vars) = paste0("var.", colnames(vars))
    results = cbind(results, vars)
}

computeResids <-  function(exp.mat, sv.mat, cov.data, scanType, mc.cores = prop$mnp$mc.cores)
{
    lfunc <- function(indexes, exp.mat, sv.mat, cov.data)
    {

        fits = fit.model(js = indexes,exp.mat = exp.mat, cov.data = cov.data, checkAnova = F, scanType=scanType, numBatch=1)
        residCols = list()
        for(jj in 1:length(fits))
        {
            if(!is.null(fits[[jj]]))
            {
                residCols[[jj]] = resid(fits[[jj]]$fit)
            }
            else
            {
                residCols[[jj]] = rep(NA, nrow(exp.mat))
                rownames(residCols[[jj]]) = rownames(exp.mat)
                print(paste0("failed on", jj))
            }
        }

        residCols = do.call(cbind, residCols)
        return(residCols)
    }

    resids = mnp.lapply.wrapper(len = ncol(exp.mat), FUN = lfunc, exp.mat = exp.mat, sv.mat = sv.mat, cov.data = cov.data, mc.cores = mc.cores)
    resids = do.call(cbind, resids)
    colnames(resids) = colnames(exp.mat)
    return(resids)
}

writeLimitedCols <- function(sigpq)
{
    limitedCols = c("variable", "pval","qval","gene_id", "gene_name", "Probe.Set.ID", "chrom", "probesetStart", "probesetEnd", "minDistToImprinted", "imprintedMGI", "Crowley_brainImprinted", "Crowley_strainEffect", "Crowley_Expressed.allele","numProbes", "hasvariant", "Est.CrossInStdCtrl", "B6.NOD.expression", "NOD.B6.expression")
    limitedTable = sigpq[,limitedCols, with=F]
    setnames(limitedTable, old = "Est.CrossInStdCtrl", new = "NOD.B6.effect.direction")
    limitedTable$NOD.B6.effect.direction= 2*((limitedTable$NOD.B6.effect.direction>0) -.5)
    
##

    setkey(limitedTable, "pval")
    write.table(file=fp(prop$mnp$output, paste0("limited_p", siglevel, ".csv")), limitedTable, row.names=FALSE,sep="\t")
    for(variablelevel in unique(limitedTable$variable))
    {
        outfile = file.path(prop$mnp$output, paste0("limited_", variablelevel, "_p", siglevel, ".csv"))
        subt = limitedTable[limitedTable$variable==variablelevel,] 
        write.table(file= outfile, subt, row.names = F, sep="\t")
    }
}

generatePermutationThresholds <- function(exp.mat, sv.mat, cov.data)
{
    permResids = formPermResid(exp.mat, sv.mat, cov.data)
    print("generated permutation residuals")
    if(prop$mnp$saveIntermediates){save(list = ls(), file=fp(prop$mnp$output, "permResid.RData"))}
    pracma::toc()
    
    pracma::tic()
    perm.ps   = generatePermsFrame(permResids)
    write.table(perm.ps$threshholds, fp(prop$mnp$output, "permutationThresholds.txt"))
    rm(permResids)
    
    gc()
    if(prop$mnp$saveIntermediates){save(list = ls(), file=fp(prop$mnp$output, "perm.pvals.RData"))}
    gc()
    return(perm.ps)
}

formPermResid <- function(exp.mat, sv.mat, cov.data, mc.cores = prop$mnp$mc.cores)
{
    permResids = list()
    for(scanType in c("nostrain", "nodiet"))##, "nostrainnodiet"))
    {
        permResids[[scanType]] = computeResids(exp.mat, sv.mat, cov.data, scanType, mc.cores)
    }
    return(permResids)
}

generatePermsFrame <- function(permResids)
{
    batchPvals <- function(Y, X, X.null)
    {
        
        RSS2 = colSums((Y - X %*%      (solve(t(X)     %*%X)      %*% (t(X)%*%Y)))^2)
        RSS1 = colSums((Y - X.null %*% (solve(t(X.null)%*%X.null) %*% (t(X.null)%*%Y)))^2)
        
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

    getPvalues <- function(Y,x)
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

        pValues = batchPvals(Y, X, X.null)
        return(pValues)
    }
    
    dfs = list()
    allMinPvals = list()
    threshholds = list()
    counter = 1
    numPerms = prop$mnp$SSVA.numperm 
    shuffles = matrix(nrow = nrow(cov.data), ncol = numPerms)
    for(i in 1:numPerms)
    {
        shuffles[,i] = sample(1:nrow(cov.data), size = nrow(cov.data), replace = F)
    }
    
    for(scanType in c("nostrain", "nodiet"))##, "nostrainnodiet"))
    {
        Y = permResids[[scanType]]

        if(scanType=="nostrain")
        {
            variable = "strain"
            x  = as.vector(cov.data[["Strain"]])
        } else if (scanType =="nodiet")
        {
            variable = "diet"
            x = as.vector(cov.data[["Diet"]])
        } else if (scanType == "nostrainnodiet")
        {
            variable = "diet:strain"
            x = as.vector(interaction(cov.data[["Diet"]], cov.data[["Strain"]]))
        } else {
            stop("error!")
        }
        
        truefit = getPvalues(Y,x)
        

        ##permfits = matrix(nrow = numPerms, ncol = ncol(Y)) 
        minPvals = rep(Inf, numPerms)
        
        for(shuffle in 1:numPerms)
        {
            xPerm = x[shuffles[,shuffle]]
            pvalsVec = getPvalues(Y,xPerm)
            ##permfits[shuffle,] = pvalsVec
            minPvals[shuffle] = min(pvalsVec)
            ## permfits[shuffle] = anova(lm(y~xPerm))["x","Pr(>F)"]
        }
        ## pctGreaterThanPerm = rep(NA, ncol(Y))
        ## for(j in 1:ncol(permfits))
        ## {
        ##     pctGreaterThanPerm[j] = sum(truefit[j]>permfits[,j])/numPerms 
        ## }
        
        
        ## dfs[[counter]] = data.frame(Probe.Set.ID = colnames(permResids[[scanType]]),
        ##                             variable = variable,
        ##                             pct.greaterperm = pctGreaterThanPerm,
        ##                             truefit = truefit)
        
        dfs[[counter]] = data.frame(Probe.Set.ID = colnames(permResids[[scanType]]),
                                    variable = variable,
                                    perm.true.pval = truefit)

        allMinPvals[[variable]]          = minPvals
        for(alpha in c(.05, .01, .001))
        {
            
            threshholds[[length(threshholds)+1]]=
                data.frame(variable = variable, alpha = alpha, threshhold = unname(quantile(allMinPvals[[variable]], alpha)))
        }

        counter = counter + 1
    }

    dfs = do.call(rbind, dfs)
    threshholds = do.call(rbind, threshholds)
    ##
    dfs$log.perm.true.pval = -log10(dfs$perm.true.pval)
    return(list(trueFits = dfs, minpvals = allMinPvals, threshholds = threshholds))
}

runBehaviorCorAnalysis <- function(exp.mat.full, probesetInfo, annot.data)
{
    pracma::tic()
    print("starting correlation with behavior")
    ##read in phenotypes
    phens = readPhens()
    ##get the behavior correlations data table
    behaviorExpCors = Calc.Cors.wrapper(phens, exp.mat.full)
    behaviorExpCors = data.table(behaviorExpCors, key="Probe.Set.ID")
    
    ##merge in probeset info
    behaviorExpCors = probesetInfo[behaviorExpCors]
    
    ##merge in annotations info TODO get rid of this if we merge probesetinfo and annot.data
    annot.descrip   = getAnnotDescriptions(annot.data)
    annot.descrip   = data.table(annot.descrip, key ="Probe.Set.ID")
    behaviorExpCors = behaviorExpCors[annot.descrip]
    write.table(behaviorExpCors, file = fp(prop$mnp$output, "behaviorExp.csv"), row.names=F)
    plotBehaviorExpCors(behaviorExpCors, prop$mnp$output)
    print("done correlating with behavior!")
    pracma::toc()
}










## print("about to go into deletable stuff")
## 

## ##TODO delete or make a function to encapsulate
## imprintFile = prop$genome$imprintExperimentFile
## litFile     = prop$genome$litFile
## ensemblhost = prop$genome$biomartHostForDNAReference

## ensembl=useMart(host = ensemblhost, "ENSEMBL_MART_ENSEMBL") #corresponds to 38.75 release, the name string might have to change for other releases
## ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
## table(annot.data$mRNA.Source)

## refseq = annot.data[ annot.data$mRNA.Source=="RefSeq" ,]
## enMart = getBM(attributes=c( "refseq_mrna","chromosome_name","strand", "start_position","end_position"), 
##                filters="refseq_mrna", values=refseq$mRNA.Accession, mart=ensembl)



##TODO: filter out the redundant meta probesets
