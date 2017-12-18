##Phenotype analysis methods
##library(WVmisc)
library(nlme)
library(data.table)
library(lme4)
library(lattice)
library(reshape)
library(ggplot2)
library(corrplot)
library(coda)
#library(car)
#library(hglm)
library(multcomp)
library(corrplot)

source("./loadParams.R")
source("./utils.R")
source("./mnp/behavior/loadData.R")
source("./lm/lm.parsing.R")
source("./lm/fitBoxCoxModels.R")
source("./mnp/general.R")

beh.analysis = new.env(hash=T)


outb <- function(...)
{
    outf("mnp/behavior", ...)
}

beh.analysis$runAll <- function(phens)
{
    df = beh.analysis$run(phens)
    df1 = beh.analysis$.adjust.pvals(df, phens)
    fname = outb("phenNew.csv")
    df2= beh.analysis$.dispSignificantPhen(df1, fname)
    
    ## merged = beh$.mergeIntoPipelines(phen=phen)
    ## ##PCA analysis of the two pipelines
    ## beh.analysis$.evalPCAphen(merged$pipel1, phen$breedLog, "pipeline_behavior_1")
    ## beh.analysis$.evalPCAphen(merged$pipel2, phen$breedLog, "pipeline_behavior_2")
    ## beh.analysis$.plotIntraPipelineCorrelations(merged)
    return(df2)
}

beh.analysis$run = function(phen, geneExp = NULL)
{
    startlechoice = "Average"
    thename = "nogene"
    if(!is.null(geneExp))
    {
        thename = geneExp$name
    }
    df = c()

    anexpType="lightdark"
    allPhenNames = c("Total.Distance",
                     "Total.Distance.Dark",
                     "Total.Distance.Light",
                     "Pct.Time.Dark",
                     "Pct.Time.Light",
                     "Total.Transitions")
    
    covariates = " ~ 1 +  as.factor(Batch)  + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID)"
    dataSet = phen$getExperiment(anexpType)

    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))
            
    allPhenNames = c("totdist", 
                     "pctctr",
                     "avgvel",
                     "jmpcts",
                     "vertcts",
                     "boli")
    
    covariates = " ~ 1 +  as.factor(Batch) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType  =  "openfield"
    dataSet = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))
##beh.analysis$.plotPhens(allPhenNames = allPhenNames, anexpType = anexpType, dataSet = dataSet, prefix = thename)

    allPhenNames = c("D3.less.D2",
        "dist_d1",
        "dist_d2",
        "dist_d3")

    covariates = " ~ 1 +  as.factor(Batch) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType  = "cocaine"

    dataSet = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c("Ten.min", "baseline", "Ten.min.less.baseline")
    covariates   = " ~ 1 +  as.factor(Batch) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType    = "cort"
    dataSet      = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c("PctFreeze.120.less.240sec")
    covariates   = " ~ 1 +  as.factor(Batch) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType    = "tail"

    dataSet      = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c("pctimmob")
    covariates   = " ~ 1 +  as.factor(Batch) + as.factor(Arena) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType    = "swim"
    dataSet      = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c("temp.1","temp.2","Difference")
    covariates   = "~ 1 +  as.factor(Batch) + Order + Diet   + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID)"
    anexpType    = "SIH"
    dataSet      = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c("PctStranger", "TRANSTOTAL")
    covariates   = " ~ 1 +  as.factor(Batch) + Box.Stranger + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1 | Dam.ID) "
    anexpType    = "sociability" 
    dataSet      = phen$getExperiment(anexpType)
    out =  beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp)
    df = rbind(df, out)

    
    allPhenNames = startlechoice
    covariates   = " ~ 1 + as.factor(Batch) + Group  + Diet + Sire.is.b6  + Group:Sire.is.b6+ Diet:Sire.is.b6 + (1|Chamber) + (1 | Dam.ID) +(1|ID) "
    anexpType    = "startle"
    dataSet      = phen$getExperiment(anexpType)
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    
    allPhenNames = c( "as50.Average_normalized",
                      "as50.Latency_normalized")
    
    covariates   = " ~ 1 + as.factor(Batch)  + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1|Chamber) + (1 | Dam.ID) "
    anexpType    = "startle"
    dataSet      =  phen$getExperiment(anexpType)
    dataSet      = dataSet[!duplicated(dataSet$ID)]
    df = rbind(df, beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp))

    for(pp in c("PP74", "PP78", "PP82", "PP86", "PP90"))
    {
        print(pp)
        allPhenNames = startlechoice
        covariates   = " ~ 1 + as.factor(Batch) + Diet + Sire.is.b6  + Diet:Sire.is.b6 + (1|Chamber) + (1 | Dam.ID) +(1|ID) "
        anexpType    = "startle"
        dataSet      = phen$getExperiment(anexpType)
        dataSet      = dataSet[dataSet$Group==pp]

        
        out = beh.analysis$.modelPhens(allPhenNames = allPhenNames, anexpType = anexpType, covariates = covariates, dataSet = dataSet, geneExp = geneExp)
        out$phenotype = paste0(out$phenotype, "ppi", substr(pp, 3,4))
    
        df = rbind(df, out)
    }


    pipel1Exp = c("lightdark", "startle", "startleByGroup", "SIH", "swim", "cocaine","weight")
    pipel2Exp = c("openfield", "sociability", "tail", "cort")              
    
    pipel1.ind = df$experiment %in% pipel1Exp
    pipel2.ind = df$experiment %in% pipel2Exp
    df$pipeline = NA
    df$pipeline[pipel1.ind]=1
    df$pipeline[pipel2.ind]=2

##    oneBigFrame = mergeIntoPipelines(phen = phen)    
    
    return(df)    
}

## fit.lambda = getBestLambda(lambdas = lambdas,
##                            pheno = pheno,
##                            covariateModelString = covariates,
##                            dataSet = dataSet,
##                            uselme = F,
##                            normalizeBeforeTransform = T,
##                            normalizeAfterTransform = T,
##                            checkAnova = T)


##
## Model every phenotype in allPhenNames using the same covariates
##
beh.analysis$.modelPhens <- function(allPhenNames, anexpType, covariates, dataSet, geneExp = NULL) 
{
    failingphen = c()
                                            #BOX COX search through these lambdas
    lambdas = c(-3, -2,-1,-.5,0,.5,1,2,3)

    if(!is.null(geneExp))
    {
        geneName = geneExp$name
        geneExp = geneExp$geneExp

        covariates = formulaWrapper$appendEffect("geneExp", covariates)$modified.string
##        names(geneExp) = sub(names(geneExp), pattern = "Mouse\\.", replacement = "")
        geneExp = data.table(ID = names(geneExp), geneExp = geneExp)
        setkey(geneExp, "ID")

        dataSet = geneExp[dataSet]

        dataSet = dataSet[!is.na(dataSet$geneExp)]
    }

    dfs = list()
    for(pheno in allPhenNames)
    {
        print("****************************")
        print(anexpType)
        print(pheno)

        nullModelString = NULL
        fit.lambda = fit.model.bc$fit(y.mat = dataSet[[pheno]],
                                      cov.data = dataSet,
                                      covariateModelString = covariates,
                                      modelParser      = fit.model.bc$getDefaultParser(covariateModelString, nullModelString),
                                      transformParams  = fit.model.bc$getDefaultTransformParams(),
                                      checkAnova       = T,
                                      strategy = fit.modelg$getDefaultModelStrategy(anovaComparison = !is.null(nullModelString),
                                                                                    prefer.lme = F))[["phen_1"]]
        
        if(is.null(fit.lambda))
        {
            failingphen = c(failingphen, paste0(anexptype, "-", pheno))
            next
        }
        
        print("just fit")
        fit.with.interaction = fit.lambda$fit
        anov = fit.lambda$anovaWrapper
        print(anov)

        if(is.null(geneExp))
        {
            dietComps = glht(fit.with.interaction, linfct = mcp(Diet="Tukey"))
            tukeyDiets =  summary(dietComps)$test$pvalues
            names(tukeyDiets) = names((summary(dietComps))$test$tstat)
            tukeyDiets = data.frame(as.list(tukeyDiets))
            
            var1 = "Diet"
            var2 = "Sire.is.b6"
            contrast.mat = lm.parsing$form.interaction.contrast.mat(fit.with.interaction, var1, var2)
            glht.out = glht(fit.with.interaction, linfct = contrast.mat)
            tukeyDietByStrains =  summary(glht.out)$test$pvalues
            names(tukeyDietByStrains) = names(((summary(glht.out))$test)$coefficients)
            tukeyDietByStrains = data.frame(as.list(tukeyDietByStrains))

            varExp = lm.parsing$varexp(fit.with.interaction)
            varExp = var(varExp$components)/var(as.vector(varExp$response))
        }
        

        modelString = as.character(formula(fit.with.interaction))[3]
        print(modelString)

        group.pval = NA
        groupByStrain.pval = NA
        groupTry = try(anov$an["Group",        anov$pvalueCol])
        if(class(groupTry)!="try-error"&!is.na(groupTry))
        {

            group.pval = groupTry
            groupByStrain.pval = anov$an["Group:Sire.is.b6", anov$pvalueCol]
        }

        tsWrapper = lm.parsing$getTsWrapper(fit.with.interaction)

        df = data.frame(experiment        = anexpType,
                        phenotype         = pheno,
                        selectedLambda    = as.character(fit.lambda$lambda),
                        model             = modelString,
                        sire.b6.effect    = unname(lm.parsing$getFixefVec(fit.with.interaction, "Sire.is.b6TRUE")),
                        sire.b6.se        = unname(tsWrapper[["ts"]][, tsWrapper[["seCol"]]]["Sire.is.b6TRUE"]),
                        strain.pval       = anov$an["Sire.is.b6",     anov$pvalueCol],
                        diet.pval         = anov$an["Diet",           anov$pvalueCol],
                        strainByDiet.pval = anov$an["Diet:Sire.is.b6",anov$pvalueCol])
        
        if(is.null(geneExp))
        {
          
            df$group.pval        = group.pval
            df$groupByStrain.pval= groupByStrain.pval
            df$strain.varexp       = varExp["Sire.is.b6", "Sire.is.b6"]
            df$diet.varexp         = varExp["Diet","Diet"]
            df$strainByDiet.varexp = varExp["Diet:Sire.is.b6", "Diet:Sire.is.b6"]
            df$batch.varexp        = varExp["as.factor(Batch)","as.factor(Batch)"]
            df$dam.varexp          = varExp["Dam.ID.(Intercept)", "Dam.ID.(Intercept)"]
            df$resid.varexp        = varExp["resid","resid"]
        }
        
        if(!is.null(geneExp))
        {
            df$gene.pval = anov$an["gene", anov$pvalueCol]
        }
        if(is.null(geneExp))
        {
            df = cbind(df,tukeyDiets)
            df = cbind(df, tukeyDietByStrains)
        }
        
        dfs = util$appendToList(dfs, df)
    }

    dfs = do.call(rbind, dfs)
    if(length(failingphen)>0)
    {
        print("FAILING PHEN:::")
        print(failingphen)
    }
    return(dfs)
}

beh.analysis$.adjust.pvals <- function(df, phen)
{
    df = copy(df)

    pvalcols = c("strain.pval","diet.pval", "strainByDiet.pval", colnames(df)[18:29])
    for(p in c(1,2))
    {
        for(pvalcol in pvalcols)
        {
            fdrcolname       = paste0(pvalcol, ".qval.fdr")
            df[[fdrcolname]][df$pipeline==p]  = p.adjust(df[[pvalcol]][df$pipeline==p],"fdr") 
        }
    }
    return(df)
}

##model every phenotype in allmodels, fitting a frequentist lmer or nlme model
beh.analysis$fitFreqModels <- function(allmodels, dataSet)
{
    outputs = list()
    for(i in 1:nrow(allmodels))
    {
        amodel = allmodels[i,]
        output = beh.analysis$fitFreqModel(amodel, dataSet)
        outputs = util$appendToList(outputs, output)
    }
    outputs = rbindlist(outputs)
    return(outputs)
}

beh.analysis$fitFreqModel <- function(amodel, dataSet)
{

    
    specificData = dataSet[[amodel$experiment]]
    if(amodel$derivedSubset!="")
    {
        specificData = specificData[[amodel$derivedSubset]]
    }
    
    afit = fit.model.bc$fit(y.mat = specificData[[amodel$phen]],
                            cov.data = specificData,
                            covariateModelString)
}


beh.analysis$.dispSignificantPhen <- function(df,outfile) 
{
    mlt = 1000
    experiment = df$experiment
    phenotype  = df$phenotype
    model      = df$model
    selectedLambda = df$selectedLambda
    pvalcols = colnames(df)[!colnames(df) %in% c("experiment", "phenotype","sire.b6.effect","sire.b6.se", "pipeline","model","selectedLambda", colnames(df)[grepl("varexp",colnames(df))])]
    
    for(colname in pvalcols)
    {
        print(colname)
        print(df[,colname])
        df[, colname] = round(df[,colname]*mlt)/mlt
        df[,colname]  = paste0(df[,colname], cut(df[,colname],c(-.0001,.001, .01,.05,.1,1),labels = c("***","**","*",".","")))
    }

    varcols = colnames(df)[grepl("varexp",colnames(df))]
    for (cname in varcols)
    {
        df[[cname]] = sprintf("%1.3f", df[[cname]])
    }

    write.table(df,file=outfile, row.names=F, sep=",")
    return(df)
}

mergeIntoPipelines <- function(phen)
{
    browser()
    prepend = function(fram, prefix)
    {

        fram = copy(fram)
        setnames(fram, setdiff(colnames(fram),"ID"), paste0(prefix, ".", setdiff(colnames(fram),"ID"))) 
        return(fram)
    }

    makepipeline = function(phenfull, outcomes)
    {
        pipel1 = phenfull[,colnames(phenfull) %in% outcomes, with=F]
        badIds = unique(pipel1[which(rowSums(is.na(pipel1))>0)]$ID)
        print("bad ids, bad records")
        
        print(badIds)
        if(length(badIds)>0)
        {
            ##
        }
        print(data.frame(pipel1[pipel1$ID %in% badIds]))
        pipel1 = na.omit(pipel1)
                                        #pipel1$ID = NULL
        setkey(pipel1,"ID")
        return(pipel1)
    }
    
    pipelineframes = list()
    framenames   = list()
    framenames[[1]] = c("lightdark", "startleByGroup", "SIH", "swim", "cocaine", "weight")
    framenames[[2]]  = c("openfield", "sociability", "tail", "cort")
    for (pipeline in c(1,2))
    {

        for(aframename in framenames[[pipeline]])
        {
            aframe = phen$getExperiment(aframename)
            aframe = aframe[,c("ID",phen$selectedoutcome[[aframename]]),with=F] ##or rawoutcome
            print(paste0("!!!!!!!", aframename))
            aframe = prepend(copy(aframe), aframename)
            
            
            if(length(pipelineframes)>=pipeline)
            {
                pipelineframes[[pipeline]] = merge(pipelineframes[[pipeline]], aframe, all = T)
            } else {
                pipelineframes[[pipeline]] = aframe
            }
        }
        pipelineframes[[pipeline]] = makepipeline(pipelineframes[[pipeline]], colnames(pipelineframes[[pipeline]])) 
    }
    
    pipel1full = phen$breedLog[pipelineframes[[1]], all=T]
    pipel2full = phen$breedLog[pipelineframes[[2]], all=T]
    
    return(list(pipel1=pipelineframes[[1]], pipel2=pipelineframes[[2]], pipel1full = pipel1full, pipel2full = pipel2full))
}






beh.analysis$nameMap <- function()
{
    old = c(paste0("lightdark.", c("Total.Distance", "Total.Distance.Dark", "Total.Distance.Light", "Pct.Time.Dark", "Pct.Time.Light", "Total.Transitions")),
            paste0("startle.", c("as50.Average_normalized", "as50.Latency_normalized", "Average_mean_PP74", "Average_mean_PP78", "Average_mean_PP82", "Average_mean_PP86", "Average_mean_PP90")),
            paste0("SIH.", c("temp.1", "temp.2", "Difference")),
            "swim.pctimmob",
            paste0("cocaine.", c("dist_d1", "dist_d2", "dist_d3","D3.less.D2")),
            "weight.Body.Weight..g.")
    
    new   = c("Total Distance", "Distance Dark", "Distance Light", "Pct Time Dark", "Pct Time Light", "Total Transitions",
              "Startle", "Latency", "PPI74", "PPI78", "PPI82", "PPI86", "PPI90",
              "SIH_T1", "SIH_T2", "SIH_\U0394",
              "Pct Immobility",
              "Day 1 Distance", "Day 2 Distance", "Day 3 Distance",  "Day 3 - Day 2",
              "Body Weight")

    names(new) = old

    old2 = c(paste0("openfield.", c("totdist", "pctctr", "avgvel","jmpcts", "vertcts", "boli")),
             paste0("sociability.",c("PctStranger", "TRANSTOTAL")),
             "tail.PctFreeze.120.less.240sec",
             paste0("cort.", c("baseline", "Ten.min",  "Ten.min.less.baseline")))
    
    new2= c("Total Distance", "Pct Center", "Avg Velocity", "Jump Cts", "Vertical Cts", "Boli",
           "Pct Time Stranger", "Transitions",
           "Pct Immobility",
           "Basal Cort",  "10 Min Cort", '\U0394')
    names(new2) = old2

    new = c(new, new2)

    return(new)
}
    

beh.analysis$.plotIntraPipelineCorrelations <- function(merged)
{
    dfm = merged$pipel1[,lapply(.SD, mean), by = ID]
    dfm = data.frame(dfm[,-1, with=F])
    newnames = gsub(colnames(dfm), pattern = "startleByGroup", replacement = "startle") 
    colnames(dfm) = newnames
    dfm = data.table(dfm)


    dfm$startle.Average_mean_No_Stim = NULL
    
    setnames(dfm, old = names(beh.analysis$nameMap), new = beh.analysis$nameMap)
                        
    colinds = 1:ncol(dfm)
    c4 = which(colnames(dfm)=="Day 3 - Day 2")
    c1 = which(colnames(dfm)=="Day 1 Distance")
    c2 = which(colnames(dfm)=="Day 2 Distance")
    c3 = which(colnames(dfm)=="Day 3 Distance")

    colinds[c4] = c1
    colinds[c1] = c2
    colinds[c2] = c3
    colinds[c3] = c4
    dfm = dfm[,colinds, with = F]

    c1 = cor(dfm)
    write.table(c1, file = outm( "pipel1cor.csv"),row.names=T, sep=",")

    cairo_pdf(outm( "pipel1cor.pdf"))
    corrplot(c1, method="ellipse", tl.col= "black")
    dev.off()


    dfm = merged$pipel2[,lapply(.SD, mean), by = ID]
    dfm = data.frame(dfm[,-1, with=F])
    
    c1  = which(colnames(dfm)== "cort.Ten.min")
    c2  = which(colnames(dfm)== "cort.baseline")
    colinds = 1:ncol(dfm)
    colinds[c2] = c1
    colinds[c1] = c2
    dfm = dfm[,colinds]
    
    dfm = data.table(dfm)
    setnames(dfm, old = names(beh.analysis$nameMap), new = beh.analysis$nameMap)
    
    c2 = cor(dfm)
    write.table(c2, file = outm( "pipel2cor.csv"),row.names=T, sep=",")

    cairo_pdf(outm( "pipel2cor.pdf"))
    corrplot(c2, method="ellipse",tl.col = "black")
    dev.off()
}

beh.analysis$.getTopHits <- function(results, vartype)
{
    tophits = results$per.variable
    tophits = tophits[variable==vartype]
    setkey(tophits, "Probe.Set.ID")
    small = results$per.probe
    small = results$per.probe[,c("Probe.Set.ID", "gene_name"), with = F]
    setkey(small, "Probe.Set.ID")
    tophits = tophits[small]
    setkey(tophits, "anova.q.value")
    return(tophits)
}




## #Generates plots of phenotype pvalues
## beh.analysis$.plotPhenPVals <- function(df) 
## {
##     pdf(outm( "phen-pvals.pdf"))
##     dfgg = df
##     dfgg$diet.pval         = -log(dfgg$diet.pval,         base = 10)
##     dfgg$strain.pval       = -log(dfgg$strain.pval,       base = 10)
##     dfgg$strainByDiet.pval = -log(dfgg$strainByDiet.pval, base = 10)
##     dfgg = melt(dfgg, id = c("experiment","phenotype"))
##     dfgg$rank = rank(-dfgg$value)
##     aplot = ggplot(dfgg,aes(x=rank, y=value, col = variable))
##     aplot = aplot + geom_point()
##     aplot = aplot + ylab("-log_10.pval")
##     aplot = aplot + xlab("pval rank")
##     aplot = aplot + geom_line(aes(y= -log(.05, base=10)), col="black")
##     print(aplot) 
##     dev.off()
    
##     pdf(outm( "pvalsPerVariable.pdf"))

##     aplot = ggplot(dfgg,aes(x=variable, y=value, col=experiment))
##     aplot = aplot + geom_point()
##     aplot = aplot + geom_line(aes(y= -log(.05, base=10), x=as.numeric(ordered(variable))), col="black")
##     aplot = aplot + ylab("-log_10.pval")
##     aplot = aplot + xlab("variable")
##     print(aplot)
##     dev.off()
## }

## beh.analysis$.plotCokeData.bydiet <- function(phen) 
## {
##     pdf(outm("cocaine123_bydiet.pdf"))
##     coke.ggpf = data.frame(phen$frame$cocaine)
##     coke.ggpf       = coke.ggpf[ ,c("Diet", "Sire.is.b6", "dist_d1","dist_d2", "dist_d3")]
    
##     coke.ggpf       = melt(coke.ggpf, c("Diet","Sire.is.b6"))
##     coke.ggpf = data.table(coke.ggpf)
    
##     coke.ggpf = coke.ggpf[,list(
##         dist_mean = mean(value), 
##         dist_sd   = sd(value)),
##         by=c("Diet", "Sire.is.b6", "variable")]
    
##     limits = aes(ymax=dist_mean+1*dist_sd, ymin=dist_mean-1*dist_sd)
    
##     coke.ggpf  = data.table(coke.ggpf)
##     aplot = ggplot(coke.ggpf, aes(x=variable, y=dist_mean, col=Sire.is.b6, group=interaction(Sire.is.b6, Diet)))
##                                         #aplot = aplot + geom_jitter(position = position_jitter(width = .2))
##     aplot = aplot + geom_point() +geom_line() +facet_wrap(~ Diet)
##     aplot = aplot + geom_errorbar(limits, width=.25)
##     print(aplot)
##     dev.off()
## }

## beh.analysis$.plotCokeData <- function(phen) 
## {
##     pdf(outm( "cocaine123.pdf"))
##     coke.ggpf = data.frame(phen$frame$cocaine)
##     coke.ggpf = coke.ggpf[ ,c( "Sire.is.b6", "dist_d1","dist_d2", "dist_d3")]
            
##     coke.ggpf = melt(coke.ggpf, c("Sire.is.b6"))
##     coke.ggpf = data.table(coke.ggpf)
    
##     coke.ggpf = coke.ggpf[,list(
##         dist_mean = mean(value), 
##         dist_sd   = sd(value)),
##         by=c("Sire.is.b6", "variable")]
    
##     limits = aes(ymax=dist_mean+1*dist_sd, ymin=dist_mean-1*dist_sd)
    
##     coke.ggpf  = data.table(coke.ggpf)
##     aplot = ggplot(coke.ggpf, aes(x=variable, y=dist_mean, col=Sire.is.b6, group=interaction(Sire.is.b6)))
##                                         #aplot = aplot + geom_jitter(position = position_jitter(width = .2))
##     aplot = aplot + geom_point() +geom_line()
##     aplot = aplot + geom_errorbar(limits, width=.25)
##     print(aplot)
##     dev.off()
## }

## beh.analysis$.plotSIH <- function(phen) 
## {
##     pdf(outm("sih.1.2.pdf"))
##     sih.ggpf = data.frame(phen$frame$SIH)
##     sih.ggpf$Batch = factor(sih.ggpf$Batch)
##     sih.ggpf       = sih.ggpf[ ,c("ID", "Batch", "Order", "Diet", "Sire.is.b6", "temp.1","temp.2")]
##                                         #sih.ggpf$interactTerm = factor(paste0(sih.ggpf$diet, sih.ggpf$Sire.is.b6))
##     sih.ggpf       = melt(sih.ggpf, c("ID","Batch", "Order", "Diet","Sire.is.b6"))
##     aplot = ggplot(sih.ggpf, aes(x=Sire.is.b6, y=value, col=variable))
##     aplot = aplot + geom_jitter(position = position_jitter(width = .2))
##                                         #sih.ggpf       = sih.ggpf[sih.ggpf$variable %in% c("temp.1","temp.2"),]
##     aplot = aplot + geom_point()
##     print(aplot)
##     dev.off()
	
##     pdf(outm( "sih.1.2.parcoord.pdf"))
##     sih.ggpf = data.frame(phen$frame$SIH)
##     sih.ggpf$Batch = factor(sih.ggpf$Batch)
##     sih.ggpf       = sih.ggpf[ ,c("ID", "Batch", "Diet", "Sire.is.b6", "temp.1","temp.2")]
## #sih.ggpf$interactTerm = factor(paste0(sih.ggpf$diet, sih.ggpf$Sire.is.b6))
##     sih.ggpf       = melt(sih.ggpf, c("ID","Batch","Diet","Sire.is.b6"))
##     sih.ggpf  = data.table(sih.ggpf)
##     setkey(sih.ggpf, ID, variable)
##     aplot = ggplot(sih.ggpf, aes(x=variable, y=value, group=ID, col=Sire.is.b6))
## #aplot = aplot + geom_jitter(position = position_jitter(width = .2))
## #sih.ggpf       = sih.ggpf[sih.ggpf$variable %in% c("temp.1","temp.2"),]
##     aplot = aplot + geom_point()+geom_line(aes(x=as.numeric(variable), y=value))
##     print(aplot)
##     dev.off()
    
##     pdf(outm("sih.delta.pdf"))
##     sih.ggpf = data.frame(phen$frame$SIH)
##     sih.ggpf$Batch = factor(sih.ggpf$Batch)
##     sih.ggpf       = sih.ggpf[ ,c("ID", "Batch", "Diet", "Sire.is.b6", "Difference")]
## #sih.ggpf$interactTerm = factor(paste0(sih.ggpf$diet, sih.ggpf$Sire.is.b6))
##     sih.ggpf       = melt(sih.ggpf, c("ID","Batch","Diet","Sire.is.b6"))
##     aplot = ggplot(sih.ggpf, aes(x=Sire.is.b6, y=value, col=variable))
##     aplot = aplot + geom_jitter(position = position_jitter(width = .2))
## #sih.ggpf       = sih.ggpf[sih.ggpf$variable %in% c("temp.1","temp.2"),]
##     aplot = aplot + geom_point()
##     print(aplot)
##     dev.off()
    
##     pdf(outm( "sih1.vs.order"))
##     plot(phen$frame$SIH$Order, phen$frame$SIH$temp.1,xlab="order", ylab="sih1")
##     dev.off()
## }

## beh.analysis$.plotVarianceExplained <- function(df1)
## {
##     colz = colnames(df1)[grepl(pattern = "varexp", colnames(df1))]
## ## colz = setdiff(colz, "var.intercept")
##     melted = melt.data.table(data.table(df1), measure.vars = colz, value.name="varexplained", variable.factor=F, variable.name="vartype")
##     aplot = ggplot(melted)
##     aplot = aplot + geom_histogram(aes(x=varexplained, y=..density..))
##     aplot = aplot + facet_grid(.~vartype )
##     pdf(width=15, height = 9, outm( "behaviorVarianceHistograms.pdf"))
##     print(aplot)
##     dev.off()

##     pdf(outm( "varianceCors.pdf"))
##     acor = cor(df1[,colz])

##     corrplot(acor, method="ellipse")
##     dev.off()
##     pdf(outm( "behviorVarianceScatter.pdf"))
##     plot(df1[,colz])
##     dev.off()
## }

## #Takes results of generating a full phenotype and evaluates them by looking at pca by strain, by diet, by strain by diet.
## beh.analysis$.evalPCAphen <- function(pipel, breedLog, prefix)
## {
##     id1 = as.numeric(phen$breedLog$ID[breedLog$Pipeline==1])
##     id2 = as.numeric(phen$breedLog$ID[breedLog$Pipeline==2])
##                                         #TODO bring back if we want to do analysis on single frame

##     ids = pipel$ID
##     pipel$ID = NULL
    
##     pcinfo = princomp(pipel, cor = T)
##     pipelall = cbind(ID = ids, pipel, pcinfo$scores)
##     setkey(pipelall, "ID")
##     setkey(breedLog, "ID")
    
##     pipelall = breedLog[pipelall]
##     pipelall$label = NA
##     pipelall$label[as.character(pipelall$Diet)=="Low Protien"]               = "p"
##     pipelall$label[as.character(pipelall$Diet)=="Control B: Vit D, Low Pro"] = "s"
##     pipelall$label[as.character(pipelall$Diet)=="Vitamin D Deficient"]       = "d"
##     pipelall$label[as.character(pipelall$Diet)=="Control A : Methyl"]        = "m"
    
        
##     aplot = ggplot(pipelall, aes(x=Comp.1, y=Comp.2, shape=as.factor(Sire.is.b6), label=label, color=label))
##     aplot = aplot + geom_text(size=5)
##     pdf(outm( paste0(prefix, "_dietpca.pdf")))
##     print(aplot)
##     dev.off()
    
##     aplot = ggplot(pipelall, aes(x=Comp.1, y=Comp.2, shape=label, label=label, color=as.factor(Sire.is.b6)))
##     aplot = aplot + geom_text(size=5)
##     pdf(outm( paste0(prefix, "_dietstrainpca.pdf")))
##     print(aplot)
##     dev.off()
## }
