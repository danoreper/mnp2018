lm.plotting = new.env(hash=T)

lm.plotting$getLabeller <- function(fdr.summary)
{
    labelfunc <- function(df)
    {
        summary.string = paste(paste0(df$numSig, " at ", df$threshhold), collapse = "\n")
    }
    
    fdr.summary = fdr.summary[,list(summaryString = labelfunc(.SD)),by = variable]
    labeller.inst = paste0(fdr.summary$variable, "\n", fdr.summary$summaryString)
    names(labeller.inst) = fdr.summary$variable
    
##    labeller.inst = labeller(labeller.inst)
    return(labeller.inst)
}

lm.plotting$getFDRSummary <- function(per.variable, threshs)
{
    fdr.summaries = list()
    for(thresh in threshs)
    {
        fdr.summary = per.variable[,list(numSig = sum(anova.q.value<thresh)) ,by=variable]
        fdr.summary$threshhold = as.character(thresh)
        fdr.summaries = util$appendToList(fdr.summaries, fdr.summary)
    }
    fdr.summary = rbindlist(fdr.summaries)    
    return(fdr.summary)
}

lm.plotting$plot.pval.histogram <- function(per.variable, threshs(c(.01,.05), postfix)
{
    fdr.summary = getFDRSummary(per.variable, threshs)
    labeller_inst = plotting$getLabeller(fdr.summary)
    
    
    pdf(fp(prop$mnp$output, paste0("p.value.hist.", postfix)))
    aplot = ggplot(data =per.variable)
    aplot = aplot + geom_histogram(aes(x=anova.p.value, y=..density..))
    aplot = aplot + facet_wrap(~variable, labeller =as_labeller(labeller_inst))
    print(aplot)
    dev.off()

    pdf(fp(prop$mnp$output, paste0("p.value.z.hist.", postfix)))
    aplot = ggplot(data = per.variable)
    aplot = aplot + geom_histogram(aes(x=qnorm(anova.p.value), y=..density..), binwidth = .2)
    aplot = aplot + facet_wrap(~variable, labeller = as_labeller(labeller_inst))
    aplot = aplot + stat_function(fun = dnorm, color = "red")
    print(aplot)
    dev.off()
}

##TODO combine with above
lm.plotting$plot.F.histogram <- function(per.variable, threshs = c(.01, .05))

{
    fdr.summary = getFDRSummary(per.variable, threshs)
    labeller_inst = plotting$getLabeller(fdr.summary)

        
    grid <- with(per.variable, seq(min(anova.F.value), max(anova.F.value), length = 100))
    
    fdens <- plyr::ddply(results$per.variable, "variable",
                         function(dframe)
                         {
                             data.frame( 
                                 anova.F.value = grid,
                                 density = df(grid,
                                              df1 = dframe$anova.numDF[1],
                                              df2 = dframe$anova.denDF[1]),
                                 numdf = dframe$anova.numDF[1],
                                 denDF = dframe$anova.denDF[1]
                             )
                         })

    aplot = ggplot(data = per.variable, aes(anova.F.value))
    aplot = aplot + geom_histogram(aes(y=..density..), binwidth = .2)
    aplot = aplot + facet_wrap(~variable, labeller=labeller)
    aplot = aplot + geom_line(aes(y=density), data = fdens, color = "red")
    print(aplot)
}



lm.plotting$plotVarianceExplained <- function(per.probe)
{
    colz = colnames(per.probe)[grepl(pattern="var\\.", colnames(per.probe))]
    melted = melt.data.table(per.probe, measure.vars = colz, value.name="varexp",
                             variable.factor=F, variable.name = "variable")

    aplot = ggplot(melted)
    aplot = aplot + geom_histogram(aes(x=varexp, y=..density..))
    aplot = aplot + facet_grid(.~variable )
    tiff(width=15, height = 9, fp(prop$mnp$output,"varianceHistograms.tiff"))
    print(aplot)
    dev.off()

    tiff(file.path(prop$mnp$output, "varianceCors.tiff"))
    acor = cor(results[,colz, with = F])

    corrplot(acor, method="ellipse")
    dev.off()
    tiff(file.path(prop$mnp$output, "varianceScatter.tiff"))
    plot(results[,colz])
    dev.off()
}
