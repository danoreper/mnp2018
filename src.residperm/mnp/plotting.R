library(corrplot)
library(ggplot2)
library(data.table)
library(grid)
library(ggrepel)
library(oligo)
source("./multipleTesting.R")

resol = 600
plotting = new.env(hash=T)

plotting$plotControls <- function(annot.data, probeTable, outdir) 
{
    mousecols = which(grepl("Mouse", colnames(annot.data)))
    controlStrings = levels(probeTable$fulltype)
    for(probetype in c(controlStrings, "regular"))
    {
        relevantProbeSets = unique(probeTable[fulltype==probetype]$Probe.Set.ID)

        try(
        {
            ii = annot.data$Probe.Set.ID %in% relevantProbeSets
            exp = annot.data[ii, mousecols]
            mouses = as.factor(colnames(exp))
            jj = rank(rowSums(as.matrix(exp)))

            ##Gene bucket is a way of grouping genes as it is infeasible to plot expression of
            ##every single probeset
            genebucket = cut(rowSums(as.matrix(exp))/ncol(exp), 100)
            genebucket = rep(genebucket, length(mouses))

            ##genebucket = cut(c(as.matrix))
            out        = c(as.matrix(exp))
            gene       = rep(jj, length(mouses))
            mouseInd   = rep(mouses, each = nrow(exp))

            df = data.frame(gene = gene, out=out, genebucket=genebucket, mouseInd=mouseInd)
            df = data.table(df)
            df2 = df[,list(mnexp.by.bucket.mouse=mean(out)),by=c("mouseInd", "genebucket")]
            
            dir.create(fp(outdir, "probeplots"), showWarnings = F)
##            
            pdf(fp(outdir, "probeplots", paste0(probetype,".pdf")))
            aplot = ggplot(aes(y=mnexp.by.bucket.mouse, x = genebucket), data=df2)
            aplot = aplot + geom_point()
            aplot = aplot + facet_wrap(~mouseInd)
            print(aplot)
##            xyplot((mnexp.by.bucket.mouse~genebucket|mouseInd), data=df2[,])
            dev.off()
        }
        )
    }
}

plotting$plot.cel.files <- function()
{

    celdir = datm(prop$mnp$cel.dir)
    cfnames =  list.celfiles(celdir)
    badcfnames = sort(c(25, 105,14,156,17,180,182,194,1,43,59,97))
    badcfnames = paste0(badcfnames, ".CEL")

    cfins  = fp(celdir, cfnames)

    celoutdir = outm("celplots")
    dir.create(celoutdir)
    cfouts = paste0(fp(celoutdir, cfnames), ".pdf")

    all = read.celfiles(cfins)
    for (i in 1:length(cfins))
    {
        print(i)
        ##cfin  = cfins[i]
        cfout = cfouts[i]
        ##x = read.celfiles(cfin)
        pdf(cfout)
        image(all, which = i, transfo = rank)
        dev.off()
    }
}

plotting$pcaExpression <- function(exp.mat, cov.data, output)
{
    pcRes = prcomp(exp.mat,scale. = T)
    pcs = (as.matrix(exp.mat)%*%(pcRes$rotation))

    
    fullthing = cbind(cov.data, pcRes$x)
    fullthing$label = as.character(fullthing$Diet)
    fullthing$label[as.character(fullthing$label)=="PD"]  = "p"
    fullthing$label[as.character(fullthing$label)=="Std"] = "s"
    fullthing$label[as.character(fullthing$label)=="VDD"] = "d"
    fullthing$label[as.character(fullthing$label)=="ME"]  = "m"
    
    aplot = ggplot(fullthing[fullthing$Pipeline==1,], aes(x=PC1, y=PC2, shape=Strain, label=label, color=label))
    aplot = aplot + geom_text(size=5)
    pdf(outm( "dietpca.pdf"))
    print(aplot)
    dev.off()
    
    aplot = ggplot(fullthing, aes(x=PC1, y=PC2, shape=Strain, label=label, color=Strain))
    aplot = aplot + geom_text(size=5)
    pdf(outm("dietstrainpca.pdf"))
    print(aplot)
    dev.off()
}

##facet plot of batches of highest genes-- [1:batchSize], [batchSize+1:2*batchSize], ... [k*batchSize:len]
plotting$plotHighestGenes <- function(stackedPQdata,
                                      output = prop$mnp$output,
                                      prefix="", len = 30,
                                      batchSize = 5,
                                      exp.mat,
                                      cov.data,
                                      logScale = T,
                                      nodietfacet = F,
                                      pooled.diet = F)
{
    createRankedDf.for.type <- function(sigpq, cov.data, vartype, rankGroup)
    {
        dfs = list()
        for (therank in rankGroup)
        {
            print(therank)
            probesetIndex = sigpq$variable==vartype & sigpq$origRank==therank
            if(sum(probesetIndex)>0)
            {
                plotdf            = cov.data
                plotdf$rank       = therank
                plotdf$gene_name  = sigpq[probesetIndex]$gene_name
                probeset          = as.character(sigpq[probesetIndex]$Probe.Set.ID)

                expressionRows = match(plotdf$ID, rownames(exp.mat))
                plotdf$expression = exp.mat[expressionRows,probeset]

                dfs = util$appendToList(dfs, plotdf)
            }
        }
        
        plotdf      = rbindlist(dfs)
        if(nrow(plotdf)>0)
        {
            plotdf$rank = convertRankToGeneName(plotdf)
        }
        return(plotdf)
    }

     ## we want to use rank to order genes according to how significant they are, but we want the rank labels in the plot to be the gene names themselves
    convertRankToGeneName <- function(plotdf)
    {
        therank = as.factor(plotdf$rank)
        oldlevs = levels(therank)
        newlevs = c()
        for(oldlev in oldlevs)
        {
            newlev = as.character(plotdf$gene_name[which(therank == oldlev)[1]])
            newlevs = c(newlevs, newlev)
        }
        levels(therank)=newlevs
        return(therank)
    }
    
    sigpq <- stackedPQdata[anova.q.value <= prop$mnp$alphalevel,]
##    
    sigpq = droplevels(sigpq)
    sigpq = data.table(sigpq)
    sigpq[,origRank:=rank(anova.q.value, ties.method="first"), by= variable]

    rankGroups = util$getIndexGroupsForLen(len, batchSize)
    for(i in 1:length(rankGroups))
    {
        rankGroup = rankGroups[[i]]

        for(vartype in levels(as.factor(sigpq$variable)))
        {
            print(vartype)
            
            plotdf = createRankedDf.for.type(sigpq, cov.data, vartype, rankGroup)
            if(!logScale)
            {
                plotdf$expression = 2^(plotdf$expression)
            }

                
            x = try({
            if(nrow(plotdf)==0)
            {
                next
            }})
            
            print(x)
            if(class(x)=="try-error")
            {
                
            }
            if(pooled.diet)
            {
                aplot = ggplot(plotdf, aes(x = "all diets/strains", color=Strain, y=expression))
                aplot = aplot + facet_grid(rank ~ .)
            } else {
                
            if(vartype!="diet" & !nodietfacet)
            {
                aplot = ggplot(plotdf, aes(x=Strain, y=expression))
                aplot = aplot + facet_grid(rank ~ Diet)
            } else {
                aplot = ggplot(plotdf, aes(x = Diet, y=expression, color = Strain))
                aplot = aplot + facet_grid(rank ~ .)
            }
            }
            ##aplot = aplot + geom_point()
            aplot = aplot + geom_jitter(width=.3)

            ##aplot = aplot + ggtitle(paste0(vartype,'_',geneName))
            afile = fp(output, paste0(prefix,"highest_",vartype,paste(rankGroup, collapse=","),"_", nodietfacet,"_", pooled.diet, ".pdf"))
            pdf(afile)
            print(aplot)
            dev.off()
        }
    }
}


plotting$prepForPlot <- function(subd, xVar, yVar, 
                                 chromz = c(as.character(1:19), "X","Y"),
                                 bp.sep = 75000000,
                                 karyo)
{
    chromz = chromz[chromz %in% subd$chrom]
    
    lens = karyo[chromz]$len + bp.sep

    
    starts = c(0, cumsum(as.numeric(lens)))
    starts = starts[-length(starts)]
    names(starts) = chromz
    mids = starts + (lens-bp.sep)/2
    names(mids) = chromz

   
    subd$chrom   = as.character(subd$chrom)
    subd         = subd[chrom %in% chromz,]
    subd$chrom   = factor(subd$chrom, levels= chromz)

    ##subd$labz    = as.character(subd$Probe.Set.ID)
    subd$labz    = as.character(subd$gene_name)
    subd$labz[subd$labz=="Lrrc16a"] = "Carmil1"
    
    ##subd$labz    = as.character(subd$Probe.Set.ID)
    as.character(subd$newname) ###subd$gene_name) ###subd$Probe.Set.ID
    ## subd$thenudge =  (nchar(str_trim(subd$labz))) + .5
    ## subd$thenudge = subd$thenudge + (subd$thenudge - 5) *.75 
    ## subd$thenudge = (ceiling(subd$thenudge/2))*20000000
    ## subd$thenudge[subd$labz == "Ndn"]           = (subd$thenudge[subd$labz == "Ndn"] + .9*25000000)
    ## subd$thenudge[subd$labz == "3830406C13Rik"] = subd$thenudge[subd$labz == "3830406C13Rik"] + 1.2*25000000

    
    subd$x       = eval(parse(text = xVar), envir = subd)
    subd$y       = eval(parse(text = yVar), envir = subd)
    subd$isinf   = is.infinite(subd$y)
    subd$y[subd$isinf] = 1.1*max(subd$y[!subd$isinf])

    ##subd$y   = subd[,yVar]
    
    subd$x.man = subd$x + starts[as.character(subd$chrom)]
    
    return(list(subd = subd, mids = mids))
}

plotting$plotPermScan <- function(full,
                                  ##per.probe,
                                  per.variable,
                                  per.level,
                                  outdir,
                                  permthresh=NULL,
                                  sz=1,
                                  point = T,
                                  filterpostfix= NULL,
                                  showGeneNames = F,
                                  alphalevel=prop$mnp$alphalevel,
                                  karyo)
{
    figwid = 8.5
    figheight = 5
    

    maxmult = 1.03

    xVar.old = "probesetStart"
    xVar = "x"
    
    ##    yVar.old = "-log10(perm.anova.p.value)"
    yVar.old = "-log10(anova.p.value)"
    yVar = "y"


    #atitle = bquote(italic(.(assay))*": microarray")


    
    ylab.str = bquote("-log"[10] *"(p-value)")

    setkey(per.level, "Probe.Set.ID")
    merged = full
    
    xVar = "x.man"
    setkey(merged, "Probe.Set.ID", "variable")
    threshTable = permthresh ##permthresh$permStatistics


    merged = merged[anova.p.value<=.05]

    subd.n.mids = plotting$prepForPlot(merged, xVar.old, yVar.old, bp.sep = 0, karyo = karyo)
    
    xmax = 1.01*max((subd.n.mids$subd$x.man))

    ##    per.level$coef.Value = per.level$coef.Value*per.level$effectMult

    per.level$coef.Value = per.level$coef.Value

    
    
    ##TODO fixme
    for(plot.type in c("Strain", "Diet.2", "Diet.4","Diet:Strain"))##, "Diet.2", "Diet:Strain")) #, "Diet.1"))
    {
        print("plotting type:")
        print(plot.type)
        
        tokens   = strsplit(plot.type, split="\\.")
        vartype  = tokens[[1]][1]
        modetype = tokens[[1]][2]
        print(vartype)

        subd        = merged[variable==vartype]
        subd.n.mids = plotting$prepForPlot(subd, xVar.old, yVar.old, bp.sep = 0, karyo = karyo)

        lower.bound = 1.3

                
        if(vartype =="Strain")
        {
           
            justvar = per.level[variable.level=="StrainNOD.B6", j=list(level1 = coef.Value, Probe.Set.ID)]
            setkey(justvar, "Probe.Set.ID")
        }
        if(vartype=="Diet")
        {
            browser()
            ##TODO extract this into a method for results writing
            v1 = per.level[variable.level=="DietME", j = list(level1 = coef.Value, Probe.Set.ID)]
            v2 = per.level[variable.level=="DietVDD",    j = list(level1 = coef.Value, Probe.Set.ID)]
            v3 = per.level[variable.level=="DietPD",     j = list(level1 = coef.Value, Probe.Set.ID)]

            justvar = data.table(Probe.Set.ID = v1$Probe.Set.ID,
                                 methyl.more = (v1$level1>0)*1 + ((v1$level1 - v2$level1)>0)*1 + ((v1$level1 - v3$level1)>0)*1,
                                 methyl.less = (v1$level1<0)*1 + ((v1$level1 - v2$level1)<0)*1 + ((v1$level1 - v3$level1)<0)*1)
                
            setkey(justvar, "Probe.Set.ID")

            
            methsuff.colz  = c("ME.vs.Std.p.value",  "PD.vs.ME.p.value",  "VDD.vs.ME.p.value")
            otherdiet.colz = c("VDD.vs.Std.p.value", "VDD.vs.PD.p.value", "PD.vs.Std.p.value")


            table(rowSums((subd[,methsuff.colz, with = F]<=alphalevel)))
            table(rowSums((subd[,otherdiet.colz, with = F]<=alphalevel)))
            
            s1 = rowSums((subd.n.mids$subd[,methsuff.colz, with = F]<=alphalevel))
            s2 = rowSums((subd.n.mids$subd[,otherdiet.colz, with = F]<=alphalevel))
            
            subd.n.mids$subd$s1 = s1
            subd.n.mids$subd$s2 = s2
            subd.n.mids$subd = justvar[subd.n.mids$subd]
            
            subd.n.mids$subd$pctMethSig = paste0(s1, "/(",s1, "+", s2,")")
            rm(pctMethSig)
        }
        
        if(vartype=="Diet:Strain")
        {
            ##TODO: implement properly when strain by diet contrast tests implemented.
            justvar = per.level[variable.level=="StrainNOD.B6", j=list(level1 = coef.Value, Probe.Set.ID)]
            setkey(justvar, "Probe.Set.ID")
        }

        allpoints        = justvar[subd.n.mids$subd]
        allpoints$anova.BY.value = p.adjust(allpoints$anova.p.value, method = "BY")

        gotThresh = F
        if(!is.null(threshTable))
        {
 
            threshTable.sub = threshTable[variable==vartype,]

            
            if(nrow(threshTable.sub)>0)
            {
                ind              = util$lookupByFloat(threshTable.sub, floatkeyCol="alpha", floatkey = alphalevel)
                thresh.alpha     = threshTable.sub[ind,]$threshhold
                thresh.alpha.gev = threshTable.sub[ind,]$threshhold.gev
                gotThresh = T
                
                allpoints$labz[allpoints$isinf] = paste0("*",allpoints$labz[allpoints$isinf])
                if(vartype != "Diet:Strain")
                {
                    aboveThresh      = allpoints[allpoints$y>-log10(thresh.alpha.gev),]
                } else {
                    aboveThresh      = allpoints[allpoints$y>-log10(thresh.alpha.gev)
                                                 |allpoints$gene_name=="Meg3",]
                }
                belowThresh      = allpoints[allpoints$y<=-log10(thresh.alpha.gev),]

            }

        }

        
        partial.inf      = allpoints[allpoints$isinf,]
        p.thresh.for.q = multipleTesting$get.empirical.p.value.for.q(allpoints$anova.p.value,
                                                                     allpoints$anova.q.value,
                                                                     alphalevel)
        p.thresh.for.BY = multipleTesting$get.empirical.p.value.for.q(allpoints$anova.p.value,
                                                                      allpoints$anova.BY.value,
                                                                      alphalevel)
        if(!gotThresh)
        {
            
            belowThresh = allpoints
            aboveThresh = allpoints[0,]
            partial.inf = allpoints[0,]
            thresh.alpha.gev = p.thresh.for.q
        }

        if(vartype =="Diet" && modetype == 4)
        {
            if(gotThresh)
            {
                ind              = util$lookupByFloat(threshTable.sub, floatkeyCol="alpha", floatkey = alphalevel)
                lower.bound = -log10(threshTable.sub[ind,]$threshhold.gev)
                lower.bound = sprintf( "%.1f",lower.bound)
                lower.bound = as.numeric(lower.bound)
            } else {
                lower.bound = -log10(p.thresh.for.q)
            }
        }


        default.y.breaks = c(lower.bound,seq(2, max(c(3,allpoints$y)), by = 2))
        default.y.labels = as.character(default.y.breaks)
        
        if(vartype=="Strain")
        {
            xoffset = 60000000/8.5 * figwid
            textmult = 5.5##3.5

            y.breaks = c(lower.bound,2, seq(6, max(max(allpoints$y),8), by = 2))
            y.labels = as.character(y.breaks)

            ##include? or comment out?
            y.breaks = c(y.breaks, -log10(p.thresh.for.q))
            y.labels = c(y.labels, "")

            y.breaks = c(y.breaks, -log10(thresh.alpha.gev))
            y.labels = c(y.labels, "")

            
            aplot = ggplot(belowThresh, aes(color = as.factor(as.integer(chrom)%%2)))
##            aplot = ggplot(belowThresh)
            aplot = aplot + geom_point(aes_string(x=xVar, y=yVar), size=sz*2.0)
            
            aplot = aplot + geom_point(data = aboveThresh, aes_string(x=xVar, y=yVar, shape="level1>0"), size=sz*2.0)
            aplot = aplot + scale_shape_manual(labels = c("B6xNOD > NODxB6"," NODxB6 > B6xNOD"), values = c(1,15), guide = guide_legend(title=""))

            
            aplot = aplot + annotation_custom(textGrob(sprintf("%.1f", -log10(p.thresh.for.q)), gp = gpar(col = "black")), 
                                              xmin=-xoffset, xmax=-xoffset,
                                              ymin=-log10(p.thresh.for.q)-.25,
                                              ymax=-log10(p.thresh.for.q)-.25)

            aplot = aplot + annotation_custom(textGrob(sprintf("%.1f", -log10(thresh.alpha.gev)), gp = gpar(col = "black")), 
                                              xmin=-xoffset, xmax=-xoffset,
                                              ymin=-log10(thresh.alpha.gev)+.25,
                                              ymax=-log10(thresh.alpha.gev)+.25)
            
            ##aplot = aplot + guides(shape =F)


            aboveThresh$offset.x = 0
            aboveThresh$offset.y = 0

            aboveThresh[Probe.Set.ID == "10553833"]$offset.y = -.05

         
            aboveThresh[Probe.Set.ID == "10547056"]$offset.y = -.1
            aboveThresh[Probe.Set.ID == "10563949"]$offset.y = 0#.3

         
            aboveThresh[Probe.Set.ID == "10398426"]$offset.y = +.257
            aboveThresh[Probe.Set.ID == "10398360"]$offset.y = -.257


            aboveThresh[Probe.Set.ID == "10563989"]$offset.y = +.2
            
            aboveThresh[Probe.Set.ID == "10563911"]$offset.x = -1.1 * 200000000
            aboveThresh[Probe.Set.ID == "10547056"]$offset.x = -1.75 * 200000000

        }
        
        if(vartype=="Diet")
        {
            y.breaks = default.y.breaks
            y.labels = default.y.labels

            y.breaks = c(y.breaks, -log10(p.thresh.for.q))
            y.labels = c(y.labels, sprintf("%.1f", -log10(p.thresh.for.q)))

            if(modetype == 2)
            {
                y.breaks = c(y.breaks, -log10(thresh.alpha.gev))
                y.labels = c(y.labels, sprintf("%.1f", -log10(thresh.alpha.gev)))
            }
            
            delt.x = 20000000
            delt.x.unitless = delt.x/diff(range(allpoints$x.man))
            delt.x.unitless = delt.x.unitless*(figwid/figheight)
            
            delt.y = delt.x.unitless*max(allpoints$y*maxmult - lower.bound)
            textmult = 3.5##3.5
            

            aboveThresh$offset.x = 0
            aboveThresh$offset.y = 0

            if(modetype == 2|| modetype == 4)
            {

                if(modetype == 2)
                {
                    aplot = ggplot(belowThresh, aes(color = as.factor(as.integer(chrom)%%2)))
                    aplot = aplot + geom_point(aes_string(x=xVar, y=yVar), size=sz*2)
                } 
                aboveThresh$methyl.more = as.character(aboveThresh$methyl.more)
                if(any(aboveThresh$methyl.more %in% c(1,2)))
                {
                    aboveThresh$methyl.more[aboveThresh$methyl.more %in% c("1","2")] = "neither"
                    aboveThresh$methyl.more = factor(aboveThresh$methyl.more, c("0", "3", "neither"))
                    shapes = c(6,17,3)
                    shapeLabels = c("Methyl < every other diet", "Methyl > every other diet", "Neither True")
                } else {
                    aboveThresh$methyl.more = factor(aboveThresh$methyl.more, c("0", "3"))
                    shapes = c(6,17)
                    shapeLabels = c("Methyl < every other diet", "Methyl > every other diet")
                }

                if(modetype == 4)
                {
                    aplot = ggplot(data = aboveThresh, aes(x=x.man, y=y, shape=methyl.more), size = sz*2.0)
                    aboveThresh[Probe.Set.ID == "10360806"]$offset.x = -.07 *200000000
                    aboveThresh[Probe.Set.ID == "10360806"]$offset.y = -.05
##                    aboveThresh[Probe.Set.ID == "10485514"]$offset.y = +.02


                    aboveThresh[grepl(gene_name, pattern = "Rif1")]$offset.x =  -2*20000000
                    aboveThresh[grepl(gene_name, pattern = "Rif1")]$offset.y = .06
                    
                    aboveThresh[grepl(gene_name, pattern = "Disp1")]$offset.x = -8.0*20000000
                    aboveThresh[grepl(gene_name, pattern = "Disp1")]$offset.y = .053
                                        
                   
                    
                    aboveThresh[Probe.Set.ID == "10550383"]$offset.y = +.01
                    
                    aboveThresh[Probe.Set.ID == "10537909"]$offset.x = -.03 * 200000000
                    aboveThresh[Probe.Set.ID == "10537909"]$offset.y = -.02

                    aboveThresh[Probe.Set.ID == "10525185"]$offset.x = -1.28 * 200000000
                    ##                    aboveThresh[Probe.Set.ID == "10525185"]$offset.y = -.06

                    aboveThresh[Probe.Set.ID == "10404069"]$offset.y = -.02
                    aboveThresh[Probe.Set.ID == "10381419"]$offset.y = +.02

                    aboveThresh[Probe.Set.ID == "10579089"]$offset.y = +.025
                    aboveThresh[Probe.Set.ID == "10454966"]$offset.x = -.4 *200000000
##                    aboveThresh[Probe.Set.ID == "10548966"]$offset.y = +.03

                    
                    ##                    aboveThresh[Probe.Set.ID == "10438726"]$offset.y = +.02
                    aboveThresh[Probe.Set.ID == "10438726"]$offset.x = -1.1 *200000000

                    aboveThresh[Probe.Set.ID == "10598422"]$offset.x = -1 *200000000

                    

                    
                    
                    aplot = aplot + geom_point()
                }

                if(modetype == 2)
                {
                    aplot = aplot + geom_point(data = aboveThresh, aes(x=x.man, y=y, shape = methyl.more), size=sz*2.0)
                }

                aplot = aplot + scale_shape_manual(values = shapes,
                                                   labels = shapeLabels,
                                                   guide = guide_legend(title ="", reverse = T))
            } 
            


        }

        if(vartype=="Diet:Strain")
        {
            y.breaks = default.y.breaks
            y.labels = default.y.labels

            y.breaks = c(y.breaks, -log10(p.thresh.for.q))
            y.labels = c(y.labels, sprintf("%.1f", -log10(p.thresh.for.q)))
            
            y.breaks = c(y.breaks, -log10(thresh.alpha.gev))
            y.labels = c(y.labels, sprintf("%.1f", -log10(thresh.alpha.gev)))

            textmult = 3.5##3.5
            aplot = ggplot(allpoints, aes(color = as.factor(as.integer(chrom)%%2)))
            aplot = aplot + geom_point(aes_string(x=xVar, y=yVar), size=2*sz)
            

            
            aboveThresh$offset.x = 0
            aboveThresh$offset.y = 0
            
            aboveThresh[Probe.Set.ID == "10515714"]$offset.y = .1

            aboveThresh[Probe.Set.ID == "10517559"]$offset.y = -.06

            ## aboveThresh[Probe.Set.ID == "10558134"]$offset.y = +.025

            ## aboveThresh[Probe.Set.ID == "10569163"]$offset.y = -.025

            ## aboveThresh[Probe.Set.ID == "10384579"]$offset.y = +.06

            ## aboveThresh[Probe.Set.ID == "10407907"]$offset.y = +.065


            ## aboveThresh[Probe.Set.ID == "10434285"]$offset.y = +.025

            ## aboveThresh[Probe.Set.ID == "10601551"]$offset.x = -3.2*20000000
            ## aboveThresh[Probe.Set.ID == "10601551"]$offset.y = +.035

             aboveThresh[Probe.Set.ID == "10398326"]$offset.y = -.08
            

            ##            aplot = aplot + scale_shape_manual(labels = c("B6xNOD > NODxB6"," NODxB6 > B6xNOD"), values = c(1,15), guide = guide_legend(title=""))

            
            ## methsuff.colz  = c("MethylSuff...StdCtrl", "MethylSuff...LowPro", "VitDDef...MethylSuff")
            ## otherdiet.colz = c("VitDDef...StdCtrl", "VitDDef...LowPro", "LowPro...StdCtrl")
            ## table(rowSums((subd[,methsuff.colz, with = F]<=alphalevel)))
            ## table(rowSums((subd[,otherdiet.colz, with = F]<=alphalevel)))
            
            ## s1 = rowSums((subd.n.mids$subd[,methsuff.colz, with = F]<=alphalevel))
            ## s2 = rowSums((subd.n.mids$subd[,otherdiet.colz, with = F]<=alphalevel))

            
        }
        
        ##aplot = aplot + geom_point(data = partial.inf, aes_string(x=xVar, y=yVar), shape = 21, size=2*sz, stroke = 2*sz, color = "black")
        
        aplot = aplot + scale_color_manual(values = c("grey50", "black"))
        aplot = aplot + guides(color=FALSE)
        
        
        
        aplot = aplot + geom_hline(yintercept=-log10(p.thresh.for.q), linetype="longdash")
##        aplot = aplot + geom_hline(yintercept=-log10(p.thresh.for.BY), linetype="longdash", color = "green")

        
        
        yoff = max(allpoints$y)/80
        
        scl = 30000000/sum(as.numeric(karyo$len))
        
        if(vartype!="Diet"|(vartype=="Diet" & modetype!=2))
        {
            aplot = aplot +  geom_text(data = aboveThresh,
                                       aes(x = x.man + offset.x + 20000000,
                                           y = y + offset.y,
                                           label = labz),
                                       fontface = "italic",
                                       size=textmult*sz, hjust = 0)
        }
        ## else {
        ##     aplot = aplot +  geom_text_repel(data = aboveThresh,
        ##                                      aes(x= x.man, y=y, label=labz),
        ##                                      nudge_x = 20000000,
        ##                                      max.iter = 80000,
        ##                                      size=textmult*sz)
            
        ## }
        
        

        ## aplot = aplot +  geom_text_repel(data = aboveThresh,
        ##                                  aes(x= x.man, y=y, label=labz),
        ##                                  force = .0001, #.2,
        ##                                  ##nudge_x = 15000000,
        ##                                  segment.color = NA,
        ##                                  segment.alpha = 1,
        ##                                  point.padding = unit(8*sz, "points"),
        ##                                  box.padding      = unit(-.05, 'lines'),
        ##                                  nudge_y = 0,#,yoff,
        ##                                  nudge_x = aboveThresh$thenudge,
        ##                                  max.iter = 80000,
        ##                                  size=textmult*sz)
##       print(aplot)
        ## browser()
    
        ## print(aplot.1)
        ## browser()
        ## x= 5
        ## pdf(fp(outdir, paste0(plot.type, "_", filterpostfix, ".pdf")),width=figwid, height=figheight);
        ## print(aplot.1);
        ## dev.off()
        
        ## browser()
        
        ## aplot = aplot +  geom_text_repel(data = aboveThresh,
        ##                                  aes_string(x=paste0(xVar,"+30000000"),
        ##                                             y=paste0(yVar,"+",yoff),
        ##                                             label="labz"),
        ##                                  size=textmult*sz)  
        
        print(vartype)
  
        if(vartype %in% c("Diet","Strain", "Diet:Strain") & (is.na(modetype)||modetype != 4))
        {
            aplot = aplot + geom_hline(yintercept=-log10(thresh.alpha.gev), color = "black")
        }
    

        aplot = aplot + scale_x_continuous(expand = c(0,0),
                                           limits = c(min(allpoints$x.man), xmax),
                                           breaks = subd.n.mids$mids,
                                           labels = names(subd.n.mids$mids))

        
        aplot = aplot + scale_y_continuous(expand = c(.001,0),
                                           limits=c(lower.bound, max(allpoints$y)*maxmult),
                                           breaks = y.breaks,
                                           labels = y.labels)
        
##        aplot = aplot + expand_limits(0,0)
        aplot = aplot + theme_bw()
        aplot = aplot + theme(panel.grid.major.y = element_blank())
        aplot = aplot + theme(panel.grid.major.x = element_blank())
        aplot = aplot + theme(panel.background = element_rect(fill = "white"))
        aplot = aplot + theme(axis.line = element_line(colour = "black"))
        aplot = aplot + theme(text = element_text(size=18))
        aplot = aplot + theme(axis.text = element_text(size=12))
        aplot = aplot + theme(legend.position = c(.005,.995))
        aplot = aplot + theme(legend.justification = c(0,1))
        aplot = aplot + theme(legend.text = element_text(size = 12))
        aplot = aplot + theme(legend.title=element_blank())


        
##        aplot = aplot + theme_bw()


        if(vartype == "Diet")
        {
            titl = "Diet"
        } else if(vartype == "Diet:Strain") {
            titl = "Diet-by-parent-of-origin"
        } else if(vartype == "Strain") {
            titl = "Parent-of-origin"
        }
        
        
        aplot = aplot + ggtitle(titl)
        aplot = aplot + xlab("chromosome")
        aplot = aplot + ylab(ylab.str)
        aplot = aplot + theme(panel.grid.minor.x = element_blank())
        aplot = aplot + theme(panel.grid.minor.y = element_blank())


        ##print(aplot)
        aplot.1 = ggplot_gtable(ggplot_build(aplot))
        aplot.1$layout$clip[aplot.1$layout$name=="panel"] <- "off"
        try(dev.off())
##        x11(height = figheight, width = figwid)
##        grid.draw(aplot.1)
        ##print(aplot)

        figh.mult = 1
        ## if(vartype=="Diet")
        ## {
        ##     figh.mult = 2
        ## }

        plotname = plot.type
        if(plotname == "Diet.2")
        {
            plotname = "Diet"
        }
        if(plotname == "Diet.4")
        {
            plotname = "Diet_fwersig"
        }

        if(filterpostfix!="")
        {
            fname = fp(outdir, paste0(plotname, "_", filterpostfix, ".pdf"))
        } else {
            fname = fp(outdir, paste0(plotname, ".pdf"))
        }

        pdf(fname,width=figwid, height=figheight*figh.mult)

        ##print(aplot)
        grid.draw(aplot.1)
        dev.off()
    }
}                        

plotting$buildManyScans <- function(results, outdir, thresh, karyo)
{

    print("building many scans!!!") 
    
    outdirpeak = fp(outdir,"micro","manhattan")
    
    dir.create(outdir, recursive = T, showWarnings =F)
    dir.create(outdirpeak, recursive = T, showWarnings = F)

    full = results$full
    fulls = list()
    fulls[[1]] = full
    nmz = c("", "_712", "7", "12")
    for(i in 1:length(fulls))
    {
        
        plotting$plotPermScan(fulls[[i]], results$per.variable, results$per.level,
                              outdirpeak,
                              permthresh = thresh,
                              filterpostfix = paste0(nmz[i]),
                              karyo = karyo)
    }
}

plotBehaviorExpCors <- function(cor.frame, output)
{
    pdf(fp(output, "behaviorExpCors.regular.pdf"))
    rankedps = cor.frame[,max(LogP.Cor),by="phen"]
    setkey(rankedps, "V1")
    cor.frame$phen = factor(cor.frame$phen, levels = rankedps$phen)

    
    aplot = ggplot(cor.frame)
    aplot = aplot + geom_histogram(aes(x=LogP.Cor, y=..density..))
    aplot = aplot + facet_wrap(~phen )
    print(aplot)
    dev.off()

    pdf(fp(output, "behaviorExpCors.spearman.pdf"))
    rankedps = cor.frame[,max(LogP.SpearCor),by="phen"]
    setkey(rankedps, "V1")
    cor.frame$phen = factor(cor.frame$phen, levels = rankedps$phen)
    
    aplot = ggplot(cor.frame)
    aplot = aplot + geom_histogram(aes(x=LogP.SpearCor, y=..density..))
    aplot = aplot + facet_wrap(~phen )
    print(aplot)
    dev.off()
}

plot.poe.data <- function(alldata, ylabel, atitle, mode, fname, legend.position = NULL)
{
    alldata$Strain = as.character(alldata$Strain)
    alldata$Strain = gsub(alldata$Strain, pattern = "\\.", replacement = "x")
    alldata$Strain = factor(alldata$Strain)
    ## alldata$Strain[alldata$Strain =="B6.NOD"] = "B6xNOD"
    ## alldata$Strain[alldata$Strain =="B6.NOD"] = "B6xNOD"
    alldata = alldata[!is.na(alldata$y),]
    factor.mult = 4
    jitter.max  = .2
    alldata$Strain.x = (as.integer(alldata$Strain)-(length(unique(alldata$Strain))+1)/2)
    factorwid = diff(range(alldata$Strain.x))
    
    if(mode==1)
    {
        meanz = alldata[j = list(themean = mean(y, na.rm=T)), by = c("Strain", "Strain.x") ]
        aplot = ggplot(alldata, aes(x = Strain, y= y, shape = Strain))
        aplot = aplot + geom_jitter(aes(colour = Diet),width = jitter.max, size = 2)
        aplot = aplot + geom_segment(
                            aes(x    =  as.integer(Strain) - factorwid/2,
                                xend =  as.integer(Strain) + factorwid/2,
                                y=themean,
                                yend = themean),
                            ##                                    linetype = Strain),
                            size = 1,
                            data=meanz)

        vals = c(1,15); names(vals) =c("B6xNOD","NODxB6")
        ##aplot = aplot + scale_x_continuous(breaks = 0, aes(labels = Strain))

        aplot = aplot + scale_shape_manual( values = vals, guide = guide_legend(title=NULL))
    }
    if(mode==2) ##split by diet
    {
        meanz = alldata[j = list(themean = mean(y,na.rm=T)), by = c("Diet", "Strain", "Strain.x")]
        aplot = ggplot(alldata, aes(x = Strain.x, y= y, shape = Strain))
        
        
        
        aplot = aplot + theme(axis.text.x = element_blank())
        
        aplot = aplot + facet_wrap(~Diet, strip.position = "bottom", nrow=1)
        aplot = aplot + theme(strip.placement = "outside",
                              strip.background = element_blank(),
                              strip.switch.pad.wrap = unit(-.4, "lines"))
        aplot = aplot + theme(panel.spacing.x = unit(0, "lines"))
        ##    aplot = aplot + guides(linetype = F)
        aplot = aplot + coord_cartesian(xlim = c(factor.mult*min(alldata$Strain.x),
                                                 factor.mult*max(alldata$Strain.x)))
        aplot = aplot + geom_jitter(width = jitter.max, size = 3, fill= "grey60")
        aplot = aplot + geom_segment(
                            aes(x    =  Strain.x - factorwid/2,
                                xend =  Strain.x + factorwid/2,
                                y=themean,
                                yend = themean),
                            ##                                    linetype = Strain),
                            size = 1.75,
                            data=meanz)

        vals = c(1,22); names(vals) =c("B6xNOD","NODxB6")
        aplot = aplot + theme(legend.position = c(1,1), legend.justification = c(1,1))
        aplot = aplot + scale_shape_manual( values = vals, guide = guide_legend(title=NULL))
    }

##    aplot = aplot + theme_bw()
    aplot = aplot + ylab(ylabel)
    aplot = aplot + theme(panel.grid.minor = element_blank())
    aplot = aplot + theme(panel.grid.major = element_blank())
    aplot = aplot + theme(axis.title.x = element_blank())
    aplot = aplot + theme(axis.ticks.x = element_blank())
    if(mode==2)
    {
        aplot = aplot + theme(axis.text.x = element_blank())
    }
    ##labels = c("B6xNOD", "NODxB6"),
    aplot = aplot + ggtitle(atitle)
    aplot = aplot + theme(panel.background = element_blank())
    aplot = aplot + theme(strip.background  = element_blank())
    ##aplot = aplot + element_rect(colour = "black", fill=NA, size=5)
    aplot = aplot + theme(axis.line.x = element_line(color="black"))
    aplot = aplot + theme(axis.line.y = element_line(color="black"))
    ##anotherplot = aplot + theme(panel.border     = element_blank())
    
    ## if(mode==2)
    ## {
    ##     browser()
    ## }
                                        #aplot = aplot + coord_fixed(ratio = .1)

    ##    aplot = aplot + guides(shape=guide_legend(title= NULL))
    if(is.null(legend.position))
    {
    aplot = aplot + guides(shape =F)
    } else {
        aplot = aplot + theme(legend.position = legend.position,
                              legend.background=element_rect(fill="white")

                              #,legend.key = element_rect(colour = "transparent", fill = "white")
                              )#, legend.text = element_text(size = 1))
    }

    aplot = aplot + theme(legend.text = element_text(size=12),
                          axis.text = element_text(size=12),
                          axis.title = element_text(size=12),
                          strip.text = element_text(size=12))
    
    return(aplot)
}

show.and.write <- function(aplot, atitle, mode, width = 4, height = 3, fname)
{
    try(dev.off())
##    x11(width = width, height = height)
##    print(aplot)
    afle = outm(paste0(fname, "_", mode,".pdf"))
    print(afle)
    pdf(afle, width = width, heigh=height)
    print(aplot)
    dev.off()
}
