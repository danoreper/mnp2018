micro.report = new.env(hash=T)

library(ggplot2)
library(officer)
library(flextable)
library(data.table)
source("./mnp/plotting.R")
source("./mnp/micro/analysis.R")

micro.report$reportAnalysis <- function(exp.mat,
                                        cov.data,
                                        originalResults,
                                        annot.data,
                                        probesetInfo,
                                        threshholds,
                                        karyo)
{
    reportDir = outm("micro")
    dir.create(reportDir, showWarnings = F)

    results = micro.report$postProcessResults(originalResults,
                                              annot.data,
                                              cov.data,
                                              probesetInfo)

    ##writing tables to file etc.
    toWrite = results$full

    
    fwrite(file=outm(fp("micro", paste0("all_results.csv"))), toWrite, sep="\t")
        
    siglevel = .05
##    toWrite = toWrite[is.na(anova.p.value)| anova.p.value<=siglevel,]

    plotting$buildManyScans(results = results, outdir=prop$mnp$output, thresh = threshholds, karyo = karyo)

    print("done plotting many scans")
    mainReport = fp(reportDir, "long_summary.txt")
#    unlink(mainReport)
#    mainReport = file(mainReport, "w")

    toReport = function(str)
    {
        print(str)
    }

    reportTable = function(df)
    {
        impTable = (df[,list(num = length(toUnique(gene_name))),by="imprinted"])
        if(sum(impTable$imprinted=="Y")==0)
        {
            impTable = rbind(impTable, data.table(data.frame(imprinted = "Y", num = 0)))
        }
        
        df.summary = data.frame(
            num.probeset.hits      = length(toUnique(df$Probe.Set.ID)),
            num.gene.hits          = length(unique(unlist(strsplit(df$gene_name, split= ",")))),
            imprinted.gene.hits    = length(unique(unlist(strsplit(df[imprinted=="Y"]$gene_name, split= ",")))))

  ##          not.imprinted.hits     = impTable[imprinted=="N"]$num)

        return(df.summary)
    }    
    
    ps=paste0

    toUnique = function(gene_name)
    {
        genez = unique(unlist(c(lapply(strsplit(gene_name,","),"[", 1),
                       setdiff(lapply(strsplit(gene_name,","),"[", 2), NA))))
        return(genez)
    }

    sink(mainReport)
    toReport("####")
    toReport(ps("Total num probesets: ", nrow(results$per.probe)))
    toReport(ps("Total num unique genes: ",              length(toUnique(results$per.probe$gene_name))))

    toReport(ps("Num imprinted probesets: ",    nrow(results$per.probe[!is.na(minDistToImprinted)&minDistToImprinted<=100])))
    toReport(ps("Num unique imprinted genes: ", length(toUnique(results$per.probe[!is.na(minDistToImprinted)&minDistToImprinted<=100]$gene_name))))
    

    df.summaries = list()
    
    for(analpha in c(.05))##unique(threshholds$permStatistics$alpha))
    {
        
        for(avar in c( "Strain", "Diet", "Diet:Strain"))#unique(threshholds$variable))
        {
            relevantThreshVal = NULL
            print(paste0("working on ", avar))
            toReport("                                          ")
            toReport("                                          ")
            toReport("                                          ")
            toReport(ps("#### variable type: ", avar))
            
            toWrite.sub=toWrite[variable==avar]
            df = micro.report$.formatTable(toWrite.sub, results$per.variable, results$per.level, avar)
            dir.create(fp(reportDir, "effect.table"), recursive=T, showWarnings=F)

            fpmap = list()
            fpmap[["Strain"]]      = "File_S27_POE_expression.csv"
            fpmap[["Diet"]]        = "File_S28_diet_expression.csv"
            fpmap[["Diet:Strain"]] = "File_S29_POE_diet_expression.csv"
            
            
            pfile = fp(reportDir, "effect.table", paste0("p_all_", avar, ".csv"))
            write.table(file=pfile, df, row.names=FALSE, sep="\t")
            
            pfile = fp(reportDir, "effect.table", fpmap[[avar]])
            write.table(file=pfile, df, row.names=FALSE, sep="\t")


            if(!is.null(threshholds))
            {
                relevantThresh = threshholds[variable==avar]
                relevantThreshVal  = util$lookupByFloat(df=relevantThresh, floatkeyCol = "alpha", floatkey = analpha, valueCol = "threshhold.gev") 
                
                if(avar=="Strain")
                {
                    toReport("##Full stats, POE imprinting enrichment by probeset")
                    
                    df$logp = df[["-log10.pval"]]
                    df$passes.FWER = as.numeric(df$logp)>-log10(relevantThreshVal)
                    df1 = df[,list(poe = any(passes.FWER), imprinted = any(imprinted=="Y")),by="gene_name"]
                    
                    df1 = df1[,list(gene_name=unlist(strsplit(as.character(gene_name), split=",")), poe = poe, imprinted = imprinted),by="gene_name"] 
                    
                    X = table(df1$imprinted==T, df1$poe)
                    dimnames(X) = list(Imprinted=c(F,T), POE=c(F,T))
                    toReport(X)
                    fish = fisher.test(X)
                    toReport("Fisher test of enrichment:")
                    toReport(fish)
                }
                
                toReport("##Permutation stats")
                toWrite.perm = toWrite.sub[anova.p.value<relevantThreshVal & anova.q.value<analpha]
                
                df = micro.report$.formatTable(toWrite.perm, results$per.variable, results$per.level, avar)
                                
                df.summary                = reportTable(df)
                
                df.prepend = data.frame(effect.type = avar,
                                        threshold.type  = "FWER", ##"Permutation",
                                        neg.log10.threshold.value = sprintf("%.2f", -log10(relevantThreshVal)))
                df.summary               = cbind(df.prepend, df.summary)

                
                df.summaries = util$appendToList(df.summaries, df.summary)
                
                if(avar=="Diet")
                {
                    toReport("methyl rank:")
                    toReport(table(df$methyl.rank))
                    
                    toReport("signifant hit types")
                    ##                    toReport(table(ps(df$significant.methyl.contrasts,
                    ##                                    "/",
                    ##                                      df$significant.non.methyl.contrasts)))
                    ##toReport(ps("Sig Methyl contrasts:", sum(df$significant.methyl.contrasts)))
                    ##toReport(ps("Non-Methyl contrasts:", sum(df$significant.non.methyl.contrasts)))
                }
                if(avar=="Strain")
                {
                    toReport("Larger expression cross")
                    
                    ##toReport(df[j=list(num=length(toUnique(gene_name))), by="larger.expression"])
                    larger1 = length(unique(unlist(strsplit(df[larger.expression == "NODxB6"]$gene_name, split= ","))))
                    larger2 = length(unique(unlist(strsplit(df[larger.expression == "B6xNOD"]$gene_name, split= ","))))
                    total   = length(unique(unlist(strsplit(df$gene_name, split= ","))))
                    largeness = data.frame(c("NODxB6", "B6xNOD"), c(larger1, larger2))
                    toReport(largeness)
                }
                
                cz  = c(setdiff(colnames(df), c("-log10.pval", "-log10.qval")), c("-log10.pval", "-log10.qval"))
                df = df[,cz, with=F]

                toReport("#p-value table")
                toReport(df)


                pfile = fp(reportDir, "effect.table", paste0("p_", analpha, "_", avar, "_fwer", ".csv"))
                write.table(file=pfile, df, row.names=FALSE, sep="\t")

                
                setnames(df, old = c("-log10.qval", "-log10.pval"), new = c("log10.qval", "log10.pval"))
                theorder = c("gene_name", "chrom", "probesetStart", "Probe.Set.ID", "imprinted", "log10.pval", "log10.qval")

                if(avar == "Diet")
                {
                    theorder = util$insertAtIndex(theorder, index = 5, elem = "methyl.rank")
                    df$methyl.rank = as.character(df$methyl.rank)
                }
                
                if(avar == "Strain")
                {
                    theorder = util$insertAtIndex(theorder, index = 5, elem = "larger.expression")
                }

                setcolorder(df, theorder)
                
                ## browser()
                df$chrom = factor(df$chrom, levels = c(as.character(1:19), "X", "Y"))
                setorder(df, "chrom", "probesetStart")
                df$chrom = as.character(df$chrom)

                
                mytab = regulartable( data = df)
                mytab = bold(mytab, part = "header")
                mytab = align( mytab, align = "center", part = "all")
                rep = list()
                rep[["x"]]             = mytab
                rep[["gene_name"]]     = "Gene"
                rep[["chrom"]]         = "Chr"
                rep[["Probe.Set.ID"]]  = "Probeset ID"
                rep[["probesetStart"]] = "Probeset Location"
                rep[["imprinted"]]     = "Imprinted"
                rep[["log10.qval"]]    = "q value"
                rep[["log10.pval"]]    = "p value"
                if(avar == "Diet")
                {
                    rep[["methyl.rank"]] = "ME Group Rank"                
                }
                if(avar == "Strain")
                {
                    rep[["larger.expression"]]="F1 with Higher Expression"
                }
                mytab = do.call(set_header_labels, rep)

                
                rep[["x"]]             = mytab
                rep[["log10.qval"]]    = "-log10"
                rep[["log10.pval"]]    = "-log10"
                mytab = do.call(add_header, rep)

                mytab = merge_h(mytab, part="header")
                mytab = merge_v(mytab, part="header")
                mytab <- italic(mytab, j = ~ gene_name, italic = TRUE)
                mytab = autofit(mytab, 0, 0)

                if(avar=="Strain")
                {
                    newlen = dim(mytab)$widths["larger.expression"]
                    mytab = flextable::width(mytab, j = ~ larger.expression, width = newlen/2)
                }
                if(avar=="Diet")
                {
                    newlen = dim(mytab)$widths["methyl.rank"]
                    mytab = flextable::width(mytab, j = ~ methyl.rank, width = newlen/2)
                }
                
                newlen = dim(mytab)$widths["probesetStart"]
                mytab = flextable::width(mytab, j = ~ probesetStart, width = newlen*.6)

                
                doc = read_docx(fp("./mnp/template.docx"))
                doc = body_add_flextable(doc, mytab)
                pfile = gsub(pfile, pattern = "csv", replacement = "docx")
                print(doc, target = pfile)
            }
            
            toReport("##FDR stats")
            toWrite.fdr   = toWrite.sub[anova.q.value <analpha]
            df = micro.report$.formatTable(toWrite.fdr, results$per.variable, results$per.level, avar, thresh = relevantThreshVal)
            pfile = fp(reportDir, "effect.table", paste0("p_", analpha, "_",  avar, "_fdr", ".csv"))

            df.summary = reportTable(df)
            df.prepend = data.frame(effect.type = avar,
                                    threshold.type = "FDR",
                                    neg.log10.threshold.value = sprintf("%.2f", min(-log10(toWrite.fdr$anova.p.value))))
            df.summary = cbind(df.prepend, df.summary)
            
            df.summaries = util$appendToList(df.summaries, df.summary)

            ## setnames(df,
            ##          old= c("perm.anova.p.value", "perm.anova.q.value"),
            ##          new = c("'-log10(p-value)'", "'-log10(q-value)'"))

            df1 = df
##            df1$Probe.Set.ID = NULL
            cz  = c(setdiff(colnames(df1), c("-log10.pval", "-log10.qval")), c("-log10.pval", "-log10.qval"))
            df1 = df1[,cz, with=F]
            if(!is.null(threshholds))
            {
                tmp = df1$passes.FWER
                df1$passes.FWER = NULL
                df1$passes.FWER = tmp
            }
            write.table(file=pfile, df1, row.names=FALSE, sep="\t")
            toReport("p-value table")
            toReport(df1)

            
            ## if(avar=="Diet:Strain")
            ## {
             
            ##     write.table(file = pfile, df1, row.names = F, sep="\t")
            ## }
        }
    }
    df.summaries = do.call(rbind, df.summaries)
    df.summaries = data.table(df.summaries, key = c("effect.type", "threshold.type" ))
    write.table(df.summaries, fp(reportDir, "short_summary.txt"), row.names = F, sep=",")
    sink()


##    write.table(threshholds$permStatistics, sep ="\t", row.names = F, file = outm( "perms.csv")
    
##    plotting$pcaExpression(exp.mat = exp.mat, cov.data = cov.data, prop$mnp$output)



##    if(T){save(file = outm("all.RData"), list=ls())}
}


micro.report$.formatTable <- function(sigpq, per.variable, per.level, variable, thresh = NULL)
{

     ## if(avar == "Strain")
     ##        {
     ##            setnames(df,
     ##                     c("gene_name","chrom", "probesetStart", "larger.expression", "imprinted"),
     ##                     c("Gene", "Chr", "Probeset Start", "F1 with Higher Expression", "Imprinted"))
     ##        }

    print(paste0("formatting ", variable))
    
    limitedCols = c(
        ##"variable",
                "gene_name",
        ##"newname",
                    "Probe.Set.ID",
                    "chrom",
                    "probesetStart",
                    "minDistToImprinted",
                 ##   "imprintedMGI",
                 ##   "Crowley_Expressed.allele",
                    "anova.p.value",
                    "anova.q.value"
                    )


    if(variable == "Strain")
    {
        justvar = per.level[variable.level=="StrainNOD.B6",
                            j=list(larger.expression =
                                       ifelse(coef.Value>0, "NODxB6", "B6xNOD"),
                                   Probe.Set.ID)]
        setkey(justvar, "Probe.Set.ID")
        sigpq = justvar[sigpq]
        limitedCols = c(limitedCols, "larger.expression")
        ##add to columns
    } else if (variable == "Diet")
    {
        print("evaluating diet")

    
        v1 = per.level[variable.level=="DietME", j = list(level1 = coef.Value, Probe.Set.ID)]
        v2 = per.level[variable.level=="DietVDD",    j = list(level1 = coef.Value, Probe.Set.ID)]
        v3 = per.level[variable.level=="DietPD",     j = list(level1 = coef.Value, Probe.Set.ID)]

        justvar = data.table(Probe.Set.ID = v1$Probe.Set.ID,
                             methyl.rank = 4 - ((v1$level1>0)*1 + ((v1$level1 - v2$level1)>0)*1 + ((v1$level1 - v3$level1)>0)*1))
                
        setkey(justvar, "Probe.Set.ID")

        methsuff.colz  = c("ME.vs.Std.p.value",  "PD.vs.ME.p.value",  "VDD.vs.ME.p.value")
        otherdiet.colz = c("VDD.vs.Std.p.value", "VDD.vs.PD.p.value", "PD.vs.Std.p.value")

##        table(rowSums((subd[,methsuff.colz, with = F]<=alphalevel)))
##        table(rowSums((subd[,otherdiet.colz, with = F]<=alphalevel)))

        alphalevel = .05
        ## sigpq$significant.methyl.contrasts     = rowSums((sigpq[,methsuff.colz, with = F]<=alphalevel))
        ## sigpq$significant.non.methyl.contrasts = rowSums((sigpq[,otherdiet.colz, with = F]<=alphalevel))
        
        sigpq = justvar[sigpq]
        limitedCols = c(limitedCols, c("methyl.rank"))#, "significant.methyl.contrasts", "significant.non.methyl.contrasts"))
    }



    limitedTable = sigpq[,limitedCols, with=F]

    imprintRange = 100
    limitedTable$imprinted = ifelse(limitedTable$minDistToImprinted<imprintRange, "Y", "N")
    limitedTable$imprinted[is.na(limitedTable$imprinted)] = "N"
    limitedTable$minDistToImprinted = NULL
    ##limitedTable$Crowley_Expressed.allele[limitedTable$imprinted=="N"] = NA
    
    setkey(limitedTable, "anova.p.value")

    if(!is.null(thresh))
    {
        limitedTable$passes.FWER = "N"
        limitedTable$passes.FWER[limitedTable$anova.p.value<thresh] = "Y"
    }
    
    setorder(limitedTable, "chrom","probesetStart")
    limitedTable$anova.q.value = sprintf( "%.1f",-log10(limitedTable$anova.q.value))
    limitedTable$anova.p.value = sprintf("%.1f",-log10(limitedTable$anova.p.value))
    
    ## limitedTable$anova.q.value = -log10(limitedTable$anova.q.value)
    ## limitedTable$anova.q.value = as.character(signif(limitedTable$anova.q.value, digits = 3))

    ## limitedTable$anova.p.value = -log10(limitedTable$anova.p.value)
    ## limitedTable$anova.p.value = as.character(signif(limitedTable$anova.p.value, digits = 3))
    
    setnames(limitedTable,
             old = c("anova.p.value","anova.q.value"),#, "imprintedMGI"),
             new = c("-log10.pval", "-log10.qval"))#, "nearestImprinted"))


    print(paste0("done formatting ", variable))

   
    return(limitedTable)
}

micro.report$postProcessResults <- function(results,
                                            annot.data,
                                            cov.data,
                                            probesetInfo)
{

    print("starting post process!!!!")
    
    ## avg level of expression of the B6.NOD vs Nod.B6 per probesetid     
    expressSummary = micro.report$.getExpressSummary(annot.data, cov.data)
    ## avg level of various diets per probeset id
    dietSummary = micro.report$.getDietSummary(annot.data, cov.data)


    ##merge the per probe data
    setkey(results$per.probe, "Probe.Set.ID")

    results$per.probe = expressSummary[results$per.probe]
    results$per.probe = dietSummary[results$per.probe]
    results$per.probe = probesetInfo[results$per.probe]

    results$per.level$lambda = results$per.probe$lambda[match(results$per.level$Probe.Set.ID, results$per.probe$Probe.Set.ID)]


    results$per.variable[,anova.q.value := p.adjust(anova.p.value , method = "fdr"), by=variable]
    
    setkey(results$per.probe,    "Probe.Set.ID")
    setkey(results$per.variable, "Probe.Set.ID")
    setkey(results$per.level,    "Probe.Set.ID")
    results$full = results$per.variable[results$per.probe]
    setkey(results$full, "Probe.Set.ID", "variable")

    setkey(results$per.variable, "Probe.Set.ID","variable")
    setkey(results$per.probe,    "Probe.Set.ID")
    return(results)
}


micro.report$.getDietSummary <- function(annot.data, cov.data)
{
    StdCtrl.id    = cov.data$ID[cov.data$Diet =="Std"]
    VitDDef.id    = cov.data$ID[cov.data$Diet =="VDD"]
    LowPro.id     = cov.data$ID[cov.data$Diet =="PD"]
    MethylSuff.id = cov.data$ID[cov.data$Diet == "ME"]

    expressSummary = data.table(data.frame(Probe.Set.ID = as.character(annot.data$Probe.Set.ID),
                                           StdCtrl.expression = apply(FUN =sum, 1, X = annot.data[,StdCtrl.id])/length(StdCtrl.id),
                                           VitDDef.expression = apply(FUN =sum, 1, X = annot.data[,VitDDef.id])/length(VitDDef.id),
                                           LowPro.expression = apply(FUN =sum, 1, X = annot.data[,LowPro.id])/length(LowPro.id),
                                           MethylSuff.expression = apply(FUN =sum, 1, X = annot.data[,MethylSuff.id])/length(MethylSuff.id)),
                                key = "Probe.Set.ID")
    return(expressSummary)
}

micro.report$.getExpressSummary <- function(annot.data, cov.data)
{
    b6id = cov.data$ID[cov.data$Strain =="B6.NOD"]
    nodid = cov.data$ID[cov.data$Strain =="NOD.B6"]
    
    expressSummary = data.table(data.frame(Probe.Set.ID = as.character(annot.data$Probe.Set.ID),
                                           B6.NOD.expression = apply(FUN =sum, 1, X = annot.data[,b6id])/length(b6id),
                                           NOD.B6.expression = apply(FUN =sum, 1, X = annot.data[,nodid])/length(nodid)), key = "Probe.Set.ID")
    return(expressSummary)
}



micro.report$.writeLimitedCols <- function(sigpq, siglevel)
{
    ##TODO include a column for the effect size... after ive got residualized working
    limitedCols = c("variable",
                    "anova.p.value",
                    "anova.q.value",
                    "gene_id",
                    "gene_name",
                    "Probe.Set.ID",

                    "chrom",
                    "probesetStart",
                    "probesetEnd",
                    "minDistToImprinted",
                    "imprintedMGI",
                    "Crowley_brainImprinted",
                    "Crowley_strainEffect",
                    "Crowley_Expressed.allele",
                    "numProbes",
                    "hasvariant",

                    "B6.NOD.expression",
                    "NOD.B6.expression"
                   ##,

                    ## "GO.Biological.Process.ID",
                    ## "GO.Biological.Process.Term",
                    ## "GO.Cellular.Component.ID",
                    ## "GO.Cellular.Component.Term",    
                    ## "GO.Molecular.Function.ID",
                    ## "GO.Molecular.Function.Term",
                    ## "Pathway.Name"
                    )

    limitedTable = sigpq[,limitedCols, with=F]

    setkey(limitedTable, "anova.p.value")
    
    write.table(file=outm( fp("micro", paste0("sigp_", siglevel, "_allvariables",".csv"))),
                limitedTable,
                row.names=FALSE,
                sep="\t")
    
    for(variablelevel in c("Strain", "Diet", "Diet:Strain"))##unique(limitedTable$variable))
    {
        outfile = outm(fp("micro", paste0("sigp_", siglevel, "_", variablelevel, ".csv")))
        subt = limitedTable[limitedTable$variable==variablelevel,] 
        write.table(file= outfile, subt, row.names = F, sep="\t")
    }
}

