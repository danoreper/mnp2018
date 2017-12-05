
##data = "../data"

to.pval <- function(logp)
{
    return(10^(-logp))
}

getSnordAccessions <- function()
{
    snordfile        = file.path(data,   "sva", "SNORD_cluster_20120524.csv")
    snord.data       <- read.csv(snordfile, stringsAsFactors=FALSE)
    ##snord.data$Location <- sub("\xca", "", snord.data$Location) ##TODO delete?
    snord.data$start <- as.integer(gsub(",", "", (sub("-.*", "", snord.data$Location))))
    snord.data$end   <- as.integer(gsub(",", "", sub(".*-", "", snord.data$Location)))
    accessions = snord.data$mRNA.Accession
    return(accessions)
}

getImprintedAccessions <- function()
{
    accessions = as.vector(read.table("../output/imprinted_transcripts.txt", header=F, sep="\t")[,1])
}


hitsWrapper <- function(fulldata, trueLabels, permLabels)
{
    true.stat = sum(fulldata$hit[trueLabels])
    perm.stat = rep(NA, ncol(permLabels))
    for(i in 1:ncol(permLabels))
    {
        perm.stat[i] = sum(fulldata$hit[permLabels[,i]])
    }
    perm.p = 1 - sum(perm.stat<true.stat)/length(perm.stat)
    return(data.frame(perm.p = perm.p))
    ##hits.perm.p = (1-sum(numHits   >numhits.perm)/length(numhits.perm)),    
}

gseaWrapper <- function(fulldata, trueLabels, permLabels)
{
    computeStat <- function(labels)
    {
        pval.vector.group = -log10(fulldata$thepval)
        pval.vector.group[!labels] =  pval.vector.group[!labels]*-1
        pval.vector.group = pval.vector.group[!is.infinite(pval.vector.group)]
        stat = max(cumsum(pval.vector.group))
        return(stat)
    }

    true.stat  = computeStat(trueLabels)
    perm.stat = rep(NA, ncol(permLabels))
    for(i in 1:ncol(permLabels))
    {
        perm.stat[i] = computeStat(permLabels[,i])
        
    }
    
    perm.p = 1-sum(perm.stat<true.stat)/length(perm.stat)
    return(data.frame(perm.p = perm.p))
}

generateStructuredLabels <- function(n, trueLabels, fulldata)
{
    if(nrow(fulldata)!=length(trueLabels))
    {
        stop("done")
        
    }
    fulldata = data.table(fulldata)
    fulldata$inGroup = trueLabels
    setkey(fulldata, "chrom", "pos.start") ##sorted by chrom, then by position
    fulldata$globalIndex = 1:nrow(fulldata) ##the index within the full data of each probe after sorting by chrom and pos
    fulldata[,pos.order:=rank(pos.start),by=chrom]


    ##probesPerChrom = fulldata[.N, by=chrom]
    
    trueSubset    = droplevels(fulldata[trueLabels])
    clusterBy = "pos.order"
    
    for(chrom in unique(trueSubset$chrom))
    {
        trueSubset.chrom = trueSubset[trueSubset$chrom==chrom]
        for(nclust in 1:2)
        {
            
            clusts = kmeans(trueSubset.chrom[[clusterBy]], nclust)
        }
    }

    ## trueSubRanges = trueSubset[,list(start.globalindex = globalIndex[which.min(pos.order)],
    ##                                  end.globalindex   = globalIndex[which.max(pos.order)],
    ##                                  start.cluster    = min(pos.order),
    ##                                  end.cluster      = max(pos.order)),by=chrom]

    ## for(i in 1:nrow(trueSubRanges))
    ## {
    ##     start.globalindex = trueSubRanges[i, "start.globalindex"]
    ##     end.globalindex   = trueSubRanges[i, "end.globalindex"]
    ##     sum(fulldata$inGroup[start.globalindex, end.globalindex])
    ## }
}
 
##dsubOut = getDataSubset(mRNA.accession  = setdiff(annot.data$mRNA.Accession, dsubIn$mRNA.subset),  annot.data=annot.data, fulldata = fulldata)
runEnrichmentAnalysis <- function(probesetGroup, annot.data, pvalType, varType, fulldataAll, probedGenes, trials = 200,nohits=F)
{
    print(varType)
    fulldata         = fulldataAll[fulldataAll$variable==varType,]
    fulldata$thepval = to.pval(fulldata[[pvalType]])
    fulldata$hit     = fulldata$thepval<.05
    if(nohits)
    {
        fulldata = fulldata[!fulldata$hit,]
    }
    print(dim(fulldata))
    fulldata$inGroup = fulldata$Probe.Set.ID %in% probesetGroup

    #TODO bring back
    ##inGroupPerm = generateStructuredLabels(n = trial, trueLabels = fulldata$inGroup, fulldata)
    inGroupPerm = util$generateShuffles(trueLabels = fulldata$inGroup, trials)
    
    metaFuncs = c(fisherCombinedWrapper, hitsWrapper, gseaWrapper)
    prefixes  = c("fisherCombined","hitscalc","gsea")
##    
    df = data.frame(groupSize = sum(fulldata$inGroup), numHits = sum(fulldata$hit[fulldata$inGroup]), varType = varType, pvalType=pvalType)

    for(j in 1:length(metaFuncs))
    {
        metaFunc = metaFuncs[[j]]
        pvalInfo = metaFunc(fulldata, trueLabels = fulldata$inGroup, permLabels = inGroupPerm)
        colnames(pvalInfo) = paste0(prefixes[j], ".", colnames(pvalInfo))
        df = cbind(df, pvalInfo)
    }
    return(df)
}

meta.inv.variance <- function(theta, theta.se)
{
    w =(1/theta.se)^2
    theta.meta = sum(w*theta)/sum(w)
    theta.meta.se = sqrt(1/sum(w))
    return(data.frame(theta.meta=theta.meta, theta.meta.se = theta.meta.se))
}


##setwd("../")
##load("../output/mnp/output_default/withvar.RData")
##stackedPQdata = stackedPQ(annot.data, results)

runEnrichmentAnalyses <- function(stackedPQdata, annot.data, probedGenes)
{
    runInverseVariance = function(fulldataAll, probesetGroup)
    {
        subd = fulldataAll[fulldataAll$Probe.Set.ID %in% probesetGroup & fulldataAll$variable=="strain",]
        subd = data.table(subd)[,j=list(theta = Est.CrossInStdCtrl, theta.se=SE.CrossInStdCtrl)] 
        df = (cbind(groupName = groupName, meta.inv.variance(subd$theta, subd$theta.se)))
    }
    
    fulldataAll   = stackedPQdata
    probesets = list()


    ##accessions$snord     = getSnordAccessions()
    imprinted = data.table(read.table(prop$genome$imprintedGenesFile, header=T, sep="\t"), key="ensembl_gene_id")
    
    
    

    setkey(probedGenes, "gene_id")
    probedAccessions = probedGenes[imprinted]

    setkey(imprinted, "ensembl_gene_id")
    missingAccessions = imprinted[J(probedAccessions[is.na(globalIndex)]$gene_id)]
    missingAccessions = buildGenomeData$buildGr(start =missingAccessions$start_position,
                                                end = missingAccessions$end_position,
                                                strand = missingAccessions$strand,
                                                chr = missingAccessions$chromosome_name)

    genez = buildGenomeData$buildGr(start = probedGenes$seq_region_start,
                            end   = probedGenes$seq_region_end,
                            strand= probedGenes$seq_region_stran,
                            chr   = probedGenes$seqname)

    overs = findOverlaps(missingAccessions, genez)
    missingAccessions = unique(probedGenes[subjectHits(overs)]$meta_probesetId)

    
    probedAccessions  = probedAccessions[!is.na(globalIndex)]

    imprinted = unique(probedAccessions$meta_probesetId)
    imprinted = c(imprinted, missingAccessions)
    imprinted = intersect(as.character(imprinted), as.character(annot.data$Probe.Set.ID))
    probesets$imprinted = imprinted
    
    
    trials=500
    dfs=list()
    dfsmeta = list()
    outfiles = c(file.path(prop$mnp$output, "withhits.csv"), file.path(prop$mnp$output, "nohits.csv"))
    nohitsVec = c(F,T)
    
    for(h in 1:length(outfiles))
    {
        nohits=nohitsVec[h]
        i = 1
        for(groupName in names(probesets))
        {
            print(groupName)
            probesetGroup = probesets[[groupName]]
##            print(runInverseVariance(fulldataAll, probesetGroup))
            
            for(varType in unique(fulldataAll$variable))
            {
                for(pvalType in c("LogP.Anova","Log.qval"))
                {
                    df = runEnrichmentAnalysis(probesetGroup=probesetGroup, annot.data, pvalType=pvalType, varType=varType, fulldataAll=fulldataAll,nohits=nohits, trials=trials)
                    df = cbind(group=groupName, df)
                    dfs[[i]] = df
                    i = i + 1
                }
            }
        }
        df = do.call(rbind, dfs)
        write.table(df,file=outfiles[h], sep="\t")

        print(df)
    }
}
