# Generate masked set of probesets and all probe info. saves it to computeProbeDataDir
# 
# Author: doreper
###############################################################################
library(affxparser)
library(data.table)
library(biomaRt)
library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(data.table)


source("./util/loadParams.R")
source("./genomerep/buildGenomeData2.R")
source("./mnp/micro/preprocess/affylibReader.R")

##either loads or recomputes probe info
getProbeInfo <- function(computedProbeDataDir, variantFile, cel.dir, recomputeApt = T, removeBadProbes = T)
{
    print("getting probe info!!")
    probeData = computeAllProbeInformation(computedProbeDataDir, variantFile, cel.dir, recomputeApt = recomputeApt, removeBadProbes = removeBadProbes)
    probeDataBundle = "probeData"

    dir.create(fp(computedProbeDataDir, probeDataBundle))
    ##fast freadable version of the data to reduce startup time; loading these from binary can take a while.
    print(paste("about to write probe info to ", fp(computedProbeDataDir, probeDataBundle)))
    fwrite(file = fp(computedProbeDataDir, probeDataBundle, "variants.csv"), probeData$variants)
    fwrite(file = fp(computedProbeDataDir, probeDataBundle, "probedGenes.csv"), probeData$probedGenes)
    fwrite(file = fp(computedProbeDataDir, probeDataBundle, "probesetInfo.csv"), probeData$probesetInfo)
    print("wrote probe info")
    save(list="probeData", file = fp(computedProbeDataDir, probeDataBundle, "probeData.RData"))
    mic = probeData$mic
    save(list="mic", file = fp(computedProbeDataDir, probeDataBundle, "mic.RData"))

    return(probeData)
}

computeAllProbeInformation <- function(computedProbeDataDir, variantFile, cel.dir,
                                       recomputeApt = T,
                                       removeBadProbes = T)
{
    ##the microarray library data
    pracma::tic()
    mic = affylibReader$readMicroArrayMetaData(pgfFile = datm(prop$mnp$pgfFile),
                                               clfFile = datm(prop$mnp$clfFile),
                                               probesetSeqLocationFile    = datm(prop$mnp$probesetSeqLocationFile),
                                               mpsFile                    = datm(prop$mnp$mpsFile))
    print("read in microarray meta data")
    pracma::toc()
  
    ##the probes and the exons/genes they are probing, with no consideration for variants
    pracma::tic()
    probedGenes = getProbedGenes(mic)
    print("got probe information, unmasked")

    if(prop$mnp$saveIntermediates) { save(list=ls(), file = outm( "probeDataIntermediate.RData"))}
##    load(outm( "probeDataIntermediate.RData"))
    pracma::toc()

    variants          =  fread(variantFile)
    genesWithVariants = unique(variants$gene_name)
    probedGenes$hasvariant = probedGenes$gene_id %in% genesWithVariants
    
    ##remove probes with variants in bad positions
    if(removeBadProbes)
    {
        probedGenes  = removeBadProbes(variants, probedGenes)
    }
    probedGenes = appendImprintingInfo(probedGenes)
    print("appended imprinting info")
    
    ##Determine what the bad probes are, write them to a kill list for apt masking
    if(removeBadProbes)
    {
        writeKillList(mic, probedGenes, computedProbeDataDir)
        print("wrote kill list")
    }

    ##probeset information summarizing the collections of probes; i.e. num probes per probeset
    probesetInfo = getProbesetInfo(probedGenes, mic)
    print("got probeset info")

    setkey(probesetInfo, "meta_probesetId")
    setkey(probedGenes,  "meta_probesetId")
        
    if(recomputeApt)
    {
        print("running apt-probeset-summarize")
        run.apt.probeset.summarize(apt.outdir = computedProbeDataDir, cel.dir = cel.dir)
    }
    
    probeData = list(probedGenes=probedGenes, probesetInfo = probesetInfo, mic=mic, variants=variants)
    
    return(probeData)    
}

## ##Purely for debugging.
## if(prop$mnp$saveIntermediates)
## {
##     save(list=ls(), file = outm("probeDataIntermediate.RData"))
## }
## ## purely for debugging.
## ##load(outm( "probeDataIntermediate.RData"))
## ##append imprinting metadata

##    probesetInfo$meta_probesetId = as.character(probeData$probesetInfo$meta_probesetId)
    

    
appendImprintingInfo <- function(probedGenes)
{
    ##get the set of all imprinted genes, in a form convenient for granges
    getImprinted <- function()
    {
        imprinted  = data.table(read.table(dat(prop$genome$imprintedGenesFile), header=T), key="ensembl_gene_id")
        setnames(imprinted, old = c("start_position", "end_position", "strainEffect"), new = c("start", "stop","Crowley_strainEffect"))
        imprinted$strand[imprinted$strand=="1"]  = "+"
        imprinted$strand[imprinted$strand=="-1"] = "-"
        return(imprinted)
    }


    getImprintedNearest <- function(probedRanges, imprintedRanges, imprinted)
    {
        imprintedNearestGene = getImprintedNearestHelper(probedRanges, imprintedRanges, imprinted)
        imprintedNearestDist = getImprintedNearestDist(probedRanges, imprintedRanges)

        imprintedJoin = imprintedNearestGene
        imprintedJoin = imprintedNearestDist[imprintedJoin]
        return(imprintedJoin)
    }
    
    getImprintedNearestHelper <- function(probedRanges, imprintedRanges, imprinted)
    {
        nearestImprinted = nearest(probedRanges, imprintedRanges)
        imprintedJoin = cbind(probeInd =  as.integer(1:length(probedRanges)),
                              imprinted[nearestImprinted])

        setnames(imprintedJoin,
                 c("start", "stop", "strand", "mgi_symbol", "chromosome_name"),
                 c("imprinted_start", "imprinted_stop", "imprinted_strand", "imprinted_mgi_symbol", "imprinted_chromosome_name"))
        ## imprintedJoin$start           = NULL
        ## imprintedJoin$stop            = NULL
        ## imprintedJoin$strand          = NULL
        ## imprintedJoin$mgi_symbol      = NULL
        ## imprintedJoin$chromosome_name = NULL

        
        imprintedJoin = data.table(imprintedJoin, key = "probeInd")
        return(imprintedJoin)
    }

    getImprintedNearestDist <-  function(probedRanges, imprintedRanges)
    {
        nearestDist      = distanceToNearest(probedRanges, imprintedRanges)
        nearestDistDf    = as.data.table(nearestDist)
        nearestDistDf$subjectHits = NULL
        setnames(nearestDistDf, "distance", "distanceToNearestImprintedGene")
        setnames(nearestDistDf, "queryHits", "probeInd")
        setkey(nearestDistDf, "probeInd")
        return(nearestDistDf)
    }


    imprinted = getImprinted()
    imprintedRanges = makeGRangesFromDataFrame(imprinted)


    probedGenes$chrom = NULL


    probedRanges    = makeGRangesFromDataFrame(probedGenes)

  
    probedGenes$probeInd    = 1:nrow(probedGenes)
    setkey(probedGenes,   "probeInd")

    ##imprintedHits    = getImprintedHits(probedRanges, imprintedRanges, imprinted)
    imprintedJoin                 = getImprintedNearest(probedRanges, imprintedRanges, imprinted)
    imprintedJoin                 = imprintedJoin[probedGenes]
##    imprintedJoin$ensembl_gene_id = NULL
    imprintedJoin$probeInd        = NULL


    for (cname in c("brainImprinted", "Crowley_strainEffect"))
    {
        imprintedJoin[[cname]][is.na(imprintedJoin[[cname]])] = F
    }
    
    return(imprintedJoin)
}    
    
getProbesetInfo <- function(probedGenes, mic)
{
    probedGenes                 = data.table(probedGenes, "meta_probesetId")
    probedGenes$meta_probesetId = as.character(probedGenes$meta_probesetId)
    probesetInfo                = probedGenes[,list(chrom = paste(unique(seqname), collapse=","),
                                                    probesetStart = min(seq_region_start),
                                                    probesetEnd = max(seq_region_end),
                                                    probeType   = paste(unique(fulltype), collapse=","),
                                                    gene_id = paste(unique(gene_id), collapse = ","),
                                                    gene_name = paste(unique(gene_name), collapse=","),
                                                    numProbes = length(unique(probeId)),
                                                    hasvariant = any(hasvariant),
                                                    minDistToImprinted = min(distanceToNearestImprintedGene),
                                                    imprintedMGI       =  imprinted_mgi_symbol[which.min(distanceToNearestImprintedGene)],
                                                    Crowley_brainImprinted  = brainImprinted[which.min(distanceToNearestImprintedGene)],
                                                    Crowley_strainEffect    = Crowley_strainEffect[which.min(distanceToNearestImprintedGene)],
                                                    Crowley_Expressed.allele = paste(unique(Expressed.allele ))),
                                              by = "meta_probesetId"]

##    
    mainProbes  = getMainProbes(mic)
    setkey(mainProbes, "meta_probesetId")
    mainProbesetInfo = mainProbes[,list(numAllProbes = length(unique(probeId))), by ="meta_probesetId"]
    mainProbesetInfo$meta_probesetId = as.character( mainProbesetInfo$meta_probesetId)
    
    setkey(mainProbesetInfo, "meta_probesetId")
    probesetInfo = mainProbesetInfo[probesetInfo]
    probesetInfo$numMaskedProbes = probesetInfo$numAllProbes - probesetInfo$numProbes
    probesetInfo$numAllProbes    = NULL

    ##Add the new names column
##    probesetInfo = updateProbesetInfo(probesetInfo)
    return(probesetInfo)
}


getProbedGenes <- function(mic)
{
    pracma::tic()
    probed_regions = getProbeAlignments(mic)
    pracma::toc()
    print("read in probe alignments")

    probed_regions.mle  = getMlePositions(probed_regions)
    print("reduced to mle positions")
    probedRanges = buildProbeRanges(probed_regions.mle)
    
    exons = buildGenomeData$getExons(dat(prop$mnp$exonReferenceFile))
    overs = findOverlaps(probedRanges,exons)

    probed_regions.mle[,hits.exon:=FALSE,]
    probed_regions.mle$hits.exon[queryHits(overs)]=T
    
    
    probedGenes = data.table(data.frame(
        transcriptname  = exons[subjectHits(overs)]$transcript_id,
        exonname        = exons[subjectHits(overs)]$exon_id,
        gene_id         = exons[subjectHits(overs)]$gene_id,
        probeId         = probed_regions.mle$probeId[queryHits(overs)],
        meta_probesetId = probed_regions.mle$meta_probesetId[queryHits(overs)]),
        key             = c("meta_probesetId","probeId"))
    

    setkey(probed_regions.mle, "meta_probesetId", "probeId")
    probedGenes = probed_regions.mle[probedGenes]


    additionalGeneInfo        = buildGenomeData$getGeneInfo()
    additionalGeneInfo$start  = NULL
    additionalGeneInfo$end    = NULL
    additionalGeneInfo$strand = NULL
    
    ## setnames(additionalGeneInfo,
    ##          old = c("start", "end", "strand"),
    ##          new= c("gene_start", "gene_end", "gene_strand"))

    setkey(additionalGeneInfo, "gene_id")
    setkey(probedGenes, "gene_id")
    probedGenes = additionalGeneInfo[probedGenes]

           
    setkey(probedGenes, "meta_probesetId")
    return(probedGenes)
}


##get the probes in the main set of the microarray that are not aligning to an unassembled contig
getMainProbes <- function(mic) 
{
    ##get the probes that are in the "main" set- the exons, and arent aligning to some unassembled contig
    probesetsToQuery = mic$probeSetTable$probeset_type=="main"
    probesetsToQuery = mic$probeSetTable[probesetsToQuery,]
    probesToQuery    = mic$probeTable[mic$probeTable$probesetId %in% probesetsToQuery$probesetId,]

    ##setidMatches = match(probesToQuery$probesetId, probesetsToQuery$probesetId)
    ##probesToQuery$seqname = probesetsToQuery[setidMatches,"seqname", with=F]
    
    ##this stuff is a little extraneuous but no point in letting it return outside the method
    probesToQuery$probeGcCount               = NULL
    probesToQuery$probeInterrogationPosition = NULL
    probesToQuery$probeLength                = NULL
    probesToQuery$probeType                  = NULL
    probesToQuery$atomId                     = NULL

    probesToQuery = data.table(probesToQuery, key="probeId")
    return(probesToQuery)
}


##reads in the 1.0 probe alignments from file as generated by connectToEnsembl.sh, discarding probes that align to too many regions, or that align to weird regions
##If grabbing alignments from a database seems unappealing, an alternate option would be to script up the ensembl alignment pipeline... seems like more trouble than its worth.
getProbeAlignments_1.0<- function()
{
    probeAligns = fread(datm(prop$mnp$probeAlignmentFile), header = T)
    limit = nrow(probeAligns)
    if(!is.na(prop$mnp$evalprobelimit))
    {
        limit = prop$mnp$evalprobelimit
        probeAligns = probeAligns[1:limit,]
    }

    probeAligns$probename = gsub(probeAligns$probename, pattern = "-.+", replacement = "")

    
    ##sometimes multiple experiments (with different analysis ids) have stored the same alignment in the database. 
    ##discard the redundant ones (easier done in R than in sql)
    setkey(probeAligns, "probename", "seqname", "seq_region_start", "seq_region_end", "seq_region_strand")
    probeAligns = unique(probeAligns)
    
    probeAligns = data.table(probeAligns, key="probename")

    ##seqname is the chromosome-- these are bad alignments as far as we are concerned
    probeAligns = probeAligns[!grepl("GL|JH", probeAligns$seqname),]
    
    ##throw out alignments with more than 1 mismatch: TODO property
    probeAligns = probeAligns[mismatches<2,]

    ##get rid of probes that align to too many positions
    countPerProbe = probeAligns[ ,list(count = length(unique(paste0(seq_region_strand, ".", seq_region_start)))),
                                by="probename"]
    
    setkey(countPerProbe, "probename")
    setkey(probeAligns,  "probename")
    
    goodProbes  = countPerProbe[countPerProbe$count<=prop$mnp$probeAlignThresh,]
    probeAligns = probeAligns[J(goodProbes)]
    
    probeAligns = data.table(probeAligns, key="probename")

    
    return(probeAligns)
}

##Gets alignments of all probes conctanated with other probe information. Discard unalligned probes
getProbeAlignments <- function(mic)
{
    ##get the main probes we want to align
    probesToQuery = getMainProbes(mic = mic)
        
    ##throw away probes aligning some place unknown
    ##probesToQuery = mainProbes[mainProbes$seqname!="chrUn"]
    
    ##get a table that converts the old r3 affy 1.0 probe ids to the r4 affy 1.1 probe ids. The alignments stored in the ensembl db are r3 probe ids
    conversionTable = getConversionTable(probesToQuery)
    print("got conversion table")
    
    ##convert the probenames, which are 1.0 r3 ids, to the 1.1 r4 probe ids
    probeAlignments_1.0 = getProbeAlignments_1.0()
    print("got1.0_alignents")

    probeAligns = conversionTable[probeAlignments_1.0]
    print("converted 1.0 probenames to 1.1 probnames")

    ##discard aligns that arent going to be queried; i.e. those that dont have 1.1 r4 probe id
    probeAligns = probeAligns[!is.na(probeAligns$probeId)]
    setkey(probeAligns, "probeId")

    
    probesToQuery$probeId = as.character(probesToQuery$probeId)
    setkey(probesToQuery, "probeId")
    probesToQuery = probeAligns[probesToQuery, allow.cartesian=T]
    
    ##discard probes that dont have a corresponding 1.0 affy id
    probesToQuery = probesToQuery[!is.na(probesToQuery$affy_1.0_r3_id)]
    ##probesToQuery[, numAlignsPerProbe:=.N, by=probeId]
    probesToQuery$probeId = as.integer(probesToQuery$probeId)

    probesToQuery = getProbedRegions(probesToQuery)
    print("got probes to query")
    return(probesToQuery)
}

#given a set of r4 1.1 probes to query, builds a table of the corresponding probe ids in the r3 affymetrix lib files
#this is necessary because the ensembl database of alignments is in terms of the old r3 probe ids.
#this method is essentially id remapping by matching the probe sequences for each probe id.
getConversionTable <- function(probesToQuery)
{
    #location of the old library files
    oldPgfFile                 = datm(prop$mnp$oldPgfFile)
    oldClfFile                 = datm(prop$mnp$oldClfFile)
    oldProbesetSeqLocationFile = datm(prop$mnp$oldProbesetSeqLocationFile)
    oldMpsFile                 = datm(prop$mnp$oldMpsFile)

    ##read in meta data about the 1.0 array in order to make correspondences from the 1.1 to 1.0 array ids
    micOld = affylibReader$readMicroArrayMetaData(pgfFile = oldPgfFile, clfFile = oldClfFile, probesetSeqLocationFile = oldProbesetSeqLocationFile,oldMpsFile)
    oldProbeTable = data.table(micOld$probeTable,key="probeSequence")
    
    uniqueOld         = oldProbeTable[!duplicated(oldProbeTable$probeSequence),]
    oldProbeId        = uniqueOld[J(probesToQuery$probeSequence),probeId]
    
    conversionTable =data.frame(probeId = probesToQuery$probeId, affy_1.0_r3_id = oldProbeId)
    conversionTable = unique(conversionTable)
    conversionTable$probeId = as.character(conversionTable$probeId)
    conversionTable$affy_1.0_r3_id = as.character(conversionTable$affy_1.0_r3_id)
    conversionTable = data.table(conversionTable, key="affy_1.0_r3_id")
    
    return(conversionTable)
}

#converts probe alignments to probe binding regions
getProbedRegions <- function(withPositions)
{
    ##for reasons that are unclear, the specified probe sequence is the reverse (NOT the reverse complement) of the sequence to which the probe BINDS (not aligns)
    withPositions$probeSequence = as.character(withPositions$probeSequence)
    withPositions$probeSequence = stringi::stri_reverse(withPositions$probeSequence)
    
    ##we are trying to specify the position of where the probe BINDS (not aligns), which is the opposite strand

    pluses = withPositions$seq_region_strand == "1"
    minuses = withPositions$seq_region_strand == "-1"
    withPositions$seq_region_strand[pluses] = "-"
    withPositions$seq_region_strand[minuses] = "+"

    print("got probed regions")
    return(withPositions)
}


##For each probe/meta_probeset combination, identifies what the intended probed region is,
## discard remaining probed genes
getMlePositions <- function(withPositions)
{
    computeDist = function(seqname1, seqname2, strand1, strand2, pos1, pos2)
    {

        out = abs(pos1-pos2)
        outmult = seqname1==seqname2 & strand1==strand2
        outmult[!outmult]=Inf
        out = out*outmult
        return(out)
    }
    
    prevBestValue = -1
    positionsToKeep = function(subdf, meta_probesetId)
    {

        if(!any(duplicated(subdf$probeId)))
        {
            return(subdf$globalIndex)
        }
        
        subdf = copy(subdf)

        prevProbe   = -1
        attempts    = 0
        maxAttempts = 3
        best.i      = NA
        best.value  = Inf
        i = 1

        while(i<nrow(subdf) & attempts<maxAttempts)
        {
            aprobe = subdf$probeId[i]
            seq1   = subdf$seqname[i]
            strand1= subdf$seq_region_strand[i]
            pos1   = subdf$seq_region_start[i]
            subdf[,dists  := (computeDist(seq1,    subdf$seqname,
                                          strand1, subdf$seq_region_strand,
                                          pos1,    subdf$seq_region_start)),]

            mindists = subdf[,list(dists=min(dists)), by=probeId]
            
            totalScore = sum((mindists$dists))
            if(is.na(totalScore))
            {
                
            }
            if(totalScore<best.value)
            {
                best.value = totalScore
                best.i = i 
            }

            ##Some sort of optimization here, TODO figure it out and document it.
            if(aprobe!=prevProbe)
            {
                attempts = attempts + 1
                prevProbe = aprobe
            }
            i = i + 1
            
        }

        if(is.infinite(best.value))
        {
            ##print("bad best value")
            return(c())
        }
        seq1    = subdf$seqname[best.i]
        pos1    = subdf$seq_region_start[best.i]
        strand1 = subdf$seq_region_strand[best.i]
        subdf[,dists  := (computeDist(seq1,    subdf$seqname,
                                      strand1, subdf$seq_region_strand,
                                      pos1,    subdf$seq_region_start)),]

        tokeep = subdf[,list(tokeep = globalIndex[which.min(dists)]), by=probeId]$tokeep
        ##print(paste(meta_probesetId, best.value))

        ##  if(best.value!=0 && best.value==prevBestValue)
        ##  {
        ##        
        ##  }
        prevBestValue <<- best.value
        
        return(tokeep) 
    }
    withPositions$globalIndex = 1:nrow(withPositions)

    setkey(withPositions, "meta_probesetId") #, "numAlignsPerProbe", "probeId")

    keepable = withPositions[,list(tokeep=positionsToKeep(.SD, meta_probesetId)), by=meta_probesetId,
                             .SDcols=c("probeId","seqname","seq_region_start","seq_region_strand","globalIndex") ]
    withPositions.mle = withPositions[withPositions$globalIndex %in% keepable$tokeep]
    
    return(withPositions.mle)
}

buildProbeRanges <- function(aligns)
{
    ranges = buildGenomeData$buildGr(start  = aligns$seq_region_start,
                                     end    = aligns$seq_region_end,
                                     strand = aligns$seq_region_strand,
                                     chr    = aligns$seqname,
                                     names  = aligns$seqname)
    return(ranges)
}
    
    
buildVariantRanges <- function(vars)
{
    varstarts = vars$pos
    varends   = varstarts + abs(max(nchar(vars$allele1) - nchar(vars$allele2),1))- 1
    ##   chr       = vars$chrom
    chr       = vars$chr
    strand    = "*"
    
    varRanges = buildGenomeData$buildGr(varstarts, varends, strand, chr, names = chr)
    return(varRanges)
}

removeBadProbes <- function(variants, probedGenes)
{
    ##shrinks probes by ignore buff
    shrinkProbeAlignmentRanges <- function(probeRanges, ignoreBuffer)
    {
        probeRanges        = copy(probeRanges)
        start(probeRanges) = start(probeRanges) + ignoreBuffer
        end(probeRanges)   = end(probeRanges)   - ignoreBuffer
        return(probeRanges)
    }
    
    variantRanges =  buildVariantRanges(variants)
    probeRanges = buildGenomeData$buildGr(probedGenes$seq_region_start,
                                          probedGenes$seq_region_end,
                                          probedGenes$seq_region_strand,
                                          probedGenes$seqname)

    
    ##TODO property for shrink amount, or pass it in.
    probeRanges.buf.3 = shrinkProbeAlignmentRanges(probeRanges, 3)
    overs = findOverlaps(variantRanges, probeRanges)
    probedGenes$badprobe = F
    
    probedGenes[subjectHits(overs), "badprobe"] = T
    probedGenes = probedGenes[!probedGenes$badprobe,]
    
    return(probedGenes)
}

writeKillList <- function(mic, probedGenes, apt.outdir)
{
    dir.create(apt.outdir, showWarnings=F)
    mainProbes  = getMainProbes(mic)
    ##any main probes not in the probedGenes structure-- i.e. which dont align to reference,
    ## or dont align to exon, or which align to too many places, or which have a variant... etc are
    ## to be masked
    tomask = !(paste0(mainProbes$probeId,  "_", mainProbes$probesetId) %in% paste0(probedGenes$probeId, "_", probedGenes$probesetId))
    tokill = mainProbes[tomask,]

    tokill = data.frame(probe_id = tokill$probeId, probeset_id=tokill$probesetId)
    
    ##TODO, pass in
    write.table(x=tokill, file =getKillFileName(apt.outdir),row.name=F,col.names=T, sep="\t",quote=F)
}

getKillFileName <- function(computedProbeDataDir)
{
    apt.killfile = file.path(computedProbeDataDir, "killList.txt")
    return(apt.killfile)
}
    
#apt.outdir must contain the killlist already. 
run.apt.probeset.summarize <- function(apt.outdir, cel.dir)
{
    ##    aptBinPath = file.path(prop$resources, "apt-1.18.0-x86_64-intel-linux/bin")
    aptBin = fp(prop$external$apt)
    aptLibPath = datm("microarray_lib_files/MoGene-1_1-st-v1.r4.analysis-lib-files")
    apt.killfile = getKillFileName(apt.outdir)


    dir.create(apt.outdir, showWarnings = F)
    apt.logfile = file.path(apt.outdir, "aptlog")

    command = paste(
        aptBin,
      "-o", apt.outdir,
      "-a rma",
      "-c ", file.path(aptLibPath, "MoGene-1_1-st-v1.r4.clf"),
      "-b ", file.path(aptLibPath, "MoGene-1_1-st-v1.r4.bgp"),
      "-m ", file.path(aptLibPath, "MoGene-1_1-st-v1.r4.mps"),
      "-p ", file.path(aptLibPath, "MoGene-1_1-st-v1.r4.pgf"),
      "--qc-probesets ", file.path(aptLibPath, "MoGene-1_1-st-v1.r4.qcc"),
      "--kill-list", apt.killfile,
      paste(cel.dir,"/*.CEL >", sep=""), apt.logfile)
    print(command)
    x = system(command)
}


##Update probeset info with the latest 38.86 names
updateProbesetInfo <- function(probesetInfo)
{
    newExons = buildGenomeData$getExons(dat("b6_reference/Mus_musculus.GRCm38.86.gtf"))
    ex = newExons[!duplicated(newExons$gene_id)]

    probesetInfo$newname = ""

    ps.clean = probesetInfo[!probesetInfo$isSnord]
##TODO generalize for case with more than two genes at a given position
    ps.clean.id1 = unlist(lapply(strsplit(ps.clean$gene_id, ","), "[",1))
    ps.clean.id2 = unlist(lapply(strsplit(ps.clean$gene_id, ","), "[",2))
    matchInd.1 = match(ps.clean.id1, ex$gene_id)
    matchInd.2 = match(ps.clean.id2, ex$gene_id)

    ps.clean$newname=""
    ps.clean$newname[!is.na(matchInd.1)] = ex[matchInd.1[!is.na(matchInd.1)]]$gene_name
    ps.clean$newname[!is.na(matchInd.2)] = paste0(ps.clean$newname[!is.na(matchInd.2)], ",",
                                                  ex[matchInd.2[!is.na(matchInd.2)]]$gene_name)
    ps.clean$newname = unlist(gsub("^,", replacement = "", ps.clean$newname))
    
    ##head(ps.clean$gene_id[is.na(matchInd)])
    probesetInfo2 = rbind(probesetInfo[probesetInfo$isSnord,], ps.clean)
    probesetInfo2$newname[probesetInfo2$isSnord]=probesetInfo2$gene_name[probesetInfo2$isSnord]
    setkey(probesetInfo2, "meta_probesetId")
    return(probesetInfo2)
}

