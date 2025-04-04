##TODO: wrap functions up in the environment
source("./util/loadParams.R")
source("./mnp/behavior/loadData.R")
source("./mnp/micro/preprocess/evalprobes2.R")
source("./mnp/loadBreeding.R")
source("./util/utils.R")
source("./mnp/general.R")

loadAllData = new.env(hash=T)

loadAllData$createAllInputs <- function()
{
    pracma::tic()

    recomputeProbeData = prop$mnp$recomputeProbeInfo
    computedProbeDataDir = datm("matnut.newmask")
    dir.create(computedProbeDataDir, recursive = T, showWarnings = F)
    
    maskedraw.expression.file <- fp(computedProbeDataDir, "/rma.summary.txt")
    print(paste0("using expression file: ", maskedraw.expression.file))

    ##Get the probe information
    ##Get the probe information
    if(!recomputeProbeData)
    {
        probeDataBundle = "probeData"
        probedGenes = suppressWarnings(fread(input = fp(computedProbeDataDir, probeDataBundle, "probedGenes.csv"), showProgress = F))
        probesetInfo = suppressWarnings(fread(input = fp(computedProbeDataDir, probeDataBundle, "probesetInfo.csv"), showProgress = F))
    } else {
        
        ##TODO different property for crowley data. OR create it on the fly? Check whether allele thing matters.
        ##Consider passing in a set of files representing each chromosome instead?
        ##Consider passing in zipped file
        ##TODO create a method which Consider passing in db folder and strain names instead?
        variantFile = datm( prop$mnp$NOD.B6.variantsFile)
        ##TODO: different property for crowley
        cel.dir = datm(prop$mnp$cel.dir)
        probeData         = getProbeInfo(computedProbeDataDir, variantFile, cel.dir = cel.dir)
        probedGenes       = probeData$probedGenes
        probesetInfo      = probeData$probesetInfo

        rm(probeData)
    }
    probesetInfo$meta_probesetId = as.character(probesetInfo$meta_probesetId)

    ##TODO move into evalProbes2
    probesetInfo$twogenes        = grepl(pattern=",", probesetInfo$gene_name)
    snords = buildGenomeData$getSnords()
    probesetInfo = adjustSnordNames(snords, probesetInfo)
##    probesetInfo = updateProbesetInfo(probesetInfo)

    probesetInfo = probesetInfo[chrom %in% c(as.character(1:19), "X", "Y")]
    

    gc()
###### Processing masked but unannotated expression data
    ## Probes were masked if they had a SNP within 3-21 inclusive
    ## Probe sets were excluded if they had less than 4 probes, or expression lt 2^x (where at time of writing, x was 5)

    ##Get the masked annotated data
    annot.data = get.annot.data(annot.old.expression.file = datm(prop$mnp$old.expression.file),
                                maskedraw.expression.file = maskedraw.expression.file,
                                probesetInfo = probesetInfo,
                                output= prop$mnp$output)

    ##the annotated data object contains two fields, which contain distinct data: the annotated data less any controls, and just the negative controls.
    ##Consider adding positive controls
    annot.data.control = annot.data$annot.data.control
    annot.data         = annot.data$annot.data
    ##TODO is this filtering needed?
##    annot.data         = annot.data[annot.data$Probe.Set.ID %in 
    
    
    ##Get meta data from annotations-- merge altogether into probesetInfo
    annot.data.description  = getAnnotDescriptions(annot.data) 
    setkey(probesetInfo, "meta_probesetId")
    setkey(annot.data.description, "Probe.Set.ID")

    probesetInfo = probesetInfo[annot.data.description]
    

    gc()

    ##for testing purposes, we sometimes want to limt the size of the data set
    if(!is.na(prop$mnp$limit))
    {
        annot.data = annot.data[1:prop$mnp$limit,]
        annot.data.control = annot.data.control[1:prop$mnp$limit,]

    }

    if(prop$mnp$tiny)
    {
        annot.data.control = annot.data.control[1:10,]
        cnamez = c("10344614", "10412701", "10408280")##, annot.data$Probe.Set.ID[1:1000] )
        cnamez = c("10412701", "10408280", "10398326")
        cnamez = c("10588903", "10408280")
        annot.data = annot.data[annot.data$Probe.Set.ID %in% cnamez,]#,"10372534", "10344624"),]
    }
    
    
    ##Get just the expression matrixes; as opposed to data frame with lots of metadata
    exp.mat.full                 = getExpressionMatrix(annot.data = annot.data)
    exp.mat.control.full         = getExpressionMatrix(annot.data = annot.data.control)


    
    ##TODO read this in the same way as for behavior files
    ##Get the covariates file
    ##only uses the mouse names from exp.mat, nothing else
    cov.data.allSamples = read.full.cov()

    cov.data.full = toExpCov(cov.data.allSamples, exp.mat.full)


    if(prop$mnp$saveIntermediates) { save(file=outm("exp.mat.RData"), list=ls())}

###Use this to check whether there are effects on control probes
    if(prop$mnp$computeEffectsOnControlProbes)
    {
        exp.mat.full = exp.mat.control.full
        prop$mnp$num.sv <<- 0
    }
    
    ##If we are only interested in one of the two pipelines, filter out the samples for the wrong pipeline
    pipelineGroup    = c("1","2") ##prop$mnp$pipelines
    pipeInd          = cov.data.full$Pipeline %in% pipelineGroup 
    exp.mat          = exp.mat.full[pipeInd,]
##    exp.mat.control  = exp.mat.control.full[pipeInd,]
##    cov.data         = cov.data.full[pipeInd,]

    phenz = loadBehavior$getPhenotypeRepository()
    ## TODO fix in evalProbes.
    setnames(probesetInfo, old = "meta_probesetId", new = "Probe.Set.ID")

    karyo = buildGenomeData$getKaryotype(dat( prop$genome$karyotype))
    
    inp = list(exp.mat            = exp.mat,
  ##             exp.mat.control    = exp.mat.control,
               exp.mat.control.full = exp.mat.control.full, ##for all pipelines, needed for sv analysis
               annot.data         = annot.data,
               annot.data.control = annot.data.control,
##               cov.data           = cov.data,
               cov.data.full      = cov.data.full, #for all pipelines
               cov.data.allSamples = cov.data.allSamples,
               pipeInd            = pipeInd,
               probedGenes        = probedGenes,
               probesetInfo       = probesetInfo,
               phens              = phenz,
               karyo              = karyo)
    pracma::toc()
    return(inp)
}

getTaqmanData <- function(toExamine)
{
    
    followup = read.table(datm("qpcr/Matnut_pilot_alltaqmanplates_datacompletev3_100716.csv"), sep=",", header=T)
    followup$Assay = as.character(followup$Assay)
    followup$Assay[followup$Assay =="Lrrc16a"] = "Carmil1"

    followup$Assay = as.factor(followup$Assay)
    followup = data.table(followup)
    followup = followup[Strain %in% c("NOD.B6", "B6.NOD")]
    
    ## followup = data.table::dcast(followup, ID + Plate + Strain ~ Assay, value.var=c("FAM.Ct", "VIC.Ct", "Delta.Ct"), fun.aggregate=mean)
    
    followup = followup[followup$Assay == toExamine]
    followup = followup[!is.na(FAM.Ct)]
    setnames(followup, old = c("FAM.Ct", "VIC.Ct"), new = c("goi.taq", "control.taq"))

    followup$ID = paste0("Mouse.", followup$ID)
    setkey(followup, "ID")
   
    return(followup)
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


toExpCov <- function(breedLog, exp.mat)
{
    exp.mice <- rownames(exp.mat)
    cov.data <- breedLog[ match(exp.mice, breedLog$ID), ]
    rownames(cov.data) = cov.data$ID
    print(all(cov.data$ID==exp.mice))
    ## if(useCurrentSamples)
    ## {
    ##     validID = getPupID.with.existingSample()
    ##     cov.data = cov.data[mouseInt %in% validID,]
    ## }
    cov.data = droplevels(cov.data)
    return(cov.data)
}


##gets just the expression matrix from the annotated data
getExpressionMatrix <- function(annot.data) 
{
    mouse.cols      <-  grep("Mouse", colnames(annot.data))
    exp.mat         <-  t(as.matrix(annot.data[ , mouse.cols]))
    colnames(exp.mat) = annot.data$Probe.Set.ID
    return(exp.mat)
}


getAnnotDescriptions <- function(annot.data)
{
    annotSub = data.table(copy(annot.data), "Probe.Set.ID")
    descrip.cols <- c(which(colnames(annot.data)=="Probe.Set.ID"), which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data))
    annotSub = annotSub[,descrip.cols, with=F]
    annotSub$Probe.Set.ID = as.character(annotSub$Probe.Set.ID)
    annotSub$Gene.Accession = NULL
    annotSub$Gene.Symbol = NULL
    annotSub$mRNA.Accession = NULL
    setkey(annotSub, "Probe.Set.ID")
    return(annotSub)
}

getPupID.with.existingSample <- function()
{
    fle = datm(prop$mnp$currentSamplesFile)
    filterFrame = fread(fle)
    filterFrame = filterFrame[filterFrame[["Left hemisphere"]] =="Yes" |
                              filterFrame[["Pulverized brain"]]=="Yes"]
    validID = filterFrame$ID
    ##added
    validID = paste0("Mouse.", validID)
    ##
    return(validID)
}

removeInvalidIds <- function(annot.data, validID)
{
    mouseCols = grepl(colnames(annot.data), pattern = "Mouse\\.[0-9]+")
    mouseCols = colnames(annot.data)[mouseCols]
    badCols = setdiff(mouseCols, validID)
    badInt = match(badCols, colnames(annot.data))
    annot.data = annot.data[,-badInt]
    return(annot.data)
}
            
get.annot.data <- function(annot.old.expression.file,
                           maskedraw.expression.file,
                           probesetInfo,
##                           probeTable,
                           output)
{
    annot.data = read.annot.data(annot.old.expression.file, outdir = output)
    
    ##remove columns from annot.data corresponding to samples that we no longer have.
    if(useCurrentSamples)
    {
        validID = getPupID.with.existingSample()
        annot.data = removeInvalidIds(annot.data, validID)
    }
    ##plotControls(annot.data = annot.data, probeTable = probeTable, outdir = output)


    annot.data = replaceAnnotExpressionWithMasked(annot.data = annot.data, maskedraw.expression.file = maskedraw.expression.file)


##save negative control probes for later use- these are the indexes into annot.data
    ##for bad probes
    ii.control = grep("^(neg_control)", perl=TRUE, annot.data$mRNA.Description)
    annot.data.control = annot.data[ ii.control, ]
    
    ##remove positive and negative controls from the analyzed expression values 
    ii <- grep("^(pos_control|neg_control)", perl=TRUE, annot.data$mRNA.Description)
    print(paste0("num control probesets: ", length(ii))) # 6546
    annot.data <- annot.data[ -ii, ]

    ## remove non-annotated
    ii <- which(is.na(annot.data$mRNA.Accession)) #the indexes into annot.data for unannotated probe
    ##at this point, ii contains indexes for which there are no annotations in existence.
        
    length(ii) # 570
    annot.data <- annot.data[ -ii, ]
    print(paste0("num non-control probesets: ", nrow(annot.data))) # 28440
    
    if(useCurrentSamples)
    {
        validID = getPupID.with.existingSample()
        annot.data = removeInvalidIds(annot.data, validID)
    }
    
    ##merge technical replicates- this MUST happen after replace call, or the work here will be merging bogus data, and will then be undone 
    ##by the replace call
    if(!useCurrentSamples)
    {
        keep  <- which(colnames(annot.data)=="Mouse.142")
        blend <- which(colnames(annot.data)=="Mouse.142_2")
        annot.data[,keep] = .5*(annot.data[,keep]+annot.data[,blend])
        annot.data        = annot.data[,-blend]
        
        keep  <- which(colnames(annot.data.control)=="Mouse.142")
        blend <- which(colnames(annot.data.control)=="Mouse.142_2")
        
        annot.data.control[,keep] = .5*(annot.data.control[,keep]+annot.data.control[,blend])
        annot.data.control        = annot.data.control[,-blend]
    }
    mouse.cols <- grep("Mouse", colnames(annot.data))
    
    mis.probes <- which( apply(as.matrix(annot.data[ , mouse.cols]), 1, function(x){ any(is.na(x)) } ))
    sum(is.na(annot.data[ mis.probes, mouse.cols ])) # 0
    ## remove single probe with missingness
    ## annot.data <- annot.data[ -mis.probes, ]
    nrow(annot.data) # 25196 (alan); 28440 (Will)
    
    
    ## filter for probe sets with duplicated expression
    ## (ie, probe sets mapping to the same gene that have *exactly* the same expression level 
    duplicates = duplicated(cbind(annot.data[,mouse.cols], annot.data$mRNA.Accession))
    annot.data = annot.data[!duplicates,]
    nrow(annot.data) # 24867 (alan); 28006 (3-21 or 0-24); 


    ##plot max expression per probeset
    maxExpressionPerProbeset = apply(annot.data[,mouse.cols], MARGIN = 1, max)
    ## pdf(outm("maxExpressionPerProbeset.pdf"))
    ## hist(maxExpressionPerProbeset)
    ## dev.off()

    ##filter out universally lowly expressed genes
    if(!is.na(prop$mnp$lowExpressionThresh))
    {
        sufficientlyExpressed = which(maxExpressionPerProbeset>=prop$mnp$lowExpressionThresh)
        annot.data = annot.data[sufficientlyExpressed,]
    }


    ##plot number of probes per probeset
    ## pdf(outm("numProbesPerProbeset.pdf"))
    ## hist(probesetInfo$numProbes)
    ## dev.off()

    ## filter out probeset with very few valid probes
    goodProbeset = probesetInfo$meta_probesetId
    if(!is.na(prop$mnp$lowValidProbeThresh))
    {
        goodProbeset = probesetInfo$meta_probesetId[probesetInfo$numProbes>=prop$mnp$lowValidProbeThresh]
    }
    
    annot.data = annot.data[as.character(annot.data$Probe.Set.ID) %in% as.character(goodProbeset),]
    

    
## length(which(apply(exp.mat, MARGIN = 2, max)<4))
    return(list(annot.data=annot.data, annot.data.control = annot.data.control))
}

read.annot.data <- function(annot.old.expression.file, outdir)
{
    annot.data <- readExpressionFile(annot.old.expression.file, old = T)
    
############# Processing annotated but unmasked expression data###
                                        # clean up irregular column names

    nmz = colnames(annot.data)
    names(nmz) = nmz
    nmz[c("mRNA..Source", "mRna...Description", "mRNA...xhyb")] =
        c("mRNA.Source", "mRNA.Description", "mRNA.xhyb")
    colnames(annot.data) = unname(nmz)
        
    ## make mice "Mouse.N" where N is mouse id
    mouse.cols <- grep("\\.rma.gene\\.default\\.Signal", colnames(annot.data))
    colnames(annot.data)[mouse.cols] <- 
        sub("X", "Mouse.", sub("\\.rma.gene\\.default\\.Signal", "", colnames(annot.data)[mouse.cols]))
    
    return(annot.data)
} 


replaceAnnotExpressionWithMasked <- function(annot.data, maskedraw.expression.file)
{
    mouse.cols <- grep("Mouse\\.[0-9]", colnames(annot.data))
    annot.data = annot.data[,setdiff(1:ncol(annot.data), mouse.cols)]
    annot.data = data.table(annot.data,key="Probe.Set.ID")
    
    mask.data <- readExpressionFile(maskedraw.expression.file, old =F)
    colnames(mask.data) <- sub("X", "Mouse.", sub("\\.CEL", "", colnames(mask.data)))
    ##This was a symptom of replacing AFTER removing technical replicate, when we should have been doing it before
    ##mask.data$Mouse.142_2 = NULL #Dont want to accidentally bring this back through the merge
    mask.data = data.table(mask.data,key="probeset_id")
    setnames(mask.data, "probeset_id", "Probe.Set.ID")
    annot.data = merge(mask.data, annot.data)
    annot.data = data.frame(annot.data)
    return(annot.data)
}

readExpressionFile <- function(fle, old)
{
    count = 1
    x = readLines(fle, 10000)
    while(substr(x[count],1,1)=="#" & count<length(x))
    {
        count = count + 1
    }
    x = suppressWarnings(data.frame(fread(fle, na.strings="---", strip.white=TRUE, skip = count-1, header = T, showProgress = F)))

    if(old)
    {
        x[["mRNA...xhyb"]] = as.character(x[["mRNA...xhyb"]])
        badcol = x[["mRNA...xhyb"]]
        badcol[badcol==0]="FALSE"
        badcol[badcol==1]="TRUE"
        x[["mRNA...xhyb"]] = as.factor(badcol)
    }
    return(x)
}


