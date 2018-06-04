##TODO: wrap functions up in the environment
source("./loadParams.R")
source("./mnp/loadAllData.R")
source("./mnp/micro/preprocess/evalprobes2.R")
source("./utils.R")
source("./mnp/general.R")

loadDataCrowley = new.env(hash=T)

loadDataCrowley$getMicroData <- function()
{
    crowley.array = fread(dat("crowley15/submission_bundle/MatrixTable/micro_expression.csv"))
    cnames = crowley.array$ID_REF

    crowley.array[,ID_REF:=NULL]
    replaces = gsub(colnames(crowley.array), pattern = "\\.CEL", replacement = "")
    setnames(crowley.array,
             old = colnames(crowley.array),
             new = replaces)
    
    crowley.exp.mat = t(as.matrix(crowley.array))
    colnames(crowley.exp.mat) = cnames
    return(crowley.exp.mat)
}

loadDataCrowley$read.full.cov <- function()
{
    ##TODO: property
    cov.data = fread(dat("crowley15/submission_bundle/Metadata/metadata.csv"))
    cov.data = data.table(data.frame(cov.data))
    colz = str_split(cov.data$Sample.name, pattern="_")
    colz = do.call(rbind, colz)
    colz = data.table(colz)
    setnames(colz, old = c("V1", "V2", "V3"), new=c("toparse", "sex", "tissue"))


    colz$cross = unlist(regmatches(colz$toparse, gregexpr(colz$toparse,  pattern = "[F|G|H][F|G|H]")))
    colz$maternal = substr(colz$cross, 1,1)
    colz$paternal = substr(colz$cross, 2,2)
    colz$cross = paste0(pmin(colz$maternal, colz$paternal),
                        pmax(colz$maternal, colz$paternal))

    colz = colz[,c(.SD, list(direction=as.integer(factor(maternal)))),by = "cross"]
    colz$direction = colz$direction-1.5
    colz[maternal==paternal, direction:=0]
        
    colz$mouse = gsub(colz$toparse,  pattern = "[F|G|H][F|G|H]", replacement = "")
    colz$ID = paste(paste0(colz$maternal, colz$paternal, colz$mouse), colz$sex, colz$tissue, sep="_")

    colz$toparse = NULL
        
    cov.data = colz
    cov.data$cross = factor(cov.data$cross, levels = c("FF", "GG", "HH", "FG","FH","GH"))
    cov.data$sex = factor(cov.data$sex)
    cov.data$tissue = factor(cov.data$tissue)
    cov.data$mouse  = factor(cov.data$mouse)

    ##TODO:refector s.t. mouse is encoded as ID instead
    setkey(cov.data, "mouse")
    return(cov.data)
}

loadDataCrowley$createAllInputs <- function()
{
    pracma::tic()

    recomputeProbeData = T#T ##prop$mnp$recomputeProbeInfo

    computedProbeDataDir = dat("crowley15", "probedata")
    dir.create(computedProbeDataDir, recursive = T, showWarnings = F)
    maskedraw.expression.file <- fp(computedProbeDataDir, "/rma.summary.txt")

    if(recomputeProbeData)
    {
        y = loadDataCrowley$getMicroData()
        z = rbind(colnames(y),y)
        rownames(z)[1] = "probeset_id"
        z = t(z)
        z = data.table(z)
        fwrite(z, file = maskedraw.expression.file)
    }

    print(paste0("using expression file: ", maskedraw.expression.file))
    ##Get the probe information
    if(!recomputeProbeData)
    {
        probeDataBundle = "probeData"
        probedGenes = suppressWarnings(fread(input = fp(computedProbeDataDir, "probeData", "probedGenes.csv"), showProgress = F))
        probesetInfo = suppressWarnings(fread(input = fp(computedProbeDataDir, "probeData", "probesetInfo.csv"), showProgress = F))
    } else {
        variantFile = dat("crowley15", "fgh_b6_variants.csv")
        cel.dir     = dat("crowley15/cel_files")
        probeData         = getProbeInfo(fp(computedProbeDataDir,"probeData"), variantFile, cel.dir = cel.dir, recomputeApt = F, removeBadProbes = T)
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
    annot.data = loadDataCrowley$get.annot.data(annot.old.expression.file = datm(prop$mnp$old.expression.file),
                                maskedraw.expression.file = maskedraw.expression.file,
                                probesetInfo = probesetInfo,
                                output= fp(prop$output, "crowley"))

    ##the annotated data object contains two fields, which contain distinct data: the annotated data less any controls, and just the negative controls.
    ##Consider adding positive controls
    annot.data.control = annot.data$annot.data.control
    annot.data         = annot.data$annot.data


    ##TODO is this filtering needed?
    annot.data         = annot.data[annot.data$Probe.Set.ID %in% probesetInfo$meta_probesetId,] 

    
    ##Get meta data from annotations-- merge altogether into probesetInfo
    annot.data.description  = getAnnotDescriptions(annot.data) 
    setkey(probesetInfo, "meta_probesetId")
    setkey(annot.data.description, "Probe.Set.ID")

    probesetInfo = probesetInfo[annot.data.description]
    

    gc()

    ##for testing purposes, we sometimes want to limt the size of the data set
    lim = prop$crowley$limit
    if(!is.na(lim))
    {
        annot.data = annot.data[1:lim,]
        annot.data.control = annot.data.control[1:lim,]
    }    
    ##Get just the expression matrixes; as opposed to data frame with lots of metadata
    exp.mat.full                 = loadDataCrowley$getExpressionMatrix(annot.data = annot.data)
    exp.mat.control.full         = loadDataCrowley$getExpressionMatrix(annot.data = annot.data.control)
    
    ##TODO read this in the same way as for behavior files
    ##Get the covariates file
    ##only uses the mouse names from exp.mat, nothing else
    cov.data.allSamples =    loadDataCrowley$read.full.cov()

    cov.data.full = toExpCov(cov.data.allSamples, exp.mat.full)

###Use this to check whether there are effects on control probes
    if(prop$mnp$computeEffectsOnControlProbes)
    {
        exp.mat.full = exp.mat.control.full
        prop$mnp$num.sv <<- 0
    }

    exp.mat          = exp.mat.full

    ## TODO fix in evalProbes.
    setnames(probesetInfo, old = "meta_probesetId", new = "Probe.Set.ID")

    karyo = buildGenomeData$getKaryotype(dat( prop$genome$karyotype))

    inp = list(exp.mat            = exp.mat,
               exp.mat.control.full = exp.mat.control.full, 
               annot.data         = annot.data,
               annot.data.control = annot.data.control,
               cov.data.full      = cov.data.full, 
               cov.data.allSamples = cov.data.allSamples,
               probedGenes        = probedGenes,
               probesetInfo       = probesetInfo,
               karyo              = karyo)
    pracma::toc()
    return(inp)
}


##TODO remove the latter 3 methods, pass in get mouse names function

loadDataCrowley$get.mouse.names <- function(namesvec)
{
    mouse.cols <- grepl("FF[0-9]|GG[0-9]|HH[0-9]|FG[0-9]|FH[0-9]|GH[0-9]|GF[0-9]|HF[0-9]|HG[0-9]",
                        namesvec)
    return(mouse.cols)
}

##gets just the expression matrix from the annotated data
loadDataCrowley$getExpressionMatrix <- function(annot.data) 
{
    mouse.cols        = loadDataCrowley$get.mouse.names(colnames(annot.data))
    exp.mat           <- t(as.matrix(annot.data[ , mouse.cols]))
    colnames(exp.mat) = annot.data$Probe.Set.ID
    return(exp.mat)
}


loadDataCrowley$get.annot.data <- function(annot.old.expression.file,
                           maskedraw.expression.file,
                           probesetInfo,
##                           probeTable,
                           output)
{
    annot.data = read.annot.data(annot.old.expression.file, outdir = output)

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
    mouse.cols = loadDataCrowley$get.mouse.names(colnames(annot.data))
    
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

    ##filter out universally lowly expressed genes
    if(!is.na(prop$mnp$lowExpressionThresh))
    {
        sufficientlyExpressed = which(maxExpressionPerProbeset>=prop$mnp$lowExpressionThresh)
        annot.data = annot.data[sufficientlyExpressed,]
    }

    ## filter out probeset with very few valid probes
    goodProbeset = probesetInfo$meta_probesetId
    if(!is.na(prop$mnp$lowValidProbeThresh))
    {
        goodProbeset = probesetInfo$meta_probesetId[probesetInfo$numProbes>=prop$mnp$lowValidProbeThresh]
    }
    annot.data = annot.data[as.character(annot.data$Probe.Set.ID) %in% as.character(goodProbeset),]

    return(list(annot.data=annot.data, annot.data.control = annot.data.control))
}
