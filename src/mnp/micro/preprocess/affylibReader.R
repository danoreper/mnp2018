# read in and parse affymetrix library files.
# 
# Author: doreper
###############################################################################
library(affxparser)
affylibReader = new.env(hash=T)

affylibReader$getCelChips <- function(basedir)
{
    celnames = dir(basedir)
    celnames = file.path(basedir, celnames)
    celChips = unlist(lapply(celnames, function(x) {readCelHeader(x)$chiptype}))
}

affylibReader$readMicroArrayMetaData <- function(pgfFile                 = datm(prop$mnp$pgfFile),
                                                 clfFile                 = datm(prop$mnp$clfFile),
                                                 probesetSeqLocationFile = datm(prop$mnp$probesetSeqLocationFile),
                                                 mpsFile                 = datm(prop$mnp$mpsFile),
                                                 qccFile                 = datm(prop$mnp$qccFile))
{
    pgf = readPgf(pgfFile)
    mps = read.table(mpsFile, header=T, sep="\t")
    mps = data.table(mps)
    mps = mps[,list(
        transcript_cluster_id = transcript_cluster_id,
        single_probeset_id = unlist(strsplit(as.character(probeset_list),split=" ")))
      , by=probeset_id]
    
    setnames(mps, c("probeset_id", "single_probeset_id"), c("meta_probeset_id", "probeset_id"))
                                        #   
    probeTable    = data.frame(probeId    = pgf$probeId, probeGcCount = pgf$probeGcCount, probeInterrogationPosition = pgf$probeInterrogationPosition, probeLength = pgf$probeLength, probeSequence = pgf$probeSequence, probeType = pgf$probeType)
    atomTable     = data.frame(atomId     = pgf$atomId,  atomExonPosition = pgf$atomExonPosition)
    probeSetTable = data.frame(probesetId = pgf$probesetId, probesetName = pgf$probesetName, probesetType = pgf$probesetType)
    
    ##fill in probeset ids for each atom
    probesetStartAtoms   = c(pgf$probesetStartAtom, length(pgf$atomId)+1)
    atomTable$probesetId = rep(probeSetTable$probesetId, diff(probesetStartAtoms))
    
    ##fill in atom ids for each probe
    atomStartProbes    = c(pgf$atomStartProbe, length(pgf$probeId)+1)
    probeTable$atomId  = rep(atomTable$atomId, diff(atomStartProbes))
    
    probeTable            = data.table(probeTable,    key="probeId")
    atomTable             = data.table(atomTable,     key="atomId")
    probeSetTable         = data.table(probeSetTable, key="probesetId")
    probeTable$probesetId = atomTable[i=J(atomId=probeTable$atomId),]$probesetId
    

	
    probeCounts = probeTable[,j=list(probe_count=.N),by=probesetId]
    setkey(probeCounts, probesetId)
	
    clf = readClf(clfFile)
    chipLocationTable = data.frame(x = clf$x, y = clf$y, probeId = clf$id)
    chipLocationTable = data.table(chipLocationTable, key="probeId")
    
    ##probesetLocations = read.delim(probesetSeqLocationFile,comment.char="#", sep=",")

    probesetLocations = fread(probesetSeqLocationFile, skip = 22, sep=",", stringsAsFactors=T)


    probesetLocations = data.table(probesetLocations, key = "probeset_id")
    probesetLocations$start = as.numeric(as.character(probesetLocations$start))
    probesetLocations$stop  = as.numeric(as.character(probesetLocations$stop))
    
    
    probeSetTable = probeSetTable[probesetLocations]




    ##probeSetTable$probesetType = NULL # remove the redundant and slightly conflicting column (in terms of the number of FLmRNA->unmapped
    
    mps$probeset_id = as.character(mps$probeset_id)
    setkey(mps, "probeset_id")
    probeTable$meta_probesetId = mps[i=J(probeset_id=as.character(probeTable$probesetId)),]$meta_probeset_id

    qcc = fread(qccFile, skip = 5)
    ## setkey(qcc, "probeset_id")
    ## qcc$probeset_name = NULL
    ## qcc$quantification_in_header = NULL
    ## mps$qcc.type = qcc$group_name[match(as.character(mps$probeset_id)),
    ##                                               as.character(qcc$probeset_id))]
    
    probeTable$qcc.type = qcc$group_name[match(as.character(probeTable$meta_probesetId),
                                           as.character(qcc$probeset_id))] 

    probeTable$qcc.type[is.na(probeTable$qcc.type)] = "main"
    ps.type  = probeSetTable$probesetType[match(as.character(probeTable$probesetId),
                                                           as.character(probeSetTable$probesetId))]

    probeTable$fulltype = factor(interaction(probeTable$qcc.type, ps.type, sep =","))


    return(list(pgf=pgf, clf=clf, mps=mps, 
                probeTable    = probeTable,
                probeSetTable = probeSetTable,
                atomTable = atomTable,
                chipLocationTable = chipLocationTable))
}

affylibReader$getControlProbesets <- function(mic)
{
    controlProbesetIds    = mic$probeSetTable$probesetId[mic$probeSetTable$probeset_type=="normgene->intron"]
    return(controlProbesetIds)
}

affylibReader$getControlProbeIds <- function(mic) 
{
    controlProbesetIds    = affylibReader$getControlProbesets(mic = mic)
    controlProbeIds       = mic$probeTable$probeId[mic$probeTable$probesetId %in% controlProbesetIds]
    return(controlProbeIds)
}

affylibReader$plotControlProbes <- function(mic) 
{
    controlProbeIds       = affylibReader$getControlProbeIds(mic = mic)
    controlLocations      = mic$chipLocationTable[J(probeId = controlProbeIds)]
    agrid                 = matrix(0, nrow= as.integer(mic$clf$header$rows), ncol = as.integer(mic$clf$header$cols))
    agrid[cbind(controlLocations$y+1, controlLocations$x+1)] = 1
    image(z=agrid, x=1:nrow(agrid), y=1:ncol(agrid))
}
