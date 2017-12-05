# TODO: Add comment
# 
# Author: doreper
###############################################################################
library("GenomicRanges")
library("IRanges")
library("data.table")

library("Biostrings")
library("rtracklayer")

source("./genomerep/buildImprintedGenes.R")
## source("./genomerep/cc_founderProbs.R")
source("./genomerep/filtervcf.R")
## source("./genomerep/variantdb/add_cc_toVariantDb.R")
source("./loadParams.R")
source("./utils.R")

buildGenomeData = new.env()
vcfParser  = "nohup python -m genomerep.variantdb.vcf_parser "

buildGenomeData$buildAllData=function()
{
    ##	source("./strainselection/hetcalc/buildIBD.R")
    ##	buildGenomeData$ibdData  = buildIBD$build(newickDir = newickDir)

    dir.create(prop$output)
    dir.create(prop$tmpdir)
    genomeData = new.env(hash=T)

    genomeData$founders = read.table(dat(prop$genome$foundersMap),header=T, sep=",")$founder
    genomeData$karyotype = buildGenomeData$getKaryotype(dat( prop$genome$dnaReferenceFile))
    if(prop$genome$rebuildImprinted) 
    {
        ##writed to imprinted file location
        buildGenomeData$rebuildImprintedFile()
    }
    print("loading imprinted genes as granges")
    genomeData$imprintedGenes = buildGenomeData$loadImprintedGranges(dat(prop$genome$imprintedGenesFile))
    print("done loading imprinted genes")

    ##TODO wrap up the slow part
    print("getting exons")
    if(prop$genome$rebuildImprintedTrancriptsFile | prop$genome$loadexons)
    {

        
        genomeData$exonFeatures               = import(dat( prop$genome$exonReferenceFile), "gtf")
        genomeData$exons                      = buildGenomeData$getExons(dat(prop$genome$exonReferenceFile))

        
        print("done getting exons")

        imprintedGenes = as.data.frame(unique(genomeData$imprintedGenes))
        genomeData$imprintedExonFeatures = genomeData$exons[genomeData$exons$gene_id %in% imprintedGenes$ensembl_gene_id]

        
        ##transcript regions rather than actual transcript names... find out from Fernando if this is what we want
        ## genomeData$imprintedExonFeatures      = genomeData$exons[subjectHits(findOverlaps(genomeData$imprintedGenes, genomeData$exons))]


        
        genomeData$imprintedGeneTranscriptIDs = unique(genomeData$imprintedExonFeatures$transcript_id)
        write.table(genomeData$imprintedGeneTranscriptIDs, file=dat(prop$genome$imprintedTranscriptsFile), row.names=F, quote=F, col.names=F)
        genomeData$exons=NULL
    }
    print("done with fast")
    ## 
    buildGenomeData$buildVariantDbs(genomeData)
    return(genomeData)
}

buildGenomeData$getSnords <- function()
{
    snords = dat(fp(prop$genome$snord))
##    
    dir(snords, "*.csv")
    dfs = list()
    for(snord in dir(snords,"*.csv"))
    {
        id  = sub("\\.csv","",snord )
        df = fread(fp(snords,snord))
        df$snord = id
        dfs = util$appendToList(dfs, df)
    }
    df = rbindlist(dfs)
    return(df)
}

buildGenomeData$rebuildImprintedFile <- function()
{
    outfile         = dat( prop$genome$imprintedGenesFile)
    fullImprintInfo = buildImprintedGenes$collateImprintInfo()
    write.table(fullImprintInfo, outfile, sep="\t", row.names = F)
}

#Returns an IRanges object containing metadata as well as all genes in imprinted genes file
buildGenomeData$loadImprintedGranges <- function(imprintedGenesFile) 
{
    geneMetaData = read.table(imprintedGenesFile, header=T, sep="\t")
    strand = as.character(geneMetaData$strand)
    strand[geneMetaData$strand=="1"]   = "+"
    strand[geneMetaData$strand=="-1"]  = "-"
    
    gr = GRanges(seqnames = geneMetaData$chromosome_name, 
                 ranges = IRanges(geneMetaData$start_position, geneMetaData$end_position, names = geneMetaData$mgi_symbol),
                 strand=strand,
                  geneMetaData[,!(colnames(geneMetaData) %in% c("chromosome_name", 
                                                                       "start_position",
                                                                       "end_position",
                                                                       "strand"))]
                 )
    return(gr)
}

buildGenomeData$buildGr <- function(start, end, strand, chr, names = chr)
{
    ir = IRanges(start = start,
                 end   = end,
                 names = names)

    gr = GRanges(ranges=ir,
                 strand   = strand,
                 seqnames = chr)

    return(gr)
}
    

buildGenomeData$getExons <- function(gtfFile = dat( prop$genome$exonReferenceFile)) 
{
    agtf=rtracklayer::import(gtfFile)
    agtf$score = NULL #this is NA anyway
    exons = agtf[agtf$type=="exon",]
    return(exons)
}

buildGenomeData$getGeneInfo <- function(gtfFile = dat( prop$genome$exonReferenceFile))
{
    print("getting exons")
    exonz = buildGenomeData$getExons(gtfFile)
    exonz = data.table(as.data.frame(exonz))

    print("breaking down by gene")
    geneRanges = exonz[ ,list(start=min(start),
                              end = max(end),
                              chrom = paste(unique(seqnames), collapse=","),
                              strand = paste(unique(strand), collapse=","),
                              gene_name = paste(unique(gene_name), collapse=",")),
                       by = gene_id]

    print("keying by gene_id")
    setkey(geneRanges, "gene_id")
    return(geneRanges)
}

                                        #create a data frame containing length per chromosome
buildGenomeData$getKaryotype <- function(referenceFile =  dat(prop$genome$dnaReferenceFile)) 
{
    refstring = readDNAStringSet(referenceFile)
    chrname   = strsplit(names(refstring), " ")
    chrname   = unlist(lapply(chrname, "[", 1))
    lengths   = width(refstring)
    rm(refstring)
    gc()
    karyotypeFrame = data.frame(chrname = chrname , len=lengths)
    rownames(karyotypeFrame) = chrname
    return(karyotypeFrame)
}

buildGenomeData$writeToBed <- function(ranges, exonBedFile)
{
    export(ranges, exonBedFile,format="bed")	
}

buildGenomeData$buildVariantDbs <- function(genomeData) 
{
    bedfiles = list()
    transcripts = list()
    
    bedfiles[["full"]] = NULL
    transcripts[["full"]] = NULL
    
    exons        = buildGenomeData$getExons(dat( prop$genome$exonReferenceFile))
    ##The location of the bed file describing where all exons are
    bedfiles[["exon"]]      = fp(prop$tmpdir, "exons.bed")
    buildGenomeData$writeToBed(exons,bedfiles[["exon"]])
    rm(exons)
    transcripts[["exon"]] = NULL

    imprintedGenes = buildGenomeData$loadImprintedGranges(dat(prop$genome$imprintedGenesFile))
    ##the output file containing imprinted gene ranges in bed format
    bedfiles[["imprinted"]] = fp(prop$tmpdir, "/imprinted.bed")
    buildGenomeData$writeToBed(imprintedGenes,bedfiles[["imprinted"]])
    rm(imprintedGenes)
    transcripts[["imprinted"]] = dat(prop$genome$imprintedTranscriptsFile)

    karyotype = buildGenomeData$getKaryotype( dat(prop$genome$dnaReferenceFile))
    for(dbtype in c("full","exon", "imprinted"))
    {      
        ##        if(dbtype=="full")
        {
            ##TODO remove
            lim = nrow(karyotype)
            ##lim = 1

            chrsToBuild = prop$variantdb$chr_range
            if(is.na(chrsToBuild))
            {
                chrsToBuild = karyotype$chrname
            }

            ##iterate over all chromosomes in karyotype that are not going to be rebuilt to get the max variant id; we want every chromosome that we subsequently build to have non-overlapping variant ids
            initialVariantId = getInitialVariantId(karyotype,chrsToBuild, dbtype)
           
            for(chr in chrsToBuild)
            {
        
                print(paste0("building db: ",dbtype, "_",chr))
                len = karyotype$len[karyotype$chrname == chr]

                gr = c(buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "+"),
                       buildGenomeData$buildGr(start=1, end = len, chr = chr, strand = "-"))

                chrbed = fp(prop$tmpdir, paste0("chr_",chr,".bed"))
                buildGenomeData$writeToBed(gr, chrbed)
                
                dbname = paste0(prop$variantdb[[dbtype]]$name, "_", chr)
                buildGenomeData$buildSingleVariantDb(dbname = dbname,
                                                     filter.bed = chrbed,
                                                     initialVariantId=initialVariantId,
                                                     transcriptsFile = NULL,
                                                     pythonlogfile   = fp(paste0("./", dbtype,".log")),
                                                     rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
                                                     rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
                                                     rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC,
                                                     limit           = prop$variantdb$var_limit)

                #dirty hack to avoid checking if db already exists, but whatever
 
                initialVariantId = try(getMaxVariantID(dbname) + 1)
                if(class(initialVariantId)=="try-error")
                {
                    initialVariantId = NA
                }
            }
        }
        ## if(dbtype!="full")
        ## {
        ##     print(paste0("building db: ",dbtype))
        ##     buildGenomeData$buildSingleVariantDb(dbname = prop$variantdb[[dbtype]]$name,
        ##                                          filter.bed = bedfiles[[dbtype]],
        ##                                          transcriptsFile = transcripts[[dbtype]],
        ##                                          pythonlogfile   = fp(paste0("./", dbtype,".log")),
        ##                                          rebuildVCF      = prop$variantdb[[dbtype]]$rebuildVCF,
        ##                                          rebuildFounder  = prop$variantdb[[dbtype]]$rebuildFounder,
        ##                                          rebuildCC       = prop$variantdb[[dbtype]]$rebuildCC,
        ##                                          limit           = prop$variantdb$var_limit)
        ## }

  
    }
}      
getInitialVariantId <- function(karyotype, chrsToBuild,dbtype)
{
    print(chrsToBuild)
    maxVariantId = 0
    for(chr in karyotype$chrname)
    {
        if(! chr %in% chrsToBuild)
        {
            dbname = paste0(prop$variantdb[[dbtype]]$name, "_", chr)
            candMax = try(getMaxVariantID(dbname))
            if(class(candMax) == "try-error")
            {
                next
            } else {
                maxVariantId = max(candMax, maxVariantId)
            }
        }
    }

    initialVariantId = maxVariantId + 1
    return(initialVariantId)
}

getMaxVariantID <- function(db,
                            dbhost = prop$variantdb$host,
                            dbuser = prop$variantdb$user,
                            dbpassword = prop$variantdb$password)
{
    con = dbConnect(MySQL(), user=dbuser, password =dbpassword, dbname=db, host=dbhost)   
    maxid = dbGetQuery(con, "select max(variant_id) from variant;")[1,1]
    dbDisconnect(con)
    return(maxid)
}
  
 buildGenomeData$buildSingleVariantDb <- function(dbname,
                                                  filter.bed = NULL,
                                                  initialVariantId = NULL,
                                                  transcriptsFile = NULL,
                                                  pythonlogfile = NULL,
                                                  rebuildVCF = T,
                                                  rebuildFounder = T,
                                                  rebuildCC = T,
                                                 limit = NA)
{
    mc.cores = 1

    founders = read.table(dat(prop$genome$foundersMap),header=T, sep=",")$founder
    
    tmpdir          = prop$tmpdir
    dbhost          = prop$variantdb$host
    dbuser          = prop$variantdb$user
    dbpassword      = prop$variantdb$password
    foundersFile    = tempfile("founder", tmpdir, ".csv")
    write.table(founders, foundersFile, col.names=F, row.names=F, quote=F)
    
    vcfExtensionIn = ".vcf.gz"
    vcftargetdir = fp(prop$tmpdir, dbname)
    print("processing vcf")
    pracma::tic()
    if(rebuildVCF)
    {
        filterVCF$filterVcfsInDir(
          vcfDirIn       = dat(prop$genome$vcfDir),
          vcfDirOut      = vcftargetdir,
          tmpdir         = tmpdir, 
          vcfExtensionIn = vcfExtensionIn, 
          founders       = founders, 
          mc.cores       = mc.cores,
          bedfile        = filter.bed)
    }
    pracma::toc()
    
    
    transcriptString = ""
    if(!is.null(transcriptsFile))
    {
        transcriptString     = paste0("-t ", transcriptsFile)
    }

    initialVariantString = ""
    if(!is.null(initialVariantId))
    {
        initialVariantString = paste0("-i ", initialVariantId)
    }

    limitString = ""
    if(!is.na(limit))
    {
        limitString = paste0("-l ", limit)
    }

    engineString = paste0("-e ", engine)
    
    if(is.null(pythonlogfile))
    {
        pythonlogfile = "./python.log"
    }

    print("parsing vcf")
    pracma::tic()
    pythonlogstring = paste0(" 2> ", pythonlogfile)
    if(rebuildFounder)
    {
##        
        command = paste(vcfParser,
                        vcftargetdir,
                        foundersFile,
                        dbname,
                        dbhost,
                        dbuser,
                        dbpassword,
                        initialVariantString,
                        transcriptString,
                        engineString,
                        limitString,
                        pythonlogstring)
        ##        
        print(command)
        system(command)##wait=T?
    }
    pracma::toc()
    
    print("done building founder db")
    
    print("began building ccdb")
    if(rebuildCC)
    {
        karyotype = buildGenomeData$getKaryotype(referenceFile = dat(prop$genome$dnaReferenceFile))
        rebuild_ccdb(dbname, dat(prop$genome$rilHaplotypeProbsDir), founders, karyotype,dbhost,dbuser, dbpassword)
        print("finished building ccdf")
    }
}

