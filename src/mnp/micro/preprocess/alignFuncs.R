
runBwa = function()
{
	acommand  = (paste(bwa, 'index', reference, paste('2>', errFile, sep="")))
	print(acommand)
	print(Sys.time());
	x = system(acommand, intern=TRUE)
	
	
	fastaAlnFile      = file.path(fastadir, 'aln.sai')
	errFile           = file.path(fastadir, 'alnOut.txt')
	acommand = paste(bwa, 'aln -t 3 -i 0 -R 5000 -n 2 -o 1', reference, fastaTargetsFile, '>', fastaAlnFile, paste('2>', errFile, sep=""), sep=" ")
	print(acommand)
	print(Sys.time());
	x = system(acommand, intern=TRUE)
	
	alignPrefix = "merged"
	errFile     = file.path(fastadir, 'samseOut.txt')
	acommand2 = paste(bwa, 'samse', reference, fastaAlnFile, fastaTargetsFile, '>', samFile, paste('2>', errFile, sep=""), sep=" ")
	print(acommand2)
	print(Sys.time());
	x = system(acommand2, intern=TRUE)
	
	
	errFile 	    = file.path(fastadir, 'bamOut.txt')
	acommand3 = paste(samtools, 'view', samFile, '-o',  bamFile, '-S', '-b', paste('2>', errFile, sep=""), sep=" ")
	print(acommand3)
	print(Sys.time());
	x = system(acommand3, intern=TRUE)
	
	x = scanBam(bamFile)
	y = x[[1]]
	z=do.call("DataFrame", y)
	z$strand=factor(as.factor(z$strand))
	
	good = z[!is.na(z$pos),]
}

tophat_pipeline <- function() 
{
	bowtieIndexFile = file.path(fastadir, "genome")
	if(rebuildIndex)
	{
		acommand        = paste(bowtie2_build, reference, bowtieIndexFile)
		x               = system(acommand, intern=TRUE)		
	}
	Sys.setenv("BOWTIE2_INDEXES"=bowtieIndexFile)
	if(!cluster)
	{
		Sys.setenv("PATH"=paste0(Sys.getenv("PATH"), ":", bowtiedir))
	}
	readFile = fastaTargetsFile
	
	algopts = "--segment-length 10 "#--b2-very-sensitive --read-edit-dist 1"
	algopts = paste0(algopts, " -G ", gtfFile, " ")
	algopts = paste0(algopts, " --output-dir ", fastadir, " ")
	
	errFile  = file.path(fastadir, 'tophatOut.txt')
	acommand = (paste(tophat, algopts, bowtieIndexFile, readFile, paste('2>', errFile, sep="")))
	print(acommand)
	x = system(acommand, intern=TRUE)
}


getProbeSequences <- function(probeSeq, probId)
{
	probeSeqs = DNAStringSet(probeSeq,use.names=TRUE)
	names(probeSeqs) = paste0("p_",as.character(probesToQuery$probeId))
	return(probeSeqs)
}

getReadFrameForPooledFile <- function(bamFile, trimFile, experimentID, startId)
{	
	print(paste("started loading frame:", Sys.time()))
	
	x=scanBam(bamFile)
	y = x[[1]]
	z=do.call("DataFrame", y)
	z$strand=factor(as.factor(z$strand))
	
	z = as.data.frame(z)
#	print(paste("ended loading frame:", Sys.time()))
	#remove all this stuff since we dont use it to save memory
	z$mrnm  = NULL
	z$rname = NULL
	z$mapq  = NULL
	z$flag  = NULL	
	z$mpos  = NULL
	z$qual  = NULL
	z$qname = as.character(z$qname)
	
	
	{
		samtools = file.path(allProps$resourcesDir, allProps$samtoolsDir, "samtools")
		if(is.null(allProps$samtoolsDir))
		{
			samtools = "samtools"
		}
		samFilename = paste(bamFile,".sam", sep="")
		command = paste(samtools," view ", bamFile, " -h -o ", samFilename,sep="")
#		print(command)
		system(command)
		cigarStrings = read.table(samFilename, skip=5)$V6
		z$cigar = as.character(cigarStrings)
		z$numMM=convertToErrorCounts(z$cigar,"X")
		z$numI=convertToErrorCounts(z$cigar,"I")
		z$numD=convertToErrorCounts(z$cigar,"D")
	}
	
	
#	print(paste("ended cigar process:", Sys.time()))
	if(allProps$turnOffCigar)
	{
		z$cigar = NULL
	}
	if(allProps$turnOffSeq)
	{
		z$seq   = NULL
		
	} else
	{
		z$seq = as.character(z$seq)
	}
	
	z$ID = rep(-1,nrow(z))
	z$readEnd = rep(NA, nrow(z))
	z$alignCounter = rep(-1,nrow(z))
	readSides = c("l","r")
	for (readSide in readSides)
	{
		splittedNames = strsplit(z$qname, readSide)
		ids = sapply(splittedNames, function (x) x[1])
		indexes= ids!=z$qname
		alignCounters = splittedNames[indexes]
		alignCounters = sapply(alignCounters, function (x) x[2])
		
		z$ID[indexes]=as.integer(ids[indexes])
		z$alignCounter[indexes]=as.integer(alignCounters)
		z$readEnd[indexes]=readSide
	}	
#	print(paste("ended readsides process:", Sys.time()))
	z$qname = NULL
	z$readEnd = as.factor(z$readEnd)
	
	d1 = z[z$readEnd=="r",]
	d2 = z[z$readEnd=="l",]
	
	z = merge(d2,d1, by=c("ID","alignCounter"), all.x=T, all.y=T)
	z$alignCounter=NULL
	
	trimData = read.table(trimFile, header=T)
	if(nrow(trimData)!=(max(z$ID)-min(z$ID)+1))
	{
		stop(paste(nrow(trimData), max(z$ID)-min(z$ID)+1))
	}
	trimData$ID = seq(min(z$ID), max(z$ID))
	z = merge(z, trimData, by="ID")
	
#	z$isize = (z$isize)	
	colnames(z) = sub(pattern=".x", replacement = ".3", x=colnames(z), fixed=T)
	colnames(z) = sub(pattern=".y", replacement = ".1", x=colnames(z), fixed=T)
#	z$isize = abs(z$isize)
	z$isize = z$isize.3 #using the 3 end to ensure we get NA's rather than 0 when we have an unpaired set
	z$isize.1 = NULL
	z$isize.3 = NULL
	z$readEnd.1=NULL
	z$readEnd.3=NULL
	z$strand = z$strand.1
	z$strand.1=NULL
	z$strand.3=NULL
	z$pos.3=NULL #THis is already taken into account by isize
	
	z$ID = z$ID + startId
	z$experimentID = experimentID
	z$jackpotIndex = allJackpotTags[z$jackpotIndex]
#	print(paste("done merging:", Sys.time()))
	
#	table(z$ID, z$numMM.1+z$numMM.3)
	colnoindel = (is.na(z$numD.1)|(z$numD.1+z$numI.1)==0) & (is.na(z$numD.3)|(z$numD.3+z$numI.3)==0)
	noindelIDs = z[colnoindel,]$ID
	remainingIDs =  setdiff(z$ID, noindelIDs)
	z = data.table(z)
	setkey(z, ID)
	neededIndel = z[J(remainingIDs),]
	z = rbind(z[colnoindel,],neededIndel)
	z$strand.1=NULL
	return(z)	
}
