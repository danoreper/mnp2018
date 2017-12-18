# TODO: Add comment
# 
# Author: doreper
###############################################################################
library(igraph)
library(IRanges)
library(GenomicRanges)
source("./strainselection/idConversion.R")

buildIBD = new.env(hash=T)

buildIBD$removeEmpties = function(vec.of.chars)
{
	toKeep = (vec.of.chars!="")
	return(vec.of.chars[toKeep])
}

# Converts an ibd newick file to a data frame containing identity by descent information per region. 
# Converts newick strings to contain strain symbols (A,B,C,D...) instead of lengthy strain names.
#
# Inputs
# fpath: the string location of a single .txt file containing IBD information for a single chromosome
# strainNameToSymbolMap: a mapping from strain name of the sort ("A/J"->"A", "C57BL/6J"->"B")
#
# Outputs
# dataframe containing
# 1)chr: chromosome integer. No notion of X or Y 
# 2)start: genomic start
# 3)end: genomic end
# 3)newick: newickstring for the start to end region.
# 4)TODO: figure out what the format of the numbers after the string is, build data structure for distances.
# 
buildIBD$getRawDataForChr <- function(fpath, postscript, lookuptable) 
{
	rawdata = read.table(file = fpath, header=F)
	rawdata = rawdata[,c(2,3,4)]
	rawdata$chr = rep(paste("chr",postscript, sep=""), nrow(rawdata))
	colnames(rawdata) = c("start", "end", "newick", "chr")
	for(pattern in names(lookuptable))
	{
		replacement = lookuptable[pattern]
		rawdata$newick = sub(rawdata$newick, pattern = pattern, replacement=replacement)
	}
	rawdata$newick = sub(rawdata$newick, pattern="\\[", replacement="")
	rawdata$newick = sub(rawdata$newick, pattern="\\]", replacement="")
	return(rawdata)
}

buildIBD$getDisjointIntervals <- function(ir) 
{
#	
	disj <- IRanges(sort(unique(c(start(ir), head(end(ir),  -1)))),
			sort(unique(c(end(ir),   tail(start(ir),-1)))))
	
	irShrunk = IRanges(pmin(start(ir)+1,end(ir)), pmax(end(ir)-1,start(ir)+2))
	disjGood = disj[disj %in% irShrunk]
	disjBad  = disj[!(disj %in% irShrunk)]
#	disj = IRanges(start(disj), end(disj))
	return(list(disjGood=disjGood, disjBad=disjBad))
}

buildIBD$getNewicksPerDisjointRegion <- function(ir,newickPerRange, disj) 
{
	overlaps          = findOverlaps(ir, disj)
	numDisjointRanges = length(disj)
	numNewickStrings  = length(newickPerRange)
	#allocate a matrix with enough space- and extra - to store
	#every newick string for every disjoint interval.
	#the rows correspond to disjoint intervals.
	#the column indexes dont mean anything, the newick string
	#are stored one after another as they are assigned to disjoint intervals
	#numHitsPerDisjoint keeps track of how many slot per disjoint were actually used,
	#and where the next free slot is.
	hits = matrix(nrow=numDisjointRanges, ncol=numNewickStrings)
	numHitsPerDisjoint = rep(0, numDisjointRanges)
	for (o in 1:length(overlaps))
	{
		newickIndex  = queryHits(overlaps)[o]
		newickString = newickPerRange[newickIndex]
		
		disjIndex     = subjectHits(overlaps)[o]
		numHitsForDisj= numHitsPerDisjoint[disjIndex]+1 #increment number of hits
		hits[disjIndex, numHitsForDisj] = newickString # put the newick string in the next empty slot
		numHitsPerDisjoint[disjIndex] = numHitsForDisj #store the incremented number of hits
	}
	newickMergedList = list()
#	
	for(disjIndex in 1:numDisjointRanges)
	{
		newickMergedList[[disjIndex]] = hits[disjIndex, 1:numHitsPerDisjoint[disjIndex]]
	}
	return(newickMergedList)
}

buildIBD$formAdjMatrix = function(newicks, lookuptable = lookuptable, leafsep=",")
{
	adjMatrix = matrix(0, nrow=length(lookuptable), ncol=length(lookuptable))
	rownames(adjMatrix) = lookuptable
	colnames(adjMatrix) = lookuptable
	splitNewicks = strsplit(newicks, split="\\|")
	splitNewicks = lapply(splitNewicks,FUN =buildIBD$removeEmpties)
#	
	for(splitNewick in splitNewicks)
	{
		leafGroups = strsplit(splitNewick,split=leafsep)
		for(leafGroup in leafGroups)
		{
			connections = expand.grid(as.character(leafGroup), as.character(leafGroup), stringsAsFactors = F)
			adjMatrix[connections[,1], connections[,2]] = 1
		}
	}
	return(adjMatrix)
}

buildIBD$formTransitiveClosure = function(adjMatrix)
{
	for(k in 1:nrow(adjMatrix))
	{
		for(i in 1:nrow(adjMatrix))
		{
			for(j in 1:nrow(adjMatrix))
			{
				adjMatrix[i,j] = adjMatrix[i,j] || (adjMatrix[i,k] && adjMatrix[k,j])
			}
		}
	}
	return(adjMatrix)
}


buildIBD$.build <- function(newickDir, lookuptable)
{
	filez = dir(path=newickDir,pattern="sef_chr[0-9XY]*\\.txt")
	

	adjMatrixesPerFile = list()
	startsPerFile      = list()
	endsPerFile        = list()
	chrPerFile         = list()
	mergedIBDPerFile   = list()
	
	for(fileCounter in 1:length(filez))
	{
		afile = filez[fileCounter]
		fpath = file.path(newickDir, afile)
		postscript = sub(x=sub(x=afile,pattern="sef_chr", replacement=""), pattern="\\.txt", replacement="")
		if(postscript=="20")
		{
			postscript = "X"
		}
#		
		print(paste("ibd processing ",postscript))
		
		rawdata = buildIBD$getRawDataForChr(fpath = fpath, postscript = postscript, lookuptable = lookuptable)
		ir = IRanges(rawdata$start, rawdata$end)
		
		disjList = buildIBD$getDisjointIntervals(ir = ir)
		
		disjGood = disjList$disjGood
		disjBad  = disjList$disjBad
		disj     = sort(c(disjGood, disjBad))
		
		#expand the bad intervals a little bit wehen looking for newick hits. need to have some sort of default
		disjForNewick = sort(c(disjGood,  IRanges(pmax(start(disjBad)-1, 0), end(disjBad)+1))) 
#		
		#shrink the iranges to avoid getting the boundary regions, makign eveything overly conservative
		irShrunk = IRanges(pmin(start(ir)+1,end(ir)), pmax(end(ir)-1,start(ir)+2))
#	newicksPerDisjGoodRegion  = getNewicksPerDisjointRegion(ir = irShrunk, disj = disjGood, newickPerRange = rawdata$newick)
#	newicksPerDisjBadRegion   = getNewicksPerDisjointRegion(ir = irShrunk, disj = disjBadExpanded, newickPerRange = rawdata$newick)
		
		newicksPerDisjointRegion  = buildIBD$getNewicksPerDisjointRegion(ir = irShrunk, disj = disjForNewick, newickPerRange = rawdata$newick)
		newicksPerDisjointRegion 
		
		founderAdjPerDisjRegion   = list() 

		mergedIBDStrings = rep("", length(newicksPerDisjointRegion))
		for(i in 1:length(newicksPerDisjointRegion))
		{
			newicks = newicksPerDisjointRegion[[i]]
			adjMatrix = buildIBD$formAdjMatrix(newicks, lookuptable = lookuptable)
			transitiveClosure = buildIBD$formTransitiveClosure(adjMatrix)
			founderAdjPerDisjRegion[[i]] = transitiveClosure
			
			grf = graph.adjacency(transitiveClosure, mode="undirected")
			clusts = clusters(grf, "strong")
			numClusts = clusts$no
		
			for(clust in 1:numClusts)
			{
				mergedIBDStrings[i] = paste(mergedIBDStrings[i], (paste(lookuptable[clusts$membership == clust], collapse="")), sep="")
				mergedIBDStrings[i] = paste(mergedIBDStrings[i], "|", sep="")
			}	
#			
		}
		
		
		chrPerFile[[fileCounter]]         = rep(paste("",postscript,sep=""), length(disj))
		startsPerFile[[fileCounter]]      = start(disj)
		endsPerFile[[fileCounter]]        = end(disj)
		adjMatrixesPerFile[[fileCounter]] = founderAdjPerDisjRegion
		mergedIBDPerFile[[fileCounter]]   = mergedIBDStrings
	}

	ibdFrame = 	data.frame(
						   chr   = do.call(c, chrPerFile), 
			   			   start = do.call(c, startsPerFile),
	   		  			   end   = do.call(c, endsPerFile),
			   			   mergedIBD = do.call(c, mergedIBDPerFile)
				   			)
	   
	write.table(ibdFrame,file="frameIBD.txt", sep=",",row.names=F,quote=F)
	
	starts      = do.call(c, startsPerFile)
	ends        = do.call(c, endsPerFile)
	chrs        = do.call(c, chrPerFile)
	ir          = GRanges(ranges = IRanges(start=starts, end=ends), seqnames = chrs, ibdID = 1:length(starts)) 
	adjMatrixes = do.call(c, adjMatrixesPerFile)
	
	return(list(ir=ir, adjMatrixes = adjMatrixes))
}

buildIBD$build=function(newickDir = newickDir)
{
	return(buildIBD$.build(newickDir = newickDir, lookuptable = idConversion$founderNameToSymbolMap))
}


