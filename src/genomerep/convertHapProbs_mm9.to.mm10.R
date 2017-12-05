# TODO: Add comment
# 
# Author: doreper
###############################################################################
library(rtracklayer)
source("./stringutils.R")

toGranges <- function(chromName, oldCoords) 
{
	oldCoords    = IRanges(oldCoords, oldCoords)
	oldCoords    = GRanges(chromName, oldCoords, strand="+")#strand is irrelevant for purposes of liftover
	return(oldCoords)
}

convertAllProbFiles = function(probDirIn, probDirOut, in.out.chain, founderSymbolFile)
{
	dir.create(probDirOut)
	unlink(probDirOut)
	dir.create(probDirOut)
	probIns = stringutils$getFilesWithExtension(probDirIn, "probs.csv")
	
	for(probIn in probIns)
	{
		probout = file.path(probDirOut, paste0(getStrainName(basename(probIn)),".csv"))
#		
		convertSingleProbFile(probIn, probout, in.out.chain, founderSymbolFile)
	}
	
}

convertToMM10Coords <- function(singleFrame, in.out.chain) 
{
	in.out.chain = import.chain(in.out.chain)

	chrNames = read.table("../data/mm9_mm10_chrnames.txt", header =T, sep=",")
	rownames(chrNames) = chrNames$mm10
	mm9Chr = chrNames[as.character(singleFrame$Chromosome), "mm9"]
	
	oldCoords = toGranges(mm9Chr, singleFrame[["Position.B37."]])
	singleFrame[["Position.B37."]] = NULL
	newCoords = liftOver(oldCoords, in.out.chain)
	newCoords = start(newCoords)
	
	lengths = unlist(lapply(newCoords, length))
	
	if(any(lengths>1))
	{
		
	}
	singleFrame[lengths==1,"Pos"] = unname(unlist(newCoords[lengths==1]))
	badInds = unname(which(lengths != 1))
	warning(paste0("could not lift over at ", paste(badInds, collapse=",")))
	singleFrame = singleFrame[lengths==1,]
	return(singleFrame)
}

convertSingleProbFile = function(probin, probout, in.out.chain, founderSymbolFile) 
{
	print(probin)
	singleFrame = read.table(probin, header = T, sep=",")
	
	#TODO get rid of this once mm10 probs come out
	singleFrame = convertToMM10Coords(singleFrame, in.out.chain)

	singleFrame$Chrom = singleFrame$Chromosome
	singleFrame$Chromosome = NULL
	#	
	founderSymbolTable = read.table(founderSymbolFile, sep=",", header=T)
	founderSymb = c("A", "B", "C", "D","E","F","G","H")
	for(founder1 in founderSymb)
	{
		fullname1 = founderSymbolTable[founderSymbolTable$abbreviation==founder1,"founder"]
		for(founder2 in founderSymb)
		{
			fullname2 = founderSymbolTable[founderSymbolTable$abbreviation==founder2,"founder"]
			cnames = colnames(singleFrame)
			cnames[colnames(singleFrame)==paste0(founder1,founder2)] = paste0("Diplo.", paste(fullname1, fullname2,sep="."))
			colnames(singleFrame) = cnames
		}
	}
#	
	write.table(singleFrame, probout,sep=",",quote=F, row.names=F)
}

#returns a vector of all strain names. 
getStrainName <- function(ccProbFilename) 
{
	strainName = strsplit(ccProbFilename,split="_")[[1]][1]
}

in.out.chain      = "../data/b6_reference/mm9ToMm10.over.chain"
probIn            = "../data/cc_probs_2015-04_mm9"
probDirOut        = "../data/cc_probs_2015-05_mm10"
founderSymbolFile = "../data/cc_founderNameMap1410.txt" 
convertAllProbFiles(probIn, probDirOut,in.out.chain, founderSymbolFile)


