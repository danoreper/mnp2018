source("./mnp/loadBreeding.R")
library(data.table)

loadBehavior = new.env(hash = T)

datp <- function(...)
{
    dat("mnp/phenotypes", ...)
##    fp(newdat, ...)
}

loadBehavior$getPhenotypeRepository <- function()
{
    phen                 = list()

    phen$experiment      = list()

##    print("reading pipeline phenotypes")
    phen$breedLog             = read.full.cov() 

    phen$pipeline[[1]] = c("lightdark", "startle", "SIH", "swim", "cocaine", "weight")
    phen$pipeline[[2]] = c("openfield", "sociability", "tail", "cort")

    
    phen$experiment$lightdark = loadBehavior$readLightDarkData()
    phen$experiment$openfield = loadBehavior$readOpenFieldData()
    phen$experiment$cocaine   = loadBehavior$readCocaineData()
    phen$experiment$tail      = loadBehavior$readTailData()
    phen$experiment$swim      = loadBehavior$readSwimData()
    phen$experiment$weight    = loadBehavior$readWeightData()
    phen$experiment$SIH       = loadBehavior$read.SIH.data()
    phen$experiment$cort           = loadBehavior$readCortData()
    phen$experiment$sociability    = loadBehavior$readSociabilityData()
    phen$experiment$startle        = loadBehavior$readStartleData(phen$experiment$weight)

##    print("done reading")
    
    for(framename in names(phen$experiment))
    {
##        print(framename)
        phen$experiment[[framename]]$ID = paste0("Mouse.", phen$experiment[[framename]]$ID)
        phen$experiment[[framename]] = data.table(phen$experiment[[framename]], key="ID")
    }

    ## phen$experiment$startleByGroup = buildStartleByGroup(phen$experiment$startle)
    ## phen$experiment$startleByGroupMean = buildStartleByGroupMean(phen$experiment$startle)
    
    phen$getExperiment <- function(experiments, all = F, useBreedLog = T)
    {
        if(length(experiments)==1)
        {
            df = phen$experiment[[experiments]]
        } else {
            
            df = lapply(FUN = function(exp){phen$experiment[[exp]]}, experiments)
            df = Reduce(f = function(x,y){merge(x,y, on="ID")}, dfs2)
            ## ##the full data set, whether or not the behaviors are recorded
            ##TODO handle repeated column names
        }
        
        ##Just the data set for which behaviors are recorded.
        if(!useBreedLog)
        {
            return(df)
        }
        if(all)
        {
            df = df[phen$breedLog]
        } else {
            df = phen$breedLog[df]
        }

        return(df)
    }



    return(phen)
}

loadBehavior$fixColumns <- function(aframe, replacements=list())
{
##    print(length(colnames(aframe)))
    dupColumns = duplicated(colnames(aframe))
    if(sum(dupColumns)>0)
    {
        print("duplicate columns: ")
        warning(paste(colnames(aframe)[dupColumns], collapse=","))
        dupColData = duplicated(cbind(colnames(aframe), t(as.matrix(aframe))))
        if(!identical(dupColData, dupColumns))
        {
            stop("some columns with duplicated column names have non-duplicated data!")
        }
    }
    aframe = aframe[,!dupColumns]
##    print(length(colnames(aframe)))
    
    pat.replace = append(list(c(" ","."), c("%", "Pct"), c("#", "num"),c("-", ".less."), c("\\(", "."), c("\\)", ".")), replacements)
    ##fix column headers that contain spaces or other symbols
    for(i in 1:length(pat.replace))
    {
        pat         = pat.replace[[i]][1]
        replacement = pat.replace[[i]][2]
        colnames(aframe) = gsub(pattern = pat, replacement=replacement, x = colnames(aframe))
    }
    
    ##remove columns that are redundant to the breeding spreadsheet.
    toRemove    = c("Dietname", "Dam.ID", "Sire.ID","Diet","Strain","DOB","Weaned","Date.on.Diet", "Age.Date","batch","Batch","Pipeline") 
    for (cname in toRemove)
    {
        aframe[[cname]] = NULL
    }
    aframe = data.table(aframe, key="ID")
    return(aframe)
}

loadBehavior$toSeconds <- function(astring)
{
    splitted = strsplit(astring,":")
    zums = 0
    zums = zums + as.numeric(lapply(splitted, "[",1))*3600
    zums = zums + as.numeric(lapply(splitted, "[",2))*60
    zums = zums + as.numeric(lapply(splitted, "[",3))*1
    return(zums)
}

loadBehavior$minutesToSeconds <- function(astring)
{
    splitted = strsplit(astring,":")
    zums = 0
    zums = zums + as.numeric(lapply(splitted, "[",1))*60
    zums = zums + as.numeric(lapply(splitted, "[",2))*1
    return(zums)
}

loadBehavior$readLightDarkData <- function()
{
    afile  <- datp("/LightDark/dietstudy_light_dark.csv")
    lightdark.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    lightdark.frame = loadBehavior$fixColumns(lightdark.frame)
    lightdark.frame$ID = (as.character(lightdark.frame$ID ))
    
    return(lightdark.frame)


}

loadBehavior$readOpenFieldData <- function()
{
    afile  <- datp("OpenField/dietstudy_openfield.csv")
    openfield.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    openfield.frame = loadBehavior$fixColumns(openfield.frame)
    openfield.frame$ID = (as.character(openfield.frame$ID ))
    
    return(openfield.frame)
}

loadBehavior$readCocaineData <- function()
{
    afile  <- datp("Cocaine/cokedata_final_1-14-16.csv")
    cocaine.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    cocaine.frame    = loadBehavior$fixColumns(cocaine.frame)
    cocaine.frame$ID = (as.character(cocaine.frame$ID ))

    return(cocaine.frame)
}

loadBehavior$readCortData <-  function()
{
    ##    afile  <- datp("RestraintStress/cortdata_1-14-16.csv")
    afile  <- datp("RestraintStress/cortdata_1-25-16.csv")

    cort.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    colnames(cort.frame)[colnames(cort.frame)=="10min"] = "Ten.min"
    colnames(cort.frame)[colnames(cort.frame)=="10min-baseline"] = "Ten.min.less.baseline"
    cort.frame = loadBehavior$fixColumns(cort.frame)
    cort.frame$ID = (as.character(cort.frame$ID ))
    return(cort.frame)
}

loadBehavior$readTailData <- function()
{
    ##Tail suspension data
    afile  <- datp("TailSuspension/dietstudy_tst.csv")
    tail.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    tail.frame = loadBehavior$fixColumns(tail.frame)
    tail.frame$ID = (as.character(tail.frame$ID ))
    goodSamples = !tail.frame$PctFreeze.120.less.240sec=="dead"
    tail.frame$PctFreeze.120.less.240sec[tail.frame$PctFreeze.120.less.240sec == "dead"] = NA
    tail.frame$PctFreeze.120.less.240sec = as.numeric(tail.frame$PctFreeze.120.less.240sec)
    ##TODO reconsider    
    tail.frame = tail.frame[goodSamples,]
    return(tail.frame)
}

loadBehavior$readSwimData <- function()
{
    afile  <- datp("ForcedSwim/newFSTdata_12-8-15.csv")
    swim.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    colnames(swim.frame)[colnames(swim.frame)=="id"] = "ID" 
    swim.frame$ID = (as.character(swim.frame$ID ))
    swim.frame    = loadBehavior$fixColumns(swim.frame)
    return(swim.frame)
}

loadBehavior$readWeightData <- function()
{
    ##weight data
    afile <- datp("BodyWeight_8.12.14.csv")
    weight.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    colnames(weight.frame)[colnames(weight.frame)=="Mouse #"] = "ID"
    weight.frame = loadBehavior$fixColumns(weight.frame)
    weight.frame$ID = (as.character(weight.frame$ID ))
    return(weight.frame)
}

loadBehavior$read.SIH.data <- function()
{
    afile  <- datp("SIH/dietstudy_allsih_1-14-16.csv")
    SIH.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    colnames(SIH.frame)[colnames(SIH.frame)=="Mouse #"] = "ID"
    colnames(SIH.frame)[colnames(SIH.frame)=="Cage #"] = "Cage.ID"
    SIH.frame = loadBehavior$fixColumns(SIH.frame)
    colnames(SIH.frame)[colnames(SIH.frame)=="Temp.1"] = "temp.1"
    colnames(SIH.frame)[colnames(SIH.frame)=="Temp.2"] = "temp.2"
    ##correcting any potential mistakes in subtraction
    SIH.frame$Difference = (SIH.frame$temp.2 - SIH.frame$temp.1)

    
    SIH.frame$ID = (as.character(SIH.frame$ID ))
    return(SIH.frame)
}

loadBehavior$readSociabilityData <- function()
{
    afile  <- datp("Sociability/SocialHabSoc.csv")
    sociability.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    sociability.frame$ID = (as.character(sociability.frame$ID ))
    sociability.frame = loadBehavior$fixColumns(sociability.frame, list(c(":",  ".to.")))
    sociability.frame$PctStranger        = as.numeric(sociability.frame$PctStranger)
    sociability.frame$Box.Stranger       = as.factor(sociability.frame$Box.Stranger)
    sociability.frame = sociability.frame[!is.na(sociability.frame$PctStranger),]
    return(sociability.frame)
}	

loadBehavior$readStartleData <- function(weight.frame)
{
    afile  <- datp("Startle_PPI/PPI_Rawdatafiles_combined.csv")
    startle.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    while(sum(colnames(startle.frame)=="")>0)
    {
        startle.frame[,which(colnames(startle.frame)=="")[1]]=NULL
    }	
    colnames(startle.frame)[colnames(startle.frame)=="Subject ID"] = "ID"
    startle.frame$ID        = as.character(startle.frame$ID )
    startle.frame           = loadBehavior$fixColumns(startle.frame)

    startle.frame$Group     = factor(startle.frame$Group)
    startle.frame$Chamber   = as.factor(startle.frame$Chamber) #TODO was this valid

    setkey(weight.frame, "ID")
    startle.frame           = weight.frame[startle.frame]
    
    ave.as50      = startle.frame[startle.frame$Group=="AS50",list(as50.Max = mean(Max),
                                                                   as50.Latency=mean(Latency),
                                                                   as50.Average=mean(Average)), by="ID"]
    startle.frame = startle.frame[as.character(startle.frame$Group)!="AS50"]
    startle.frame = ave.as50[startle.frame]

    startle.frame$as50.Max_normalized     = startle.frame$as50.Max/startle.frame$Body.Weight..g. 
    startle.frame$as50.Average_normalized = startle.frame$as50.Average/startle.frame$Body.Weight..g.
    startle.frame$as50.Latency_normalized = startle.frame$as50.Latency/startle.frame$Body.Weight..g.
    
    startle.frame$Max     = (1 - startle.frame$Max/startle.frame$as50.Max)
    startle.frame$Latency = (1 - startle.frame$Latency/startle.frame$as50.Latency)
    startle.frame$Average = (1 - startle.frame$Average/startle.frame$as50.Average)
    startle.frame$GroupInt = substr(startle.frame$Group,3,5)
    startle.frame$GroupInt = as.numeric(startle.frame$GroupInt)
    return(startle.frame)
}


##                                        #merges everything into a single frame.
## processBehavior$mergeIntoPipelines <- function(phen)
## {
##     prepend = function(fram, prefix)
##     {

##         fram = copy(fram)
##         setnames(fram, setdiff(colnames(fram),"ID"), paste0(prefix, ".", setdiff(colnames(fram),"ID"))) 
##         return(fram)
##     }

##     makepipeline = function(phenfull, outcomes)
##     {
##         pipel1 = phenfull[,colnames(phenfull) %in% outcomes, with=F]
##         badIds = unique(pipel1[which(rowSums(is.na(pipel1))>0)]$ID)
##         print("bad ids, bad records")
        
##         print(badIds)
##         if(length(badIds)>0)
##         {
##             ##
##         }
##         print(data.frame(pipel1[pipel1$ID %in% badIds]))
##         pipel1 = na.omit(pipel1)
##                                         #pipel1$ID = NULL
##         setkey(pipel1,"ID")
##         return(pipel1)
##     }
    
##     pipelineframes = list()
##     framenames   = list()
##     framenames[[1]] = c("lightdark", "startleByGroup", "SIH", "swim", "cocaine", "weight")
##     framenames[[2]]  = c("openfield", "sociability", "tail", "cort")
##     for (pipeline in c(1,2))
##     {

##         for(aframename in framenames[[pipeline]])
##         {
##             ## if(aframename=="sociability")
##             ## {
##             ##     
##             ## }
##             aframe = phen$frame[[aframename]]
##             aframe = aframe[,c("ID",phen$selectedoutcome[[aframename]]),with=F] ##or rawoutcome
##             print(paste0("!!!!!!!", aframename))
##             aframe = prepend(copy(aframe), aframename)
            
            
##             if(length(pipelineframes)>=pipeline)
##             {
##                 pipelineframes[[pipeline]] = merge(pipelineframes[[pipeline]], aframe, all = T)
##             } else {
##                 pipelineframes[[pipeline]] = aframe
##             }
##         }
##         pipelineframes[[pipeline]] = makepipeline(pipelineframes[[pipeline]], colnames(pipelineframes[[pipeline]])) 
##     }
    
##     pipel1full = phen$breedLog[pipelineframes[[1]], all=T]
##     pipel2full = phen$breedLog[pipelineframes[[2]], all=T]
    
##     return(list(pipel1=pipelineframes[[1]], pipel2=pipelineframes[[2]], pipel1full = pipel1full, pipel2full = pipel2full))
## }
