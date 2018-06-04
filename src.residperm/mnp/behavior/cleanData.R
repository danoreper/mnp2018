source("./mnp/general.R")
source("./mnp/loadBreeding.R")
library(data.table)

loadBehavior = new.env(hash = T)
newdat = "../output/newdat/mnp/phenotypes"
dir.create(newdat, recursive = T, showWarnings=F)
datp <- function(...)
{
    dat("mnp/phenotypes", ...)
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


loadBehavior$readLightDarkData <- function()
{
    afile  <- datp("/LightDark/dietstudy_light_dark.csv")

    lightdark.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)

    lightdark.frame = lightdark.frame[,c("ID",
                                         "Total Distance",
                                         "Total Distance Dark",
                                         "Total Distance Light",
                                         "% Time Dark",
                                         "% Time Light",
                                         "Total Transitions")]
    dir.create(fp(newdat, "LightDark"))
    fwrite(lightdark.frame,   fp(newdat, "LightDark/dietstudy_light_dark.csv"))
}

loadBehavior$readOpenFieldData <- function()
{
    afile  <- datp("OpenField/dietstudy_openfield.csv")

    openfield.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)

    openfield.frame = openfield.frame[,c("ID",
                                         "totdist", 
                                         "pctctr",
                                         "avgvel",
                                         "jmpcts",
                                         "vertcts",
                                         "boli")]
    dir.create(fp(newdat, "OpenField"))
    fwrite(openfield.frame, fp(newdat, "OpenField/dietstudy_openfield.csv"))
}

loadBehavior$readCocaineData <- function()
{
    afile  <- datp("Cocaine/cokedata_final_1-14-16.csv")

    cocaine.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)

    cocaine.frame = cocaine.frame[,c("ID",
                                     "D3-D2",
                                     "dist_d1",
                                     "dist_d2",
                                     "dist_d3")]
    
    dir.create(fp(newdat, "Cocaine"))
    fwrite(cocaine.frame, fp(newdat, "Cocaine/cokedata_final_1-14-16.csv"))
}

loadBehavior$readCortData <-  function()
{
    ##    afile  <- datp("RestraintStress/cortdata_1-14-16.csv")
    afile  <- datp("RestraintStress/cortdata_1-25-16.csv")
    cort.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    cort.frame = cort.frame[,c("ID",
                               "baseline",
                               "10min",
                               "10min-baseline")]
                        
    dir.create(fp(newdat, "RestraintStress"))
    ##    fwrite(cort.frame, fp(newdat, "RestraintStress/cortdata_1-14-16.csv"))
    fwrite(cort.frame, fp(newdat, "RestraintStress/cortdata_1-25-16.csv"))
}

loadBehavior$readTailData <- function()
{
    ##Tail suspension data
    afile  <- datp("TailSuspension/dietstudy_tst.csv")
    tail.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    tail.frame = tail.frame[,c("ID", "%Freeze 120-240sec")]

    dir.create(fp(newdat, "TailSuspension"))
    fwrite(tail.frame, fp(newdat, "TailSuspension/dietstudy_tst.csv"))
}

loadBehavior$readSwimData <- function()
{
    afile  <- datp("ForcedSwim/newFSTdata_12-8-15.csv")
    swim.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    swim.frame = swim.frame[,c("id", "Arena","pctimmob")]
    dir.create(fp(newdat,"ForcedSwim"))
    fwrite(swim.frame, fp(newdat, "ForcedSwim/newFSTdata_12-8-15.csv"))
}

loadBehavior$readWeightData <- function()
{
    ##weight data
    afile <- datp("BodyWeight_8.12.14.csv")
    weight.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    weight.frame = weight.frame[,c("Mouse #","Body Weight (g)")]
    fwrite(weight.frame, fp(newdat, "BodyWeight_8.12.14.csv"))
}

loadBehavior$read.SIH.data <- function()
{
    afile  <- datp("SIH/dietstudy_allsih_1-14-16.csv")
    SIH.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    SIH.frame = SIH.frame[,c("Mouse #", "Order", "Temp 1", "Temp 2")]
    SIH.frame$deltaT = (SIH.frame[["Temp 2"]] - SIH.frame[["Temp 1"]])
    dir.create(fp(newdat, "SIH"))
    fwrite(SIH.frame, fp(newdat, "SIH/dietstudy_allsih_1-14-16.csv"))
}

    
loadBehavior$readSociabilityData <- function()
{

    afile  <- datp("Sociability/SocialHabSoc.csv")
    sociability.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)
    sociability.frame = sociability.frame[,c("ID", "Box Stranger", "PctStranger", "TRANSTOTAL")]

    dir.create(fp(newdat, "Sociability"))
    fwrite(sociability.frame, fp(newdat, "Sociability/SocialHabSoc.csv"))
}	

loadBehavior$readStartleData <- function(weight.frame)
{
    afile  <- datp("Startle_PPI/PPI_Rawdatafiles_combined.csv")
    startle.frame <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=F)

    startle.frame = startle.frame[,c("Subject ID", "Group", "Chamber",  "Max", "Latency", "Average")]
    dir.create(fp(newdat, "Startle_PPI"))
    fwrite(startle.frame, fp(newdat, "Startle_PPI/PPI_Rawdatafiles_combined.csv"))
}

