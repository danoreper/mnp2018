source("./loadParams.R")

useCurrentSamples = prop$mnp$useCurrentSamples

read.full.cov <- function()
{
    afile               <- datm(prop$mnp$breedinglog)
    breedLog            <- read.csv(afile, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE, check.names=T)
    breedLog$Run        = NULL
    breedLog$ID         = as.character(breedLog$ID)
    breedLog$ID         = factor(paste0("Mouse.", breedLog$ID))
    breedLog            = data.table(breedLog, key="ID")
    breedLog$Sire.is.b6 = grepl("B6", breedLog$Sire.ID)

    
    ##TODO figure out what order names to use here.
    dietMap = c("Ctl", "ME", "PD", "VDD")
    names(dietMap) = c( "Control B: Vit D, Low Pro", "Control A : Methyl", "Low Protien", "Vitamin D Deficient")
    breedLog$Diet = dietMap[breedLog$Diet]
    breedLog$Diet       = factor(breedLog$Diet, levels = unname(dietMap))

    #breedLog$DietBySire.is.b6 = interaction(breedLog$Diet, breedLog$Sire.is.b6)
    breedLog = breedLog[!duplicated(breedLog),]
    breedLog$Strain = "B6.NOD"
    breedLog$Strain[breedLog$Sire.is.b6] ="NOD.B6"
    breedLog$Strain = factor(breedLog$Strain)

    breedLog$Pipeline = factor(breedLog$Pipeline)
    breedLog$Dam.ID   = factor(breedLog$Dam.ID)
    breedLog$Dam.ID   = factor(paste0(breedLog$Strain, breedLog$Dam.ID))
    breedLog$Batch    = factor(breedLog$Batch)

    if(useCurrentSamples)
    {
        validID = getPupID.with.existingSample()
        breedLog = breedLog[breedLog$ID %in% validID,]
    } 
    
    setkey(breedLog, "ID")
    
    return(breedLog)
}
