processBehavior = new.env(hash=T)

##startlechoice = c("Average", "Latency")
startlechoice = c("Average")


processBehavior$buildStartleByGroup <- function(startle.frame, breedLog, all)
{

    sf = copy(startle.frame)
    sf[,trialCount := rank(Trial), by = c("ID", "Group")]
    setnames(sf,old=colnames(sf), new = sub(colnames(sf),pattern="_mean", replacement=""))
    setkey(sf, "ID")
    sf = data.table::dcast(sf, ID + trialCount + Chamber + as50.Average_normalized + as50.Latency_normalized ~ Group, value.var = startlechoice, fun=mean) 

    if(length(startlechoice)==1)
    {
        phenz = c("No_Stim", "PP74", "PP78", "PP82", "PP86", "PP90")
        setnames(sf, phenz, paste0(startlechoice, "_",  phenz))
    }
    setkey(sf, "ID")

    if(is.null(breedLog))
    {
        return(sf)
    }
       
    if(all)
    {
        sf = sf[breedLog]
    } else {
        sf = breedLog[sf]
    }

    return(sf)
}

processBehavior$buildStartleByGroupMean <- function(startle.frame, breedLog, all)
{

    sf = copy(startle.frame)
    sf[,trialCount := rank(Trial), by = c("ID", "Group")]
    setnames(sf,old=colnames(sf), new = sub(colnames(sf),pattern="_mean", replacement=""))
    setkey(sf, "ID")
    sf = data.table::dcast(sf, ID  ~ Group, value.var = startlechoice, fun=mean) 

    ##    
    if(length(startlechoice)==1)
    {
        phenz = c("No_Stim", "PP74", "PP78", "PP82", "PP86", "PP90")
        setnames(sf, phenz, paste0("mean_", startlechoice, "_",  phenz))
    }
    
    setkey(sf, "ID")
    if(is.null(breedLog))
    {
        return(sf)
    }
    
    if(all)
    {
        sf = sf[breedLog]
    } else {
        sf = breedLog[sf]
    }
        
    return(sf)
}   

processBehavior$getAS50 <- function(startle.frame, breedLog, all)
{
    sf = startle.frame[!duplicated(startle.frame$ID)]

    if(is.null(breedLog))
    {
        return(sf)
    }
    if(all)
    {
        derivedSubset = sf[breedLog]
    } else {
        derivedSubset = breedLog[sf]
    }

    return(derivedSubset)
}

processBehavior$getPcs <- function(df, relevPhen)
{
    relevPhen = intersect(colnames(df), relevPhen)
    print(relevPhen)
    
    
    pcFrame   = try(data.frame(df)[,relevPhen,drop=F])
    goodRows = try(rowSums(is.na(pcFrame))==0)
    if(class(goodRows) == "try-error")
    {
        browser()
    }
    pcFrame = pcFrame[goodRows,,drop=F]
    if(class(pcFrame)=="try-error")
    {
        print("ERROR!")
        print(relevPhen)
        return(pcFrame)
    }
    
    
    pcRes = prcomp(formula = ~., data=pcFrame, na.action="na.omit")
    plot(pcRes)
    allpcs   = matrix(NA, nrow(df), ncol = ncol(pcRes$x))
    colnames(allpcs) = colnames(pcRes$x)
    rownames(allpcs) = rownames(df)
    allpcs[rownames(pcRes$x),] = pcRes$x #keeping NAs where appropriate

    return(list(pcs= pcRes, df = cbind(df, allpcs)))    
}

 
