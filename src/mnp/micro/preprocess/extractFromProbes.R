## function to parse probesetInfo metadata about microarray probesets
getProbesetId <- function(probesetInfo, geneName, like = F)
{
    f = function(avec)
    {
        return (geneName %in% avec)
    }

    splitnames = strsplit(as.character(probesetInfo$gene_name), ",")
    psid = probesetInfo[unlist(lapply(splitnames,f))]$Probe.Set.ID
    if(like == T)
    {
        psid = probesetInfo[grepl(probesetInfo$gene_name, pattern=geneName)]$meta_probesetId
    }
    return(psid)
}


getProbesetIds <- function(probesetInfo, goi)
{
    psids = rep("", length(goi))
    for(i in 1:length(goi))
    {
        psid = getProbesetId(probesetInfo, goi[i])
        if(length(psid)!=1)
        {
            browser()
            stop(paste0("decide which probeset for ", goi[i], " to use."))
        } else {
            psids[i] = psid
        }
    }
    
    names(psids) = goi
    return(psids)
}

getGeneName <- function(probesetInfo, psid)
{
    geneName = probesetInfo[probesetInfo$meta_probesetID == psid]$gene_name
    return(geneName)
}


getMicroarrayData <- function(probesetInfo, goi, exp.mat)
{
    psids = getProbesetIds(probesetInfo, goi) ##checks that only one probeset corresponds to each gene.
    expression.mic    = exp.mat[, psids]
    colnames(expression.mic) = paste0(goi, ".mic")
    expression.mic    = data.frame(expression.mic)
    expression.mic$ID = rownames(exp.mat)
#    expression.mic$ID = gsub(expression.mic$ID, pattern="Mouse\\.", replacement = "")
    expression.mic    = data.table(expression.mic, key = "ID")
    return(expression.mic)
}
