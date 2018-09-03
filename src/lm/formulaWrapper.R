library(stringr)
source("./util/utils.R")

formulaWrapper = new.env(hash=T)
##expects a covariate string (without the outcome part) of lmer form, e.g., "~ 1 + X1 +X2 + X1:X2 + (Z|X1)"
formulaWrapper$parseCovariateString <- function(covariateString)
{
    covariateInfo = list()
    covariateInfo$original.string = covariateString

    tokens = gsub(pattern = "~", replacement = "", covariateInfo$original.string)
    ##pull out the random effects terms
    ranpat = "\\([^\\(]+\\|.+?\\)"

    ranefs = gregexpr(tokens, pattern=ranpat, perl=T)
    starts = try(as.vector(ranefs[[1]]))
    if(class(starts)=="try-error")
    {
        browser()
    }
    
    ends = starts + attr(ranefs[[1]],"match.length")-1

    covariateInfo$ranef = list()
    if(any(starts>0))
    {
        for(i in seq(1, length.out = length(starts)))
        {
            ranef = substr(tokens, starts[i], ends[i])
            ranef = gsub(ranef, pattern = "\\(|\\)", replacement = "")
            lhs.rhs = strsplit(ranef, split="\\|")[[1]]
            lhs.rhs = str_trim(lhs.rhs)
            
            covariateInfo$ranef[[i]] = list()
            
            covariateInfo$ranef[[i]]$group = str_trim(strsplit(lhs.rhs[[2]], ":")[[1]])
            covariateInfo$ranef[[i]]$components = list()
            components = strsplit(lhs.rhs[[1]], split = "\\+")[[1]]
            components = str_trim(components)


            for (j in 1:length(components))
            {
                if(components[[j]]=="-1"|components[[j]]=="0")
                {
                    next()
                }
                covariateInfo$ranef[[i]]$components[[j]] = str_trim(strsplit(components[[j]], split=":")[[1]])
            }
        }
    }
    


    ##remove the random effects before getting the fixed effects
    tokens = gsub(tokens, pattern=ranpat, replacement = "", perl=T)
    tokens = stringr::str_trim(strsplit(tokens, "\\+")[[1]])
    tokens = tokens[tokens!=""]
    covariateInfo$fixef  = list()
    for(i in 1:length(tokens))
    {
        covariateInfo$fixef[[i]] = str_trim(strsplit(tokens[[i]], ":")[[1]])
    }


    covariateInfo = formulaWrapper$.rebuildString(covariateInfo)

    return(covariateInfo)
}


formulaWrapper$mergeTerms <- function(terms)
{
    terms = lapply(terms, paste, collapse = ":")
    terms = do.call(paste, c(terms, sep=" + "))
    return(terms)
}

formulaWrapper$.rebuildString <- function(covariateInfo)
{
    fixedterms = formulaWrapper$mergeTerms(covariateInfo$fixef)
    modifiedString = paste("~", fixedterms)

    for(aranef in covariateInfo$ranef)
    {
        ranterms.left  = formulaWrapper$mergeTerms(aranef$components)
        ranterms.right = formulaWrapper$mergeTerms(aranef$group)
        ranterms = paste("(",ranterms.left, "|", ranterms.right, ")", sep="")
        modifiedString = paste0(modifiedString, " + ", ranterms)
        
    }
    covariateInfo$modified.string = modifiedString

    return(covariateInfo)
}

formulaWrapper$removeEffect <- function(effect.string, covariateString)
{

    effect.string = formulaWrapper$.rebuildString(formulaWrapper$parseCovariateString(effect.string))$modified.string
    effect.string = gsub(effect.string, pattern = "~ ", replacement = "") 
    
    sinfo = formulaWrapper$parseCovariateString(covariateString)


    fixedsnew = list()
    for(i in 1:length(sinfo$fixef))
    {
        afix = formulaWrapper$mergeTerms(list(sinfo$fixef[[i]]))
        if(effect.string != afix)
        {
            fixedsnew = util$appendToList(fixedsnew, sinfo$fixef[[i]])
        }
    }

    ranefnew = list()
    for(aranef in sinfo$ranef)
    {
        ranterms.left  = formulaWrapper$mergeTerms(aranef$components)
        ranterms.right = formulaWrapper$mergeTerms(aranef$group)
        ranterms = paste("(",ranterms.left, "|", ranterms.right, ")")
        if(ranterms!=effect.string)
        {
            ranefnew = util$appendToList(ranefnew, aranef)
        }
    }

    covariateInfo = list()
    covariateInfo$fixef = fixedsnew
    covariateInfo$ranef = ranefnew
    newstring = formulaWrapper$.rebuildString(covariateInfo)$modified.string
    return(newstring)
    
}

formulaWrapper$onlyKeepEffectAndInteractions <- function(effect.string, covariateString)
{
##    effect.string = c("1", effect.string)
    return(formulaWrapper$.setOperationHelper(effect.string, covariateString, keepMatchingTerms = T))
}


##TODO rewrite to allow removal of just interactions?
#effect.string is one (or several) of the main effect terms (fixed or random) in covariateInfo
formulaWrapper$removeEffectAndInteractions <- function(effect.string, covariateString)
{
    return(formulaWrapper$.setOperationHelper(effect.string, covariateString, keepMatchingTerms = F))
}

formulaWrapper$.setOperationHelper <- function(effect.strings, covariateString, keepMatchingTerms)
{

    matches = function(effect.strings, covTerm)
    {
        foundMatch = F
        for(effect.string in effect.strings)
        {
            effectTerm  = formulaWrapper$parseCovariateString(effect.string)$fixef[[1]]
            if(all(effectTerm %in% covTerm))
            {
                foundMatch = T
                break
            }
        }
        return(foundMatch)
    }
    
    covariateInfo = formulaWrapper$parseCovariateString(covariateString)
    newCovariateInfo       = list()
    newCovariateInfo$fixef = list()
    newCovariateInfo$ranef = list()

    
    for(covTerm in covariateInfo$fixef)
    {
        foundMatch = matches(effect.strings, covTerm)
        if(foundMatch == keepMatchingTerms)
        {
            newCovariateInfo$fixef = util$appendToList(newCovariateInfo$fixef, covTerm)
        }
    }


    for(aranef in covariateInfo$ranef)
    {
        agroup = aranef$group

        
        newranef = list()
        newranef$group = agroup
        newranef$components = list()
        
        for(component in aranef$components)
        {
            if(length(component) ==1 && component == "1")
            {
                covTerm = agroup
            } else {
                covTerm = c(component, agroup)
            }
                
            foundMatch = matches(effect.strings, covTerm)
            if(foundMatch == keepMatchingTerms)
            {
                newranef$components = util$appendToList(newranef$components, component) 
            }
        }
        if(length(newranef$components)>0)
        {
            newCovariateInfo$ranef = util$appendToList(newCovariateInfo$ranef, newranef)
        }
    }

    
    newCovariateInfo = formulaWrapper$.rebuildString(newCovariateInfo)
    return(newCovariateInfo)
}


formulaWrapper$fixefstring <- function(covariateInfo)
{
    
}

formulaWrapper$ranefstring <- function(covariateInfo)
{
    
}



formulaWrapper$appendEffect <- function(effectString, covariateString)
{
    covInfo = formulaWrapper$parseCovariateString(covariateString)
    fixefs = covInfo$fixef
    ranefs = covInfo$ranef

    effectInfo = formulaWrapper$parseCovariateString(effectString)
    newInfo = list()
    newInfo$ranef = c(ranefs, effectInfo$ranef)
    newInfo$fixef = c(fixefs, effectInfo$fixef)
    newInfo = formulaWrapper$.rebuildString(newInfo)

    return(newInfo)
}

formulaWrapper$insertEffect <- function(effectString, index, covariateString)
{

    covInfo = formulaWrapper$parseCovariateString(covariateString)
    fixefs = covInfo$fixef
    ranefs = covInfo$ranef

    effectInfo = formulaWrapper$parseCovariateString(effectString)
    newInfo = list()
    
    newInfo$ranef = util$insertAtIndex(effectInfo$ranef, index, ranefs)
    newInfo$fixef = util$insertAtIndex(effectInfo$fixef, index, fixefs)
    newInfo = formulaWrapper$.rebuildString(newInfo)
    return(newInfo)
} 


