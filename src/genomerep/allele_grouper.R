source("./loadParams.R")
allele_grouper = new.env(hash=T)

allele_grouper$alleleGroups = fread(dat(prop$genome$merged_alleles_apm))
allele_grouper$alleleGroups$newtarget = unlist(lapply(strsplit(as.character(allele_grouper$alleleGroups$target), ":"), "[", 2))
setkey(allele_grouper$alleleGroups, "newtarget", "allele")

allele_grouper$getGroup <- function(target, allele)
{

    inp = data.table(newtarget = as.character(target), allele = as.character(allele))
   ## setkey(inp, "newtarget", "allele")
    groups = allele_grouper$alleleGroups[inp]$grouping
    return(groups)
}

allele_grouper$groupDiplotypes <- function(target, diplotypes)
{
    allele1 = unlist(lapply(strsplit(diplotypes, split="/"), "[", 1))
    allele2 = unlist(lapply(strsplit(diplotypes, split="/"), "[", 2))
    allele1 = allele_grouper$getGroup(allele = allele1, target = target)
    allele2 = allele_grouper$getGroup(allele = allele2, target = target)
    return(list(allele1 = allele1, allele2=allele2))
}

