library(tools)
library(yaml)

tmpdir = prop$tmp
dir.create(tmpdir,showWarnings=F)

source("./loadParams.R") 

limits          = prop$mnp$limit #NA             ##c("1:200") ##NA
pipelines       = list(c(1,2))  #, c(1), c(2))   ##")##, "c(2)", "c(1,2)")
svas            = c("SSVA")      ##, "'SVA'")
num.svs         = c(7)              ##c(0:15,33)
dietLeaveOneOut = c(NA) #c("StdCtrl", "LowPro", "VitDDef", "MethylSuff")
residualizes   = T
medianCorrects = c(T)


numCores   = prop$mnp$mc.cores
outputroot = outm("multi5") 
dir.create(outputroot, showWarnings = F, recursive = T)

getOutputFile <- function(pipel, sva, num.sv, limit, dietOut, medianCorrect, residualize)
{
    pipel.string = paste(pipel, collapse=".")
    
    limit.string = paste(limit, collapse=".")
    
    outfile = paste0("p",".",pipel.string, ".", sva, ".", num.sv, ".",limit.string, ".", dietOut,
                     ".", "medcorrect_",medianCorrect, ".","residualize_",residualize)
    outfile = fp(outputroot, outfile)
    ##
    return(outfile)
}

runcommand = function(overfile)
{
    ##memlim = max(10, ceiling(numCores*3.8))
    memlim = 22
    argz = paste0(" '--args ", "../config/defaultCluster.yaml ", overfile,"' ")
    command = paste0("bsub -R 'span[hosts=1]' -n ", numCores, " -M ", memlim ," -q week R CMD BATCH --no-save --no-restore", argz, "./mnp/runAll.R ",
                     fp(gsub(overfile, pattern = "override.yaml", replacement = "out.ROut")))
    print(command)
    
    system(command)
}

for(limit in limits)
{
    for(residualize in residualizes)
    {
        for(medianCorrect in medianCorrects)
        {
            for (pipel in pipelines)
            {
                for(sva in svas)
                {
                    for(num.sv in num.svs)
                    {
                        for(dietOut in dietLeaveOneOut)
                        {
                            output = getOutputFile(pipel,sva, num.sv, limit, dietOut, medianCorrect, residualize)
                            dir.create(output, recursive = T, showWarnings=F)
                            overfile = file.path(output, "override.yaml")
                            print(overfile)
                            
                            mnpover = list()
                            mnpover$mnp = list()
                            mnpover$mnp[["output"]] = output
                            mnpover$mnp[["pipelines"]] = pipel
                            mnpover$mnp[["surrogat"]]  = sva
                            mnpover$mnp[["limit"]]     = limit
                            mnpover$mnp[["mc.cores"]]  = numCores 
                            mnpover$mnp[["num.sv"]] = num.sv
                            mnpover$mnp[["groupedDiets"]] = NA
                            mnpover$mnp[["medianAdjust.p.values"]]= medianCorrect
                            mnpover$mnp[["residualize.expression.pre.straindiet"]] = residualize
                            
                            if(!is.na(dietOut))
                            {
                                mnpover$mnp$groupedDiets = list()
                                otherGroup = setdiff(dietLeaveOneOut, dietOut)
                                mnpover$mnp$groupedDiets[["other"]] = otherGroup
                                mnpover$mnp$groupedDiets[[dietOut]] = dietOut
                            }
                            
                            write(yaml::as.yaml(mnpover), overfile,append=F)
                            runcommand(overfile)
                        }
                    }
                }
            }
        }
    }
}
