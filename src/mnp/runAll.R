print("starting")
#################### 
# Author: doreper
###############################################################################
##
## The main entry point for running all analyses
##
##

try(rm("prop"))
sources <- function()
{
    if(!exists("prop"))
    {
        source("./util/loadParams.R")
    }    
    source("./mnp/loadAllData.R")
    source("./mnp/behavior/analysis.R")
    source("./mnp/micro/analysis.R")
    source("./mnp/micro/report.R")
    source("./mnp/qpcr/analysis.R")
}

sources()
dir.create(prop$mnp$output, showWarnings = F, recursive = T)
dir.create(fp(prop$mnp$output,"micro"), showWarnings = F, recursive = T)

##load inputs to all analysis; microarray, qpcr, behavior. 
inp  = loadAllData$createAllInputs()
alphas = c(.05)

#Run microarray expression analysis. If fromFile, load it up from file.
fromFile = F
if(!fromFile)
{
    print("running no perm")
    out = micro.analysis$run.noperm(inp)
    print("done running no perm")
    
    util$save(list =ls(), file = outm(fp("micro", "noperm.RData")))
    
    if(prop$mnp$SSVA.numperm>0)
    {
        print("running all perms")
        out = micro.analysis$runallPerms(inp)
        print("done all perms")
        util$save(file = outm(fp("micro", "perm.RData")), list = ls())
    }
} else {
    load(outm(fp("micro","perm.RData")))
    sources()
    inp  = loadAllData$createAllInputs()
    save(file = outm(fp("micro","permz.RData")), list = ls())
}

#Generate plots for microarray analysis
micro.report$reportAnalysis(inp$exp.mat,
                            inp$cov.data.full,
                            out$ident.full$results,
                            ##                            out$results$results,
                            inp$annot.data,
                            inp$probesetInfo,
                            out$threshholds,
                            inp$karyo)

##Run qpcr analysis, and also generate plots for Carmil1 and meg3, (on microarray data too)
qpcr.analysis$run(inp)

##Run behavior analysis
beh.analysis$runAll(inp$phens)

print("done with behavior")

##Run mediation analysis.
source("./mnp/mediation/runMedBayes3.R")

##Generate plots of mediation analysis results.
source("./mnp/mediation/mediationPlots2.R")


## make a few more behavior plots
df=inp$phens$getExperiment("openfield")
df$y=df$pctctr
aplot = plot.poe.data(df, ylabel = "% Center Time", "", mode = 2)
show.and.write(aplot, "x", width = 4, height = 3, fname = "pctctr", mode = 2)

df=inp$phens$getExperiment("SIH")
df$y=df$temp.2
aplot = plot.poe.data(df, ylabel = "Stress-induced Temp", "", mode = 2)
show.and.write(aplot, "SIH_T2", width = 4, height = 3, fname = "SIH", mode = 2)

df=inp$phens$getExperiment("ppi")
df$y=df$ppi.82
aplot = plot.poe.data(df, ylabel = "startle amplitude", "", mode = 2)
show.and.write(aplot, "ppi_82", width = 4, height = 3, fname = "startle", mode = 2)
