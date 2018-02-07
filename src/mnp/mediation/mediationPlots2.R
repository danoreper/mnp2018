library(data.table)
library(Cairo)
library(flextable)
library(officer)
##library(ReporteRs)
source("./enrichmentTesting.R")
source("./mnp/loadAllData.R")
source("./mnp/mediation/mediationBayes3.R")
source("./mnp/mediation/BDmodel2.R")
source("./utils.R")
source("./mnp/micro/preprocess/extractFromProbes.R")

## discardMissing = T
## mergeQPCR      = F

getFile = function( mediator, outcome, discardMissing, mergeQPCR, inp=NULL)
{
    if(is.null(inp))
    {
        inp  = loadAllData$createAllInputs()
    }

    froot = fp(prop$mnp$output, "mediation")
    if(outcome == "behavior")
    {
        fle = fp(froot, paste0(discardMissing, "_", mergeQPCR, "_", mediator,"_",outcome),"mediationAll.csv")
    } else {
        fle = fp(froot, paste0(discardMissing, "_", mergeQPCR, "_", mediator,"_",outcome),"mediationAll.csv")
    }

    
    df = fread(fle)
    df[mediator_name=="Lrrc16a", mediator_name:="Carmil1"]
    df$coef.a = NULL
    df$coef.b = NULL
    df$moderators = NULL
    df$mediator.id = as.character(df$mediator.id)
    df = inp$probesetInfo[df, on = c(Probe.Set.ID="mediator.id")]
    setnames(df, old="Probe.Set.ID", new="mediator.id")
    df[,imprinted:= !is.na(minDistToImprinted)& minDistToImprinted <=100]
    df$p.value.a = NULL
    df$p.value.b = NULL
    df$lb.ab     = NULL
    df$ub.ab     = NULL
    df$outcome = as.character(df$outcome.id)
    if(outcome =="micro")
    {
        df$outcome ="Carmil1"
    }
    df$outcome.id = NULL
    setnames(df, c("coef.c", "p.value.c"), c("coef.cprime", "p.value.cprime"))
    ##Remove lrrc16 as a mediator for itself
    df = df[df$outcome!=mediator_name]


    
    ##merge the mediator name and id, in such a manner that the id is appended if the name isnt unique.
    singleOut = df$outcome[1]
    x = df[outcome==singleOut]
    x = table(x$mediator_name)
    bads = names(x[x>1])

    df$mediator = df$mediator_name
    bads = df$mediator %in% bads
    df$mediator[bads] = paste0(df$mediator[bads],"_",df$mediator.id[bads])
    df$mediator.id = NULL
    df$mediator_name = NULL

    


    ##TODO merge id and name, get rid of both subscripts.
##    df$p.strain.on.mediator = NULL
    df$suppressor = (df$coef.cprime*df$coef.ab<0)
    colz = c("p.value.ab", "coef.ab", "coef.cprime", "p.value.cprime", "imprinted",  "suppressor",  "mediator","outcome")
    df = df[,colz, with = F]
    setcolorder(df, colz)
    return(df)
}



toAlias <- function(behaviorCol)
{
    behaviorCol = sub(behaviorCol, pattern="_", replacement = ".")
    nameMap = getNameMap()
    nameMap2 = paste0("PPI", c("74", "78", "82", "86", "90"))
    
    names(nameMap2) =  c(paste0("startle.mean_Average_PP", c("74", "78", "82", "86", "90")))
    nameMap = c(nameMap, nameMap2)
    behaviorCol = nameMap[behaviorCol]
    return(behaviorCol)
}

getNameMap <- function()
{

    nameMap = beh.analysis$nameMap()
    nameMap["openfield.totdist"]             = "OF Total Distance"
    nameMap["lightdark.Total.Distance"]      = "LD Total Distance"
    nameMap["tail.PctFreeze.120.less.240sec"] = "TailSusp Pct Immobility"
    nameMap["swim.pctimmob"]                 = "Swim Pct Immobility"
    return(nameMap)

}



plotEff <- function()
{
    aplot = ggplot(full, aes(x=micro.obs, y = qpcr.obs, color = Strain))
    aplot = aplot + geom_point()
    pdf(fp(outdir, "qpcr_micro_scatter.pdf"))
    print(aplot)
    dev.off()
}


plot.lrrc.airn <- function()
{
    adf = full[!is.na(Strain)]
##    aplot = ggplot(adf, aes(x=micro.obs, y = airn.obs, color = Strain))
    aplot = ggplot(adf, aes(y=micro.obs, x = airn.obs, color = Strain))
    
    aplot = aplot + geom_point()
    pdf(fp(outdir, "airn_micro_scatter.pdf"))
    print(aplot)
    dev.off()
}


plot.p.comp <- function()
{
    df.m = getFile("micro", "behavior", discardMissing=T, mergeQPCR = F)
    df.q = getFile("qpcr", "behavior", discardMissing = T, mergeQPCR =T)
    df.m.lrrc = df.m[mediator =="Lrrc16a"]
    setorder(df.q, outcome)
    setorder(df.m.lrrc, outcome)
    df.compare = data.frame(outcome = df.q$outcome,
                            qpcr.p.value = df.q$p.value.ab,
                            micro.p.value = df.m.lrrc$p.value.ab)
    aplot = ggplot(df.compare, aes(x = micro.p.value, y = qpcr.p.value))
    aplot = aplot + geom_point()
    atitle = cor(df.compare$qpcr.p.value, df.compare$micro.p.value)
    atitle = paste0("cor=", atitle)
    aplot = aplot + ggtitle(atitle)
    pdf(fp(outdir, "qpcr_micro_pvalue_scatter.pdf"))
    print(aplot)
    dev.off()
}



## gg_color_hue <- function(n) {
##   hues = seq(15, 375, length = n + 1)
##   hcl(h = hues, l = 65, c = 100)[1:n]
## }


get.dec.points <- function(acol)
{
    decpoints = (min(as.integer(unlist(strsplit(format(acol, scientific = T), "e"))))*-1)+1
    return(decpoints)
}

round.to.dec.point <- function(acol,pts=NULL)
{
    if(is.null(pts))
    {
        pts = get.dec.points(acol)
    }
    rounded = sprintf(acol, fmt = paste0("%.",pts,"f"))
    return(rounded)
}


writeSig <- function(df.m, outcome.type)
{
    cnames = c("mediator", "imprinted", "coef.ab", "p.value.ab",  "coef.cprime", "p.value.cprime", "suppressor")
    if(outcome.type =="behavior")
    {
        cnames = c("outcome", cnames)
        df.m$outcome = as.character(df.m$outcome)
        df.m$outcome[df.m$outcome == "Δ"] = "Δ Cort"
    }
    if(outcome.type == "micro")
    {
        df.m = df.m[! mediator %in% "Carmil1"]
    }
    z = (df.m[p.value.ab<.05,cnames, with = F])
    cpind = z[["p.value.cprime"]]==0
    

    z[["p.value.ab"]]     = formatC(z[["p.value.ab"]], digits=3)
    z[["p.value.cprime"]] = formatC(z[["p.value.cprime"]], digits =3)
    z[["coef.ab"]]     = formatC(z[["coef.ab"]], digits=3)
    z[["coef.cprime"]] = formatC(z[["coef.cprime"]], digits =3)
    z[["p.value.cprime"]][cpind] = paste0("<",formatC(1/(16000*.8), digits=3)) 
    
    
    if(outcome.type!="behavior")
    {
        setorder(z, p.value.ab)
    }
    print(z)
    fwrite(z, file = fp(outdir, paste0(outcome,"_sig.pvalues.csv")), sep = "\t")


    mytab = regulartable(data = z)
    mytab = bold(mytab, part = "header")
    mytab = align( mytab, align = "center", part = "all")
    rep = list()
    
    rep[["x"]] = mytab
    rep[["mediator"]]="Mediator Gene"
    rep[["imprinted"]]="Imprinted"
    rep[["coef.ab"]]="ab"
    rep[["p.value.ab"]]="CTP"
    rep[["coef.cprime"]]="c\'"
    rep[["p.value.cprime"]]="CTP"
    rep[["suppressor"]]="Suppressor"
    if(outcome.type=="behavior"){rep[["outcome"]] = "Behavior"}

    mytab = do.call(set_header_labels, rep)
    mytab = autofit(mytab, 0, 0)

    

    rep[["x"]] = mytab
    rep[["coef.ab"]]        = "Mediation Effect"
    rep[["p.value.ab"]]     = "Mediation Effect"
    rep[["coef.cprime"]]    = "Direct Effect"
    rep[["p.value.cprime"]] = "Direct Effect"
    
    mytab = do.call(add_header, rep)
    mytab = merge_h(mytab, part="header")
    mytab = merge_v(mytab, part="header")
    mytab <- italic(mytab, j = ~ mediator, italic = TRUE)
    mytab = autofit(mytab, 0, 0)

   
    newlen = dim(mytab)$widths["p.value.ab"]
    mytab = flextable::width(mytab, j = ~ coef.ab,     width = newlen*1.1)
    mytab = flextable::width(mytab, j = ~ coef.cprime, width = newlen*1.1)
        
    doc = read_docx(fp("./mnp/template.docx"))
    doc = body_add_flextable(doc, mytab)
    print(doc, target = fp(outdir, paste0(outcome,"_sig.pvalues.docx")))
}

getBehLevels <- function()
{
    nameMap = unname(getNameMap())
    return(nameMap)
}


##strain.results = fread(outm("limited_Strain_p_0.05.csv"))    ##fread(datm("2017-05_strainResults.csv"))
## beh = fread(datm("2016-05-02_behavior_pvals.csv"))
## ##limitedPhen = beh[strain.pval.qval.fdr =="o" | grepl(pattern = "\\*", strain.pval.qval.fdr)]$phen
## limitedPhen = beh[strain.pval =="o" | grepl(pattern = "\\*", strain.pval)]
## limitedPhen = paste0(limitedPhen$experiment, "_", limitedPhen$phenotype)
##limitedPhen = setdiff(limitedPhen, c("PC1", "PC2"))

outdir = fp(outm("mediation", "plots"))
raw.data = loadAllData$createAllInputs()
lrrc16a  = getProbesetId(raw.data$probesetInfo, "Lrrc16a")
airn     = "10441787"
taq.data = getRawTaqData(merge.qpcr.plate = T)
micro.data.orig  = mnp.med$get.sv.corrected.genes(raw.data, T)$exp.mat
##micro.data.orig  = mnp.med$get.sv.corrected.genes(raw.data, F)$exp.mat


##micro.data = data.frame(raw.data$exp.mat, check.names = F)

micro.data = micro.data.orig[ ,c(lrrc16a, airn),drop = F]
colnames(micro.data) = c("micro.obs", "airn.obs")
micro.data = data.table(ID = rownames(micro.data), micro.data, key = "ID")
cov.data    = raw.data$phens$breedLog

full = cov.data[taq.data, on = c(ID="qpcr.ID")]
full = full[micro.data, on="ID"]


##df.m = getFile("micro", "behavior", discardMissing, mergeQPCR)

##setorder(df.m, p.value.ab)

dir.create(outdir, showWarnings = F, recursive = T)
##plotEff()
plot.lrrc.airn()
##plot.p.comp()

##setorder(df.m, p.value.ab)

for(outcome in c( "behavior", "micro"))
{
    df.m = getFile("micro", outcome, discardMissing=T, mergeQPCR = F)
    df.m$outcome2 = df.m$outcome
    df.m$xoffset = 0
    df.m$yoffset = 0

    if(outcome == "behavior")
    {
        df.m$outcome = toAlias(df.m$outcome)
        df.m$outcome = factor(df.m$outcome, levels = getBehLevels())
    }

    if(outcome=="micro")
    {
        df.m$yoffset[df.m$mediator=="Pcdhb2"] = 1000
##        df.m$yoffset[df.m$mediator=="Mir485,Mirg"] = 1000
    }
    
    if(outcome == "behavior")
    {
        df.m$yoffset[df.m$outcome=="PPI78" & df.m$mediator=="Carmil1"] = 1500
        
        df.m$yoffset[df.m$outcome=="Pct Time Stranger" & df.m$mediator=="s113_10398354"] = 2700
##        df.m$yoffset[df.m$outcome=="Pct Time Stranger" & df.m$mediator=="s116_10564209"] = 800
        df.m$yoffset[df.m$outcome=="Pct Time Stranger" & df.m$mediator=="s116_10564209"] = 1200
        
        df.m$yoffset[df.m$outcome=="Basal Cort"        & df.m$mediator=="s115_10563915"] = 2700
        df.m$yoffset[df.m$outcome=="Basal Cort"       & df.m$mediator=="s113_10398354"] = 800
        
        df.m$yoffset[df.m$outcome=="10 Min Cort"       & df.m$mediator=="s115_10563949"] = 800


        df.m$yoffset[df.m$mediator=="3830406C13Rik"] = -5
    } else {
##        df.m$yoffset[df.m$mediator == "Irak1bp1"] = 150
    }

    writeSig(df.m, outcome.type = outcome)
    
    ##y = df.m[therank<10]
    df.m[,therank:=frank(p.value.ab), by = "outcome"]
   
    ##df.m$sigstrain = df.m$outcome %in% limitedPhen

    countpoint = 3
    top.med = df.m[,.(Imprinted = imprinted[1],
                      top3.count=sum(therank<=countpoint),
                      Suppressed = sum(suppressor),
                      
                      CTP = enrichmentTesting$fisherCombined(p.value.ab)$chisq.analytic.p),
                   by = "mediator"]
    ##setkey(top.med, "count")
    setorder(top.med, -top3.count)

    if(outcome=="behavior")
    {
        topwrite = top.med[CTP<.05]
        topwrite[["CTP"]] = as.character(formatC(topwrite[["CTP"]], digits =3))
        fwrite(topwrite, file = fp(outdir, paste0(outcome,"_top3.pvalues.csv")), sep = "\t")

        print("top3")
        mytab = regulartable(topwrite)
        mytab = bold(mytab, part = "header")
        mytab = align( mytab, align = "center", part = "all")

        
        rep = list()
        rep[["x"]] = mytab
        rep[["mediator"]]="Mediator Gene"
        rep[["Imprinted"]]="Imprinted"
        rep[["top3.count"]]="strongly (top-3) mediated"
        rep[["Suppressed"]]="suppressed"
        mytab = do.call(set_header_labels, rep)

        rep = list()
        rep[["x"]]         = mytab
        rep[["top"]]       = T
        rep[["mediator"]]  ="Mediator Gene"
        rep[["Imprinted"]] ="Imprinted"
        rep[["top3.count"]]="# Behavior POEs"
        rep[["Suppressed"]]="# Behavior POEs"
        rep[["CTP"]]       ="CTP"
        mytab = do.call(add_header, rep)

        mytab = merge_h(mytab, part="header")
        mytab = merge_v(mytab, part="header")
        mytab <- italic(mytab, j = ~ mediator, italic = TRUE)
        mytab = autofit(mytab, 0, 0)
        doc = read_docx(fp("./mnp/template.docx"))
        doc = body_add_flextable(doc, mytab)
        print(doc, target = fp(outdir, paste0(outcome,"_top3.pvalues.docx")))
    }
    
    ##df.m[,therank:=frank(p.value.ab, ties.method ="min"), by = "imprinted"]
    


    levs = as.character(levels(df.m$outcome))
    df.m$outcome = as.character(df.m$outcome)
    df.m$outcome[df.m$outcome == "Δ"] = "Δ Cort"
    levs[levs == "Δ"] = "Δ Cort"
    df.m$outcome = factor(df.m$outcome, levels = levs)
    
    aplot = ggplot(df.m, aes(x=-log10(p.value.ab), fill = imprinted))
    aplot = aplot + geom_histogram(alpha =.2, position = "identity", bins = 100)
    aplot = aplot + scale_y_sqrt()

    
    xlab.str = bquote("-log"[10] *"(CTP"[mediation(ab)]*")")

    aplot = aplot + xlab(xlab.str)
    aplot = aplot + geom_vline(xintercept = -log10(.05), linetype = "dashed")
    if(outcome=="behavior")
    {
    aplot = aplot + facet_wrap(~outcome, ncol =7)
##    aplot = aplot + theme(legend.position="top")
    }
    
    aplot = aplot + theme_bw()
    aplot = aplot + theme(panel.grid.minor = element_blank())
    aplot = aplot + theme(panel.grid.major = element_blank())
    aplot = aplot + theme(panel.background = element_rect(fill = "white"))
    aplot = aplot + theme(axis.line = element_line(colour = "black"))
    aplot = aplot + scale_color_manual(values = c("red", "blue"), guide=F)
    aplot = aplot + scale_fill_manual(values = c("red", "blue"),
                                      ## breaks =c("red","blue"),
                                      labels = c( "not imprinted", "imprinted"))

    hite = ifelse(outcome=="behavior", .04, .92)
    aplot = aplot  + theme(legend.position = c(.9, hite))
    aplot = aplot  + theme(legend.title=element_blank())
    aplot = aplot  + theme(legend.text = element_text(size = 20))

##    aplot = aplot + theme(strip.background = element_blank())

    if(outcome !="behavior")
    {
        aplot = aplot + theme(strip.text = element_text(size=20))
        aplot = aplot + theme(axis.text = element_text(size=20))
        aplot = aplot + theme(axis.title = element_text(size=20))

    }
    
    df.high = df.m[p.value.ab<.05 | (mediator %in% c("Airn_10441787", "Carmil1") & therank<=countpoint)]
    
    df.high$mediator = gsub(df.high$mediator, pattern = "_.*", replacement = "")

    sz = ifelse(outcome =="behavior", 3.8, 8)
    aplot = aplot + geom_text(fontface = "italic", hjust = 0, size = sz, data = df.high, aes( angle = 90, x = -log10(p.value.ab), y =5 + yoffset, label = mediator, color = imprinted))##paste0(mediator.id,"_", mediator_name)))

    

##    aplot = aplot + theme(strip.background = element_rect(fill = outcome2, colour = NA))
    

    ggsave(fp(outdir, paste0(outcome, "_facet.pdf")), width = 12, height = 8, device = cairo_pdf)
}

