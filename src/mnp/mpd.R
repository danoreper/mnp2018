##df = fread(dat("mpd/Chesler3.csv")) ##no locomotor activity assesed for B6 or NOD
##df = fread(dat("mpd/Crowley1.csv")) ## all male mice, only 10
##df = fread(dat("mpd/Donahue3.csv")) ## only b6 mice
##df = fread(dat("mpd/Gershenfield1.csv")) ##no NOD mice
##df = fread(dat("mpd/Golani1.csv")) ##No NOD mice
##df = fread(dat("mpd/Koide3.csv")) ##No NOD mice
##df = fread(dat("mpd/Metten1.csv")) ##No NOD mice
##df = fread(dat("mpd/Richfield1.csv")) ## No Nod mice
##df = fread(dat("mpd/Brown1.csv")) ##No Nod mice
##df = fread(dat("mpd/Koide1.csv")) ##No Nod mice
##df = fread(dat("mpd/Loos2.csv")) ##No distance data
##df = fread(dat("mpd/Palmer4.csv")) ##Hybrids, no b6xnod
##df = fread(dat("mpd/Schalkwyk1.csv")) ##No NOD
##SPijker2 No individual level animals
##df = fread(dat("mpd/Thomsen1.csv")) ##No NOD

getTarantino1 = function()
{
    df = fread(dat("mpd/Tarantino1.csv")) ##No NOD
    df = df[strain %in% c("C57BL/6J", "NOD/ShiLtJ")]# & !is.na(distance10_OFT)]
    df = df[,c("id", "strain", "sex", "batch", "dist_ctrl", "dist_cocaine")]
    setnames(df, "id", "mouse_ID")
    df = melt.data.table(df, id.vars = c("mouse_ID","strain", "sex", "batch"),
                         measure.vars = c("dist_ctrl",
                                          "dist_cocaine"),
                         variable.name = "varname")
    
    summary(lm(data = df[varname == "dist_cocaine"], formula = value~batch + strain))
    return(df)
}






library(data.table)
library(flextable)
library(officer)

source("./lm/fitBoxCoxModels.R")
source("./loadParams.R")
source("./mnp/general.R")


df0 = fread(dat("mpd/Pletcher1.csv")) ##No Female mice
df0 = df0[strain %in% c("C57BL/6J", "NOD/ShiLtJ")]# & !is.na(distance10_OFT)]
df0 = df0[,c("id", "strain", "sex", "OFT_activity")]
df0$varname = "OFT_activity"
df0$value = df0[["OFT_activity"]]
df0$OFT_activity = NULL


df1 = fread(dat("mpd/Chesler4.csv"),skip = 6)
df1 = df1[strain %in% c("C57BL/6J", "NOD/ShiLtJ") & OF_distance_first10!="="]
df1 = melt.data.table(df1, id.vars = c("mouse_ID","strain", "sex"),
                      measure.vars = c("OF_distance_first10",
                                       "LD_distance_light"),
                      variable.name = "varname")

df2 = fread(dat("mpd/Wiltshire1.csv"))
df2 = df2[strain %in% c("C57BL/6J", "NOD/ShiLtJ")]
df2$projsym = NULL
df2$stocknum = NULL
df2$zscore   = NULL
df2$measnum  = NULL
df2$strainid = NULL
setnames(df2, "animal_id", "mouse_ID")
setcolorder(df2, colnames(df1))


toeval = 
list(Chesler4   = list(df = df1, phens = c("OF_distance_first10", "LD_distance_light")),
     Wiltshire1 = list(df = df2, phens = c("act_OFT")),
     Pletcher1  = list(df = df0, phens = c("OFT_activity")),
     Tarantino1 = list(df = getTarantino1(), phens = c("dist_ctrl", "dist_cocaine")))


tabls = list()
i = 1
for(dfname in names(toeval))
{
    dataset = toeval[[dfname]]
    df = dataset$df
    for(phen in dataset$phens)
    {
##        browser()
        df.phen = df[varname == phen]
        
        sexes = unique(df.phen$sex)

        if(length(sexes)==1)
        {
            sexiters = list(sexes)
        } else {
            sexiters = list(c("m","f"), "m", "f")
        }
        
        for(j in 1:length(sexiters))
        {
            thesex = sexiters[[j]]
            covariateModelString = ifelse(paste(thesex, collapse = "")=="mf", "~ sex + strain", "~ strain")
            df.phen.sex = df.phen[sex %in% thesex]
            
            
            out = fit.model.bc$fit(y.mat = as.numeric(df.phen.sex$value),
                                   cov.data = df.phen.sex,
                                   covariateModelString = covariateModelString)$phen_1
            
            df.phen.sex$transformed.value = out$y.transformed
            fname = paste(dfname, phen, paste(thesex, collapse=""), "txt", sep=".")
            fname = outm("mpd", fname)
            print(fname)
            fwrite(df.phen.sex, fname)
            
            anov = data.frame(out$anovaWrapper$an)
            ##browser()
            pval = anov["strain", "Pr..F."]
            acoef = coef(out$fit)["strainNOD/ShiLtJ"]
            
            tabl = data.frame(dataset = dfname,
                              phenotype = phen,
                              sex = paste(thesex, collapse = ","),
                              n       = nrow(df.phen.sex),
                              NOD.further  = ifelse(acoef>0, "NOD", "B6"),
                              anova.pval = pval)
            tabls[[i]] = tabl
            print(i)
            i = i + 1
        }
    }
}
tabls = rbindlist(tabls)

strs = util$stars.pval(tabls$anova.pval, cutpoints = c(0, .001, .01, .05, .2,1))
tabls$anova.pval =  as.character(signif(tabls$anova.pval, digits = 2))
tabls$anova.pval = paste(tabls$anova.pval, strs, sep="")

library(flextable)
mytab = regulartable(data = tabls)
mytab = theme_box(mytab)
mytab = bold(mytab, part = "header")
#mytab = bold(mytab,    j = ~ pipeline)  
mytab = align( mytab, align = "center", part = "all")
rep = list()
rep[["x"]] = mytab
rep[["dataset"]]="Dataset"
rep[["sex"]]="Sexes"
rep[["n"]] = "n"
rep[["phenotype"]] = "Phenotype"
rep[["anova.pval"]]="p-value"
rep[["NOD.further"]] = "Larger Distance Traveled"
mytab = do.call(set_header_labels, rep)
mytab = autofit(mytab, 0, 0)

wid = max(strwidth(tabls[["anova.pval"]], font = 10, units = 'in'))
mytab = flextable::width(mytab, j = as.formula("~anova.pval"),  width = wid)

tmplate = fp("./mnp/template2.docx")
doc = read_docx(tmplate)
##browser()
doc = body_add_flextable(doc, mytab)
target = outm("mpd.docx")


print(doc, target)





