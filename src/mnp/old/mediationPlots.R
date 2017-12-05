library("data.table")
source("./enrichmentTesting.R")



outputs = fread("../output.mediation/mediation/behavior_TRUE_TRUE/mediationAll.csv")
outputs$suppressor = outputs$coef.ab*outputs$coef.c < 0
tokeep = c("outcome.id", "mediator_name","mediator.id",  "imprinted", "suppressor", "coef.c",  "coef.ab", "coef.a", "coef.b", "p.value.c", "p.value.ab")
x  = outputs[,tokeep, with = F]
##x = fread("~/Desktop/mediation.csv")
## x = x[moderators == "Diet=Ave"]
## x$moderators = NULL
x[,therank:=frank(p.value.ab), by = "outcome.id"]
setorder(x, p.value.ab)


##y = x[therank<10]
top5.med = x[,.(count=sum(therank<5),
                suppressor = sum(suppressor),
                enhancer   = sum(!suppressor),
                 p.fisher = enrichmentTesting$fisherCombined(p.value.ab)$chisq.analytic.p,
                 mediator_name = mediator_name[1]), by = "mediator.id"]
setorder(top5.med, -count)
print(top5.med[p.fisher<.05])
x$therank = NULL
x$mediator.id = NULL
print(x[p.value.ab<.05])


## beh = fread("~/Desktop/2016-05-02_behavior_pvals.csv")
## limitedPhen = beh[strain.pval.qval.fdr =="o" | grepl(pattern = "\\*", strain.pval.qval.fdr)]$phen

## y = x[outcome.id %in% limitedPhen]
## print(y[p.value<.05])
## top10.med = y[,
##               .(count=sum(therank<10),
##                 p.fisher = enrichmentTesting$fisherCombined(p.value)$chisq.analytic.p,
##                 mediator_name = mediator_name[1]), by = "mediator.id"]

## setorder(top10.med, p.fisher)
## print(top10.med)


load("~/Desktop/out.big2.RData")
df.all = output2$mediation
df.all[,imprinted := !is.na(minDistToImprinted)&minDistToImprinted<5000]


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for(moderator in unique(df.all$moderators))
{
    print(moderator)
    df = df.all[moderators == moderator]
    df[,therank:=frank(p.value, ties.method ="min"), by = "imprinted"]
    
    aplot = ggplot(df, aes(x=-log10(p.value), fill = imprinted))
    aplot = aplot + ggtitle(moderator)
    aplot = aplot + geom_histogram(alpha =.2, position = "identity", bins = 100)
    aplot = aplot + scale_y_sqrt()
    aplot = aplot + geom_vline(xintercept = -log10(.05))


    if(moderator == "Diet=Ave")
    {
        df.high = df[p.value<.05|mediator_name %in% c("Airn")]
        aplot = aplot + geom_text(hjust = 0, size =3, data = df.high, aes(angle = 90, x = -log10(p.value), y =5, label = mediator_name, color = imprinted), show.legend = F)##paste0(mediator.id,"_", mediator_name)))

    } else {
        df.high = df[mediator_name =="Airn"|(p.value<.05 & therank<5 & imprinted)]
        aplot = aplot + geom_text(hjust = 0, size =3, data = df.high, aes(angle = 90, x = -log10(p.value), y =5, label = mediator_name), color = gg_color_hue(2)[2], show.legend = F)
    }
    

    fle = outm("mediation.p.values", paste0(moderator, ".mediation.pdf"))
    print(fle)
    pdf(file = fle, width = 7, height = 3)
    print(aplot)
    dev.off()
}

