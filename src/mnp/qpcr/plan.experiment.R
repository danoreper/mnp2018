library(data.table)
followup = fread("../data/mnp/MatNut_validation_samplelist_072516.csv")
## followup$original = paste0(followup$ID) %in% original$ID


original = loadAllData$createAllInputs()$cov.data.full
followup$original = paste0("Mouse.",followup$ID) %in% original$ID

browser()

b1   = followup[original==T&Batch==1]$ID
b2   = followup[original==T&Batch==2]$ID
b3   = followup[original==T&Batch==3]$ID
bnew = followup[original==F]$ID


plate2 = c(followup[ID %in% bnew & Diet =="LowPro"]$ID,
           b2,
           bnew,
           rep(followup[(ID %in% bnew) & Diet =="MethylSuff" & Strain=="B6.NOD"]$ID, 3))
           

plate1 = c(b1, b3, plate2[1:(80-length(b1)-length(b3))])

plate3 = c(plate2[(41-(80-length(plate2))):length(plate2)], plate1[1:40])

browser()
sink("./plan.txt")
cat("plate1\n")
cat(plate1)
cat("\n\nplate2\n")
cat(plate2)
cat("\n\nplate3\n")
cat(plate3)
sink()
