library(WVmisc)
source("./datacleaning.R")
source("./loadAllData.R")
dat <- read.csv(file.path(data, "allcovariates_2011_12_08.csv"))[,1:6]
dat <- trim.factor.levels(dat)
dat$Diet <- merge.factor.levels(lookup=list(stdctrl="StdCtrl"), dat$Diet)
dat$DamID <- paste(dat$Strain, sep=".", dat$DamID)

output = "../output"
dir.create(output, showWarnings = F)
#6 animals per (2) strain per diet (4) per pipeline (2), balanced over dams


rbalbin <- function(n, prob=0.5){
  # balanced binomial draw
  if (0==n) return (numeric(0))
  n1 <- ifelse(rbinom(1,prob=0.5, size=1)==1,
      ceiling(prob*n),
      floor(prob*n)
      )
  x <- c(rep(0,n-n1), rep(1, n1))
  sample(x, size=n, replace=FALSE)     
}

chooseMice <- function(subDat, u, nWant=6){
  # Aim: choose nWant mice for expression analysis in a way that is 
  # 1) randomized
  # 2) involves minimal repeat sampling of dams
  # 3) relatively balanced over batches
  # Approach:
  # a) identify all numComb non-empty unique combinations of Batch/DamID
  # b) sample one mouse randomly without replacement from each combination
  # c) repeat (b) until nWant have been chosen
  # Caveats/Assumptions:
  # i) the selection will typically be imbalanced unless numComb == nWant
  # ii) balancing over batches is prioritized equally with balancing over
  #    dams. However, if nWant < numComb, that may produce undesirable
  #    results whereby a new batch/used_dam combination is randomly
  #    randomly prioritized over a batch/new_dam combination.
  # iii) the order in which the combinations are selected from is randomized
  #    between (a) and (b), but, for simplicity, is constant between
  #    (b) and (c).
  
  # get unique combinations
  combs <- paste(subDat$DamID, subDat$Batch, sep=".")
  # make each combination a factor level, randomizing seniority
  facCombs <- factor(combs, levels=sample(unique(combs), replace=FALSE))
  # label combinations by seniority
  subDat$facNum <- as.integer(facCombs)
  # randomize order within factor levels (by simply randomizing generally)
  subDat <- subDat[ order(runif(nrow(subDat))), ]
  # reorder data so that factors occur in runs, eg, 1 2 3 4 1 2 3 4 etc
  # do this by making a factor string as "<#th occurance>_<fac numb>"
  facNumCount <- integer(nlevels(facCombs))
  subDat$facNumString <- NA
  for (i in 1:nrow(subDat)) {
    fn <- subDat$facNum[i]
    facNumCount[fn] <- facNumCount[fn] + 1 
    subDat$facNumString[i] <- paste(facNumCount[fn], sep="_", fn)
  }
  subDat <- subDat[order(subDat$facNumString),]
  subDat$Expression <- FALSE
  # choose first nWant
  subDat$Expression[1:nWant] <- TRUE
  # clean up
  subDat$facNum <- NULL
  subDat$facNumString <- NULL
  subDat
}

set.seed(1)
decided <- dfapply(dat, dat[,c("Pipeline", "Strain","Diet")], 
    FUN=chooseMice,
    results=NULL,
    results.add.FUN=rbind)
decided <- decided[order(!decided$Expression, decided$ID),]
write.csv(decided, file=file.path(output, "expression_choice_20120411.csv"), row.names=FALSE,quote=F)


sep=function(){cat("\n")}
catSec=function(...){cat("-------------------------------\n",...,"\n\n")}
sink(file.path(output, "expression_choice_stats.txt"))
ok <- decided$Expression
catSec("Choices by strain (balanced by design):")
table(decided[,c("Strain", "Expression")]); sep()
cat("\nAll mice:\n")
table(decided[,c("Strain", "Pipeline")]);
cat("\nExpression mice:\n")
table(decided[ok,c("Strain", "Pipeline")]);
catSec("Choices by diet (balanced by design):")
table(decided[,c("Diet", "Expression")]);  sep()
cat("\nAll mice:\n")
table(decided[,c("Diet", "Pipeline")])
cat("\nExpression mice:\n")
table(decided[ok,c("Diet", "Pipeline")])
catSec("Choices by batch (randomized in a balanced way):")
table(decided[,c("Batch", "Expression")]);  sep()
cat("\nAll mice:\n")
table(decided[,c("Batch", "Pipeline")])
cat("\nExpression mice:\n")
table(decided[ok,c("Batch", "Pipeline")])
cat("Choices by dam (randomized):")
table(decided[,c("DamID", "Expression")]);  sep()
cat("\nAll mice:\n")
table(decided[,c("DamID", "Pipeline")])
cat("\nExpression mice:\n")
table(decided[ok,c("DamID", "Pipeline")])
sink()



#====Earlier mistaken objective to assign a pipeline randomly====

cols <- c("Strain", "Diet", "Batch", "DamID")
out <- dfapply(dat, dat[,cols], 
    FUN=function(x, ...){
      rbalbin(nrow(x))
    },
    results=c())
dat$Pipeline <- unlist(out)+1

rdat <- dfapply(dat, dat[,cols], 
    FUN=function(x, ...){
      x$Pipeline <- 1+rbalbin(nrow(x))
      x
    },
    results=NULL,
    result.add.FUN=rbind)

rdat$Pipeline <- dfapply(dat, dat[,cols], 
    FUN=function(x, ...){1+rbalbin(nrow(x))},
    matched.vector=TRUE)


rand.summary <- dfapply(rdat, rdat[,cols], 
    FUN=function(x, u, ...){
      u$NumMice <- nrow(x)
      u$NumInPipeline1 <- sum(x$Pipeline==1)
      u$NumInPipeline2 <- sum(x$Pipeline==2)
      u
      },
    results = NULL,
    result.add.FUN = rbind
    )

table(rdat$Pipeline)
table(rdat$Strain, dat$Pipeline)
table(rdat$Diet, dat$Pipeline)
table(rdat$Batch, dat$Pipeline)
table(rdat$DamID, dat$Pipeline)


write.csv(dat, file=file.path(output, "randomizedpipeline_2012_12_21.csv"))
write.csv(rand.summary, file="randomizedpipeline_summary_2012_12_12.csv")

tapply(dat$RandPipeline, list(dat$Strain,dat$Diet), mean)

out <- by(dat, list(dat$Diet, dat$Strain, dat$Batch, dat$DamID),
  FUN = function(x){x$RandPipeline<-rbalbin(nrow(x)); x})
