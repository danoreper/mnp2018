#The expression processing of microarray data before I refactored all the calls to lmer, as well as the inclusino of multicore, and 
#properly computing the box cox transform

library(WVmisc)
library(nlme)
library(MASS)
library(data.table)
library(lme4)
library(lattice)
library(reshape)
library(ggplot2)
#library(lmerTest) #caTools, bitOps, and possibly others need to be specifically installed. really annoying.
options(warn=1)
#setwd("../../mnp/WillPilotAnalysis/src/")

data                      = "../data"
output                    = "../output/mnp/output_newkill"
surrogat                  = "SSVA"
lambdasToTry              = 1#c(-2,-1,-.5,0,.5,1,2)
mergeDiets                = F
mc.cores                  = 4
dir.create(output, showWarnings = F)

allPhenNames = c("OpenField.PctCenterTime","OpenField.TotalDistance",
		"LightDark.TotalDistance","LightDark.TimeInLight",
		"LightDark.TotalTransitions")

get.annot.data = function(annot.old.expression.file,
		maskedraw.expression.file,
		probeset.info.file,
		output)
{
	annot.data = read.annot.data(annot.old.expression.file, outdir = output)
	
	
# remove positive and negative controls
	
	ii <- grep("^(pos_control|neg_control)", perl=TRUE, annot.data$mRNA.Description)
	
	plotControls(annot.data = annot.data, outdir = outdir)
	
	print(length(ii)) # 6546
	annot.data <- annot.data[ -ii, ]
	
# remove non-annotated
	ii <- which(is.na(annot.data$mRNA.Accession))
	badIds = annot.data$Probe.Set.ID[ii]
	
	probeset.info <- read.csv(probeset.info.file)
	probeset.info$Num.Pure.Probes <- probeset.info$Num.Probes - probeset.info$Num.Probes.With.SNP
	
#	probes = read.csv(file = file.path("../output/mnp/output", "probe.info"))
#	probes = data.table(probes, key="Probe.Set.ID")
	
#	library("biomaRt")
#	ensembl=useMart("ensembl") 
#	listDatasets(ensembl)
#	ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
#	filters = listFilters(ensembl)
	
#	
	#Figure out which locations have annotations in biomart
	counter = 0
	recoveredIs = NULL
	for(ind in 1:length(ii))
	{
		badId = badIds[ind]
		single.i = ii[ind]
#		subprobes = probes[J(badId),]
#		print(badId)
		
#		results <- try(
#				getBM(attributes = c("ensembl_gene_id","refseq_mrna", "refseq_ncrna", "mgi_symbol", "chromosome_name", "start_position", "end_position"),
#					  filters = c("chromosome_name", "start", "end"), 
#					  values = list(as.character(subprobes$chr), subprobes$start, subprobes$stop), mart = ensembl))
#		#TODO put this back
#	    if(!caught.error(results))
#		{
#			annot.data[annot.data,"mgi_symbol"]
#			print(results)
#			if(nrow(results)>0)
#			{
#				recoveredIs = c(recoveredIs,single.i)
#				counter = counter +1
#			}
#		}
	}
	ii = setdiff(ii, recoveredIs)
	#at this point, ii contains indexes for which there are no annotations in existence.
	
	length(ii) # 570
	annot.data <- annot.data[ -ii, ]
	print(paste0("remaining: ", nrow(annot.data))) # 28440
	
	annot.data = replaceAnnotExpressionWithMasked(annot.data = annot.data, 
			maskedraw.expression.file = maskedraw.expression.file)
	
	mouse.cols <- grep("Mouse", colnames(annot.data))
	
	mis.probes <- which( apply(as.matrix(annot.data[ , mouse.cols]), 1, function(x){ any(is.na(x)) } ))
	sum(is.na(annot.data[ mis.probes, mouse.cols ])) # 0
# remove single probe with missingness
#annot.data <- annot.data[ -mis.probes, ]
	nrow(annot.data) # 25196 (alan); 28440 (Will)
	
	
# filter for probe sets with duplicated expression
# (ie, probe sets mapping to the same gene that have *exactly* the same expression level 
	duplicates = duplicated(cbind(annot.data[,mouse.cols], annot.data$mRNA.Accession))
	annot.data = annot.data[!duplicates,]
	nrow(annot.data) # 24867 (alan); 28006 (3-21 or 0-24); 
	
	
#probeset.info = data.table(probeset.info,key="Probe.Set.ID")
	
	ii <- match(annot.data$Probe.Set.ID, probeset.info$Probe.Set.ID)
	annot.data$Num.Pure.Probes <- probeset.info$Num.Pure.Probes[ii]
	annot.data$Num.SNP.Probes <- probeset.info$Num.Probes.With.SNP[ii]
	return(annot.data)
}


read.annot.data=function(annot.old.expression.file, outdir)
{
	annot.data <- read.delim(annot.old.expression.file, na.string="---", strip.white=TRUE)
	
	############# Processing annotated but unmasked expression data###
# clean up irregular column names
	colnames(annot.data) <- map.eq(colnames(annot.data),
			list(
					"mRNA..Source"       = "mRNA.Source",
					"mRna...Description" = "mRNA.Description",
					"mRNA...xhyb"        = "mRNA.xhyb"))
	
# make mice "Mouse.N" where N is mouse id
	mouse.cols <- grep("\\.rma.gene\\.default\\.Signal", colnames(annot.data))
	colnames(annot.data)[mouse.cols] <- 
			sub("X", "Mouse.", sub("\\.rma.gene\\.default\\.Signal", "", colnames(annot.data)[mouse.cols]))
	
	return(annot.data)
} 

plotControls <- function(annot.data, outdir) 
{
	mousecols = which(grepl("Mouse", colnames(annot.data)))
	for(controlString in c("(neg_control)","(pos_control)" ))
	{
		ii <- grep(pattern=controlString, perl=TRUE, annot.data$mRNA.Description)
		exp = annot.data[ii, mousecols]
		mouses = as.factor(colnames(exp))
#	fit = lm(exp2~mouses+1)
		jj = rank(rowSums(as.matrix(exp)))
		genebucket = cut(rowSums(as.matrix(exp))/ncol(exp), 100)
		out = c(as.matrix(exp))
		gene   = rep(jj,length(mouses))
		mouseInd = rep(mouses, each = nrow(exp))
		genebucket = rep(genebucket, length(mouses))
		df = data.frame(gene = gene, out=out, genebucket=genebucket, mouseInd=mouseInd)
		df = data.table(df)
		df2 = df[,list(mnexp=mean(out)),by=c("mouseInd", "genebucket")]
		
		pdf(file.path(outdir, paste0(controlString,".pdf")))
		#xyplot((out~x|mouseInd), data=df[1:(95*5522),])
		xyplot((mnexp~genebucket|mouseInd), data=df2[,])
		dev.off()
	}
}

replaceAnnotExpressionWithMasked <- function(annot.data, maskedraw.expression.file)
{
	mouse.cols <- grep("Mouse\\.[0-9]", colnames(annot.data))
	annot.data = annot.data[,setdiff(1:ncol(annot.data), mouse.cols)]
	annot.data = data.table(annot.data,key="Probe.Set.ID")
	
	mask.data <- read.delim(maskedraw.expression.file, na.string="---", strip.white=TRUE, comment.char="#")
	colnames(mask.data) <- sub("X", "Mouse.", sub("\\.CEL", "", colnames(mask.data)))
	mask.data = data.table(mask.data,key="probeset_id")
	setnames(mask.data, "probeset_id", "Probe.Set.ID")
	annot.data = merge(mask.data, annot.data)
	annot.data = data.frame(annot.data)
	return(annot.data)
}

form.Cov.Data <- function(exp.mat, design.file)
{
	exp.mice <- rownames(exp.mat)
# covariate data
	design.data <- read.csv(design.file)
# Prepare mouse-specific covariates
	design.data$ID <- paste("Mouse", sep=".", as.character(design.data$ID))
	print(setdiff(exp.mice, design.data$ID)) # identifies an accidental technical replicate 142
	
# allow for accidental technical replicate in design file
	design.data <- rbind(design.data[ design.data$ID=="Mouse.142", ], design.data)
	design.data$ID[1] <- "Mouse.142_2"
	
# make covariate file from design data
	cov.data <- design.data[ match(exp.mice, design.data$ID), ]
	
	if(mergeDiets)
	{
		cov.data$Diet = as.character(cov.data$Diet)
		cov.data$Diet[cov.data$Diet %in% c("StdCtrl", "MethylSuff")] = "StdCtrl"
		cov.data$Diet[cov.data$Diet != "StdCtrl"] = "Deficient"
		cov.data$Diet= as.factor(cov.data$Diet) 
	}
	print(all(cov.data$ID==exp.mice))
	
	return(cov.data)
}

###########################
# Make cluster dendrogram
makeClusterDendogram <- function(exp.mat, cov.data, outdir, maskedraw.expression.file) 
{
	mouse.dist  <- as.dist(1 - abs(cor(t(exp.mat))))
	mouse.clust <- hclust(mouse.dist, method="average")
	
	mother.col  <- rainbow(nlevels(cov.data$DamID))
	
	pdf(file.path(outdir, "expression_dendrogram.pdf"), width=14, height=5)
	par(mar=c(left=6, right=0.5, bottom=5, top=3))
#	par(mar=sides(left=6, right=0.5, bottom=5, top=3))
	plot(mouse.clust, xlab="", las=1, labels=sub("Mouse.", "", mouse.clust$labels),
			main=paste(maskedraw.expression.file, "\n", nrow(exp.mat), "x", ncol(exp.mat)))
	
	mtext(side=1, at=0, adj=1, line=0, "Cross (col=dam)")
	mtext(side=1, at=1:nrow(cov.data), line=0, c("B","N")[ cov.data$Strain[mouse.clust$order]],
			col = mother.col[ cov.data$DamID ] )
	
	mtext(side=1, at=0, adj=1, line=1, "Diet")
	mtext(side=1, at=1:nrow(cov.data), line=1, c("L","M","S","V")[ cov.data$Diet[ mouse.clust$order ] ])
	
	mtext(side=1, at=0, adj=1, line=2, "Pipeline")
	mtext(side=1, at=1:nrow(cov.data), line=2, c(1,2)[ cov.data$Pipeline[ mouse.clust$order ] ])
	
	mtext(side=1, at=0, adj=1, line=3, "Batch")
	mtext(side=1, at=1:nrow(cov.data), line=3, c(1,2,3)[ cov.data$Batch[ mouse.clust$order ] ])
	dev.off()
}

draw.hists <- function(d, prefix="", suffix="", logP=FALSE, hits.at.fdr=NULL)
{
	effects <- colnames(d)[-1]
	for (effect in effects) 
	{
		y <- d[, effect]
		if (logP) { y <- 10^-y }
		effect.suffix <- suffix
		if (logP & !is.null(hits.at.fdr))
		{
			qvals <- p.adjust(y, method="fdr")
			effect <- sub("LogP.", "", effect)
			effect.suffix <- paste0(suffix, 
					"\n(", 
					sum(qvals <= hits.at.fdr, na.rm=TRUE), 
					"/", nrow(d), 
					" at FDR<=", hits.at.fdr, ")")
		}
		if (!logP) {
			effect.suffix <- paste0(suffix, "\n(", 
					signif(min(na.rm=TRUE, y), 3), ", ",
					signif(max(na.rm=TRUE, y), 3), ")")
		}
		if(length(unique(y))==1 && is.na(unique(y)))
		{
			next
		} else {
			(truehist(y, main=paste(sep="", prefix, effect, effect.suffix), xlab=effect))
		}
		
	}
}

#return the scaled residuals: step 1.1
getResids <- function(exp.mat, cov.data)
{
#	lambdas =c(-2,-1,-.5,0,.5,1,2)
	lambdas =c(1)
	
	exp.resid.mat <- matrix(NA, nrow=nrow(exp.mat), ncol=ncol(exp.mat)) 
#why do we represent cross this way?
	cross <- ifelse(as.character(cov.data$Strain)=="B6.NOD", -0.5, 0.5)
	
	modelPipeline = length(unique(cov.data$Pipeline))>1
	bestLambdas = rep(NA, ncol(exp.mat))
	pheno ="y"
	for (jj in 1:ncol(exp.mat))
	{
		if (jj %% 100 == 0) 
		{ cat(sep="","[",jj,"]") } 
		
		y <- exp.mat[ , jj]
		cov.data$y=y
		
		fit   <- 
				try( 
						{
							if(modelPipeline)
							{
#								
								covariates                                          = "~ 1 + as.factor(Batch) + Pipeline + Diet + Strain +  Strain:Diet + (1|DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data) 
								bestLambdas[jj] = bestLambda
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data =cov.data, as.formula(paste0(lhs,            "~ 1 + as.factor(Batch) + Pipeline + Diet + Strain +  Strain:Diet")), random = ~ 1 | DamID)
#								      lme(y ~ 1 + as.factor(Batch) + Pipeline + Diet + Strain + Strain:Diet, data=cov.data, random = ~ 1 | DamID)
							} else {
								
								
								covariates                                          = "~ 1 + as.factor(Batch)             + Diet + Strain +  Strain:Diet + (1|DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data)
								bestLambdas[jj] = bestLambda
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data=cov.data, as.formula(paste0(lhs,            "~ 1 + as.factor(Batch)             + Diet + Strain +  Strain:Diet")), random = ~ 1 | DamID)
								
#								      lme(y ~ 1 + as.factor(Batch)            + Diet + Strain + Strain:Diet, data=cov.data, random = ~ 1 | DamID)
								#strain is the same as cross
							}
						}
				)
		if (caught.error(fit))
		{
			
			fit = try(
					{
						if(modelPipeline)
						{
							covariates                                          = "~ 1 + as.factor(Batch) + Pipeline + Diet + Strain + Strain:Diet+ (1|DamID)"
							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data)
							bestLambdas[jj] = bestLambda
							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, "~ 1 + as.factor(Batch) + Pipeline + Diet + Strain + Strain:Diet+ (1|DamID)")))
							
#							lmerTest::lmer(y ~ 1 + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet+ (1|DamID) , data=cov.data)
							
							
							
						} else {
							
							covariates                                          = "~ 1 + as.factor(Batch) +           + Diet + Strain + Strain:Diet+ (1|DamID)"
							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data)
							bestLambdas[jj] = bestLambda
							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, "~ 1 + as.factor(Batch) +           + Diet + cross + cross:Diet+ (1|DamID)")))
							
#							lmerTest::lmer(y ~ 1 + as.factor(Batch)            + Diet + cross + cross:Diet+ (1|DamID) , data=cov.data)
						}
					}
			)
			
			if (caught.error(fit))
			{
				next
			}
		}
		
		
		exp.resid.mat[ , jj] <- resid(fit)
	}
	
	date()
	na.cols <- which(apply(exp.resid.mat, 2, function(x){all(is.na(x))}))
	if(length(na.cols)>0)
	{
		exp.resid.mat <- (exp.resid.mat[ , -na.cols])
	}
	exp.resid.mat <- scale(exp.resid.mat) # columns should be 0-centered anyway, so this is optional
	
	return(list(exp.resid.mat=exp.resid.mat, bestLambdas=bestLambdas))
}

refine.sv = function(exp.mat, exp.resid.mat)
{
	#Step 1.2
	exp.resid.pca <- prcomp(exp.resid.mat)
	#Modified form of step 1.3-1.10. N.B. we are fitting a mixed model, so permutation doesnt make
	#sense as a way to pick the number of principal components. Instead, we just eyeball the 
	#eigenvalues and see where the inflection point is.
	plot(exp.resid.pca$sdev^2, pch=as.character(c(1:9,0))) # suggests 7 pcs
	num.mice <- nrow(exp.resid.mat)
	num.sv <- 7 # number of surrogate variables
	sv.mat <- matrix(nrow=num.mice, ncol=num.sv) # holds refined sv's
	for (ii in 1:num.sv) #Step 2.3 is the result of 1.0-1.10 
	{
		#step 2.2 (step 2.1 is cached from 1.1)
		rough.sv      <- exp.resid.pca$x[,ii]
		#getting the highest correlated genes to the eigengene seems to be equivalent,
		#but just in case, following the exact text and computing a p value for every linear model.
#		gene.cors     <- apply(exp.mat, 2, function(y){ abs(cor(rough.sv, y)) })
		#step 2.4
		pval = rep(0, ncol(exp.mat))
		for(m in 1:ncol(exp.mat))
		{
			mthGeneExpresion =  exp.mat[,m]
			fit = lm(rough.sv ~ mthGeneExpresion )
			pval[m] = coefficients(summary(fit))["mthGeneExpresion","Pr(>|t|)"]
		}
		#modified form of step 2.5 TODO look into this
		top.thou.jj   <- order(pval, decreasing=F)[1:1000]
#		top.thou.jj   <- order(gene.cors, decreasing=TRUE)[1:1000]
#		
		#step 2.6
		refined.pca   <- prcomp(exp.mat[ , top.thou.jj])
		
		#we want to find the eigengene of the reduced gene set, that is most correlated with 
		# the original eigengene
		#step 2.7 Note that we use the absolute value of correlation, since flipping a principal component positive or nega5tive
		#should account for the same amount of variance, it will just change the regression coefficient 
		#for that component. i.e. corr = -.95 is a better pc than corr = .02 
		pca.cors      <- apply(refined.pca$x, 2, function(y){ abs(cor(rough.sv, y)) })
		indx          = which.max(pca.cors)
		print(indx)
#		sv.mat[ , ii] <- refined.pca$x[ , 1]
		sv.mat[ , ii] <- refined.pca$x[ , indx]
	}
	
	return(sv.mat)
}


formulaWrapper.boxCoxString = function(phenotype, bestLambda)
{
	if(bestLambda == 0)
	{
		return(paste0("log(", phenotype, ")"))
	}
	else
	{
		return (paste0("((((", phenotype, "^", bestLambda, ") -1))/",bestLambda,")"));
#		return (paste0("(",  phenotype, "^", bestLambda, ")"));
	}
}

formulaWrapper.yjString = function(phenotype, bestLambda)
{
	ifelsestring = paste("(", phenotype, ">=0)")
#	
	
	yjPhen = paste("(", phenotype, " + 1 )", sep="");
	yjLambda = bestLambda;
	fullString1 = paste(ifelsestring, "*",formulaWrapper.boxCoxString(yjPhen, yjLambda), sep="");	
	
	
	yjPhen = paste("( abs(", phenotype, ") + 1 )", sep="");
	yjLambda = 2 - bestLambda;
	fullString2 = paste(paste("(1 - ", ifelsestring, ")", sep=""), "*",
			formulaWrapper.boxCoxString(yjPhen, yjLambda), sep="");
	
	combinedString = paste(fullString1, " + ", fullString2);
#	print(paste("combined string is:", combinedString));
	return(combinedString);
}

getFormla <- function(lambda, pheno, covariates)
{
#	if(lambda==0)
#	{
#		formla = as.formula(paste0("log(", pheno, ")", covariates))
#	} else {
#		formla = as.formula(paste0("((", pheno, "^", lambda, ")", "-1)/", lambda, covariates))
#	}
	
#	formla = (paste0(formulaWrapper.boxCoxString(pheno, lambda), covariates))
	formla = (paste0(formulaWrapper.yjString(pheno, lambda), covariates))
	return(formla)
}


getBestLambda <- function(lambdas, pheno, covariates, dataSet, sv.mat) 
{
	pvals = rep(NA, length(lambdas))
	for(l in 1:length(lambdas))
	{
		lambda = lambdas [l]
		formla = as.formula(getFormla(lambda = lambda, pheno = pheno, covariates = covariates))
#		print(formla)
#		print(formla)
		fit.with.interaction = NULL
		fit.with.interaction = try(lme4::lmer(formla, data=dataSet))
		if(class(fit.with.interaction)=="try-error")
		{
			next
		}
#		print(fit.with.interaction)
		
#		
		resid.fit = NULL
		resid.fit            = try(resid(fit.with.interaction))
		if(class(resid.fit)=="try-error")
		{
			
		}
#		
		pval                 = NULL
		pval                 = try(shapiro.test(resid.fit)$p.value)
		if(class(pval)=="try-error")
		{
			
		}
		
		pvals[l] = pval
	}
	selectedLambda  = lambdas[which.max(pvals)]
#	
	print(selectedLambda)
	return(selectedLambda)
}


# Test for effects of cross at any diet
#Step 2.8
crossEffectsAtDiet <- function(cov.data, exp.mat, sv.mat, annot.data, bestLambdas) 
{
	cross <- ifelse(as.character(cov.data$Strain)=="B6.NOD", -0.5, 0.5)	
	results.pipe <- NULL
	cov.data$Diet <- relevel(cov.data$Diet, ref="StdCtrl")
	cov.data$cross = cross
	
	modelPipeline = length(unique(cov.data$Pipeline))>1
	date()
	
#	lambdas = c(-2,-1,-.5,0,1,.5,1)
	
	pheno = "y"
	
	for (jj in 1:ncol(exp.mat)) 
	{
		lambdas = bestLambdas[jj]
#		if (jj %% 100 == 0) { cat(sep="","[",jj,"]") } 
		
		
		y <- exp.mat[ , jj]
		cov.data$y = y
		
		
		fit <- 
				try({
							if(modelPipeline)
							{
#				
								covariates =                              " ~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet +( 1 | DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet")), random = ~ 1 | DamID)
							} else {
								
								covariates =                               "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + cross:Diet +( 1 | DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + cross:Diet")), random = ~ 1 | DamID)
							}
						})
		if (caught.error(fit)) 
		{
			print(paste0("recovering fit:", jj)) 
			fit <-
					try({
								if(modelPipeline)
								{
									covariates =                                          "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet +( 1 | DamID)"
									bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
									lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
									lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
								} else {
									covariates =                                          "~ 1 + sv.mat + as.factor(Batch)            + Diet + cross + cross:Diet +( 1 | DamID)"
									bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
									lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
									
									lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
								}
							})
		}
		
		
		
		fit.main.only <- 
				try({
							if(modelPipeline)
							{
#				
								covariates =                              " ~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + ( 1 | DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross ")), random = ~ 1 | DamID)
							} else {
								
								covariates =                               "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + ( 1 | DamID)"
								bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
								lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
								lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross ")), random = ~ 1 | DamID)
							}
						})
		if (caught.error(fit.main.only)) 
		{
			print(paste0("recovering main fit:", jj)) 
			fit.main.only <-
					try({
								if(modelPipeline)
								{
									covariates =                                          "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross +( 1 | DamID)"
									bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
									lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
									lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
								} else {
									covariates =                                          "~ 1 + sv.mat + as.factor(Batch)            + Diet + cross +( 1 | DamID)"
									bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
									lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
									
									lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
								}
							})
		}
		
		
		if (caught.error(fit) | caught.error(fit.main.only)) 
		{
			next
		}
		
		# parse results
		res <- try(CalcSignif(fit, fit.main.only))
		
		if (!caught.error(res)) 
		{
			if (is.null(results.pipe)) 
			{ # initiate results data frame
				tmp.mat <- matrix(nrow=ncol(exp.mat), ncol=ncol(res)+1)
				colnames(tmp.mat) <- c("Probe.Set.ID", colnames(res))
				results.pipe <- as.data.frame(tmp.mat)
				results.pipe$Probe.Set.ID <- annot.data$Probe.Set.ID
			}
			results.pipe[ jj, 2:(ncol(res)+1) ] <- res
		}
	} # about 1hr
	date()
	results <- results.pipe
	return(results)
}

############ Test for effects of cross at any diet
##Step 2.8
#crossEffectsAtDiet <- function(cov.data, exp.mat, sv.mat, annot.data) 
#{
#	cross <- ifelse(as.character(cov.data$Strain)=="B6.NOD", -0.5, 0.5)	
#	results.pipe <- NULL
#	cov.data$Diet <- relevel(cov.data$Diet, ref="StdCtrl")
#	cov.data$cross = cross
#	
#	modelPipeline = length(unique(cov.data$Pipeline))>1
#	date()
#	
##	lambdas = c(-2,-1,-.5,0,1,.5,1)
#	lambdas = c(1)
#	pheno = "y"
#	
#	for (jj in 1:ncol(exp.mat)) 
#	{
#
##		if (jj %% 100 == 0) { cat(sep="","[",jj,"]") } 
#		
#		
#		y <- exp.mat[ , jj]
#		cov.data$y = y
# 
#		
#		fit <- 
#		try({
#		    if(modelPipeline)
#			{
##				
#				covariates =                              " ~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet +( 1 | DamID)"
#				bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#				lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#				lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet")), random = ~ 1 | DamID)
#			} else {
#				
#				covariates =                               "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + cross:Diet +( 1 | DamID)"
#				bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#				lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#				lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + cross:Diet")), random = ~ 1 | DamID)
#			}
#		})
#		if (caught.error(fit)) 
#		{
#			print(paste0("recovering fit:", jj)) 
#			fit <-
#			try({
#						if(modelPipeline)
#						{
#							covariates =                                          "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + cross:Diet +( 1 | DamID)"
#							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
#						} else {
#							covariates =                                          "~ 1 + sv.mat + as.factor(Batch)            + Diet + cross + cross:Diet +( 1 | DamID)"
#							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#							
#							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
#						}
#				})
#		}
#		
#		
#		
#		fit.main.only <- 
#		try({
#					if(modelPipeline)
#					{
##				
#						covariates =                              " ~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross + ( 1 | DamID)"
#						bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#						lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#						lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross ")), random = ~ 1 | DamID)
#					} else {
#						
#						covariates =                               "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross + ( 1 | DamID)"
#						bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#						lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#						lme (data=cov.data, as.formula(paste0(lhs, "~ 1 + sv.mat + as.factor(Batch)             + Diet + cross ")), random = ~ 1 | DamID)
#					}
#		})
#		if (caught.error(fit.main.only)) 
#		{
#			print(paste0("recovering main fit:", jj)) 
#			fit.main.only <-
#			try({
#						if(modelPipeline)
#						{
#							covariates =                                          "~ 1 + sv.mat + as.factor(Batch) + Pipeline + Diet + cross +( 1 | DamID)"
#							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
#						} else {
#							covariates =                                          "~ 1 + sv.mat + as.factor(Batch)            + Diet + cross +( 1 | DamID)"
#							bestLambda = getBestLambda(lambdas, pheno, covariates, cov.data, sv.mat) 
#							lhs = formulaWrapper.boxCoxString("y",bestLambda = bestLambda)
#							
#							lmerTest::lmer (data=cov.data, as.formula(paste0(lhs, covariates)))
#						}
#			})
#		}
#
#
#		if (caught.error(fit) | caught.error(fit.main.only)) 
#		{
#			next
#		}
#		
#		# parse results
#		res <- try(CalcSignif(fit, fit.main.only))
#			  
#		if (!caught.error(res)) 
#		{
#			if (is.null(results.pipe)) 
#			{ # initiate results data frame
#				tmp.mat <- matrix(nrow=ncol(exp.mat), ncol=ncol(res)+1)
#				colnames(tmp.mat) <- c("Probe.Set.ID", colnames(res))
#				results.pipe <- as.data.frame(tmp.mat)
#				results.pipe$Probe.Set.ID <- annot.data$Probe.Set.ID
#			}
#			results.pipe[ jj, 2:(ncol(res)+1) ] <- res
#		}
#	} # about 1hr
# 	date()
#	results <- results.pipe
#	return(results)
#}


CalcSignif <- function(fit, fit.main.only) 
{
	aw = getAnovaWrapper(fit)
	includesPipe = "Pipeline" %in% rownames(aw$an)
	pipeP = NA 
	if(includesPipe)
	{
		pipeP = -log10(aw$an["Pipeline",         aw$pvalueCol])
	}
	out <- data.frame(
			LogP.Anova.SV          = -log10(aw$an["sv.mat",           aw$pvalueCol]),
			LogP.Anova.Batch       = -log10(aw$an["as.factor(Batch)", aw$pvalueCol]),
			LogP.Anova.Pipeline    =  pipeP,
			LogP.Anova.Diet        = -log10(aw$an["Diet",             aw$pvalueCol]),
			LogP.Anova.Cross       = -log10(aw$an["cross",            aw$pvalueCol]),
			LogP.Anova.CrossByDiet = -log10(aw$an["Diet:cross",       aw$pvalueCol]))
	
	if(is.na(out$LogP.Anova.Cross))
	{
		print("calc signif failed at some point!!!")
	}
	w = getTsWrapper(fit)
#	
	out[ , paste0("Est.CrossInStdCtrl") ]  <- w$ts[ "cross", w$estCol ]
	out[ , paste0("Tval.CrossInStdCtrl") ] <- w$ts[ "cross", w$tvalueCol ]
	out[ , paste0("LogP.CrossInStdCtrl") ] <- -log10(w$ts[ "cross", w$pvalueCol ])
	
	if(!mergeDiets)
	{
		dietz =  c("LowPro", "MethylSuff", "VitDDef")
	} else {
		dietz = c("Deficient")
	}
	for (diet in dietz) 
	{
		pred.string <- paste("Diet", diet, ":cross",sep="")
		out[ , paste0("Est.CrossBy", diet, ".vs.StdCtrl") ]  <-  w$ts[ pred.string, w$estCol ]
		out[ , paste0("Tval.CrossBy", diet, ".vs.StdCtrl") ] <-  w$ts[ pred.string, w$tvalueCol ]
		out[ , paste0("LogP.CrossBy", diet, ".vs.StdCtrl") ] <- -log10(w$ts[ pred.string, w$pvalueCol ])
	}
	
	
	aw = getAnovaWrapper(fit.main.only)
	out$LogP.Anova.MainOnly.Cross = -log10(aw$an["cross", aw$pvalueCol])
	
	w <- getTsWrapper(fit.main.only)
	out[ , paste0("Est.MainOnly.Cross") ]  <- w$ts[ "cross", w$estCol ]
	out[ , paste0("SE.MainOnly.Cross") ]   <- w$ts[ "cross", w$seCol ]
	out[ , paste0("DF.MainOnly.Cross") ]   <- w$ts[ "cross", w$dfCol ]
	out[ , paste0("Tval.MainOnly.Cross") ] <- w$ts[ "cross", w$tvalueCol ]
	out[ , paste0("LogP.MainOnly.Cross") ] <- -log10(w$ts[ "cross", w$pvalueCol ])
	if( sum(unlist(lapply(out, is.na)>0)) )#| sum(unlist(lapply( out, is.infinite)))>0 )#|)
	{
#		
	}
	return(out)
}

stackedPQ = function(annot.data, results)
{
	pres.data <- results[ , c("Probe.Set.ID", "LogP.Anova.Cross")]#, "Est.MainOnly.Cross", "SE.MainOnly.Cross")]
	colnames(pres.data) =  c("Probe.Set.ID", "LogP.Anova")
	pres.data$variable  = "strain" 
	if(sum(pres.data$Probe.Set.ID!=annot.data$Probe.Set.ID)>0)
	{
		stop("id mismatch")
	}
	descrip.cols <- which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data)
	pres.data                         <- cbind(pres.data, annot.data[ , descrip.cols])
	alldata = pres.data
	
	
	pres.data <- results[ , c("Probe.Set.ID", "LogP.Anova.Diet")]#, "Est.MainOnly.Cross", "SE.MainOnly.Cross")]
	colnames(pres.data) =  c("Probe.Set.ID", "LogP.Anova")
	pres.data$variable  = "diet" 
	if(sum(pres.data$Probe.Set.ID!=annot.data$Probe.Set.ID)>0)
	{
		stop("id mismatch")
	}
	descrip.cols <- which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data)
	pres.data                         <- cbind(pres.data, annot.data[ , descrip.cols])
	alldata = rbind(alldata, pres.data)
	
	pres.data <- results[ , c("Probe.Set.ID", "LogP.Anova.CrossByDiet")]#, "Est.MainOnly.Cross", "SE.MainOnly.Cross")]
	colnames(pres.data) =  c("Probe.Set.ID", "LogP.Anova")
	pres.data$variable  = "diet:strain" 
	if(sum(pres.data$Probe.Set.ID!=annot.data$Probe.Set.ID)>0)
	{
		stop("id mismatch")
	}
	descrip.cols <- which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data)
	pres.data                         <- cbind(pres.data, annot.data[ , descrip.cols])
	alldata = rbind(alldata, pres.data)
	
	pres.data = alldata
	pres.data$pval                    <- 10^-pres.data$LogP.Anova
	pres.data$qval                    <- p.adjust(pres.data$pval, method="fdr")
	pres.data$Log.qval                <- -log10(pres.data$qval)
#	
	
#	
	pres.data <- pres.data[ order(pres.data$pval, pres.data$variable),]
	return(pres.data)
	
}

getAnovaWrapper = function(fit)
{
	fit.lmer           = (class(fit) =="merModLmerTest")
	if(!fit.lmer)
	{
		return(list(
						an = anova(fit),
						pvalueCol = "p-value"))
	}
	else
	{
		return(list(
						an = lmerTest::anova(fit),
						pvalueCol = "Pr(>F)"))
	}
}

getTsWrapper=function(fit)
{
	fit.lmer = (class(fit) =="merModLmerTest")
	
	if(!fit.lmer)
	{
		return(list(
						ts        = summary(fit)$tTable,
						estCol    = "Value",
						tvalueCol = "t-value",
						pvalueCol = "p-value",
						seCol    = "Std.Error",
						dfCol    = "DF"
				))
	}
	else
	{
		return(list(
						ts        = coefficients(summary(fit)),
						estCol    = "Estimate",
						tvalueCol = "t value",
						pvalueCol = "Pr(>|t|)",
						seCol     = "Std. Error",
						dfCol     = "df"))
	}
}



pqplots = function(annot.data, results)
{
	descrip.cols <- which(colnames(annot.data)=="Gene.Accession"):ncol(annot.data)
	pres.data <- results[ , c("Probe.Set.ID", "LogP.Anova.Cross", "Est.MainOnly.Cross", "SE.MainOnly.Cross")]
	pres.data$pval                    <- 10^-pres.data$LogP.Anova.Cross
	pres.data$qval                    <- p.adjust(pres.data$pval, method="fdr")
	pres.data$Log.qval                <- -log10(pres.data$qval)
	pres.data                         <- cbind(pres.data, annot.data[ , descrip.cols])
	ii.signif                         <- pres.data$qval <= 0.05 & is.finite(pres.data$qval)
	pres.data <- pres.data[ ii.signif, ]
	pres.data <- pres.data[ order(pres.data$pval, -abs(pres.data$Est.MainOnly.Cross)), ]
	write.csv(file=file.path(outdir, "SignifCrossGenes_fdr05.csv"), pres.data, row.names=FALSE)
	
	
	
	###############
	##TODO redo as data.table
	library(lattice)
	pdf(file=file.path(outdir, "SignifCrossGenes_plots_fdr05.pdf"), height=7, width=7)
	labels <- rep(NA, nrow(cov.data))
	Pipeline.fact <- as.factor(paste("Pipeline", as.character(cov.data$Pipeline)))
	groups <- tapply(cov.data$ID, list(cov.data$Strain, cov.data$Diet, Pipeline.fact))
	
	for (g in unique(sort(groups))) 
	{
		ii <- which(g==groups)
		labels[ii] <- 1:length(ii)
	}
	labels <- as.character(labels)
	cex <- 1
	jitter.x <- TRUE
	for (ii in 1:nrow(pres.data)) 
	{
		id <- pres.data$Probe.Set.ID[ii]
		ai <- which(id==annot.data$Probe.Set.ID)
		gid <- as.character(pres.data$mRNA.Accession[ii])
		y <- exp.mat[ , ai]
		xyp <- xyplot( y ~ Strain | Diet * Pipeline.fact, data=cov.data, 
				ylab="Expression",
				main=paste0(
						gid, " (1/", sum(gid==pres.data$mRNA.Accession), " signif non-dup probe sets)\n",
						" probe set: ", id, " (1/", length(unlist(strsplit("daf_sdfs", "_"))), " duplicate probe sets)\n",
						"p-val=", signif(pres.data$pval[ii], 5), ", ",
						"q-val=", signif(pres.data$qval[ii], 5), "", "\n",
						"Effect of B6.NOD->NOD.B6 = ", round(pres.data$Est.MainOnly.Cross[ii],4), ""),
				col="black", jitter.x=jitter.x,
				panel = function(x, subscripts, ...) {
					panel.xyplot(pch=labels[subscripts], x=x, subscripts=subscripts, cex=cex, ...)
				})
		print(xyp, newpage=TRUE)
	}
	dev.off()
	return(pres.data)	
}

readPipelinePhen = function(cov.data)	
{
	pipe1.file <- file.path(data, "phenotypes/pipeline1_20120528.csv")
	pipe2.file <- file.path(data, "phenotypes/pipeline2_20120528.csv")
	
	phen.pipe1 <- read.csv(pipe1.file, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE)
	phen.pipe1$Pipeline <- 1
	phen.pipe2 <- read.csv(pipe2.file, stringsAsFactors=FALSE, strip.white=TRUE, blank.lines.skip=TRUE)
	phen.pipe2$Pipeline <- 2
	phen.data <- merge(phen.pipe1, phen.pipe2, by=c("Strain", "Diet", "ID", "Pipeline"), all.x=TRUE, all.y=TRUE)
#	phen.data <- merge(phen.pipe1, phen.pipe2, by=c("Diet", "ID", "Pipeline"), all.x=TRUE, all.y=TRUE)
	phen.data$ID <- paste("Mouse", sep=".", phen.data$ID)
	phen.data$Diet <- map.eq(phen.data$Diet, 
			list(
					methylsuff = "MethylSuff",
					vitddef    = "VitDDef",
					stdctrl    = "StdCtrl",
					lowpro     = "LowPro"))
	
	
	z = table(merge(cov.data, phen.data, by=c("Strain", "ID", "Pipeline", "Diet"), all.x=TRUE, all.y=TRUE, suffixes=c("",".phen"))$ID)
	print(z[z>1])
	phen.data <- merge(cov.data, phen.data, by=c("Strain", "ID", "Pipeline", "Diet"), all.x=TRUE, all.y=TRUE, suffixes=c("",".phen"))
	
	#take the mean of mouse 142 observations, drop redundant obs.
	for(colnam in allPhenNames)
	{
		phen.data[phen.data$ID=="Mouse.142",colnam] = mean(phen.data[phen.data$ID=="Mouse.142",colnam])
	}
	phen.data = phen.data[-which(phen.data$ID=="Mouse.142")[2],]
	if(length(table(table(phen.data$ID)))>1)
	{
		print(table(phen.data$ID)[table(phen.data$ID)>1])
		stop("duplicate data for at least one id, see printout above")
	}
	
	return(phen.data)
}



Calc.Cors <- function(phen, phen.data, exp.mat) 
{
	d <- data.frame(Probe.Set.ID = annot.data$Probe.Set.ID)
	d$Est.Cor       <- NA
	d$LogP.Cor      <- NA
	d$Est.SpearCor  <- NA
	d$LogP.SpearCor <- NA
	d$Est.KendCor   <- NA
	d$LogP.KendCor  <- NA
	x <- phen.data[, phen]
	ok <- !is.na(x)
	x <- x[ok]
	cat(sep="", "\n", phen, ":\n")
	for (i in 1:ncol(exp.mat))
	{
		if (i %% 1000 == 0) { cat(sep="","[",i,"]") } 
		y <- exp.mat[ok, i]
		
		ct <- cor.test(x,y, method=c("pearson"))
		d$Est.Cor[i] <- ct$estimate
		d$LogP.Cor[i] <- -log10(ct$p.value)
		
		ct <- suppressWarnings(cor.test(x,y, method=c("spearman"), exact=TRUE))
		d$Est.SpearCor[i] <- ct$estimate
		d$LogP.SpearCor[i] <- -log10(ct$p.value)
		
		ct <- suppressWarnings(cor.test(x,y, method=c("kendall"), exact=TRUE))
		d$Est.KendCor[i] <- ct$estimate
		d$LogP.KendCor[i] <- -log10(ct$p.value)
	}
	return(d)
}

plotPhenCor <- function(phen.cols, pres.data, cor.list, outdir)
{
	names(cor.list) <- phen.cols
	pdf(file.path(outdir, "phenotype_express_correlations_20120528.pdf"), height=10, width=12)
	par(mfrow=c(5,6))
	for (phen in phen.cols)
	{
		d <- cor.list[[phen]]
		draw.hists(d[, c(1, grep("LogP", colnames(d)))],             prefix=paste0(phen, "\n"), logP=TRUE, hits.at.fdr=0.05)
		draw.hists(d[,      grep("LogP", colnames(d), invert=TRUE)], prefix=paste0(phen, "\n"))
	}
	par(mfrow=c(5,6))
	for (phen in phen.cols) 
	{
		d <- cor.list[[phen]]
		d <- d[ d$Probe.Set.ID %in% pres.data$Probe.Set.ID[pres.data$qval<=0.05], ]
		draw.hists(d[, c(1, grep("LogP", colnames(d)))],             prefix=paste0(phen, "\n"), logP=TRUE, hits.at.fdr=0.05)
		draw.hists(d[,      grep("LogP", colnames(d), invert=TRUE)], prefix=paste0(phen, "\n"))
	}
	dev.off()
}

getSnordResults <- function(results, snordfile) 
{
	snord.data       <- read.csv(snordfile, stringsAsFactors=FALSE)
#snord.data$Location <- sub("\xca", "", snord.data$Location)
	snord.data$start <- as.integer(gsub(",", "", (sub("-.*", "", snord.data$Location))))
	snord.data$end   <- as.integer(gsub(",", "", sub(".*-", "", snord.data$Location)))
	
	
	snord.results <- dfapply(snord.data, snord.data$mRNA.Accession, results.add.FUN = rbind,
			FUN = function(d, ...){
				ai <- which(as.character(d$mRNA.Accession)==annot.data$mRNA.Accession)
				
				lp <- NA
				beta <- NA
				n <- NA
				se <- NA
				degf <- NA
				if (1<=length(ai)) {
					lp <- results$LogP.Anova.MainOnly.Cross[ai]
					beta <- results$Est.MainOnly.Cross[ai]
					se   <- results$SE.MainOnly.Cross[ai]
					degf <- results$DF.MainOnly.Cross[ai] 
					n <- annot.data$Num.Pure.Probes[ai]      
				}
				d <- cbind(d, data.frame(LogP.Cross=lp, Est.Cross=beta, SE.Cross=se, DF.Cross=degf, Num.Probes=n))
				
				# d$Anova.LogP.Cross <- results$Anova.LogP.Cross[ai]
				# d$Est.MainOnly.Cross <- results$Est.MainOnly.Cross[ai]
				return (d)
			})
	snord.results <- snord.results[!is.na(snord.results$LogP.Cross),]
	
	snord.results$Lower95 <- snord.results$Est.Cross - qt(0.975, df=snord.results$DF.Cross)*snord.results$SE.Cross
	snord.results$Upper95 <- snord.results$Est.Cross + qt(0.975, df=snord.results$DF.Cross)*snord.results$SE.Cross
	#the midpoint in megabases of the relevant mRNA
	snord.results$Mid.Mb <- 0.5*(snord.results$start + snord.results$end)/10^6
	return(snord.results)
}

#also writes to csv
plot.snord.results <- function(outdir, snord.results) 
{
	snord.basename <- "snordeffects_20140610"
	pdf(file.path(outdir, paste0(snord.basename, ".pdf")), height=2.5, width=6)
	par(mar=c(top=0.5, bottom=3.5, left=4, right=0.2))
	
#z <- snord.results$Num.Probes/max(snord.results$Num.Probes)
	z <- (snord.results$Num.Probes-min(snord.results$Num.Probes))/diff(range(snord.results$Num.Probes))
	cols <- gray(1-z)
	
#par(mfrow=c(1,1))
	plot(snord.results$Mid.Mb, snord.results$Est.Cross,
			main="", #paste0("Estimated effects and 95% CIs for ", nrow(snord.results), " probesets in SNORD cluster"), 
			type="n",
			ylim=range(c(snord.results[, c("Lower95", "Upper95")])),
			las=1,
			xlab="",
			ylab="Effect of B6.NOD -> NOD.B6            ")
	title(xlab="Position on Chr 7 (Mb)", line=2)
	abline(h=0, col="gray")
	segments(snord.results$Mid.Mb, snord.results$Lower95, snord.results$Mid.Mb, snord.results$Upper95)
	
	
	points.discs(snord.results$Mid.Mb, snord.results$Est.Cross, col.fill=cols)    
	
	legend("topright", pch=c(1,19), ncol=2, legend=c(paste0(min(snord.results$Num.Probes), " probes"),
					paste0(max(snord.results$Num.Probes), " probes")), 
			title="Shading shows no. of probes contributed")    
	dev.off()
	write.csv(row.names=FALSE, file=file.path(outdir, paste0(snord.basename, ".csv")), snord.results)
}

points.discs <- function(x, y, pch=19, pch.fill=pch, pch.outline=1,col="black", col.fill=col, col.outline="black",...)
{
	col.fill=rep(col.fill, length.out=length(x))
	col.outline=rep(col.outline, length.out=length(x))
	for (i in 1:length(x))
	{
		points(x[i],y[i], pch=pch.fill, col=col.fill[i], ...)
		points(x[i],y[i], pch=pch.outline, col=col.outline[i], ...)
	}
}

plotsigpq <- function(outdir, sigpq) 
{
	pdf(file=file.path(outdir, "pval_005_scatter.pdf"), height=10, width=15)
	aplot = ggplot(sigpq,aes(x=variable, y=LogP.Anova))
	aplot = aplot + geom_point()
	aplot = aplot + geom_line(aes(y= -log(.05, base=10), x=as.numeric(ordered(variable))), col="black")
	aplot = aplot + ylab("-log_10.pval")
	aplot = aplot + xlab("variable")
	aplot = aplot + geom_jitter(position = position_jitter(width = .2))
	print(aplot)
	dev.off()
	
	pdf(file=file.path(outdir, "pval_005_grouphist.pdf"), height=10, width=15)
	aplot = ggplot(sigpq,aes(x=LogP.Anova,))
	s1 = sigpq[sigpq$variable=="diet",]
	s2 = sigpq[sigpq$variable=="strain",]
	s3 = sigpq[sigpq$variable=="diet:strain",]
	aplot = aplot + geom_histogram(data = s1, alpha = .2,  fill="red",binwidth = .2)
	aplot = aplot + geom_histogram(data = s2, alpha = .2,  fill="blue",binwidth = .2)
	aplot = aplot + geom_histogram(data = s3, alpha = .2,  fill="yellow",binwidth = .2)
#	aplot = aplot + geom_line(aes(y= -log(.05, base=10), x=as.numeric(ordered(variable))), col="black")
	aplot = aplot + ylab("-log_10.pval")
	aplot = aplot + xlab("variable")
#	aplot = aplot + geom_jitter(position = position_jitter(width = .2))
	print(aplot)
	dev.off()
}

##### GLOBAL OPTIONS #####

# mask0-24 and mask3-21 give exactly the same rma-sketch summary
#MASK.TYPE                <- "mask0-24"
MASK.TYPE                 <- "rel1410_3-21"
outdir                    <- file.path(output, paste0("expression_analysis_", MASK.TYPE))
dir.create(outdir, showWarnings=F)
apt.outdir                <- file.path(data, "sva", paste0("apt_output_mask_", MASK.TYPE))

maskedraw.expression.file <- file.path(apt.outdir, "/rma.summary.txt")
#maskedraw.expression.file <- file.path(apt.outdir, "/rma-sketch.summary.txt")
#TODO fix this to reside in output
design.file               <- file.path(data,   "sva", "expression_choice_20120411.csv")
annot.old.expression.file <- file.path(data,   "sva", "log2RMA-GENE-DEFAULT-Group1.txt")
probeset.info.file        <- file.path(output, "probeset_info_summary.csv")
snordfile                  = file.path(data,   "sva", "SNORD_cluster_20120524.csv")

#expression.file <- "linearRMA-GENE-DEFAULT-Group1.txt"
#maskedraw.expression.file <- "log2RMA_masked_raw_20120510.txt"
# annotated unmasked data

###### Processing masked but unannotated expression data
# Probes were masked if they had a SNP within 3-21 inclusive
# Probe sets were excluded only if all probes were masked
# masked unannotated data
annot.data = get.annot.data(annot.old.expression.file = annot.old.expression.file,
		maskedraw.expression.file = maskedraw.expression.file,
		probeset.info.file = probeset.info.file, 
		output= outdir)

write.delim(file=file.path(outdir, "log2RMA_maskedannot.txt"), annot.data)

mouse.cols <- grep("Mouse", colnames(annot.data))
exp.mat  <- t(as.matrix(annot.data[ , mouse.cols]))
cov.data = form.Cov.Data(exp.mat = exp.mat, design.file = design.file)

makeClusterDendogram(exp.mat = exp.mat, 
		cov.data = cov.data, 
		outdir = outdir, 
		maskedraw.expression.file = maskedraw.expression.file)

#############################
# Analysis Pre-Step 

# collapse technical replicate

keep  <- which(cov.data$ID=="Mouse.142")
blend <- which(cov.data$ID=="Mouse.142_2")

cov.data       <- cov.data[-blend,]
exp.mat[keep,] <- 0.5 * (exp.mat[keep,] + exp.mat[blend,])
exp.mat        <- exp.mat[-blend,]

pcRes = prcomp(exp.mat,scale. = T)
pcs = (as.matrix(exp.mat)%*%(pcRes$rotation))

fullthing = cbind(cov.data, pcRes$x)
fullthing$label = as.character(fullthing$Diet)
fullthing$label[as.character(fullthing$label)=="LowPro"] = "p"
fullthing$label[as.character(fullthing$label)=="StdCtrl"] = "s"
fullthing$label[as.character(fullthing$label)=="VitDDef"]= "d"
fullthing$label[as.character(fullthing$label)=="MethylSuff"]= "m"

#df = 
aplot = ggplot(fullthing[fullthing$Pipeline==1,], aes(x=PC1, y=PC2, shape=Strain, label=label, color=label))
aplot = aplot + geom_text(size=5)
pdf(file.path(output, "dietpca.pdf"))
print(aplot)
dev.off()

aplot = ggplot(fullthing, aes(x=PC1, y=PC2, shape=Strain, label=label, color=Strain))
aplot = aplot + geom_text(size=5)
pdf(file.path(output, "dietstrainpca.pdf"))
print(aplot)
dev.off()

### Weird p-value distributions found, so...

##############################
# surrogate variable analysis

#save(file="exp.mat.object", exp.mat)
pipelineGroups = list("p1.2.pool"=c(1,2), "p1"=1, "p2"=2)
#pipelineGroups = list("p1"=1, "p2"=2)
#pipelineGroups = list("p1.2.pool"=c(1,2))
exp.mat.full = exp.mat
cov.data.full = cov.data

outdir.full = outdir
# Code to do this separately for each pipeline

for(pipename in names(pipelineGroups))
{
	print(pipename)
	pipelineGroup    = pipelineGroups[[pipename]]
	pipeInd          = cov.data.full$Pipeline %in% pipelineGroup 
	exp.mat          = exp.mat.full[pipeInd,]
	cov.data         = cov.data.full[pipeInd,]
	#Step 1.1 from paper
	exp.resid.mat_bestLamdas = getResids(exp.mat = exp.mat, cov.data = cov.data)
	exp.resid.mat            = exp.resid.mat_bestLamdas$exp.resid.mat 
	bestLambdas              = exp.resid.mat_bestLamdas$bestLambdas
	#Steps 1 and 2 from paper
	sv.mat           = refine.sv(exp.mat, exp.resid.mat)
	#load(file=file.path(outdir, "sva_exp.resid.mat.object"))
	
	outdir = file.path(outdir.full, pipename)
	dir.create(file.path(outdir), showWarnings = F)
	saveddir = file.path(outdir, "savedobj")
	dir.create(saveddir, recursive=T)
	print(saveddir)
	save(file=file.path(saveddir, "post.sv2.object"), list=ls())
#	load("../output/expression_analysis_3-21/p1.2.pool/savedobj/post.sv2.object")
#===================================================================
	
	#step 2.8
	results <- crossEffectsAtDiet(cov.data = cov.data, exp.mat = exp.mat, sv.mat = sv.mat, annot.data = annot.data, bestLambdas = bestLambdas)
	write.csv(file=file.path(outdir, "CrossDietEffects_BothPipes.txt"), results, row.names=FALSE)
	save(file=file.path(saveddir, "post.results2.object"), list=ls())
#	load("../output/expression_analysis_3-21/p1.2.pool/savedobj/post.results2.object")
	
	pdf(file=file.path(outdir, "results_histograms.pdf"), height=10, width=15)
	par(mfcol=c(4,6))
	draw.hists(results[ , c(1, grep("LogP", colnames(results))[-1])], logP=TRUE, hits.at.fdr=0.05)
	draw.hists(results[ , c(1, grep("LogP", colnames(results))[-1])], logP=TRUE, hits.at.fdr=0.01)
	dev.off()
	
	#######################################
# Comparison of significant results from Alan's (wrong kill file)
# vs current (Will kill file)
	
#	new.q <- p.adjust(10^-results$LogP.Anova.Cross, method="fdr")
#	new.hits <- results$Probe.Set.ID[new.q <= 0.05 & is.finite(new.q)]
	
	############################
# Presenting the new results. as a side effect returns results in a handy form as pres.data
	
#	load("../output/expression_analysis_3-21/p1/savedobj/post.results2.object")
	pres.data = pqplots(annot.data, results)
	
	stackedPQdata = stackedPQ(annot.data,results)
	siglevel = .05
	stackedPQdata = stackedPQdata[is.na(stackedPQdata$pval)| stackedPQdata$pval<=siglevel,]
	write.table(file=file.path(outdir, paste0("SignifGenes_p", siglevel, ".csv")), stackedPQdata, row.names=FALSE,sep="\t")
	
	ii.signif  <- stackedPQdata$pval <= 0.005 & is.finite(stackedPQdata$pval)
	sigpq <- stackedPQdata[ ii.signif, ]
	
	plotsigpq(outdir = outdir, sigpq = sigpq)
	
	ii.signif  <- stackedPQdata$qval <= 0.05 & is.finite(stackedPQdata$pval)
	sigpq <- stackedPQdata[ ii.signif, ]
	write.table(file=file.path(outdir, "SignifGenes_fdr05.csv"), sigpq, row.names=FALSE,sep="\t")
	
	################################################
	## Looking at SNORD cluster
	snord.results <- getSnordResults(results,snordfile)
	plot.snord.results(outdir = outdir, snord.results = snord.results)
	
	########################
	# Phenotype correlations
	phen.data = readPipelinePhen(cov.data)	
	phen.data = phen.data[pipeInd,]
	phen.cols <- grep("(OpenField|LightDark)", colnames(phen.data), value=TRUE)
	library(parallel)
	date()
	#this is actually not that slow, so just use lapply if mclapply gives problems. 
#	cor.list <- lapply(phen.cols, Calc.Cors, phen.data = phen.data, exp.mat=exp.mat)
	cor.list <- mclapply(phen.cols, Calc.Cors, phen.data = phen.data, exp.mat=exp.mat, mc.cores=mc.cores)
	
	date() # 25 minutes if mc.cores >= length(phen.cols)
	
	plotPhenCor(phen.cols = phen.cols, cor.list=cor.list, pres.data = pres.data, outdir = outdir)
	save(file = file.path(saveddir, "all.RData"), list=ls())
}




expressedCov = read.table("../data/expression_choice_20120411.csv", header=T, sep=",")
expressedCov$ID = paste0("Mouse.", expressedCov$ID)

