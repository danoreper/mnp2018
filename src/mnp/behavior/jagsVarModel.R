library(rjags)

n.chains = 1
n.iter   = 100000

#Used to evaluate JAGS model of variance
runJags <- function(allData,phencol) 
{

	allData$Batch.j            = match(as.character(allData$Batch), unique(as.character(allData$Batch)))
	allData$Order.j            = match(as.character(allData$Order), unique(as.character(allData$Order)))
	allData$Diet.j             = match(as.character(allData$Diet), unique(as.character(allData$Diet))) 
	allData$Sire.is.b6.j       = match(as.character(allData$Sire.is.b6),  unique(as.character(allData$Sire.is.b6)))
 	allData$Dam.j              = match(as.character(allData$Dam.ID),  unique(as.character(allData$Dam.ID))) 
	
	data= list()
	data$d.batch.j      = allData$Batch.j
	data$d.order.j      = allData$Order.j
	data$d.order        = allData$Order
	data$d.diet.j       = allData$Diet.j
	data$d.sire.is.b6.j = allData$Sire.is.b6.j
	data$d.dam.j        = allData$Dam.j
	data$d.dietBySire.j = allData$Diet.j* allData$Sire.is.b6.j
	data$d.phen         = allData[[phencol]]
	
	data$n.sire.is.b6   = length(unique(allData$Sire.is.b6))
	data$n.order        = length(unique(allData$Order))
	data$n.diet         = length(unique(allData$Diet))
	data$n.batch        = length(unique(allData$Batch))
	data$n.dam          = length(unique(allData$Dam.j))
	data$n.dietBySire   = max(data$d.dietBySire.j)
	data$n              = nrow(allData)
	
	strainname         <- allData$Order
	strainindex        <- allData$Order.j
	mappings           = returnMappings(strainname = strainname, strainindex = strainindex)
	orderToIndex       = mappings$strainToIndex
	indexToOrder       = mappings$indexToStrain
	
	strainname         <- allData$Sire.is.b6
	strainindex        <- allData$Sire.is.b6.j
	mappings           = returnMappings(strainname = strainname, strainindex = strainindex)
	Sire.is.b6.ToIndex = mappings$strainToIndex
	indexToSire.is.b6  = mappings$indexToStrain
	

	
	strainname         <- allData$Batch
	strainindex        <- allData$Batch.j
	mappings           = returnMappings(strainname = strainname, strainindex = strainindex)
	batchToIndex       = mappings$strainToIndex
	indexToBatch       = mappings$indexToStrain
	
	strainname         <- allData$Diet
	strainindex        <- allData$Diet.j
	mappings           =  returnMappings(strainname = strainname, strainindex = strainindex)
	dietToIndex        = mappings$strainToIndex
	indexToDiet        = mappings$indexToStrain
	
	strainname         <- allData$Dam.ID
	strainindex        <- allData$Dam.j
	mappings           =  returnMappings(strainname = strainname, strainindex = strainindex)
	damToIndex        = mappings$strainToIndex
	indexToDam        = mappings$indexToStrain
	
	startTime = proc.time()
#Jags fit
#
	vals = c()	
	for(bugname in c("sih.bug"))
	{
		print(bugname)
		jags <- jags.model(bugname,data=data, n.chains = n.chains, n.adapt = n.iter)
		update( jags , n.iter=.5*n.iter )
		
		#uncomment variables of interest.
		vals = coda.samples(jags, c(
								   'intercept',
#									"batch",
#									"order",
								   "intercept.v",
#								   "diet",
#								   "diet.v",
								   "sire.is.b6",
								   "sire.is.b6.v"
#								   "dietBySire",
#								   "dietBySire.v"
#								   "dam",
#								   "dam.v",
#								   "sigma.inv.batch",
#								   "sigma.inv.diet",
#								   "sigma.inv.sire.is.b6",
#									"sigma.inv.order",
#									"sigma.inv.batch.v",
#									"sigma.inv.diet.v",
#									"sigma.inv.sire.is.b6.v",
#									"sigma.inv.order.v"
											),
									n.iter,thin=1)
		
	}
	
	totalTime = proc.time() - startTime
	print("total time for jags:")
	print(totalTime)
	return(list(vals=vals,
					indexToOrder       = indexToOrder,
					indexToSire.is.b6  = indexToSire.is.b6,
					indexToBatch       = indexToBatch,
					indexToDiet        = indexToDiet,
					indexToDam         = indexToDam,
					data               = data,
					allData = allData))
}


#Helper function for Jags variance model
returnMappings <- function(strainname, strainindex)
{
#	
	strainMapping = data.frame(momstrain=strainname, momstrain.j= strainindex)
	strainMapping$momstrain = as.character(strainMapping$momstrain)
	#TODO make sure this works
	strainMapping = na.omit(unique(strainMapping))
	
	indexToStrain = rep(0, max(strainMapping$momstrain.j))
	indexToStrain[strainMapping$momstrain.j] = strainMapping$momstrain
	strainToIndex = 1:nrow(strainMapping)
	names(strainToIndex) = indexToStrain
	return(list(strainToIndex = strainToIndex, indexToStrain=indexToStrain))
}


applyJagsToData = function(phen)
{
	valz=runJags(phen$SIH.full, "Difference")
	hpd = HPDinterval(as.mcmc(valz$vals))
	hpd = hpd[hpd[,1]*hpd[,2]>0,]
	print(hpd)
	return()
}

simulateJags = function()
{
#Runs jags simulation
	numNodB6 = 50
	nodMean = 1.5
	b6Mean  = 1.5 
	
	nodVar =.4
	b6Var  =.2
	
	Sire.is.b6 = c(rep(F,numNodB6), rep(T, numNodB6))
	mu      = rep(0, length(Sire.is.b6))
	mu[Sire.is.b6==F] = nodMean
	mu[Sire.is.b6==T] = b6Mean
	
	sigma2             = rep(0, length(Sire.is.b6))
	sigma2[Sire.is.b6==F] = nodVar
	sigma2[Sire.is.b6==T] = b6Var
	phen2 = data.frame(phen = rnorm(n=length(mu), mu, sqrt(sigma2)), Sire.is.b6=Sire.is.b6)
	phen2$Batch = 0
	phen2$Order = 1
	phen2$Diet  = 1
	phen2$ID = 1:nrow(phen2)
	phen2$Dam.ID= phen2$ID#1
	valz=runJags(phen2, "phen")
	print(summary(as.mcmc(valz$vals)))
    }



runHGLM <- function(bootnum, aframe, phenName)
{
	ests = rep(0,bootnum)
	Dam.ID.j = match(aframe$Dam.ID, unique(aframe$Dam.ID))
	for(i in 1:bootnum)
	{
		print(i)
		afit = hglm(fixed = as.formula(paste0(phenName, " ~ Sire.is.b6")), disp = ~Sire.is.b6, random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
#	afit2 = hglm(fixed = as.formula(paste0(phenName, " ~ Sire.is.b6")),  random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
		asum = summary(afit)
		
		intercept   = asum$FixCoefMat["(Intercept)","Estimate"]
		b6True      = asum$FixCoefMat["Sire.is.b6TRUE","Estimate"]
#	print(b6True)
		
		Dam.ID.var  = asum$varRanef
		Dam.ID = rnorm(n=length(unique(aframe$Dam.ID)), mean=0, sd = sqrt(Dam.ID.var))
		
		mu = intercept + b6True*aframe$Sire.is.b6 + Dam.ID[Dam.ID.j]
#	print(range(mu))
		
		intercept.v = asum$SummVC1["(Intercept)","Estimate"]
		b6True.v    = asum$SummVC1["Sire.is.b6TRUE","Estimate"]
		vars = exp(intercept.v+b6True.v*aframe$Sire.is.b6)
		
		aframe[[phenName]] = rnorm(n=length(mu), mean = mu, sqrt(vars))
		ests[i] = b6True.v
		if(i==1)
		{
			theEst = b6True.v
		}
	}
	
	phenName = "Difference"
	aframe = phen$SIH.full
	
	afit1 = hglm(fixed = as.formula(paste0(phenName, " ~ as.factor(Batch) + Diet + Sire.is.b6")), random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	afit2 = hglm(fixed = as.formula(paste0(phenName, " ~ as.factor(Batch) + Diet")),              random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	afit3 = hglm(fixed = as.formula(paste0(phenName, " ~ as.factor(Batch)")),  random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	
	
	afit1 = hglm(fixed = as.formula(paste0(phenName, " ~ Batch + Diet + Sire.is.b6")), disp = ~Batch + Diet +Sire.is.b6, random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	afit2 = hglm(fixed = as.formula(paste0(phenName, " ~ Batch + Diet + Sire.is.b6")), disp = ~Batch + Diet, random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	
	
	afit = hglm(fixed = as.formula(paste0(phenName, " ~ Sire.is.b6")), disp = ~Sire.is.b6, random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
	afit2 = hglm(fixed = as.formula(paste0(phenName, " ~ Sire.is.b6")),  random = ~1 | Dam.ID, data = aframe, calc.like=T, method="HL11")
}

#Runs jags analsysis of variance.
#applyJagsToData(phen)

#examine the effect of various parameter choices on the ability of jags to reconstruct anything.
#simulateJags()

#hglm modeling of variance
#aframe = phen$SIH.full
#phenName = "Difference"
#bootnum = 10000
#runHGLM(bootnum = bootnum, aframe = aframe, phenName = phenName)

#levene test of variance
#leveneTest(phen$SIH.full$Difference,group=as.factor(phen$SIH.full$Sire.is.b6), center="median")
