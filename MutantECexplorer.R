# R CMD BATCH MutantECexplorer.R
# Execute from R by typing "source("MutantECexplorer.R")"
# We will use BoolNet
if(!require(BoolNet)){
	install.packages("BoolNet")
}
if(!require(foreach)){
	install.packages("foreach")
}
if(!require(doParallel)){
	install.packages("doParallel")
}
library(BoolNet)
library(foreach)
library(doParallel)
source("matrixPlotter.R")
source("ECexplorerTools.R")


####  The parameters that define a microenvironment ####
# Calcium ions enter the cell through the membrane, however this is controlled by the activity of several ion channels and CCE.
# during angiogenesis the ER releases calcium ion when needed -> Calcium is not a direct input.
# Similarly, NO acts mostly autocrinally in ECs.
### Main parameters
# 1. ANG1 - The main blood vessel stabilizer
# 2. Oxygen - Lack of oxygen triggers the secretion of angiogenic signals
# 3. AMPATP - A high AMP/ATP ratio triggers the secretion of angiogenic signals (Energy)
# 4. ShearStress - mechanical forces
### The WNT ligands
# 5. WNT5a
# 6. WNT11
# 7. WNT7a
### FGF & IGF
# 8. FGF
# 9. IGF
### The NOTCH ligands
# 10. JAGp
# 11. DLL4p
### The VEGF ligands
# 12. VEGFAxxxP
# 13. VEGFC_D
# 14. VEGFC_Dp
# 15. PlGF
### The TGF ligands
# 16. BMP9
# 17. BMP10
# 18. TGFB1


exploreECbehavior <- function(net, microenvironment_variables, filename = "none"){
	Attractors <- getAttractors(net,
				type = "synchronous",
				method = "sat.exhaustive",
				maxAttractorLength = 36)
	classed <- classifyAttractors(Attractors, net, microenvironment_variables)
	sortedEnvironments <- sortEnvironments(classed$types)
	Tip <- length(sortedEnvironments$Tip)
	Stalk <- length(sortedEnvironments$Stalk)
	Phalanx <- length(sortedEnvironments$Phalanx)
	Atypical <- length(sortedEnvironments$Atypical)
	envcount <- c("Tip" = Tip, "Stalk" = Stalk, "Phalanx" = Phalanx, "Atypical" = Atypical)
	behavior <- list("Attractors" = Attractors,"Classed" = classed, "Sorted" = sortedEnvironments, "EnvCount" = envcount)
	if (filename != "none"){
		save(behavior, file = filename)
	}
	return(envcount)
}

countEnvironmentTypes <- function(classed){
	sortedEnvironments <- sortEnvironments(classed$types)
	Tip <- length(sortedEnvironments$Tip)
	Stalk <- length(sortedEnvironments$Stalk)
	Phalanx <- length(sortedEnvironments$Phalanx)
	Atypical <- length(sortedEnvironments$Atypical)
	envcount <- c("Tip" = Tip, "Stalk" = Stalk, "Phalanx" = Phalanx, "Atypical" = Atypical)
	return(envcount)
}

summarizeMutantResults <- function(mutant_folder_path){
	file_paths <- list.files(path=mutant_folder_path, pattern="*.Rdata", full.names=TRUE, recursive=FALSE)
	file_names <- list.files(path=mutant_folder_path, pattern="*.Rdata", full.names=FALSE, recursive=FALSE)
	df <- data.frame()
	for(f in file_paths){
		print(f)
		load(f)
		ce <- countEnvironmentTypes(behavior$Classed)
		df <- rbind(df, ce)
		print(ce)
	}
	rownames(df) <- file_names
	colnames(df) <- names(ce)
	nTip <- subset(df, Tip == 0)
	lTip <- subset(df, Tip > 0 & Tip < behavior$EnvCount["Tip"])
	eTip <- subset(df, Tip == behavior$EnvCount["Tip"])
	mTip <- subset(df, Tip > behavior$EnvCount["Tip"])
	TipMutants <- list(none = nTip, less = lTip, equal = eTip, more = mTip)
	nStalk <- subset(df, Stalk == 0)
	lStalk <- subset(df, Stalk > 0 & Stalk < behavior$EnvCount["Stalk"])
	eStalk <- subset(df, Stalk == behavior$EnvCount["Stalk"])
	mStalk <- subset(df, Stalk > behavior$EnvCount["Stalk"])
	StalkMutants <- list(none = nStalk, less = lStalk, equal = eStalk, more = mStalk)
	nPhalanx <- subset(df, Phalanx == 0)
	lPhalanx <- subset(df, Phalanx > 0 & Phalanx < behavior$EnvCount["Phalanx"])
	ePhalanx <- subset(df, Phalanx == behavior$EnvCount["Phalanx"])
	mPhalanx <- subset(df, Phalanx > behavior$EnvCount["Phalanx"])
	PhalanxMutants <- list(none = nPhalanx, less = lPhalanx, equal = ePhalanx, more = mPhalanx)
	nAtypical <- subset(df, Atypical == 0)
	lAtypical <- subset(df, Atypical > 0 & Atypical < behavior$EnvCount["Atypical"])
	eAtypical <- subset(df, Atypical == behavior$EnvCount["Atypical"])
	mAtypical <- subset(df, Atypical > behavior$EnvCount["Atypical"])
	AtypicalMutants <- list(none = nAtypical, less = lAtypical, equal = eAtypical, more = mAtypical)
	mutant_behavior <- list(all = df, Tip = TipMutants, Stalk = StalkMutants, Phalanx = PhalanxMutants, Atypical = AtypicalMutants)
	return(mutant_behavior)
}

# For each gene in range [first, last] within the gene list in net, 
# the analisis of each gene is saved in a genelf.Rdata or genegf.Rata
analizeMutantBehavior <- function(net, first, last, microenvironment_variables, plot_folder = "none", verbose = FALSE){
	genes <- net$genes
	testvector <- c()
	rnames <- c()
	analisisgenes <- genes[c(first:last)]
	for (gene in analisisgenes){
		netl <- fixGenes(net, gene, 0)
		namel <- paste(gene, "lf", sep = " ")
		namefl <- paste(gene, "lf.Rdata", sep = " ")
		analisisl <- exploreECbehavior(netl, microenvironment_variables, filename = namefl)
		if(verbose){
			print(namel)
			print(analisisl)
		}
		rnames <- append(rnames, namel)
		testvector <- append(testvector,analisisl)
		netg <- fixGenes(net, gene, 1)
		nameg <- paste(gene, "gf", sep = " ")
		namefg <- paste(gene, "gf.Rdata", sep = " ")
		analisisg <- exploreECbehavior(netg, microenvironment_variables, filename = namefg)
		if(verbose){
			print(nameg)
			print(analisisg)
		}
		rnames <- append(rnames, nameg)
		testvector <- append(testvector,analisisg)
	}
	mutantMatrix <- matrix(testvector, nrow = 2*length(analisisgenes), byrow = TRUE)
	colnames(mutantMatrix) <- analisisgenes
	rownames(mutantMatrix) <- rnames
	# Check if the plot_folder exists, if not create it
	if(plot_folder != "none"){
		if(file.exists(plot_folder)){
		}else{
			dir.create(plot_folder)
		}
		dim1 <- dim(mutantMatrix)
		r1 <- dim1[1]
		r2 <- dim1[2]
		w <- r2 
		h <- r1 / 4
		plot_file <- paste(plot_folder, "mutantEffect", sep = "")
		plot_file <- paste(plot_file, toString(first), sep = "")
		plot_file <- paste(plot_file, toString(last), sep = "")
		png(file = plot_file, width = w, height = h, units = "in", res = 100)
		MatrixPlot(mutantMatrix, "marlist" = c(2,6,6,2), "aboveaxis" = TRUE, "title" = "Mutation effects on EC behavior")
		dev.off()
	}
	return(mutantMatrix)
}

# Analize all genes in parallele 
# The analisis of each gene is saved in a genelf.Rdata or genegf.Rata
analizeMBparallel <- function(net, microenvironment_variables, plot_folder = "none", verbose = FALSE){
	no_cores <- detectCores() - 1
	print("Number of processes:") 
	print(no_cores)
	registerDoParallel(cores=no_cores)  
	cl <- makeCluster(no_cores, type="FORK")
	genes <- net$genes
	testvector <- c()
	rnames <- c()
	foreach(i = 2:(2*length(genes) + 1)) %dopar% {
		print(i)
		gene <- genes[i%/%2]
		if(i%%2 == 0){
			netl <- fixGenes(net, gene, 0)
			namel <- paste(gene, "lf", sep = " ")
			namefl <- paste(gene, "lf.Rdata", sep = " ")
			analisisl <- exploreECbehavior(netl, microenvironment_variables, filename = namefl)
			if(verbose){
				print(namel)
				print(analisisl)
			}
			rnames <- append(rnames, namel)
			testvector <- append(testvector,analisisl)
		}else{
			netg <- fixGenes(net, gene, 1)
			nameg <- paste(gene, "gf", sep = " ")
			namefg <- paste(gene, "gf.Rdata", sep = " ")
			analisisg <- exploreECbehavior(netg, microenvironment_variables, filename = namefg)
			if(verbose){
				print(nameg)
				print(analisisg)
			}
			rnames <- append(rnames, nameg)
			testvector <- append(testvector,analisisg)
		}
	}
	mutantMatrix <- matrix(testvector, nrow = 2*length(analisisgenes), byrow = TRUE)
	colnames(mutantMatrix) <- analisisgenes
	rownames(mutantMatrix) <- rnames
	# Check if the plot_folder exists, if not create it
	if(plot_folder != "none"){
		if(file.exists(plot_folder)){
		}else{
			dir.create(plot_folder)
		}
		dim1 <- dim(mutantMatrix)
		r1 <- dim1[1]
		r2 <- dim1[2]
		w <- r2 
		h <- r1 / 4
		plot_file <- paste(plot_folder, "mutantEffect", sep = "")
		plot_file <- paste(plot_file, toString(first), sep = "")
		plot_file <- paste(plot_file, toString(last), sep = "")
		png(file = plot_file, width = w, height = h, units = "in", res = 100)
		MatrixPlot(mutantMatrix, "marlist" = c(2,6,6,2), "aboveaxis" = TRUE, "title" = "Mutation effects on EC behavior")
		dev.off()
	}
	stopCluster(cl)
	return(mutantMatrix)
}

# A paralelized function to analize the robustness of the update rules of a network
getUpdateRuleRobustnessParalel <- function(net, nSamples = 1000){
	no_cores <- detectCores() - 1
	print("Number of processes:") 
	print(no_cores)
	registerDoParallel(cores = no_cores)  
	cl <- makeCluster(no_cores, type="FORK")
	genes <- net$genes
	foreach(i = 1:length(genes)) %dopar% {
		r <- perturbTrajectories(net, measure = "sensitivity", numSamples = nSamples, flipBits = 1, gene = genes[i])
		print(genes[i])
		print(r$value)
	}
	stopCluster(cl)
	return(1)
}

# A function to analize the robustness of the update rules of a network
getUpdateRuleRobustness <- function(net, nSamples = 1000){
	genes <- net$genes
	sensitivities <- c()
	for(g in genes){
		r <- perturbTrajectories(net, measure = "sensitivity", numSamples = nSamples, flipBits = 1, gene = g)
		sensitivities <- append(sensitivities,r$value)
		#print(genes[i])
		#print(r$value)
	}
	names(sensitivities) <- genes
	return(sensitivities)
}

# A function to calculate the robustness of Phalanx, Stalk and Tip ECs.
getECrobustness <- function(net, numSamples = 10000, flipBits = 1){
	tip <- 0
	stalk <- 0
	phalanx <- 0
	ctip <- 0
	cstalk <- 0
	cphalanx <- 0
	ctotal <- 0
	for(i in 1:numSamples){
		# Generate a random state
		si <- floor(runif(length(net$genes), 0, 1.999))
		names(si) <- net$genes
		# Generate a perturbed version of the random state
		sii <- si
		for(j in 1:flipBits){
			pj <- floor(runif(1, 1, length(net$genes) + 0.999))
			if(sii[pj] == 0){
				sii[pj] <- 1
			}else{
				sii[pj] <- 0
			}
		}
		# Get and sort the atractor of the random state and the perturbed state
		attr <- getAttractor(net, si)
		kind <- sortAtractor(analyzeAttractor(net,  attr))
		attrp <- getAttractor(net, sii)
		kindp <- sortAtractor(analyzeAttractor(net,  attrp))
		# Store the result
		if(kind == kindp){
			ctotal <- ctotal + 1
		}
		if(kind == "Tip"){
			tip <- tip + 1
			if(kindp == "Tip"){
				ctip <- ctip + 1
			}
		}else if(kind == "Stalk"){
			stalk <- stalk + 1
			if(kindp == "Stalk"){
				cstalk <- cstalk + 1
			}
		}else if(kind == "Phalanx"){
			phalanx <- phalanx + 1
			if(kindp == "Phalanx"){
				cphalanx <- cphalanx + 1
			}
		}
	}
	rtip <- ctip/tip
	rstalk <- cstalk/stalk
	rphalanx <- cphalanx/phalanx
	rtotal <- ctotal/numSamples
	return(c(rtip = rtip, rstalk = rstalk, rphalanx = rphalanx, rtotal = rtotal))
}

## First verify that the models are in the working directory
WD <- getwd()
simple_model <- paste(WD, "angiogenesis.net", sep = "/")
if(file.exists(simple_model)){
}else{
	print("The file containing the simplified model is not in the working directory")
}
# Load the simplified network from a file
net <- loadNetwork(simple_model)
fixed <- net$fixed
microenvironment_variables <- c("WNT5a", "WNT7a", "FGF","IGF", "BMP9", "BMP10", "TGFB1", "VEGFC_D", "AMPATP", "ANG1", "Oxygen", "ShearStress", "JAGp", "DLL4p", "VEGFC_Dp", "VEGFAxxxP")
mutant_plot_folder <- paste(WD, "MutantPlotFolder/", sep = "/")
mutants <- analizeMBparallel(net, microenvironment_variables, plot_folder = mutant_plot_folder)










