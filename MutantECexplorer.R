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
	sortedEnvironments <- sortEnvironments(classed)
	Tip <- length(sortedEnvironments$Tip)
	Stalk <- length(sortedEnvironments$Stalk)
	Phalanx <- length(sortedEnvironments$Phalanx)
	EC <- length(sortedEnvironments$EC)
	Other <- length(sortedEnvironments$Other)
	envcount <- c("Tip" = Tip, "Stalk" = Stalk, "Phalanx" = Phalanx, "EC" = EC, "Other" = Other)
	behavior <- list("Attractors" = Attractors,"Classed" = classed, "Sorted" = sortedEnvironments, "EnvCount" = envcount)
	if (filename != "none"){
		save(behavior, file = filename)
	}
	return(envcount)
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

# Load the simplified network from a file
net <- loadNetwork("~/Desktop/BoolNet/angiogenesis.txt")
fixed <- net$fixed
microenvironment_variables <- c("WNT5a", "WNT7a", "FGF","IGF", "BMP9", "BMP10", "TGFB1", "VEGFC_D", "AMPATP", "ANG1", "Oxygen", "ShearStress", "JAGp", "DLL4p", "VEGFC_Dp", "VEGFAxxxP")

mutants <- analizeMBparallel(net, microenvironment_variables, plot_folder = "~/Desktop/BoolNet/")

# A cada una de estas lineas hay que correrlas en otro procesador
#mutants1to16 <- analizeMutantBehavior(net, 1, 16, microenvironment_variables, plot_folder = "~/Desktop/BoolNet/")
#mutants17to32 <- analizeMutantBehavior(net, 17, 32, microenvironment_variables, plot_folder = "~/Desktop/BoolNet/")
#mutants33to48 <- analizeMutantBehavior(net, 33, 48, microenvironment_variables, plot_folder = "~/Desktop/BoolNet/")
#mutants49to64 <- analizeMutantBehavior(net, 49, 64, microenvironment_variables, plot_folder = "~/Desktop/BoolNet/")









