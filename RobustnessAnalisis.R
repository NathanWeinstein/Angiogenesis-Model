# R CMD BATCH MutantECexplorer.R
# Execute from R by typing "source("MutantECexplorer.R")"
# We will use BoolNet
if(!require(BoolNet)){
	install.packages("BoolNet")
}
library(BoolNet)
source("matrixPlotter.R")
source("ECexplorerTools.R")

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
simple_model <- paste(WD, "angiogenesis.txt", sep = "/")
if(file.exists(simple_model)){
}else{
	print("The file containing the simplified model is not in the working directory")
}
# Load the simplified network from a file
net <- loadNetwork(simple_model)

# 1) Robustnes of attractor detemination to molecular activation noise
attractor_robustness <- perturbTrajectories(net, measure = "attractor",numSamples = 1000, flipBits = 1)

# 2) Robustness of the Tip, Stalk and Phalanx EC behavior
ECbehavior_robustness <- getECrobustness(net, numSamples = 1000, flipBits = 1)

# 3) Sensitivity of the components of the update function
update_rule_sensitivity <- getUpdateRuleRobustness(net, nSamples = 1000)

# Save the results
robustness <- list("attractor" = attractor_robustness, "ECbehavior" = ECbehavior_robustness, "sensitivities" = update_rule_sensitivity)
save(robustness, file = "robustness.Rdata")








