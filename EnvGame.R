#Use after runing ECexplorer.R only




analyzePhalanxBehavior <- function(net, Attractors, sortedEnvironments, classed, microenvironment_variables){
	phalanx_environments <- sortedEnvironments$Phalanx
	classed_indexes <- classed$indexes
	phalanx_attractor_indexes <- c()
	for(env in phalanx_environments){
		env_attr_indexes <- unlist(classed_indexes[env])
		phalanx_attractor_indexes <- append(phalanx_attractor_indexes, env_attr_indexes)
	}
	phalanx_attractors <- matrix(NA, nrow = length(net$genes), ncol = length(phalanx_attractor_indexes))
	rownames(phalanx_attractors) <- net$genes
	for(i in 1:length(phalanx_attractor_indexes)){
		atri <- phalanx_attractor_indexes[i]
		atr <- getAttractorSequence(Attractors, atri)
		st <- generateState(net, atr[1,])
		phalanx_attractors[,i] <- unlist(st)
		panl <- analizePhalanx(net, st1, verbose = FALSE, plot_folder = "none")
	}
	dim1 <- dim(phalanx_attractors)
	r1 <- dim1[1]
	r2 <- dim1[2]
	w <- (r1 / 4) 
	h <- r2 / 6
	png(file = "PhalanxAttractors.png", width = w, height = h, units = "in", res = 200)
	MatrixPlot(phalanx_attractors, "marlist" = c(2,6,6,2), "aboveaxis" = TRUE, "title" = "Phalanx Attractors")
	dev.off()
	allactive <- c()
	allinactive <- c()
	numatr <- length(phalanx_attractor_indexes)
	for(i in 1:length(net$genes)){
		print(net$genes[i])
		f <- factor(phalanx_attractors[i,], c("0", "1"))
		sf <- summary(f)
		if(sf["1"] == numatr){
			allactive <- append(allactive, net$genes[i])
		}else if(sf["0"] == numatr){
			allinactive <- append(allinactive, net$genes[i])
		}
	}
	np <- length(phalanx_environments)
	mel <- length(microenvironment_variables)
	matphalanx <- matrix(NA, nrow = np, ncol = mel)
	for(i in 1:np){
		pbin <- decimal_to_binary(phalanx_environments[i], microenvironment_variables)
		matphalanx[i,] <- pbin
	}
	colnames(matphalanx) <- microenvironment_variables
	rownames(matphalanx) <- phalanx_environments
	for(i in 1:length(microenvironment_variables)){
		print(microenvironment_variables[i])
		f <- factor(matphalanx[,i])
		print(summary(f))
	}
	results <- list("environment_matrix" = matphalanx, "attractor_matrix" = phalanx_attractors, 
			"allactive" = allactive, "allinactive" = allinactive)
	return(results)
}

analyzeECBehavior <- function(net, Attractors, sortedEnvironments, classed, microenvironment_variables, kind){
	# First get the right environments
	if(kind == "Phalanx"){
		environments <- sortedEnvironments$Phalanx
	}else if(kind == "Stalk"){
		environments <- sortedEnvironments$Stalk
	}else if(kind == "Tip"){
		environments <- sortedEnvironments$Tip
	}else if(kind == "EC"){
		environments <- sortedEnvironments$EC
	}else{
		environments <- sortedEnvironments$Other
	}
	#print(environments)
	# Build the environment matrix
	np <- length(environments)
	mel <- length(microenvironment_variables)
	environment_matrix <- matrix(NA, nrow = np, ncol = mel)
	for(i in 1:np){
		pbin <- decimal_to_binary(environments[i], microenvironment_variables)
		environment_matrix[i,] <- pbin
	}
	colnames(environment_matrix) <- microenvironment_variables
	rownames(environment_matrix) <- environments
	# Chack which molecules have fixed values
	for(i in 1:length(microenvironment_variables)){
		#print(microenvironment_variables[i])
		f <- factor(environment_matrix[,i])
		#print(summary(f))
	}
	return(environment_matrix)
}

#phalanx_behavior <- analyzePhalanxBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables)
otherm <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Other")
failedtip <- otherm[which(otherm[, "VEGFAxxxP"] == 1 | otherm[, "VEGFC_Dp"] == 1),]
	dim1 <- dim(failedtip)
	r1 <- dim1[1]
	r2 <- dim1[2]
	w <- (r2 / 2) + 4 
	h <- r1 / 6
	png(file = "failedtip.png", width = w, height = h, units = "in", res = 200)
	MatrixPlot(failedtip, "marlist" = c(2,6,6,2), "aboveaxis" = TRUE, "title" = "VEGF but no Tip")
	dev.off()
	allactiveft <- c()
	allinactiveft <- c()
	numft <- nrow(failedtip)
	for(i in 1:length(microenvironment_variables)){
		f <- factor(failedtip[,i], c("0", "1"))
		sf <- summary(f)
		print(microenvironment_variables[i])
		print(sf)
		if(sf["1"] == numft){
			allactiveft <- append(allactiveft, microenvironment_variables[i])
		}else if(sf["0"] == numft){
			allinactiveft <- append(allinactiveft, microenvironment_variables[i])
		}
	}
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 1, "ANG1" = *, "Oxygen" = *, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 0, "WNT7a" = 0, "FGF" = 0,"IGF" = *, "BMP9" = 0, 
#   "BMP10" = 0, "TGFB1" = 0, "VEGFC_D" = *, "AMPATP" = *) 2^7
# c("VEGFC_Dp" = 1, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = *, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 0, "WNT7a" = 0, "FGF" = 0,"IGF" = *, "BMP9" = 0, 
#   "BMP10" = 0, "TGFB1" = 0, "VEGFC_D" = *, "AMPATP" = *) 2^7
# c("VEGFC_Dp" = 1, "VEGFAxxxP" = 1, "ANG1" = *, "Oxygen" = *, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 0, "WNT7a" = 0, "FGF" = 0,"IGF" = *, "BMP9" = 0, 
#   "BMP10" = 0, "TGFB1" = 0, "VEGFC_D" = *, "AMPATP" = *) 2^7, 2^7 = 128, 128*3 = 384

# Now lets analyze tip cell behavior

tip_env_matrix <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Tip")
# Get the list of VEGF active environments that cause Tip cell behavior
vegf_tip_matrix <- tip_env_matrix[which(tip_env_matrix[, "VEGFAxxxP"] == 1 | tip_env_matrix[, "VEGFC_Dp"] == 1),]
# The next step is to get all the attractors that correspond to those micro-environments
vegf_tip_envs <- rownames(vegf_tip_matrix)
# Then get a list of the indexes of vegf_tip attractors
classed_indexes <- classed$indexes
vegf_tip_attractor_indexes <- c()
for(env in vegf_tip_envs){
	env_attr_indexes <- unlist(classed_indexes[env])
	vegf_tip_attractor_indexes <- append(vegf_tip_attractor_indexes, env_attr_indexes)
}
i1 <- vegf_tip_attractor_indexes[1]
vegf_tip_frame <- getAttractorSequence(Attractors, i1)
atr0 <- vegf_tip_frame[1,]
pdgfb <- 1
cyclind <- 1
vegf_tip_dividing <- c()
vegf_tip_pdgfb <- c()
if(atr0["FOXO1"] == 0 && atr0["SMAD2"] == 0){
	vegf_tip_pdgfb <- append(vegf_tip_pdgfb, atri)
}
if(atr0["Bcatenin"] == 0 || atr0["LEF1"] == 0){
	vegf_tip_dividing <- append(vegf_tip_dividing, atri)
}
fixed <- rep(1, length(net$genes))
names(fixed) <- net$genes
# Here for each attractor I want to know in addition to checking the fixed genes if it divides and secretes PDGFB
for(i in 2:length(vegf_tip_attractor_indexes)){
	atri <- vegf_tip_attractor_indexes[i]
	atr_frame <- getAttractorSequence(Attractors, atri)
	l <- nrow(atr_frame)
	pdgfb <- 1
	cyclind <- 1
	for(s in 1:l){
		st <- atr_frame[s, ]
		if(st["FOXO1"] == 0 && st["SMAD2"] == 0){
			pdgfb <- 0
		}
		if(st["Bcatenin"] == 0 || st["LEF1"] == 0){
			cyclind <- 0
		}
		for(j in 1:length(atr0)){
			if(st[j] != atr0[j]){
				fixed[j] <- 0
			}
		}
	}
	if(pdgfb == 1){
		vegf_tip_pdgfb <- append(vegf_tip_pdgfb, atri)
	}
	if(cyclind == 1){
		vegf_tip_dividing <- append(vegf_tip_dividing, atri)
	}
}
constantvt <- fixed[fixed ==1]
nconstvt <- names(constant)
const_valsvt <- atr0[nconst]

# Now lets analyze the environments without VEGFAxxxP or VEGFC_Dp that cause tip cell behavior
no_vegf_tip_matrix <- tip_env_matrix[which(tip_env_matrix[, "VEGFAxxxP"] == 0 & tip_env_matrix[, "VEGFC_Dp"] == 0),]
no_vegf_tip_summary <- summary(no_vegf_tip_matrix)
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = 1, "ShearStress" = *, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = *, "WNT7a" = *, "FGF" = *,"IGF" = 1, "BMP9" = *, 
#   "BMP10" = *, "TGFB1" = *, "VEGFC_D" = *, "AMPATP" = 0) 2^11 = 2048 > 1804
# It appears that ShearStress WNT5a and WNT7a may account for 1024 of the micro-environments (mean 0.5676 1024/1804)
# lets choose the first
no_vegfss_tip_matrix <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 1),]
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = 1, "ShearStress" = 1, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = *, "WNT7a" = *, "FGF" = *,"IGF" = 1, "BMP9" = *, 
#   "BMP10" = *, "TGFB1" = *, "VEGFC_D" = *, "AMPATP" = 0) 2^10 = 1024
# Lets explore the remaining 780 environments
no_vegfnss_tip_matrix <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 0),]
summary_no_vegfnss_tip <- summary(no_vegfnss_tip_matrix)
# Now WNT5a or WNT7a appear to account for another 512 (median 0.6564 512/780)
n_v_nsswnt5_tip_matrix <- no_vegfnss_tip_matrix[which(no_vegfnss_tip_matrix[, "WNT5a"] == 1),]
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = 1, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 1, "WNT7a" = *, "FGF" = *,"IGF" = 1, "BMP9" = *, 
#   "BMP10" = *, "TGFB1" = *, "VEGFC_D" = *, "AMPATP" = 0) 2^9 = 512
# Lets explore the remaining 268 environments
n_v_nssnw5_tip_matrix <- no_vegfnss_tip_matrix[which(no_vegfnss_tip_matrix[, "WNT5a"] == 0),]
summary_n_v_nssnw5 <- summary(n_v_nssnw5_tip_matrix)
# Now WN7a appears to account for 256/268 with a mean of 0.9552
n_v_nssnw5w7_tip_matrix <- n_v_nssnw5_tip_matrix[which(n_v_nssnw5_tip_matrix[, "WNT7a"] == 1),]
# If Oxygen, no AMPATP and (ShearStress WNT5a or WNT7a) is sufficient, the following will have 2048*(7/8) = 1792 envs
no_vegfss_orw5_or_w7_tip_matrix <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 1 | no_vegf_tip_matrix[, "WNT5a"] == 1 | no_vegf_tip_matrix[, "WNT7a"] == 1),]
# Indeed it does, now for the other 12
o12_tip_matrix <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 0 & no_vegf_tip_matrix[, "WNT5a"] == 0 & no_vegf_tip_matrix[, "WNT7a"] == 0),]
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = 1, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 0, "WNT7a" = 0, "FGF" = 1,"IGF" = 1, "BMP9" = 0, 
#   "BMP10" = 0, "TGFB1" = 0, "VEGFC_D" = *, "AMPATP" = 0) 2^4 = 16 however add ("JAGp" == 1 or "DLL4p" == 0) and we get the 12






