# Execute from R by typing "source("angiofull.R")"
# We will use BoolNet
if(!require(BoolNet)){
	install.packages("BoolNet")
}
library(BoolNet)
source("matrixPlotter.R")


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



# A function to classify a molecular activation pattern
sortingHat <- function(state) {
	if (state["AKT"] == 1) {
		stable <- "TRUE"
	}else if (state["FOXO1"] == 1 || state["SMAD2"] == 1){
		stable <- "PDGF-B"
	}else{
		stable <- "FALSE"
	}
	if (state["AKT"] == 0 && state["NRP1"] == 1 && state["DLL4a"] == 1) {
		fate <- "Tip"
	}else if (state["NRP1"] == 0 && state["JAGa"] == 1){ 
		fate <- "Stalk"
	}else if (state["AKT"] == 1 && state["NRP1"] == 0 && state["JAGa"] == 0){
		fate <- "Phalanx"
	}else{
		fate <- "Atypical"
	}
	celltype <- c(stable, fate)
	names(celltype) <- c("Stable", "Fate")
	return(celltype)
}

# A function to get only the rows that correspond to the attractor reached from a state
getAttractor <- function(net, state, plot_path_file = "none") {
	path <- getPathToAttractor(net, state, includeAttractorStates = "all")
	dim1 <- dim(path)
	r1 <- dim1[1]
	r2 <- dim1[2]
	path0 <- getPathToAttractor(net, state, includeAttractorStates = "none")
	if(plot_path_file != "none"){
		w <- (r1 / 2) + 5
		h <- r2 / 4
		png(file = plot_path_file, width = w, height = h, units = "in", res = 100)
		plotSequence("sequence" = path)
		dev.off()
	}
	dim0 <- dim(path0)
	r0 <- dim0[1]
	return(path[(r0 + 1):r1,])
}

# A function to analyze an atractor
analyzeAttractor <- function(net,  attr){
	analisis <- data.frame(Stable = character(nrow(attr)), 
			       Fate = character(nrow(attr)),
			       stringsAsFactors = FALSE)
	for (i in 1 : nrow(attr)){
		description <- sortingHat(attr[i,])
		analisis[i, "Stable"] <- description["Stable"]
		analisis[i, "Fate"] <- description["Fate"]
	}
	return(analisis)
}

# A function to decide the type of an attractor
sortAtractor <- function(atractorAnalisis){
	nr <- nrow(atractorAnalisis)
	Fate <- atractorAnalisis[,"Fate"]
	fateFactor <- factor(Fate, c("Tip", "Stalk", "Phalanx", "Atypical"))
	s <- summary(fateFactor)
	if(s["Tip"] == nr){
		kind <- "Tip"
	}else if(s["Stalk"] == nr){
		kind <- "Stalk"
	}else if(s["Phalanx"] == nr){
		kind <- "Phalanx"
	}else{
		kind <- "Atypical"
	}
	return(kind)
}

# A function to analyze a path
analizePath <- function(net, state){
	attr <- getAttractor(net, state)
	analisis <- analyzeAttractor(net,  attr)
	return(analisis)
}

# A function to analize the main transitions in a certain context
analizeTransitions <- function(net, context, verbose = FALSE, plot_folder = "none"){
	analisis <- c("Phalanx->Tip" = 0, "Phalanx->TipCD" = 0, "Phalanx->Stalk" = 0, "Tip->Phalanx" = 0,"Tip->Stalk" = 0, "Stalk->Phalanx" = 0, "Stalk->Tip" = 0, "Stalk->TipD4" = 0)
	fixed <- net$fixed
	if(length(context) > 0){
		namelist <- names(context)
		for(gene in namelist){
			if((fixed[gene] != -1) & (fixed[gene] != context[gene])){
				print("Context incompatible with fixed gene value:")
				print(gene)
				print(fixed[gene])
				return(-1)
			}
		}
	}else{
		context <- append(context, c("VEGFAxxxP" = 0))
	}
	# Check if the plot_folder exists, if not create it
	if(plot_folder != "none"){
		if(file.exists(plot_folder)){
		}else{
			dir.create(plot_folder)
		}
	}
	# Check if the context by itself leads to an EC.
	EC <- generateState(net, context)
	if(plot_folder == "none"){
		attrEC <- getAttractor(net, EC)
	}else{
		plot_file <- paste(plot_folder, "pathToEC", sep = "")
		attrEC <- getAttractor(net = net, state = EC, plot_path_file = plot_file)
	}
	ECanalisis <- analyzeAttractor(net,  attrEC)
	initialState <- generateState(net, unlist(attrEC[1,]))
	# EC to Tip
	# Now add some VEGFAxxxP and we should get a tip cell. 
	initTip <- initialState
	if(fixed["VEGFAxxxP"] != 0){
		initTip["VEGFAxxxP"] <- 1
	}
	if(fixed["DLL4p"] != 0){ 
		initTip["DLL4p"] <- 1
	}
	if(fixed["JAGp"] != 0){
		initTip["JAGp"] <- 1
	}
	if(plot_folder == "none"){
		attrTip <- getAttractor(net, initTip)
	}else{
		plot_file <- paste(plot_folder, "ECtoTip", sep = "")
		attrTip <- getAttractor(net = net, state = initTip, plot_path_file = plot_file)
	}
	TipAnalisis <- analyzeAttractor(net,  attrTip)
	if(sortAtractor(TipAnalisis) == "Tip"){
		analisis["Phalanx->Tip"] <- 1
	}
	# Now add some VEGFC_Dp and we should get a tip cell. 
	initTip1 <- initialState
	if(fixed["VEGFC_Dp"] != 0){
		initTip1["VEGFC_Dp"] <- 1
	}
	if(fixed["DLL4p"] != 0){ 
		initTip1["DLL4p"] <- 1
	}
	if(fixed["JAGp"] != 0){
		initTip1["JAGp"] <- 1
	}
	if(plot_folder == "none"){
		attrTip1 <- getAttractor(net, initTip1)
	}else{
		plot_file <- paste(plot_folder, "ECtoTipVEGFC_Dp", sep = "")
		attrTip1 <- getAttractor(net = net, state = initTip1, plot_path_file = plot_file)
	}
	Tip1Analisis <- analyzeAttractor(net,  attrTip1)
	if(sortAtractor(Tip1Analisis) == "Tip"){
		analisis["Phalanx->TipCD"] <- 1
	}
	# EC to Stalk
	# With NOTCH WNT and TGF signalling
	initStalk <- initialState
	if(fixed["DLL4p"] != 0){ 
		initStalk["DLL4p"] <- 1
	}
	if(fixed["WNT5a"] != 0){ 
		initStalk["WNT5a"] <- 1 # WNT11 miight also work but not WNT7a
	}
	if(fixed["TGFB1"] != 0){ 
		initStalk["TGFB1"] <- 1
	}
	if(plot_folder == "none"){
		attrStalk <- getAttractor(net, initStalk)
	}else{
		plot_file <- paste(plot_folder, "ECtoStalk", sep = "")
		attrStalk <- getAttractor(net = net, state = initStalk, plot_path_file = plot_file)
	}
	StalkAnalisis <- analyzeAttractor(net,  attrStalk)
	if(sortAtractor(StalkAnalisis) == "Stalk"){
		analisis["Phalanx->Stalk"] <- 1
	}
	# Tip from initial state
	pretipcontext <- context
	if(fixed["VEGFAxxxP"] != 0){
		pretipcontext <- append(pretipcontext, c("VEGFAxxxP" = 1))
	}
	if(fixed["DLL4p"] != 0){
		pretipcontext <- append(pretipcontext, c("DLL4p" = 1))
	}
	if(fixed["JAGp"] != 0){
		pretipcontext <- append(pretipcontext, c("JAGp" = 1))
	}
	preTip <- generateState(net, pretipcontext)
	attrTip <- getAttractor(net, preTip)
	initTip <- generateState(net, unlist(attrTip[1, ]))
	# Is tip cell behavior reveersible?
	# Now from tipical tip cell to EC.
	tipToStable <- initTip
	if(fixed["VEGFAxxxP"] != 1){ 
		tipToStable["VEGFAxxxP"] <- 0
	}
	if(fixed["JAGp"] != 1){ 
		tipToStable["JAGp"] <- 0
	}
	if(fixed["DLL4p"] != 1){ 
		tipToStable["DLL4p"] <- 0
	}
	if(plot_folder == "none"){
		attrTip2 <- getAttractor(net, generateState(net, tipToStable))
	}else{
		plot_file <- paste(plot_folder, "TipToEC", sep = "")
		attrTip2 <- getAttractor(net = net, state = generateState(net, tipToStable), plot_path_file = plot_file)
	} 
	Tip2Analisis <- analyzeAttractor(net,  attrTip2)
	if(sortAtractor(Tip2Analisis) == "Phalanx"){ # sortAtractor(Tip2Analisis) == "EC" || sortAtractor(Tip2Analisis) == "Phalanx")
		analisis["Tip->Phalanx"] <- 1
	}
	# Now from tipical tip cell to stalk cell.
	tipToStalk <- initTip
	if(fixed["VEGFAxxxP"] != 1){ 
		tipToStalk ["VEGFAxxxP"] <- 0
	}
	if(fixed["DLL4p"] != 0){ 
		tipToStalk ["DLL4p"] <- 1
	}
	if(fixed["JAGp"] != 1){ 
		tipToStalk ["JAGp"] <- 0
	}
	if(fixed["WNT5a"] != 0){ 
		tipToStalk ["WNT5a"] <- 1 # WNT11 would activate the same pathways.
	}
	if(fixed["TGFB1"] != 0){
		tipToStalk ["TGFB1"] <- 1
	}
	if(plot_folder == "none"){
		attrTip3 <- getAttractor(net, tipToStalk)
	}else{
		plot_file <- paste(plot_folder, "TipToStalk", sep = "")
		attrTip3 <- getAttractor(net = net, state = tipToStalk, plot_path_file = plot_file)
	}
	Tip3Analisis <- analyzeAttractor(net,  attrTip3)
	if(sortAtractor(Tip3Analisis) == "Stalk"){
		analisis["Tip->Stalk"] <- 1
	}
	# Stalk from initial state
	prestalkcontext <- context
	if(fixed["DLL4p"] != 0){
		prestalkcontext <- append(prestalkcontext, c("DLL4p" = 1))
	}
	if(fixed["WNT5a"] != 0){
		prestalkcontext <- append(prestalkcontext, c("WNT5a" = 1))
	}
	if(fixed["TGFB1"] != 0){
		prestalkcontext <- append(prestalkcontext, c("TGFB1" = 1))
	}
	preStalk <- generateState(net, prestalkcontext)
	attrStalk1 <- getAttractor(net, preStalk)
	initStalk <- generateState(net, unlist(attrStalk1[1, ]))
	# Is the secondary fate reveersible?
	# Now from tipical tip cell to Phalanx.
	stalkToStable <- initStalk
	if(fixed["DLL4p"] != 1){
		stalkToStable["DLL4p"] <- 0
	}
	if(fixed["WNT5a"] != 1){
		stalkToStable["WNT5a"] <- 0
	}
	if(fixed["TGFB1"] != 1){
		stalkToStable["TGFB1"] <- 0
	}
	if(plot_folder == "none"){
		attrStalk2 <- getAttractor(net, stalkToStable)
	}else{
		plot_file <- paste(plot_folder, "StalkToEC", sep = "")
		attrStalk2 <- getAttractor(net = net, state = stalkToStable, plot_path_file = plot_file)
	}
	Stalk2Analisis <- analyzeAttractor(net,  attrStalk2)
	if(sortAtractor(Stalk2Analisis) == "Phalanx"){#sortAtractor(Stalk2Analisis) == "EC" || sortAtractor(Stalk2Analisis) == "Phalanx")
		analisis["Stalk->Phalanx"] <- 1
	}
	# Now from tipical stalk cell to tip cell.
	stalkToTip <- initStalk
	if(fixed["VEGFAxxxP"] != 0){
		stalkToTip["VEGFAxxxP"] <- 1
	}
	if(fixed["JAGp"] != 0){
		stalkToTip["JAGp"] <- 1
	}
	if(fixed["WNT5a"] != 1){
		stalkToTip["WNT5a"] <- 0
	}
	if(fixed["TGFB1"] != 1){
		stalkToTip["TGFB1"] <- 0
	}
	if(plot_folder == "none"){
		attrStalk3 <- getAttractor(net, stalkToTip)
	}else{
		plot_file <- paste(plot_folder, "StalkToTip", sep = "")
		attrStalk3 <- getAttractor(net = net, state = stalkToTip, plot_path_file = plot_file)
	}
	Stalk3Analisis <- analyzeAttractor(net,  attrStalk3)
	if(sortAtractor(Stalk3Analisis) == "Tip"){
		analisis["Stalk->Tip"] <- 1
	}
	# Now from tipical stalk cell to tip cell by loosing DLL4p.
	stalkToTip1 <- initStalk
	if(fixed["VEGFAxxxP"] != 0){
		stalkToTip1["VEGFAxxxP"] <- 1
	}
	if(fixed["DLL4p"] != 1){
		stalkToTip1["DLL4p"] <- 0
	}
	if(plot_folder == "none"){
		attrStalk4 <- getAttractor(net, stalkToTip1)
	}else{
		plot_file <- paste(plot_folder, "StalkToEC", sep = "")
		attrStalk4 <- getAttractor(net = net, state = stalkToTip1, plot_path_file = plot_file)
	}
	Stalk4Analisis <- analyzeAttractor(net,  attrStalk4)
	if(sortAtractor(Stalk4Analisis) == "Tip"){
		analisis["Stalk->TipD4"] <- 1
	}
	if(verbose){
		print("Endothelial cell")
		print(ECanalisis)
		print("Tip cell")
		print(TipAnalisis)
		print("Tip cell with VEGFC_Dp")
		print(Tip1Analisis)
		print("DLL4 + TGF + WNT -> Stalk cell")
		print(StalkAnalisis)
		print("Tip -> Stable EC")
		print(Tip2Analisis)
		print("Tip -> Stalk")
		print(Tip3Analisis)
		print("Stalk -> Stable EC")
		print(Stalk2Analisis)
		print("Stalk -> Tip JAGp=1")
		print(Stalk3Analisis)
		print("Stalk -> Tip DLL4=0")
		print(Stalk4Analisis)
	}
	return(analisis)
}

# A function to analize the main transitions from a phalanx cell state (fixed point attractor)
analizePhalanx <- function(net, phalanx, verbose = FALSE, plot_folder = "none"){
	analisis <- c("Phalanx->Tip" = 0, "Phalanx->TipCD" = 0, "Phalanx->Stalk" = 0)
	fixed <- net$fixed
	# Check if the plot_folder exists, if not create it
	if(plot_folder != "none"){
		if(file.exists(plot_folder)){
		}else{
			dir.create(plot_folder)
		}
	}
	# EC to Tip
	# Now add some VEGFAxxxP and we should get a tip cell. 
	initTip <- phalanx
	if(fixed["VEGFAxxxP"] != 0){
		initTip["VEGFAxxxP"] <- 1
	}
	if(fixed["DLL4p"] != 0){ 
		initTip["DLL4p"] <- 1
	}
	if(fixed["JAGp"] != 0){
		initTip["JAGp"] <- 1
	}
	attrTip <- getAttractor(net, initTip)
	if(plot_folder == "none"){
		attrTip <- getAttractor(net, initTip)
	}else{
		plot_file <- paste(plot_folder, "PhalanxtoTip", sep = "")
		attrTip <- getAttractor(net = net, state = initTip, plot_path_file = plot_file)
	}
	TipAnalisis <- analyzeAttractor(net,  attrTip)
	if(sortAtractor(TipAnalisis) == "Tip"){
		analisis["Phalanx->Tip"] <- 1
	}
	# Now add some VEGFC_Dp and we should get a tip cell. 
	initTip1 <- phalanx
	if(fixed["VEGFC_Dp"] != 0){
		initTip1["VEGFC_Dp"] <- 1
	}
	if(fixed["DLL4p"] != 0){ 
		initTip1["DLL4p"] <- 1
	}
	if(fixed["JAGp"] != 0){
		initTip1["JAGp"] <- 1
	}
	if(plot_folder == "none"){
		attrTip1 <- getAttractor(net, initTip1)
	}else{
		plot_file <- paste(plot_folder, "PhalanxtoTipVEGFC_Dp", sep = "")
		attrTip1 <- getAttractor(net = net, state = initTip1, plot_path_file = plot_file)
	}
	Tip1Analisis <- analyzeAttractor(net,  attrTip1)
	if(sortAtractor(Tip1Analisis) == "Tip"){
		analisis["Phalanx->TipCD"] <- 1
	}
	# EC to Stalk
	# With NOTCH WNT and TGF signalling
	initStalk <- phalanx
	if(fixed["DLL4p"] != 0){ 
		initStalk["DLL4p"] <- 1
	}
	if(fixed["WNT5a"] != 0){ 
		initStalk["WNT5a"] <- 1 # WNT11 miight also work but not WNT7a
	}
	if(fixed["TGFB1"] != 0){ 
		initStalk["TGFB1"] <- 1
	}
	if(plot_folder == "none"){
		attrStalk <- getAttractor(net, initStalk)
	}else{
		plot_file <- paste(plot_folder, "PhalanxtoStalk", sep = "")
		attrStalk <- getAttractor(net = net, state = initStalk, plot_path_file = plot_file)
	}
	StalkAnalisis <- analyzeAttractor(net,  attrStalk)
	if(sortAtractor(StalkAnalisis) == "Stalk"){
		analisis["Phalanx->Stalk"] <- 1
	}

	if(verbose){
		print("Tip cell")
		print(TipAnalisis)
		print("Tip cell with VEGFC_Dp")
		print(Tip1Analisis)
		print("DLL4 + TGF + WNT -> Stalk cell")
		print(StalkAnalisis)
	}
	celltypes <- list("Phalanx->Tip" = attrTip, "Phalanx->TipCD" = attrTip1, "Phalanx->Stalk" = attrStalk)
	results = list("celltypes"=celltypes, "analysis" = analisis)
	return(results)
}

# A function to test the effect of all single gain and loss of function mutations on EC behavior
testMutants <- function(net, verbose = FALSE, plot_folder = "none"){
	genes <- net$genes
	testvector <- c()
	rnames <- c()
	for (gene in genes){
		netl <- fixGenes(net, gene, 0)
		fixedl <- netl$fixed
		if(length(fixedl[fixedl == 1]) == 0){
			context <- c()
		}else{
			context <- fixedl[fixedl == 1]
		}
		# context1101 <- append(mincontext, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
		if(fixedl["ANG1"] != 0){
			context <- append(context, c("ANG1" = 1))
		}
		if(fixedl["Oxygen"] != 0){
			context <- append(context, c("Oxygen" = 1))
		}
		if(fixedl["ShearStress"] != 0){
			context <- append(context, c("ShearStress" = 1))
		}
		analisisl <- analizeTransitions(netl, context)
		namel <- paste(gene, "lf", sep = " ")
		if(verbose){
			print(namel)
			print(analisisl)
		}
		rnames <- append(rnames, namel)
		testvector <- append(testvector,analisisl)
		netg <- fixGenes(net, gene, 1)
		fixedg <- netg$fixed
		if(length(fixedg[fixedg == 1]) == 0){
			contextg <- c()
		}else{
			contextg <- fixedg[fixedg == 1]
		}
		# context1101 <- append(mincontext, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
		if(fixedg["ANG1"] != 0){
			contextg <- append(contextg, c("ANG1" = 1))
		}
		if(fixedg["Oxygen"] != 0){
			contextg <- append(contextg, c("Oxygen" = 1))
		}
		if(fixedg["ShearStress"] != 0){
			contextg <- append(contextg, c("ShearStress" = 1))
		}
		analisisg <- analizeTransitions(netg, contextg)
		nameg <- paste(gene, "gf", sep = " ")
		if(verbose){
			print(nameg)
			print(analisisg)
		}
		rnames <- append(rnames, nameg)
		testvector <- append(testvector,analisisg)
	}
	mutantMatrix <- matrix(testvector, nrow = 2*length(genes), byrow = TRUE)
	colnames(mutantMatrix) <- names(analisisg)
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
		w <- r2 + 1 
		h <- r1 / 4
		plot_file <- paste(plot_folder, "mutantEffect", sep = "")
		png(file = plot_file, width = w, height = h, units = "in", res = 100)
		MatrixPlot(mutantMatrix, "marlist" = c(2,6,6,2), "aboveaxis" = TRUE, "title" = "Mutation effects on EC behavior")
		dev.off()
	}
	return(mutantMatrix)
}

summarizeMutants <- function(mutantMatrix){
	resume <- list()
	cnames <- colnames(mutantMatrix)
	mutations <- rownames(mutantMatrix)
	for(i in c(1:ncol(mutantMatrix))){
		effectiveMutations <- c()
		for(j in c(1:nrow(mutantMatrix))){
			if(mutantMatrix[j,i] == 0){
				effectiveMutations <- append(effectiveMutations, mutations[j])
			}
		}
		resume[[i]] <- effectiveMutations
	}
	names(resume) <- cnames
	return(resume)
}
 
compareRules <- function(){
	net <- loadNetwork("~/Desktop/BoolNet/angiofull.txt")
	nets <- loadNetwork("~/Desktop/BoolNet/angiogenesis.txt")
	geness <- nets$genes
	interactions <- net$interactions
	interactionss <- nets$interactions
	for(gene in geness){
		print(gene)
		print("Large network")
		print(interactions[gene])
		print("Small network")
		print(interactionss[gene])
	}
}

exploreECmicroenvironment <- function(net, verbose = FALSE, plot_folder = "none"){
	fixed <- net$fixed
	mincontext <- fixed[fixed == 1]
	results <- c()
	for(i in c(0:15)){
 		bin <- (as.integer(intToBits(i)))
 		bins <- bin[1:4]
		names(bins) <- c("ANG1","Oxygen", "AMPATP", "ShearStress")
		context <- append(mincontext, bins)
		analisis <- append(bins, analizeTransitions(net, context))
		results <- append(results, analisis)
 	}
	resultmatrix <- matrix(results, nrow = length(results)/length(analisis), byrow = TRUE)
	colnames(resultmatrix) <- names(analisis)
	if(verbose){
		print(resultmatrix)
	}
	# Check if the plot_folder exists, if not create it
	if(plot_folder != "none"){
		if(file.exists(plot_folder)){
		}else{
			dir.create(plot_folder)
		}
		dim1 <- dim(resultmatrix)
		r1 <- dim1[1]
		r2 <- dim1[2]
		w <- r2 
		h <- r1 / 4
		plot_file <- paste(plot_folder, "microenvironmentEffect", sep = "")
		png(file = plot_file, width = w, height = h, units = "in", res = 100)
		MatrixPlot(resultmatrix, "title" = "Effect of the microenvironment on EC behavior")
		dev.off()
	}
}


binary_to_decimal <- function(x){
	if (!(all(x %in% c(0,1)))){stop("Not a binary sequence")}
	decimal <- sum(x*2^((length(x):1) -1))
	return(decimal)
}

decimal_to_binary <- function(y, microenvironment_variables){
	len <- length(microenvironment_variables)
	bin <- (as.integer(intToBits(y)))
 	bins <- bin[1:len]
	bins <- rev(bins)
	names(bins) <- microenvironment_variables
	return(bins)
}

# This only considers the inputs for the smaller EC network!
getEnvironmentNumber <- function(state, microenvironment_variables){
	environment_number <- 0
	microenvironment <- state[microenvironment_variables]
	#print(microenvironment)
	environment_number <- binary_to_decimal(microenvironment)
	return(environment_number)
}

classifyAttractors <- function(attractor_info, net, microenvironment_variables){
	types_by_environment <- vector(mode = "list", length = 2^length(microenvironment_variables)) # Creates a list of vectors
	indexes_by_environment <- vector(mode = "list", length = 2^length(microenvironment_variables)) # Creates a list of vectors
	attractor_lengths_by_env <- vector(mode = "list", length = 2^length(microenvironment_variables)) # Creates a list of vectors
	envnames <- c(0:((2^length(microenvironment_variables))-1))
	names(types_by_environment) <- envnames
	names(indexes_by_environment) <- envnames
	attractors <- attractor_info$attractors
	for(i in 1:length(attractors)){
		attractor <- getAttractorSequence(attractor_info, i)
		analisis <- analyzeAttractor(net,  attractor)
		sort <- sortAtractor(analisis)
		state <- attractor[1,]
		envnumber <- getEnvironmentNumber(state, microenvironment_variables)
		types_by_environment[[envnumber + 1]] <- append(types_by_environment[[envnumber + 1]], sort)
		indexes_by_environment[[envnumber + 1]] <- append(indexes_by_environment[[envnumber + 1]], i)
		attractor_lengths_by_env[[envnumber + 1]] <- append(attractor_lengths_by_env[[envnumber + 1]], nrow(attractor))
	}
	classed <- list("types" = types_by_environment, "indexes" = indexes_by_environment, "attractor_lengths" = attractor_lengths_by_env)
	return(classed)
}

sortEnvironments <- function(attractors_by_environment){
	tip <- c()
	stalk <- c()
	phalanx <- c()
	atypical <- c() 
	aenvnames <- names(attractors_by_environment)
	for(i in 1:length(attractors_by_environment)){
		aenv <- attractors_by_environment[[i]]
		nenv <- aenvnames[i]
		laenv <- length(aenv)
		if(length(aenv[aenv == "Tip"]) == laenv){
			tip <- append(tip, nenv)
		}else if(length(aenv[aenv == "Stalk"]) == laenv){
			stalk <- append(stalk, nenv)
		}else if(length(aenv[aenv == "Phalanx"]) == laenv){
			phalanx <- append(phalanx, nenv)
		}else{
			atypical <- append(atypical, nenv)
		}
	}
	sortedEnvironments <- list(Tip = tip, Stalk = stalk, Phalanx = phalanx, Atypical = atypical)
	return(sortedEnvironments)
}

exploreEC <- function(simple_model, detailed_model, plot_folder, detailed_plot_folder){
	### This part is to test that the model behaves as expected. 
	### We where willing to fith the model to preserve the expected EC behavior transitions.
	### However the boolean update rules as written based on the molecular basis and 
	### the model simplification procedures allready behaved as expected.
	# Load the simplified network from a file
	net <- loadNetwork(simple_model)
	fixed <- net$fixed
	# Explore the role of the microenvironment
	mincontext <- fixed[fixed == 1]
	context1101 <- append(mincontext, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
	# analizeTransitons simulates the EC behavior transitions that have been reported in the literature.
	analisis1101 <- analizeTransitions(net, context1101, verbose = FALSE, plot_folder = plot_folder)
	exploreECmicroenvironment(net, plot_folder = plot_folder)
	# Test the effect of all gain and loss of function mutants on the most relevant EC behavior transitions.
	mutantMatrix <- testMutants(net, plot_folder = plot_folder)
	mutantSummary <- summarizeMutants(mutantMatrix)
	microenvironment_variables <- c("VEGFC_Dp", "VEGFAxxxP", "ANG1", "Oxygen", "ShearStress", "JAGp", "DLL4p", "WNT5a", "WNT7a", "FGF","IGF", "BMP9", "BMP10", "TGFB1", "VEGFC_D", "AMPATP")
	### Now we analyze the dynamic behavior of our model in depth.
	# First, we obtain all the attractors. This takes about 20 hours to run
print(1)
	Attractors <- getAttractors(net, type = "synchronous", method = "sat.exhaustive", maxAttractorLength = 36)
print(2)
	# classifyAttractors builds three lists of length 2^16, one slot for each microenvironment
	# classed$types contains for each micro-environment the kind of EC behavior represented
	# classed$indexes contains the corresponding index of each attractor in the list of attractors
	# classed$attractor_lengths contains the sizes of the corresponding attractors
	classed <- classifyAttractors(Attractors, net, microenvironment_variables)
	# sortEnvironments sorts the environments according to the EC behavior represented by the attractors that 
	# are reachable in that micro-environment, if some of the attractors in a micro-environment
	# represent atypical behaviors or the behaviors are heterogeneous, the environment is classified as Atypical.
	sortedEnvironments <- sortEnvironments(classed$types)
	Tip <- length(sortedEnvironments$Tip)
	Stalk <- length(sortedEnvironments$Stalk)
	Phalanx <- length(sortedEnvironments$Phalanx)
	Atypical <- length(sortedEnvironments$Atypical)
	envcount <- c("Tip" = Tip, "Stalk" = Stalk, "Phalanx" = Phalanx, "Atypical" = Atypical)
	behavior <- list("Network" = net, "Attractors" = Attractors,"Classed" = classed, "Sorted" = sortedEnvironments, "EnvCount" = envcount)
	# Load the full network from a file
	netf <- loadNetwork(detailed_model)
	fixedf <- netf$fixed
	# Explore the role of the microenvironment
	mincontextf <- fixedf[fixedf == 1]
	context1101f <- append(mincontextf, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
	analisis1101f <- analizeTransitions(netf, context1101f, verbose = FALSE, plot_folder = detailed_plot_folder)
	exploreECmicroenvironment(netf, plot_folder = detailed_plot_folder)
	# Test the effect of all gain and loss of function mutants
	mutantMatrixf <- testMutants(netf, plot_folder = detailed_plot_folder)
	mutantSummaryf <- summarizeMutants(mutantMatrixf)
	return(behavior)
}








