# Execute from R by typing "source("angiofull.R")"
# We will use BoolNet
#joshej04
library(BoolNet)

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
	if (state["TIE2"] == 1) {
		stable <- "TRUE"
	}else if (state["PDGFB"] == 1){
		stable <- "PDGF-B"
	}else{
		stable <- "FALSE"
	}
	if (state["NPR1"] == 1 && state["STAT3"] == 1 && state["DLL4a"] == 1 && state["HEY1"] == 0) {
		fate <- "Tip"
	}else if(state["HEY1"] == 1 && state["JAGa"] == 1){
		fate <- "Stalk"
	}else{
		fate <- "EC"
	}
	celltype <- c(stable, fate)
	names(celltype) <- c("Stable", "Fate")
	return(celltype)
}

# A function to get only the rows that correspond to the attractor reached from a state
getAttractor <- function(net, state) {
	path <- getPathToAttractor(net, state, includeAttractorStates = "all")
	dim1 <- dim(path)
	r1 <- dim1[1]
	r2 <- dim1[2]
	path0 <- getPathToAttractor(net, state, includeAttractorStates = "none")
	dim0 <- dim(path0)
	r0 <- dim0[1]
	path[(r0 + 1):r1,]
}

# A function to analyze an atractor
analyzeAttractor <- function(net,  attr){
	analisis <- data.frame(Stable = logical(nrow(attr)), 
			       Fate = character(nrow(attr)),
			       stringsAsFactors = FALSE)
	for (i in 1 : nrow(attr)){
		description <- sortingHat(attr[i,])
		analisis[i, "Stable"] <- description["Stable"]
		analisis[i, "Fate"] <- description["Fate"]
	}
	return(analisis)
}

# A function to analyze a path
analizePath <- function(net, state){
	attr <- getAttractor(net, state)
	analisis <- analyzeAttractor(net,  attr)
	return(analisis)
}

# Load the network from a file
net <- loadNetwork("~/Desktop/BoolNet/angiofull.txt")

# Generate a state that represents an EC with no external stimuli.
bEC <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "Oxygen" = 1))
# Get the path from there to an attractor
pathbEC <- getPathToAttractor(net, bEC, includeAttractorStates = "all")
pathbEC0 <- getPathToAttractor(net, bEC, includeAttractorStates = "none")
attrbEC <- getAttractor(net, bEC)
bECanalisis <- analyzeAttractor(net,  attrbEC)
print("EC with normal oxygen")
print(bECanalisis)
png(file = "EC_normal_oxygen.png", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathbEC)
dev.off()


# Generate a state that represents a stable EC with ANG1.
sEC <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "Oxygen" = 1, "ANG1" = 1))
# Get the path from there to an attractor
pathsEC <- getPathToAttractor(net, sEC, includeAttractorStates = "all")
pathsEC0 <- getPathToAttractor(net, sEC, includeAttractorStates = "none")
attrsEC <- getAttractor(net, sEC)
sECanalisis <- analyzeAttractor(net,  attrsEC)
print("Stable EC with ANG1")
print(sECanalisis)
png(file = "Stable_EC_with_ANG1", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathsEC)
dev.off()

# Generate a state that represents an EC in a hypoxic environment (not realistic in a blood vessel).
hEC <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "ANG1" = 1))
# Get the path from there to an attractor
pathhEC <- getPathToAttractor(net, hEC, includeAttractorStates = "all")
pathhEC0 <- getPathToAttractor(net, hEC, includeAttractorStates = "none")
attrhEC <- getAttractor(net, hEC)
hECanalisis <- analyzeAttractor(net,  attrhEC)
print("EC in a hypoxic environment")
print(hECanalisis)
png(file = "EC_in_a_hypoxic_environment", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathhEC)
dev.off()

# Generate a state that represents an EC in a low energy environment (not realistic in a blood vessel).
leEC <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "Oxygen" = 1, "AMPATP" = 1, "ANG1" = 1))
# Get the path from there to an attractor
pathleEC <- getPathToAttractor(net, leEC, includeAttractorStates = "all")
pathleEC0 <- getPathToAttractor(net, leEC, includeAttractorStates = "none")
attrleEC <- getAttractor(net, leEC)
leECanalisis <- analyzeAttractor(net,  attrleEC)
print("EC in a low energy environment")
print(leECanalisis)
png(file = "EC_in_a_low_energy_environment", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathleEC)
dev.off()

# Generate a relevant initial state for a quiescent EC affected by shear stress and ANG1.
qEC <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "ShearStress" = 1, "Oxygen" = 1, "ANG1" = 1))
# Get the path from there to an attractor
pathqEC <- getPathToAttractor(net, qEC, includeAttractorStates = "all")
pathqEC0 <- getPathToAttractor(net, qEC, includeAttractorStates = "none")
attrqEC <- getAttractor(net, qEC)
qECanalisis <- analyzeAttractor(net,  attrqEC)
print("EC induced by shear stress and ANG1")
print(qECanalisis)
png(file = "EC_induced_by_shear_stress_and_ANG1", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathqEC)
dev.off()

# The attractor reached by a quiescent endothelial cell includes inactive ANG2-mediated vessel 
# destabilization.
initialState <- attrqEC[1,]

# Now add some VEGFAxxxP and we should get a tip cell. 
initTip <- initialState
initTip["VEGFAxxxP"] <- 1
pathTip <- getPathToAttractor(net, initTip, includeAttractorStates = "all")
attrTip <- getAttractor(net, initTip)
TipAnalisis <- analyzeAttractor(net,  attrTip)
print("Tip cell")
print(TipAnalisis)
png(file = "EC_to_Tip", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip)
dev.off()

# Now add some VEGFC_Dp and we should get a tip cell. 
initTipCD <- initialState
initTipCD["VEGFC_Dp"] <- 1
pathTipCD <- getPathToAttractor(net, initTipCD, includeAttractorStates = "all")
attrTipCD <- getAttractor(net, initTipCD)
TipAnalisisCD <- analyzeAttractor(net,  attrTipCD)
print("Tip cell VEGFC_Dp")
print(TipAnalisisCD)
png(file = "EC_to_Tip_VEGFC_Dp", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTipCD)
dev.off()

# In the tip cell TGF activity is inhibited.
initTip0 <- initialState
initTip0["VEGFAxxxP"] <- 1
initTip0["TGFB1"] <- 1
pathTip0 <- getPathToAttractor(net, initTip0, includeAttractorStates = "all")
attrTip0 <- getAttractor(net, initTip0)
Tip0Analisis <- analyzeAttractor(net,  attrTip0)
print("Tip cell with TGF ligands")
print(Tip0Analisis)
png(file = "Tip_Inhibits_TGF", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip0)
dev.off()


# Now what happens when VEGF and Notch active.
initTip1 <- initialState
initTip1["VEGFAxxxP"] <- 1
initTip1["DLL4p"] <- 1
pathTip1 <- getPathToAttractor(net, initTip1, includeAttractorStates = "all")
attrTip1 <- getAttractor(net, initTip1)
Tip1Analisis <- analyzeAttractor(net,  attrTip1)
print("VEGFA and DLL4")
print(Tip1Analisis)
png(file = "VEGFA_and_DLL4", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip1)
dev.off()

# Now what happens when VEGF, DLL4, and JAG1p are active.
initTip2 <- initialState
initTip2["VEGFAxxxP"] <- 1
initTip2["DLL4p"] <- 1
initTip2["JAGp"] <- 1
pathTip2 <- getPathToAttractor(net, initTip2, includeAttractorStates = "all")
attrTip2 <- getAttractor(net, initTip2)
Tip2Analisis <- analyzeAttractor(net,  attrTip2)
print("VEGFA JAG and DLL4, tipical tip cell specification")
print(Tip2Analisis)
png(file = "VEGFAp_JAG_and_DLL4", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip2)
dev.off()

# Now do we get a stalk cell?. No, it may require Wnt and TGF signalling
initStalk <- initialState
initStalk["DLL4p"] <- 1
pathStalk <- getPathToAttractor(net, initStalk, includeAttractorStates = "all")
attrStalk <- getAttractor(net, initStalk)
StalkAnalisis <- analyzeAttractor(net,  attrStalk)
print("DLL4")
print(StalkAnalisis)
png(file = "DLL4", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathStalk)
dev.off()

# With Wnt and TGF signalling
initStalk1 <- initialState
initStalk1["DLL4p"] <- 1
initStalk1["WNT5a"] <- 1 # WNT11 miight also work but not WNT7a
initStalk1["TGFB1"] <- 1
pathStalk1 <- getPathToAttractor(net, initStalk1, includeAttractorStates = "all")
attrStalk1 <- getAttractor(net, initStalk1)
Stalk1Analisis <- analyzeAttractor(net,  attrStalk1)
print("DLL4 + TGF + WNT -> Stalk cell")
print(Stalk1Analisis)
png(file = "DLL4_+TGF_+_WNT", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathStalk1)
dev.off()

print("Checking the reversibility of tip cell behavior")
# Tip from initial state
preTip <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "ShearStress" = 1, "Oxygen" = 1, "ANG1" = 1, "VEGFAxxxP" = 1, "DLL4p" = 1, "JAGp" = 1))
attrTip3 <- getAttractor(net, preTip)
checktip <- (attrTip3 == attrTip2)
fchecktip <- factor(checktip)
summarycheck <- summary(fchecktip) # We do get the same attractor as attrTip2
initTip3 <- attrTip3[1, ]

# Is the primary fate reveersible?
# Now from tipical tip cell to stalk cell.
print("Tip -> Stable EC")
tipToStable <- initTip3
tipToStable["VEGFAxxxP"] <- 0
tipToStable["JAGp"] <- 0
tipToStable["DLL4p"] <- 0 
#print(initTip3)
pathTip3 <- getPathToAttractor(net, tipToStable, includeAttractorStates = "all")
attrTip3 <- getAttractor(net, tipToStable)
Tip3Analisis <- analyzeAttractor(net,  attrTip3)
print(Tip3Analisis)
png(file = "Tip_to_Stable_EC", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip3)
dev.off()


# Now from tipical tip cell to stalk cell.
print("Tip -> Stalk")
tipToStalk <- initTip3
tipToStalk ["VEGFAxxxP"] <- 0
tipToStalk ["DLL4p"] <- 1
tipToStalk ["WNT5a"] <- 1 # WNT11 miight also work but not WNT7a
tipToStalk ["TGFB1"] <- 1
#print(initTip3)
pathTip4 <- getPathToAttractor(net, tipToStalk, includeAttractorStates = "all")
attrTip4 <- getAttractor(net, tipToStalk)
Tip4Analisis <- analyzeAttractor(net,  attrTip4)
print(Tip4Analisis)
png(file = "Tip_to_Stalk", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathTip4)
dev.off()

print("Checking the reversibility of stalk cell behavior")
# Stalk from initial state
preStalk <- generateState(net, c("ACVR2A" = 1, "BMPRII" = 1, "TGFBRII" = 1, "SMAD4" = 1, "sGC" = 1, "ADAM10" = 1, "gSecretase" = 1, "ShearStress" = 1, "Oxygen" = 1, "ANG1" = 1, "DLL4p" = 1, "WNT5a" = 1, "TGFB1" = 1))
attrStalk2 <- getAttractor(net, preStalk)
checkstalk <- (attrStalk2 == attrStalk1)
fcheckstalk <- factor(checkstalk)
summarycheck1 <- summary(fcheckstalk) # We do get the same attractor as attrTip2
# print(summarycheck1) this one is slightly different perhaps only the start state of the atractor?
initStalk2 <- attrStalk2[1, ]

# Is the primary fate reveersible?
# Now from tipical tip cell to stalk cell.
print("Stalk -> Stable EC")
stalkToStable <- initStalk2
stalkToStable["DLL4p"] <- 0
stalkToStable["WNT5a"] <- 0
stalkToStable["TGFB1"] <- 0
pathStalk3 <- getPathToAttractor(net, stalkToStable, includeAttractorStates = "all")
attrStalk3 <- getAttractor(net, stalkToStable)
Stalk3Analisis <- analyzeAttractor(net,  attrStalk3)
print(Stalk3Analisis)
png(file = "Stalk_to_EC", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathStalk3)
dev.off()


# Now from tipical tip cell to stalk cell.
print("Stalk -> Tip")
stalkToTip <- initStalk2
stalkToTip["VEGFAxxxP"] <- 1
stalkToTip["JAGp"] <- 1
stalkToTip["WNT5a"] <- 0
stalkToTip["TGFB1"] <- 0
#print(initTip3)
pathStalk4 <- getPathToAttractor(net, stalkToTip, includeAttractorStates = "all")
attrStalk4 <- getAttractor(net, stalkToTip)
Stalk4Analisis <- analyzeAttractor(net,  attrStalk4)
print(Stalk4Analisis)
png(file = "Stalk_to_Tip", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathStalk4)
dev.off()

# Now from tipical tip cell to stalk cell.
print("Stalk -> Tip DLL4")
stalkToTip1 <- initStalk2
stalkToTip1["VEGFAxxxP"] <- 1
stalkToTip1["DLL4p"] <- 0
stalkToTip1["WNT5a"] <- 0
stalkToTip1["TGFB1"] <- 0
#print(initTip3)
pathStalk5 <- getPathToAttractor(net, stalkToTip1, includeAttractorStates = "all")
attrStalk5 <- getAttractor(net, stalkToTip1)
Stalk5Analisis <- analyzeAttractor(net,  attrStalk5)
print(Stalk5Analisis)
png(file = "Stalk_to_TipDLL4", width = 12, height = 30, units = "in", res = 300)
plotSequence("sequence" = pathStalk5)
dev.off()

# getAttractors takes too long on this network!
# attractorinfo0 <- getAttractors(net, 
#				type = "synchronous", 
#				method = "sat.exhaustive", 
#				maxAttractorLength = 36,
#				genesON = c("ANG1", "Oxygen", "ShearStress", "AMPATP"), 
#				genesOFF = c())






