# Execute from R by typing "source("angiofull.R")"
# We will use BoolNet
library(BoolNet)
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


# Load the simplified network from a file
net <- loadNetwork("~/Desktop/BoolNet/angiogenesis.txt")
fixed <- net$fixed
# Explore the role of the microenvironment
mincontext <- fixed[fixed == 1]
context1101 <- append(mincontext, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
analisis1101 <- analizeTransitions(net, context1101, verbose = FALSE, plot_folder = "~/Desktop/BoolNet/PlotFolder/")
exploreECmicroenvironment(net, plot_folder = "~/Desktop/BoolNet/PlotFolder/")
# Test the effect of all gain and loss of function mutants
mutantMatrix <- testMutants(net, plot_folder = "~/Desktop/BoolNet/PlotFolder/")
mutantSummary <- summarizeMutants(mutantMatrix)
microenvironment_variables <- c("VEGFC_Dp", "VEGFAxxxP", "ANG1", "Oxygen", "ShearStress", "JAGp", "DLL4p", "WNT5a", "WNT7a", "FGF","IGF", "BMP9", "BMP10", "TGFB1", "VEGFC_D", "AMPATP")
#Attractors <- getAttractors(net,
#				type = "synchronous",
#				method = "sat.exhaustive",
#				maxAttractorLength = 36)
#classed <- classifyAttractors(Attractors, net, microenvironment_variables)
sortedEnvironments <- sortEnvironments(classed$types)
Tip <- length(sortedEnvironments$Tip)
Stalk <- length(sortedEnvironments$Stalk)
Phalanx <- length(sortedEnvironments$Phalanx)
EC <- length(sortedEnvironments$EC)
Other <- length(sortedEnvironments$Other)
envcount <- c("Tip" = Tip, "Stalk" = Stalk, "Phalanx" = Phalanx, "EC" = EC, "Other" = Other)
behavior <- list("Attractors" = Attractors,"Classed" = classed, "Sorted" = sortedEnvironments, "EnvCount" = envcount)
#classedAttractors <- getAttractorsForMicroenvironments(net, microenvironment_variables)

# Load the full network from a file
netf <- loadNetwork("~/Desktop/BoolNet/angiofull.txt")
fixedf <- netf$fixed
# Explore the role of the microenvironment
mincontextf <- fixedf[fixedf == 1]
context1101f <- append(mincontextf, c("ANG1" = 1,"Oxygen" = 1, "ShearStress" = 1))
analisis1101f <- analizeTransitions(netf, context1101f, verbose = FALSE, plot_folder = "~/Desktop/BoolNet/PlotFolderFull/")
exploreECmicroenvironment(netf, plot_folder = "~/Desktop/BoolNet/PlotFolderFull/")
# Test the effect of all gain and loss of function mutants
mutantMatrixf <- testMutants(netf, plot_folder = "~/Desktop/BoolNet/PlotFolderFull/")
mutantSummaryf <- summarizeMutants(mutantMatrixf)








