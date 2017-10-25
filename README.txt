## This project contains the information needed to run the analisis of our model of EC behavior during angiogenesis.
# In order to analyze the model, go to the directory that contains this file. Then run R 
#############################################################################################################################
# Then, run the script that analyzes wild type behavior by typing:
source("ECexplorer.R")

#The script runs for about 20 hours plots the main transitions, finds the attractors, sorts the attractors and stores the results in "wild_type.Rdata"
#This data can be accesed by typing:
load("wild_type.Rdata")
# The results are stored in the list "behavior"
summary(behavior)
#           Length Class          Mode    
#Network    3      BooleanNetwork list      #The network as a BooleanNetwork class used by BoolNet 
#Attractors 2      AttractorInfo  pairlist  #The list of Attractors returned by BoolNet
#Classed    3      -none-         list      #The attractors classified by micro-environment
#Sorted     4      -none-         list      #The microenvironments sorted according to the EC behavior they induce
#EnvCount   4      -none-         numeric   #The number of microenvironments that induce each EC behavior

summary(behavior$Classed)
#                  Length Class  Mode
#types             65536  -none- list # For each microenvironment, contains the behavior that corresponds to each attractor
#indexes           65536  -none- list # For each microenvironment, contains the index of each attractor in Attractors$atractors
#attractor_lengths 65536  -none- list # For each microenvironment, contains the size of each attractor
################################################################################################################################
# After that, we group the microenvironments and summarize the characteristics of each group
source("ECbehaviorExplorer.R") # Loads the data from "wild_type.Rdata", takes about 4 hours to run
# Saves the results in "grouped_environments.Rdata" in the list "grouped_environments"
# Load the data
load("grouped_environments.Rdata")
# Contains two lists
summary(grouped_environments)
#                   Length Class  Mode
#grouped_by_behavior 4      -none- list # The micro-environments grouped by the EC behavior they induce
#grouped_by_division 3      -none- list # The micro-environments grouped by the EC division
> summary(grouped_environments$grouped_by_behavior)
#         Length Class  Mode
#Phalanx  2      -none- list
#Stalk    3      -none- list
#Tip      3      -none- list
#Atypical 9      -none- list
> summary(grouped_environments$grouped_by_behavior$Tip)
#       Length Class  Mode The name of the groups correspond to those on Table 1 in the article
#TipI   2      -none- list
#TipII  2      -none- list
#TipIII 2      -none- list
> summary(grouped_environments$grouped_by_behavior$Phalanx)
#                   Length Class  Mode   
#Environment Matrix 1536   -none- numeric # A matrix containing the micro-environments as rows
#Analisis              9   -none- list # Some information about the behavior induced by the group 
> summary(grouped_environments$grouped_by_behavior$Phalanx$Analisis)
#                     Length Class      Mode     
#all_divide            0     -none-     NULL     
#some_divide           0     -none-     NULL     
#none_divide          96     -none-     character
#all_pdgfb             0     -none-     NULL     
#some_pdgfb            0     -none-     NULL     
#none_pdgfb           96     -none-     character
#number_of_attractors  1     -none-     numeric  
#attractor_lengths    96     -none-     numeric  
#fixed                47     data.frame list  
##############################################################################################################
# To analize the robustness to molecular activation noise:
>source("RobustnessAnalisis.R") # Our current script uses 1000 initial random states to run faster
# The results are stored in robustness.Rdata, which contains a list "robustness".
>load("robustness.Rdata")
##############################################################################################################
# Last, to simulate all single gain and loss of function mutations: It will take a long time to run 
# (About 128 days/Available threads), it is paralelized
source("MutantECexplorer.R")
# Basically the same as running ECexplorer.r for each mutant, stores the result in (mutantdescription).Rdata

