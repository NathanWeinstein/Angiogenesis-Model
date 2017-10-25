# Execute from R by typing "source("ECexplorer.R")"

source("ECexplorerTools.R")

## First verify that the models are in the working directory
WD <- getwd()
simple_model <- paste(WD, "angiogenesis.txt", sep = "/")
if(file.exists(simple_model)){
}else{
	print("The file containing the simplified model is not in the working directory")
}
detailed_model <- paste(WD, "angiofull.net", sep = "/")
if(file.exists(simple_model)){
}else{
	print("The file containing the detailed model is not in the working directory")
} 
## Set the names of the plot folders, if they do not exist, thay will be created by the program.
plot_folder <- paste(WD, "PlotFolder/", sep = "/") 
detailed_plot_folder <- paste(WD, "DetailedPlotFolder/", sep = "/")

## This function verifies and plots the transitions, finds the attractors, classifies them and sumarizes them.
# It takes anout 20 hours to run, spends most of the time finding the attractors
behavior <- exploreEC(simple_model, detailed_model, plot_folder, detailed_plot_folder)
save(behavior, file = "wild_type.Rdata")






