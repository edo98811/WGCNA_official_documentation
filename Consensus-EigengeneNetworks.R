#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Basic settings: we work with two data sets
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "Consensus-NetworkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Create a variable weight that will hold just the body weight of mice in both sets
weight = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
  weight[[set]] = list(data = as.data.frame(Traits[[set]]$data$weight_g));
  names(weight[[set]]$data) = "weight"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, weight));


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(8,10);
pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off();


