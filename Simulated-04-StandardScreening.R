#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory.  On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the previously saved data
load("Simulated-dataSimulation.RData");
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


GS1= as.numeric(cor(y, datExpr, use="p"))
# Network terminology: GS1 will be referred to as signed gene significance measure
p.Standard=corPvalueFisher(GS1, nSamples =length(y) )
# since the q-value function has problems with missing data, we use the following trick
p.Standard2=p.Standard
p.Standard2[is.na(p.Standard)]=1
q.Standard=qvalue(p.Standard2)$qvalues
# Form a data frame to hold the results
StandardGeneScreeningResults=data.frame(GeneName,PearsonCorrelation=GS1, p.Standard, q.Standard)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


NoiseGeneIndicator=is.element( truemodule, c("turquoise", "blue", "yellow", "grey"))+.0
SignalGeneIndicator=1-NoiseGeneIndicator


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


table(q.Standard<.20)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


mean(NoiseGeneIndicator[q.Standard<=0.20]) 


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


save.image(file = "Simulated-StandardScreening.RData")


