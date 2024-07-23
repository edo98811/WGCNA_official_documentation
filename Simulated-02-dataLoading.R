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


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datGeneSummary=read.csv("GeneSummaryTutorial.csv")
datTraits=read.csv("TraitsTutorial.csv")
datMicroarrays=read.csv("MicroarrayDataTutorial.csv")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# This vector contains the microarray sample names
ArrayName= names(data.frame(datMicroarrays[,-1]))
 # This vector contains the gene names
GeneName= datMicroarrays$GeneName
# We transpose the data so that the rows correspond to samples and the columns correspond to genes
# Since the first column contains the gene names, we exclude it.
datExpr=data.frame(t(datMicroarrays[,-1]))
names(datExpr)=datMicroarrays[,1]
dimnames(datExpr)[[1]]=names(data.frame(datMicroarrays[,-1]))
#Also, since we simulated the data, we know the true module color:
truemodule= datGeneSummary$truemodule
rm(datMicroarrays)
collectGarbage()


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# First, make sure that the array names in the file datTraits line up with those in the microarray data 
table( dimnames(datExpr)[[1]]==datTraits$ArrayName)
# Next, define the microarray sample trait 
y=datTraits$y


