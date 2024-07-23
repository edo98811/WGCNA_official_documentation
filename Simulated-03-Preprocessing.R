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


meanExpressionByArray=apply( datExpr,1,mean, na.rm=T)  
NumberMissingByArray=apply( is.na(data.frame(datExpr)),1, sum)  


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(9, 5)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:50), cex.names = 0.7)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Keep only arrays containing less than 500 missing entries
KeepArray= NumberMissingByArray<500
table(KeepArray)
datExpr=datExpr[KeepArray,]
y=y[KeepArray]
ArrayName[KeepArray]


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


NumberMissingByGene =apply( is.na(data.frame(datExpr)),2, sum)
# One could do a barplot(NumberMissingByGene), but the barplot is empty in this case.
# It may be better to look at the numbers of missing samples using the summary method:
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)
# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)
datExpr=datExpr[, KeepGenes]
GeneName=GeneName[KeepGenes]


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


sizeGrWindow(9, 5)
plotClusterTreeSamples(datExpr=datExpr, y=y)


