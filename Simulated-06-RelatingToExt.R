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
load("Simulated-NetworkConstruction.RData"); 
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME, use="p"), 2)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


sizeGrWindow(8,9)
plotMEpairs(datME,y=y)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


signif(cor(datME, ModuleEigengeneNetwork1[,-1]),2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise"; 
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
         clabels=T,rcols=which.module,
         title=which.module )
# for the second (blue) module we use
which.module="blue";  
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
         clabels=T,rcols=which.module,
         title=which.module )
which.module="brown"; 
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
         clabels=T,rcols=which.module,
         title=which.module )


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


sizeGrWindow(8,7);
which.module="green"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
         nrgcols=30,rlabels=F,rcols=which.module,
         main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


signif(cor(y,datME, use="p"),2)


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


cor.test(y, datME$MEbrown)


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


p.values = corPvalueStudent(cor(y,datME, use="p"), nSamples = length(y))


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, colorh1, mean, na.rm=T)


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1)


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


collectGarbage()
save.image("Simulated-RelatingToExt.RData")


