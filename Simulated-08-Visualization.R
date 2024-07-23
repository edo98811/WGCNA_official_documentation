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
# Load WGCNA package
library(WGCNA)
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the previously saved data
load("Simulated-RelatingToExt.RData"); 
load("Simulated-Screening.RData")


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1),  main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


power=6
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "TOM heatmap plot, module genes" )


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


power=6
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-adjacency( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "Adjacency heatmap plot, module genes" )


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sizeGrWindow(7,7)
topList=rank(NS1$p.Weighted,ties.method="first")<=30
gene.names= names(datExpr)[topList]
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=FALSE,
                   power=1, main="signed correlations")


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


sizeGrWindow(7,7)
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


sizeGrWindow(7,7)
# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=TRUE,
                   power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="D. TOM in an unsigned network")


