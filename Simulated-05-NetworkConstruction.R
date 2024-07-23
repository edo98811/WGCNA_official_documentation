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
# Load additional necessary packages
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the previously saved data
load("Simulated-StandardScreening.RData"); 
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^6
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6) 
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


datExpr=datExpr[, rank(-k,ties.method="first" )<=3600]


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


dissTOM=TOMdist(ADJ1)
collectGarbage()


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pam4$clustering, truemodule)
table(pam5$clustering, truemodule)
table(pam6$clustering, truemodule)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


pamTOM4=pam(as.dist(dissTOM), 4)
pamTOM5=pam(as.dist(dissTOM), 5)
pamTOM6=pam(as.dist(dissTOM), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pamTOM4$clustering, truemodule)
table(pamTOM5$clustering, truemodule)
table(pamTOM6$clustering, truemodule)


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


hierADJ=hclust(as.dist(dissADJ), method="average" )
# Plot the resulting clustering tree together with the true color assignment
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03, 
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule, colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ, 
                              cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

# Plot results of all module detection methods together:
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ, 
                    colors=data.frame(truemodule, colorStaticADJ, 
                                     colorDynamicADJ, colorDynamicHybridADJ), 
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1 
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                       deepSplit=2, pamRespectsDendro = FALSE))
# Now we plot the results
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM, 
               colors=data.frame(truemodule, colorStaticTOM, 
                                 colorDynamicTOM, colorDynamicHybridTOM), 
               dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
               main = "Gene dendrogram and module colors, TOM dissimilarity")


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


tabStaticADJ=table(colorStaticADJ,truemodule)
tabStaticTOM=table(colorStaticTOM,truemodule)
tabDynamicADJ=table(colorDynamicADJ, truemodule)
tabDynamicTOM=table(colorDynamicTOM,truemodule)
tabDynamicHybridADJ =table(colorDynamicHybridADJ,truemodule)
tabDynamicHybridTOM =table(colorDynamicHybridTOM,truemodule)


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


randIndex(tabStaticADJ,adjust=F)
randIndex(tabStaticTOM,adjust=F)
randIndex(tabDynamicADJ,adjust=F)
randIndex(tabDynamicTOM,adjust=F)
randIndex(tabDynamicHybridADJ ,adjust=F)
randIndex(tabDynamicHybridTOM ,adjust=F)


#=====================================================================================
#
#  Code chunk 15
#
#=====================================================================================


colorh1= colorDynamicHybridTOM
# remove the dissimilarities, adjacency matrices etc to free up space
rm(ADJ1); rm(dissADJ);              
collectGarbage()
save.image("Simulated-NetworkConstruction.RData")


