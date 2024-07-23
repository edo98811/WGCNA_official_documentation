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
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the previously saved data
load("Simulated-RelatingToExt.RData"); 


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


colorlevels=unique(colorh1)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1], 
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     main=whichmodule, 
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
table(FilterGenes)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


dimnames(data.frame(datExpr))[[2]][FilterGenes]


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


sizeGrWindow(8,6)
par(mfrow=c(2,2))
# We choose 4 modules to plot: turquoise, blue, brown, green. 
# For simplicity we write the code out explicitly for each module.
which.color="turquoise"; 
restrictGenes=colorh1==which.color 
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes], 
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color, 
                   xlab="Intramodular Connectivity", 
                   ylab="(Module Membership)^6")

which.color="blue"; 
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

which.color="brown"; 
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

which.color="green";
restrictGenes=colorh1==which.color 
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes], 
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color, 
                   xlab="Intramodular Connectivity", 
                   ylab="(Module Membership)^6")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


NS1=networkScreening(y=y, datME=datME, datExpr=datExpr,
         oddPower=3, blockSize=1000, minimumSampleSize=4,
         addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# network screening analysis
mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=100])
# standard analysis based on the correlation p-values (or Student T test)
mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=100]) 


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


topNumbers=c(10,20,50,100)
for (i in c(1:length(topNumbers)) ) 
{
  print(paste("Proportion of noise genes in the top", topNumbers[i], "list"))
  WGCNApropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=topNumbers[i]])
  StandardpropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=topNumbers[i]])
  print(paste("WGCNA, proportion of noise=", WGCNApropNoise, 
        ", Standard, prop. noise=", StandardpropNoise))
  if (WGCNApropNoise< StandardpropNoise) print("WGCNA wins")
  if (WGCNApropNoise==StandardpropNoise) print("both methods tie")
  if (WGCNApropNoise>StandardpropNoise) print("standard screening wins")
} 


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


rm(dissTOM); collectGarbage()


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


#Form a data frame containing standard and network screening results
CorPrediction1=data.frame(GS1,NS1$cor.Weighted)
cor.Weighted=NS1$cor.Weighted
# Plot the comparison
sizeGrWindow(8, 6)
verboseScatterplot(cor.Weighted, GS1,
         main="Network-based weighted correlation versus Pearson correlation\n",
         col=truemodule, cex.main = 1.2)
abline(0,1)


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


set.seed(2)
nSamples2=2000
MEgreen=rnorm(nSamples2)
scaledy2=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(nSamples2)
y2=ifelse( scaledy2>median(scaledy2),2,1)
MEturquoise= ESturquoise*scaledy2+sqrt(1-ESturquoise^2)*rnorm(nSamples2)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(nSamples2)
MEbrown= ESbrown*scaledy2+sqrt(1-ESbrown^2)*rnorm(nSamples2)
MEyellow= ESyellow*scaledy2+sqrt(1-ESyellow^2)*rnorm(nSamples2)
# Put together a data frame of eigengenes
ModuleEigengeneNetwork2=data.frame(y=y2,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)
# Simulate the expression data
dat2=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork2$MEturquoise,
   MEblue=ModuleEigengeneNetwork2$MEblue,MEbrown=ModuleEigengeneNetwork2$MEbrown,
   MEyellow=ModuleEigengeneNetwork2$MEyellow,
   MEgreen=ModuleEigengeneNetwork2$MEgreen,simulateProportions=simulateProportions1, 
   nGenes=nGenes1)
# recall that this is the signed gene significance in the training data
GS1= as.numeric(cor(y, datExpr, use="p"))
# the following is the signed gene significance in the test data
GS2=as.numeric( cor(ModuleEigengeneNetwork2$y, dat2$datExpr, use="p"))


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


sizeGrWindow(8,6)
par(mfrow=c(1,1))
verboseScatterplot(GS1,GS2,
       main="Trait-based gene significance in test set vs. training set\n",
       xlab = "Training set gene significance",
       ylab = "Test set gene significance",
       col=truemodule, cex.main = 1.4)


#=====================================================================================
#
#  Code chunk 15
#
#=====================================================================================


EvaluationGeneScreening1 = corPredictionSuccess(
           corPrediction = CorPrediction1, 
           corTestSet=GS2,
           topNumber=seq(from=20, to=500, length=30) )
par(mfrow=c(2,2))
listcomp = EvaluationGeneScreening1$meancorTestSetOverall
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main="Predicting positive and negative correlations",
        ylab="mean cor, test data", 
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreening1$meancorTestSetPositive
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main="Predicting positive correlations",
        ylab="mean cor, test data", 
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreening1$meancorTestSetNegative
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main = "Predicting negative correlations",
        ylab = "mean cor, test data", 
        xlab = "top number of genes in the training data")


#=====================================================================================
#
#  Code chunk 16
#
#=====================================================================================


relativeCorPredictionSuccess(corPredictionNew = NS1$cor.Weighted,
                             corPredictionStandard = GS1, 
                             corTestSet=GS2,
                             topNumber=c(10,20,50,100,200,500) )


#=====================================================================================
#
#  Code chunk 17
#
#=====================================================================================


# Create a data frame holding the results of gene screening
GeneResultsNetworkScreening=data.frame(GeneName=row.names(NS1), NS1)
# Write the data frame into a file
write.table(GeneResultsNetworkScreening, file="GeneResultsNetworkScreening.csv",
row.names=F,sep=",")
# Output of eigengene information:
datMEy = data.frame(y, datME)
eigengeneSignificance = cor(datMEy, y);
eigengeneSignificance[1,1] = (1+max(eigengeneSignificance[-1, 1]))/2
eigengeneSignificance.pvalue = corPvalueStudent(eigengeneSignificance, nSamples = length(y))
namesME=names(datMEy)
# Form a summary data frame
out1=data.frame(t(data.frame(eigengeneSignificance,
eigengeneSignificance.pvalue, namesME, t(datMEy))))
# Set appropriate row names
dimnames(out1)[[1]][1]="EigengeneSignificance"
dimnames(out1)[[1]][2]="EigengeneSignificancePvalue"
dimnames(out1)[[1]][3]="ModuleEigengeneName"
dimnames(out1)[[1]][-c(1:3)]=dimnames(datExpr)[[1]]
# Write the data frame into a file
write.table(out1, file="MEResultsNetworkScreening.csv", row.names=TRUE, col.names = TRUE, sep=",")
# Display the first few rows:
head(out1)


#=====================================================================================
#
#  Code chunk 18
#
#=====================================================================================


# Write out gene information
GeneName=dimnames(datExpr)[[2]]
GeneSummary=data.frame(GeneName, truemodule, SignalGeneIndicator,  NS1)
write.table(GeneSummary, file="GeneSummaryTutorial.csv", row.names=F,sep=",")
# here we output the module eigengenes and trait y without eigengene significances
datTraits=data.frame(ArrayName, datMEy)
dimnames(datTraits)[[2]][2:length(namesME)]=paste("Trait",  
                                             dimnames(datTraits)[[2]][2:length(namesME)], 
                                             sep=".")
write.table(datTraits, file="TraitsTutorial.csv", row.names=F,sep=",")
rm(datTraits)
# here we output the simulated gene expression data
MicroarrayData=data.frame(GeneName, t(datExpr))
names(MicroarrayData)[-1]=ArrayName
write.table(MicroarrayData, file="MicroarrayDataTutorial.csv", row.names=F,sep=",")
rm(MicroarrayData)


#=====================================================================================
#
#  Code chunk 19
#
#=====================================================================================


# Perform network screening
NS1GS=networkScreeningGS(datExpr=datExpr, datME = datME, GS=GS1)
# Organize its results for easier plotting
GSprediction1=data.frame(GS1,NS1GS$GS.Weighted)
GS.Weighted=NS1GS$GS.Weighted
# Plot a comparison between standard gene significance and network-weighted gene significance
sizeGrWindow(8, 6)
par(mfrow=c(1,1))
verboseScatterplot(GS1, GS.Weighted, 
                   main="Weighted gene significance vs. the standard GS\n",
                   col=truemodule)
abline(0,1)


#=====================================================================================
#
#  Code chunk 20
#
#=====================================================================================


EvaluationGeneScreeningGS = corPredictionSuccess(corPrediction=GSprediction1, corTestSet=GS2,
                                    topNumber=seq(from=20, to=500, length=30) )
sizeGrWindow(8, 6)
par(mfrow=c(2,2))
listcomp= EvaluationGeneScreeningGS$meancorTestSetOverall
matplot(x=listcomp$topNumber,
        y=listcomp[,-1],
        main="Predicting positive and negative correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreeningGS$meancorTestSetPositive
matplot(x=listcomp$topNumber, 
        y=listcomp[,-1], 
        main="Predicting positive correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreeningGS$meancorTestSetNegative
matplot(x=listcomp$topNumber,
        y=listcomp[,-1],
        main="Predicting negative correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")


#=====================================================================================
#
#  Code chunk 21
#
#=====================================================================================


collectGarbage()
save.image("Simulated-Screening.RData")


