# DEseq2 filter
library(DESeq2)
library(dplyr)
options(stringsAsFactors=F)

setwd("<YOUR_WORK_DIR>")
group <- read.delim('group.txt', header = TRUE, sep = '\t', row.names = 1,check.names = FALSE)
group$group <- factor(group$time)
data <- read.delim('counts.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
data <- data[rowSums(data) > 10, ]
data1 <- data.frame(gene = rownames(data), data)
dds <-DESeqDataSetFromMatrix(countData=data1, 
                             colData=group, 
                             design=~1,
                             tidy=TRUE)
Mat_normal <- vst(dds, blind = FALSE) %>% assay()
write.csv(Mat_normal,file= "WGCNA_input.csv", row.names=T)

# WGCNA
# step 1: filtration
library(WGCNA)
library(reshape2)
library(stringr)
library(GO.db)
library(doParallel)
options(stringsAsFactors = FALSE)
setwd("<YOUR_WORK_DIR>")
dataExpr<-read.csv("WGCNA_input.csv", header = T,row.names = 1) 
dim(dataExpr)
head(dataExpr)[,1:8]
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))

# step 2: detect missing value
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)

# step 3: sample clustering and detect outliers
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

# step 4: calculate soft threshold
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5,
                          networkType = "signed")
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
power = sft$powerEstimate

# step 5: Step-by-step construction of the gene network and identification of modules
enableWGCNAThreads(10)
adjacency = adjacency(dataExpr, power = power)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize = 50
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
table(mergedColors)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# step 6: Associate with phenotype
traitData<-read.csv("traits.csv", row.names=1,header=T,comment.char = "",check.names=F)
allTraits = traitData
dim(allTraits)
names(allTraits)

fpkmSamples = rownames(dataExpr)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits)
collectGarbage()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrixWithStars <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
textMatrixWithStars[moduleTraitPvalue <= 0.0001] <- "****"
textMatrixWithStars[moduleTraitPvalue <= 0.001 & moduleTraitPvalue > 0.0001] <- "***"
textMatrixWithStars[moduleTraitPvalue <= 0.01 & moduleTraitPvalue > 0.001] <- "**"
textMatrixWithStars[moduleTraitPvalue <= 0.05 & moduleTraitPvalue > 0.01] <- "*"

textMatrix <- matrix(nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
for (i in 1:nrow(moduleTraitCor)) {
  for (j in 1:ncol(moduleTraitCor)) {
    textMatrix[i, j] <- paste(signif(moduleTraitCor[i, j], 2), "\n", textMatrixWithStars[i, j], sep = "")
  }
}

par(mar = c(5, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(c("#0da9ce", "white", "#e74a32"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))

# step 7: calculate ModuleMembership (MM) and geneTraitSignificance (GS)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1)
      #sizeGrWindow(7, 7)
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

names(dataExpr)
probes = names(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}

geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "Rsults_GS_and_MM.xls",sep="\t",row.names=F)