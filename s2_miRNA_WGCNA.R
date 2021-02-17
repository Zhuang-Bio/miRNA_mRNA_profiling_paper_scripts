library(WGCNA)
library(doParallel)
library(multtest)
library(flashClust)
library(tidyverse)
allowWGCNAThreads(10)
enableWGCNAThreads(10)
rm(list = ls())

setwd("C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/s1_wb_miRNA_WGCNA")
load("../DEresults_wb_miRNAExpression_new.RData")
load("Results_wb_miRNA_WGCNA_analysis_NEW.RData")

####Step1: read data####
#read miRNA expression
WBmiRNA <- allmiRNAdata_wb$Ori_miRNA_readCount_TPM %>% column_to_rownames(var = "miRNA")
#miRNA expressed in at least 10 samples
keep=which(rowSums(WBmiRNA[,1:20] > 5) >= 10)
temp=WBmiRNA[keep, ] %>% select(21:40)
colnames(temp) <- gsub("_TPM","",colnames(temp))
colnames(temp) <- gsub("SK","Skin",colnames(temp))
colnames(temp)[-1] <- gsub("W1","Wound1",colnames(temp)[-1])
colnames(temp) <- gsub("W7","Wound7",colnames(temp))
##Sample Info
subname <- gsub("[0-9]$","",colnames(temp))
datTraits <- data.frame(samples=names(temp), subtype=subname)
rownames(datTraits) <- datTraits[,1]
datTraits$subtype <- parse_factor(datTraits$subtype)
str(datTraits)

datExpr <- t(temp)
dim(datExpr);rm(temp)

####Step2: optimize and choosea a set of soft-thresholding powers####
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.8,
                         corFnc="bicor",corOptions = list(use = 'p', maxPOutliers = 0.1),
                         verbose = 5, blockSize = 5000, networkType = "signed")
#beta value
best_beta <- sft$powerEstimate;best_beta #> best_beta 
# Plot the results:
pdf(file = "wb_allmiRNA_TPM_SoftThreshold.pdf")
par(mfrow = c(1,2))
cex1 <-0.85
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, Signed R^2",type="n",
     main = paste("Scale independence"),axes=FALSE,ylim=c(-0.5,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
axis(side = 2, at=seq(-0.5,1.2,by=0.3))
axis(side=1, at=seq(0,20,by=5))
box()
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

####Step3: One-step to construct co-expression netwrok and modle detection####
net <- blockwiseModules(datExpr, power = 18,
                        maxBlockSize = 20000,networkType = "signed",
                        corFnc="bicor",corOptions = list(use = 'p', maxPOutliers = 0.1), pearsonFallback="individual",
                        minModuleSize = 10, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,checkMissingData = FALSE,
                        verbose = 5)
##module detection and visualization
table(net$colors) #check the module number
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
table(mergedColors)
rm(mergedColors)

#load("temp.RData")
##calculate the weighted adjacency matrix
adjacency <- adjacency(datExpr,corFnc="bicor",corOptions = list(use = 'p', pearsonFallback="individual", maxPOutliers = 0.1, robustX=FALSE),
                       type = "signed",power = 18)
dim(adjacency)
rownames(adjacency)=colnames(adjacency)=colnames(datExpr)
#adjacency1 <- adjacency %>% as_tibble()
#rownames(adjacency1) <- colnames(adjacency1)
#data.table::fwrite(adjacency1,file = "wb_allmiRNA_TPM_adjacency.txt",row.names = T,quote = F,sep = "\t")
##relabel modules: grey means no mapped module genes/circRNAs
modules=as.data.frame(table(net$colors))
colnames(modules)=c("Label", "N") 
modules$Label=paste("M", modules$Label, sep="")
modules$Color=c("grey",labels2colors(modules$Label[-1]))
write.table(modules, file="wb_allmiRNA_TPM_Color_Instruction_for_module.txt",row.names = F,quote = F,sep = "\t")
moduleLabel=paste("M",net$colors, sep="") 
moduleColor=modules$Color[match(moduleLabel, modules$Label)] 

##add the traits
str(datTraits)
design= model.matrix( ~ -1 + datTraits$subtype) #We obtain an identity matrix if we remove the intercept from the model
colnames(design) = levels(datTraits$subtype)
rownames(design) <- datTraits$samples 
#Try with individual samples
design= model.matrix( ~ -1 + datTraits$samples)
colnames(design) = datTraits$samples
rownames(design) <- datTraits$samples

Skin <- as.data.frame(design[,2]);names(Skin) <- "Skin"
Wound1 <- as.data.frame(design[,3]);names(Wound1) <- "Wound1"
Wound7 <- as.data.frame(design[,4]);names(Wound7) <- "Wound7"
CW <- as.data.frame(design[,1]);names(CW) <- "CW"
GS.Skin <- as.numeric(cor(datExpr, Skin, use = "p"))
GS.Skincolor <- numbers2colors(GS.Skin,signed = T)
GS.Wound1 <- as.numeric(cor(datExpr,Wound1, use = "p"))
GS.Wound1color <- numbers2colors(GS.Wound1,signed = T)
GS.Wound7 <- as.numeric(cor(datExpr,Wound7, use = "p"))
GS.Wound7color <- numbers2colors(GS.Wound7,signed = T)
GS.CW <- as.numeric(cor(datExpr,CW, use = "p"))
GS.CWcolor <- numbers2colors(GS.CW,signed = T)
datColor <- data.frame(moduleColor ,GS.Skincolor, GS.Wound1color, GS.Wound7color, GS.CWcolor)[net$blockGenes[[1]],]
# Plot the dendrogram and the module colors underneath
pdf("wb_allmiRNA_TPM_dendrogram_modulesAddTraits.pdf",height = 7,width = 8,useDingbats = F)
plotDendroAndColors(net$dendrograms[[1]], colors = datColor,
                    groupLabels = c("Module colors","Skin","Wound1","Wound7","CW"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs_try <- moduleEigengenes(datExpr, moduleColor)

####Step4: Modules related to Traits####
#define the miRNA number and sample number
nGenes <- ncol(datExpr);nGenes
nSamples <- nrow(datExpr);nSamples
me <- data.frame(rownames(datExpr), net$MEs,stringsAsFactors = F)
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me)[1]="Sample"
#write the moduleEigengene value for whole biopsy
write.table(me, "wb_allmiRNA_TPM_ME_ModuleEigenes.txt",sep = "\t",quote = F,row.names = F)
#MEs0 = moduleEigengenes(datExpr, moduleColor)$eigengenes
MEs = orderMEs(me[,-1])
##draw the plot of module eigengene
which.module <- "M9"
module.color <- modules[modules$Label==which.module,]$Color
MEplot <- MEs[, which.module]
names(MEplot) <- rownames(MEs)
MEplot <- MEplot[c(6:20,1:5)]
max(MEplot);min(MEplot)
barplot(MEplot, col=module.color,ylim = c(-0.6,0.4))

pdf(file = paste0("wb_allmiRNA_Modules_EigengenePlot_",which.module,".pdf"),useDingbats = F,width = 8,height = 5)
barplot(MEplot, col=module.color,ylim = c(-0.6,0.4))
dev.off()

moduleTraitCor = cor(MEs, design , use = "p")#, robustY=FALSE, maxPOutliers = 0.1)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue2 = mt.rawp2adjp(corPvalueStudent(moduleTraitCor, nSamples), proc="BH")
temp <- moduleTraitPvalue2$adjp[order(moduleTraitPvalue2$index),]
moduleTraitPvalue_1 <- matrix(unlist(temp[,2]),ncol=4,byrow=FALSE)
colnames(moduleTraitPvalue_1) <- colnames(moduleTraitCor)
rownames(moduleTraitPvalue_1) <- rownames(moduleTraitCor)
modulTrait <- cbind(moduleTraitCor,moduleTraitPvalue,moduleTraitPvalue_1)
colnames(modulTrait)[5:8] <- paste0(colnames(modulTrait)[5:8],"pvalue")
colnames(modulTrait)[9:12] <- paste0(colnames(modulTrait)[9:12],"padj")
write.table(modulTrait,"wb_allmiRNA_TPM_correlation_Modules_Traits_values.txt",sep = "\t",quote = F,row.names = T)
rm(moduleTraitPvalue2,temp)
sizeGrWindow(10,6)
#moduleTraitCor_modi <- apply(moduleTraitCor, 2, FUN = function(x){ifelse(x<0,0,x)})
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue_1, 2), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
colnames(textMatrix) = colnames(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
library(RColorBrewer)

textMatrix <- textMatrix[-5,c(2:4,1)] ##Be careful need to delete the M0: grey module
moduleTraitCor_modi <- moduleTraitCor[-5,c(2:4,1)]
addrowcolor <- modules[match(names(MEs),modules$Label),]

rowLabels = paste("ME", addrowcolor$Color, sep="")
pdf("wb_allmiRNA_TPM_correlation_Modules_Traits.pdf",useDingbats = F,width = 9,height = 10)
labeledHeatmap(Matrix = moduleTraitCor_modi,
               xLabels = colnames(design)[c(2:4,1)],
               yLabels = rowLabels[-5],
               ySymbols = names(MEs)[-5],
               xColorLabels = FALSE,
               yColorLabels = TRUE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.76,cex.lab.y=1.2,cex.lab.x = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

##Relationship to trait and important modules: Gene significance and membership
#calculate the correlation between genes and each module
#modNames <- substring(names(MEs),3) ##different to substr, which need the start and end position
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("M", names(MEs), sep="")
names(MMPvalue) = paste("p.M", names(MEs), sep="")
MMPadj = mt.rawp2adjp(corPvalueStudent(as.matrix(geneModuleMembership), nSamples), proc="BH")

temp2 <- MMPadj$adjp[order(MMPadj$index),]
MMPadj <- matrix(unlist(temp2[,2]),ncol=14,byrow=FALSE)
colnames(MMPadj) <- paste0(colnames(MMPvalue),"_BH") 

kme=data.frame(moduleColor,moduleLabel, geneModuleMembership,MMPvalue,MMPadj,stringsAsFactors = F)
kme_export <- data.frame(ID=rownames(kme),kme)
write.table(kme_export,"wb_allmiRNA_TPM_SignedKME.txt",col.names = T,row.names = F,quote = F,sep = "\t")


###export data to cytoscape
#filter signkMEs less than 0.5 AND filter adjacency connectivity less than 0.02
export_module_cyto <- function(x,y,z){
  export_module <- rownames(kme)[which(kme$moduleLabel==x)]
  kme_module <- kme[export_module,]
  kme_module_sort <- kme_module[order(kme_module[, y],decreasing = TRUE),]
  kme_module_sort <- kme_module_sort[which(kme_module_sort[, y] > 0.5),]
  #adj <- adjacency[rownames(kme_module_sort), rownames(kme_module_sort)
  #rownames(adj)=colnames(adj)=rownames(kme_module_sort)
  #only show the top20 miRNAs
  adj <- adjacency[rownames(kme_module_sort)[1:20], rownames(kme_module_sort)[1:20]]
  rownames(adj)=colnames(adj)=rownames(kme_module_sort)[1:20]
  # Export the network into edge and node list files for Cytoscape
  cyt = exportNetworkToCytoscape(adj,
                                 edgeFile=paste0("WB_allmiRNA_Cyto_",z,"_",y,"edge.txt"),
                                 nodeFile=paste0("WB_allmiRNA_Cyto_",z,"_",y,"node.txt"),
                                 weighted = TRUE, threshold = 0.02,nodeNames=rownames(adj),
                                 altNodeNames = rownames(adj), nodeAttr = rep(x, nrow(adj)))
}
#only export each trait the most two correlated modules
export_module_cyto("M8","MM8","CW_upTop20")
export_module_cyto("M12","MM12","CWup_SkindownTop20")
export_module_cyto("M7","MM7","CWdown_SkinupTop20")
export_module_cyto("M3","MM3","CWdown_SkinupTop20")
export_module_cyto("M9","MM9","CWdown_SkinupTop20")

export_module_cyto("M5","MM5","Wound7_upTop20")
export_module_cyto("M6","MM6","Wound7_upTop20")

export_module_cyto("M2","MM2","Wound1_upTop20")
export_module_cyto("M10","MM10","Wound1_upTop20")
export_module_cyto("M11","MM11","Wound1_upTop20")


#####Step5: Modules preservation between mRNA and circRNA####
setLabels = c("test", "others")

test_preservation <- datExpr
others_preservation <- datExpr

# We now set up the multi-set expression data
multiExpr=list(test=list(data=test_preservation),
               others=list(data=others_preservation))
names(multiExpr)
labels = kme$moduleColor
names(labels)=rownames(kme)
multiColor=list(test=t(labels),others=t(labels))

system.time({
  mp = modulePreservation(multiExpr, multiColor,networkType = "signed",
                          referenceNetworks = c(1:2), nPermutations = 200,
                          checkData = FALSE,
                          randomSeed = 1, quickCor = 0, verbose = 5)
})
# Save the results of the module preservation analysis
# specify the reference and the test networks
ref=1; test = 2
Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
# Z statistics from the permutation test analysis
sum(Z.PreservationStats$moduleSize)
# Let us now visualize the data.
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]
point.label2 <- data.frame(ID=point.label) %>% 
  left_join(.,addrowcolor,by=c("ID"="Color")) %>% 
  select(2)
'point.label2' <- as.character(point.label2$Label)
#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres

pdf("wb_mRNA_modulePreservation.pdf",useDingbats = F,width = 11,height = 7)
par(mfrow=c(1,2))
#plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x",ylim=c(0,14))
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label2,cex=1,offs=0.03)
#plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Preservation of modules defined in \n CW samples only in other samples",
     cex=2,ylab ="Preservation Z-summary", xlab = "Module size", ylim=c(0,52))#, xlim=c(0,350),ylim=c(-1,15)
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label2,cex=1,offs=0.1)
# Add threshold lines for Zsummary
abline(h=2, col = "red", lty = 2)
abline(h=10, col = "red", lty = 2)

dev.off()

