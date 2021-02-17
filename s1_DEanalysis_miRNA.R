library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(ggforce)
library(concaveman)
library(sva) #ComBat function to remove batch effect

setwd("C:/Users/zhuliu/Desktop/miRNA/s3DEanalysis/DEresults_wb_mRNA/")
rm(list = ls())
#read the phenotype of all samples
pheData <- read_tsv("../readCounts_TPM_expression_Phenotypes/allSample_PheData.txt")

pheData_wb <- pheData %>% filter(Batch_seq=="Third") %>% select(2,3,6,8,9)
str(pheData_wb)

#####----whole biopsy----####
#wb_miRNA <- read_tsv("../readCounts_TPM_expression_Phenotypes/all_wb_miRNA_readCounts_tpm_re_someNovel.txt")
#data1 <- wb_miRNA %>% select(1:21) %>% column_to_rownames(.,var = "miRNA")
#name_data1 <- c(paste0("CW",1:5),paste0("Skin",1:5),paste0("Wound1",1:5),paste0("Wound7",1:5))
#colnames(data1) <- name_data1
#str(data1);dim(data1)
#keep=which(rowSums(data1[,] > 5) >= 10) ##mRNA expressed in at least half of all samples and average 10 read counts
#data1 <- as.matrix(data1[keep, ])
##data1 <- as.matrix(data1)
#dim(data1)
#
#tpm <- function(counts) {
#  rate <- counts / sum(counts)
#  rate * 1e6
#}
#tpms <- apply(data1, 2, function(x) tpm(x))
#colnames(tpms) <- paste0(colnames(tpms),"_TPM")
#all_wb_readCounts_tpm <- data.frame(miRNA=rownames(data1),data1,tpms)
#data.table::fwrite(all_wb_readCounts_tpm,"all_wb_miRNA_readCounts_tpm_re_someNovel_filter.txt",sep = "\t")

wb_miRNA <- read_tsv("../readCounts_TPM_expression_Phenotypes/all_wb_miRNA_readCounts_tpm_re_someNovel_filter.txt")
data1 <- wb_miRNA %>% select(1:21) %>% column_to_rownames(.,var = "miRNA")
##DE results for the Wound7,Wound1 and skin using the paired sample test
##Compared to CW, using unpaired sample test

DE_analysis <- function(x,y){
  data_x <- data1[,grep(x,colnames(data1))]
  data_y <- data1[,grep(y,colnames(data1))]
  data <- cbind(data_x,data_y)
  colData_x <- pheData_wb %>% filter(Groups==x)
  colData_y <- pheData_wb %>% filter(Groups==y)
  colData1 <- bind_rows(colData_x, colData_y)
  colData <- data.frame(sampleID=colnames(data),colData1)
  colData$Patient <- parse_factor(colData$Patient)
  colData$Sex <- parse_factor(colData$Sex)
  colData$Age <- parse_factor(colData$Age)
  colData$Groups <- parse_factor(colData$Groups)
  print(colData)
  str(colData)
  #For the unpaired DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = data,colData = colData, design = ~  Groups) 
  #For the paired DESeq2 analysis
  #dds <- DESeqDataSetFromMatrix(countData = data,colData = colData, design = ~ Patient + Groups) 
  dds
  dds$Groups <- relevel(dds$Groups, ref = y)
  dds <- DESeq(dds)
  res <- results(dds, pAdjustMethod = "BH", cooksCutoff = FALSE, contrast = c("Groups", x, y)) 
  print(res)
  mcols(res, use.names=TRUE)
  print(sum(res$pvalue < 0.05, na.rm = TRUE))
  print(sum(res$padj < 0.05, na.rm = TRUE))
  print(resultsNames(dds))
  # obtain mean value of groupA
  baseA <- counts(dds, normalized=TRUE)[, colData(dds)$Groups == x]
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <- x
  # obtain mean value of groupB
  baseB <- counts(dds, normalized=TRUE)[, colData(dds)$Groups == y]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- y
  ##combine the results
  res_dataframe <- cbind(ID=rownames(data), baseMeanA, baseMeanB, as.data.frame(res))
  res_dataframe$padj[is.na(res_dataframe$padj)] <- 1
  res_dataframe$pvalue[is.na(res_dataframe$pvalue)] <- 1
  res_dataframe$significant <- as.factor(ifelse(res_dataframe$padj < 0.05 & abs(res_dataframe$log2FoldChange) >=1,ifelse(res_dataframe$log2FoldChange >= 1 ,'Up','Down'),'Not'))
  res_dataframe <- res_dataframe %>% left_join(.,wb_miRNA[,c(1,22:41)],by=c("ID"="miRNA"))
  write.table(res_dataframe, paste0("DE_wb_miRNA10Exp_",x,y,".txt"), row.names = F,col.names = T, quote=F, sep = "\t")
  significantCirc <- res_dataframe %>% 
    filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% 
    filter(pvalue < 0.05)
  #write.table(significantCirc, paste0("DE_wb_miRNA10Exp_",x,y,"sig.txt"), row.names = F,col.names = T, quote=F, sep = "\t")
  significantCirc2 <- res_dataframe %>% 
    filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% 
    filter(padj < 0.05)
  write.table(significantCirc2, paste0("DE_wb_miRNA10Exp_",x,y,"sigPadj.txt"), row.names = F,col.names = T, quote=F, sep = "\t")
  #volcano plot
  p1 <- ggplot(data=res_dataframe, aes(x=log2FoldChange, y=-log10(padj), 
                                       color=significant,fill=significant,shape =significant)) +
    geom_point(size=1.5)+
    scale_color_manual(values=c("blue","#B3B3B3", "red"))+
    coord_fixed()+
    geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+ 
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
    #scale_x_continuous(limits = c(-4,6.2),breaks = c(-4,-2,0,2,4,6)) +
    scale_y_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10)) +
    scale_x_continuous(limits = c(-5,4),breaks = c(-4,-2,0,2,4)) +
    #scale_y_continuous(limits = c(0,50),breaks = c(0,10,20,30,40,50)) +
    theme(
      plot.background = element_rect(fill = NA), #color="grey50",size=2
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "black",fill = NA,size = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = "right",#top,bottom,left,right
      legend.spacing = unit(0,"cm"),
      aspect.ratio = 1)
  ggsave(filename = paste0("DE_wb_miRNA10Exp_",x,y,"volcano.pdf"), plot = p1,width = 8,height = 8)
}
DE_analysis("CW","Skin")
DE_analysis("CW","Wound1")
DE_analysis("CW","Wound7")
#using the paired DE analysis
DE_analysis("Wound1","Skin")
DE_analysis("Wound7","Skin")
DE_analysis("Wound7","Wound1")

##combined the DE analysis results
path <- getwd()
filename <- dir(path,pattern = "*.txt")
filepath <- sapply(filename, function(x){ 
  paste(path,x,sep='/')}) 
allmRNAdatas <- lapply(filepath, function(x){
  read_tsv(x)
})
rm(path,filepath)
names(allmRNAdatas) <- gsub(".txt","",names(allmRNAdatas))

allmRNAdatas[["Ori_miRNA_readCount_TPM"]] <- wb_miRNA
allmRNAdatas[["allData_miRNA_phenoInfo"]] <- pheData
allmiRNAdata_wb <- allmRNAdatas
save(allmiRNAdata_wb,file = "DEresults_wb_miRNAExpression_new.RData")

##using miRNA together do the PCA analysis
data2 <- wb_allgene_NewDE_mRNAresults$readCount_wb
data1 <- as.matrix(data2[,2:21])
Groups <- factor(c(rep("CW",5), rep("Skin",5), rep("Wound1",5), rep("Wound7",5)))
colData <- data.frame(sampleID=colnames(data1),Groups=Groups)
colData
keep=which(rowSums(data1[,] > 5) >= 10) ##mRNA expressed in at least half of all samples and average 10 read counts
data1 <- as.matrix(data1[keep, ])
data1 <- as.matrix(data1)
dim(data1)
#first need to construct the DEseq matrix
dds <- DESeqDataSetFromMatrix(countData = data1,colData = colData, design = ~ Groups);dds
dds <- DESeq(dds)
vst1 <- varianceStabilizingTransformation(dds,blind = FALSE)
pcaData <- plotPCA(vst1, intgroup=c( "Groups"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
max(pcaData$PC1);min(pcaData$PC1)
max(pcaData$PC2);min(pcaData$PC2)
nudge <- position_nudge(y = 1)
p1<-ggplot(pcaData, aes(PC1, PC2, color=Groups, shape =Groups)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #scale_x_continuous(limits = c(-20,21),breaks = c(-20,-15,-10,-5,0,5,10,15,20)) +
  #scale_y_continuous(limits = c(-15,30),breaks = c(-15,-10,-5,0,5,10,15,20,25,30)) +
  coord_fixed()+
  theme(legend.position = "bottom") +
  scale_colour_manual(values = brewer.pal(4,'Dark2')) +
  theme(
    plot.background = element_rect(fill = NA), #color="grey50",size=2
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(colour = "black",fill = NA,size = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",#top,bottom,left,right
    legend.spacing = unit(0,"cm"),
    aspect.ratio = 1)
p1
p2 <- p1 + geom_mark_hull(expand = unit(2.7, "mm")) 
p2
pdf("wb_mRNA_PCA.pdf",height = 8,width = 8,useDingbats = F)
p2 + geom_text(aes(label = name), hjust = -0.2, nudge_x = -0.2, vjust = -0.2, nudge_y = -0.2,show.legend = F,check_overlap = F,stat = "identity")
dev.off()

##Also consider the age sex and patient information
#check the batch effect
dist_mat <- dist(t(data1))
clustering <- hclust(dist_mat, method = "complete")
pdf("wb_miRNA501_clustering.pdf",height = 14,width = 8,useDingbats = F)
par(mfrow = c(4, 1), mar = c(3,5,3,1))
plot(clustering, labels = pheData_wb$Patient)
plot(clustering, labels = pheData_wb$Sex)
plot(clustering, labels = pheData_wb$Age)
plot(clustering, labels = pheData_wb$Groups)
dev.off()
