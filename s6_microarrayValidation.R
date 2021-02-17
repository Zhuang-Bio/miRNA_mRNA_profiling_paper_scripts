library(tidyverse)
library(GEOquery)
library(limma)
library(data.table)
rm(list = ls())

setwd("C:/Users/zhuliu/Desktop/miRNA/s7miRNA_validation/microarray_miRNA_oridata/")
####----inhibitor 149-5p microarray data----####
rm(list = ls())
si149_5p <- data.table::fread("149-5p inhibitor FB sample signal .txt", skip = 4)
colnames(si149_5p) <- gsub("-", "_", colnames(si149_5p))
si149_5p_meta <- data.table::fread("149-5p inhibitor FB.txt", skip = 4)
##add the gene annotation to expression data and modify the geneSymbol columns
si149_5p_anno <- si149_5p %>% select(1:3,7,4:6) %>% left_join(., si149_5p_meta[,c(1,17,13)], by=c("ID"="ID")) %>% 
  separate(`Gene Symbol`, into = c("GeneSymbol","Others"),sep = ";") %>% filter(Group == "Coding" | Group == "Multiple_Complex") %>% select(1:7,9)
##transform the signal data into expression data
#si149_5p_anno2 <- apply(si149_5p_anno[,2:7], 2, FUN = function(x) 2^x)
#si149_5p_anno2 <- cbind(si149_5p_anno,si149_5p_anno2) %>% select(-2:-7)
##filter the probes which shared the same gene symbol, select the largestly expressed one
si149_5p_anno$avg <- apply(si149_5p_anno[,2:7],1,mean)
si149_5p_anno2 <- si149_5p_anno %>% arrange(desc(avg)) %>% 
  distinct(GeneSymbol,.keep_all = TRUE) %>% select(8,2:7) %>% column_to_rownames(var = "GeneSymbol")

boxplot(si149_5p[,-1])
#make the design matrix and check if the Test and Control Group are correct
##Shoud pay attention to the factor levels (The order)
Group <- factor(c(rep("NC",3), rep("si149_5p",3)), levels = c("NC", "si149_5p"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- c("Control","Test")
rownames(design) <- colnames(si149_5p_anno2)
design
#Fit linear model for each gene given a series of arrays
fit <- lmFit(si149_5p_anno2, design)
summary(fit$Amean)
keep <- fit$Amean > 5
fit <- fit[keep,]

#Construct the contrast matrix corresponding to specified contrasts
#This step is very important
contrast.matrix <- makeContrasts(Test-Control, levels = design)
#Given a linear model fit to microarray data,
#compute estimated coefficients and standard errors for a given set of contrasts
fit2 <- contrasts.fit(fit, contrast.matrix); fit2
fit2 <- eBayes(fit2) 
#Extract a table of the top-ranked genes from a linear model fit.
si149_5p_anno2_f <- si149_5p_anno2 %>% rownames_to_column(var="GeneSymbol")
si149_5p_DEG <- topTable(fit2, n=Inf) %>%
  rownames_to_column(var = "probeID") %>% 
  left_join(.,si149_5p_anno2_f,by=c("probeID"= "GeneSymbol")) %>% 
  add_column(Source = "siRNA_149_5p") %>% arrange(desc(logFC))
#write all the DE results into file
fwrite(si149_5p_DEG, file = "../DEG_all_149_5p_InhibitorFB.txt", sep="\t")

####----mimics 149-5p microarray data----####
rm(list = ls())
mic149_5p <- data.table::fread("149-5p mimics FB sample signal.txt", skip = 4)
colnames(mic149_5p) <- gsub("-", "_", colnames(mic149_5p))
mic149_5p_meta <- data.table::fread("149-5p mimics FB.txt", skip = 4)
mic149_5p_anno <- mic149_5p %>% left_join(., mic149_5p_meta[,c(1,17,13)], by=c("ID"="ID")) %>% 
  separate(`Gene Symbol`, into = c("GeneSymbol","Others"),sep = ";") %>% filter(Group == "Coding" | Group == "Multiple_Complex") %>% select(1:7,9)
mic149_5p_anno$avg <- apply(mic149_5p_anno[,2:7],1,mean)
mic149_5p_anno2 <- mic149_5p_anno %>% arrange(desc(avg)) %>% 
  distinct(GeneSymbol,.keep_all = TRUE) %>% select(8,2:7) %>% column_to_rownames(var = "GeneSymbol")
boxplot(mic149_5p[,-1])
Group <- factor(c(rep("NC",3), rep("mic149_5p",3)), levels = c("NC", "mic149_5p"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- c("Control","Test")
rownames(design) <- colnames(mic149_5p_anno2)
design
fit <- lmFit(mic149_5p_anno2, design)
summary(fit$Amean)
keep <- fit$Amean > 5
fit <- fit[keep,]
contrast.matrix <- makeContrasts(Test-Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix); fit2
fit2 <- eBayes(fit2) 
mic149_5p_anno2_f <- mic149_5p_anno2 %>% rownames_to_column(var="GeneSymbol")
mic149_5p_DEG <- topTable(fit2, n=Inf) %>%
  rownames_to_column(var = "probeID") %>% 
  left_join(.,mic149_5p_anno2_f,by=c("probeID"= "GeneSymbol")) %>% 
  add_column(Source = "micRNA_149_5p") %>% arrange(desc(logFC))
fwrite(mic149_5p_DEG, file = "../DEG_all_149_5p_mimicsFB.txt", sep="\t")

####----mimics 424-5p microarray data----####
rm(list = ls())
mic424_5p <- data.table::fread("424-5p mimics  FB sample signal .txt", skip = 4)
colnames(mic424_5p) <- gsub("-", "_", colnames(mic424_5p))
mic424_5p_meta <- data.table::fread("424-5p mimics  FB.txt", skip = 4)
mic424_5p_anno <- mic424_5p %>% select(1:3,7,4:6) %>% left_join(., mic424_5p_meta[,c(1,17,13)], by=c("ID"="ID")) %>% 
  separate(`Gene Symbol`, into = c("GeneSymbol","Others"),sep = ";") %>% filter(Group == "Coding" | Group == "Multiple_Complex") %>% select(1:7,9)
mic424_5p_anno$avg <- apply(mic424_5p_anno[,2:7],1,mean)
mic424_5p_anno2 <- mic424_5p_anno %>% arrange(desc(avg)) %>% 
  distinct(GeneSymbol,.keep_all = TRUE) %>% select(8,2:7) %>% column_to_rownames(var = "GeneSymbol")
Group <- factor(c(rep("NC",3), rep("mic424_5p",3)), levels = c("NC", "mic424_5p"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- c("Control","Test")
rownames(design) <- colnames(mic424_5p_anno2)
design
fit <- lmFit(mic424_5p_anno2, design)
summary(fit$Amean)
keep <- fit$Amean > 5
fit <- fit[keep,]
contrast.matrix <- makeContrasts(Test-Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix); fit2
fit2 <- eBayes(fit2) 
mic424_5p_anno2_f <- mic424_5p_anno2 %>% rownames_to_column(var="GeneSymbol")
mic424_5p_DEG <- topTable(fit2, n=Inf) %>%
  rownames_to_column(var = "probeID") %>% 
  left_join(.,mic424_5p_anno2_f,by=c("probeID"= "GeneSymbol")) %>% 
  add_column(Source = "micRNA_424_5p") %>% arrange(desc(logFC))
fwrite(mic424_5p_DEG, file = "../DEG_all_424_5p_mimicsFB.txt", sep="\t")

####----mimics 7704_FB microarray data----####
rm(list = ls())
mic7704_FB <- data.table::fread("7704 mimics FB sample signal.txt", skip = 4)
colnames(mic7704_FB) <- gsub("-", "_", colnames(mic7704_FB))
mic7704_FB_meta <- data.table::fread("7704 mimics FB.txt", skip = 4)
mic7704_FB_anno <- mic7704_FB %>% left_join(., mic7704_FB_meta[,c(1,17,13)], by=c("ID"="ID")) %>% 
  separate(`Gene Symbol`, into = c("GeneSymbol","Others"),sep = ";") %>% filter(Group == "Coding" | Group == "Multiple_Complex") %>% select(1:7,9)
mic7704_FB_anno$avg <- apply(mic7704_FB_anno[,2:7],1,mean)
mic7704_FB_anno2 <- mic7704_FB_anno %>% arrange(desc(avg)) %>% 
  distinct(GeneSymbol,.keep_all = TRUE) %>% select(8,2:7) %>% column_to_rownames(var = "GeneSymbol")
Group <- factor(c(rep("NC",3), rep("mic7704_FB",3)), levels = c("NC", "mic7704_FB"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- c("Control","Test")
rownames(design) <- colnames(mic7704_FB_anno2)
design
fit <- lmFit(mic7704_FB_anno2, design)
summary(fit$Amean)
keep <- fit$Amean > 5
fit <- fit[keep,]
contrast.matrix <- makeContrasts(Test-Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix); fit2
fit2 <- eBayes(fit2) 
mic7704_FB_anno2_f <- mic7704_FB_anno2 %>% rownames_to_column(var="GeneSymbol")
mic7704_FB_DEG <- topTable(fit2, n=Inf) %>%
  rownames_to_column(var = "probeID") %>% 
  left_join(.,mic7704_FB_anno2_f,by=c("probeID"= "GeneSymbol")) %>% 
  add_column(Source = "micRNA_7704_FB") %>% arrange(desc(logFC))
fwrite(mic7704_FB_DEG, file = "../DEG_all_7704_mimicsFB.txt", sep="\t")

####----mimics 7704_KC microarray data----####
rm(list = ls())
mic7704_KC <- data.table::fread("7704 mimics KC sample signal.txt", skip = 4)
colnames(mic7704_KC) <- gsub("-", "_", colnames(mic7704_KC))
mic7704_KC_meta <- data.table::fread("7704 mimics KC.txt", skip = 4)
mic7704_KC_anno <- mic7704_KC %>% left_join(., mic7704_KC_meta[,c(1,17,13)], by=c("ID"="ID")) %>% 
  separate(`Gene Symbol`, into = c("GeneSymbol","Others"),sep = ";") %>% filter(Group == "Coding" | Group == "Multiple_Complex") %>% select(1:7,9)
mic7704_KC_anno$avg <- apply(mic7704_KC_anno[,2:7],1,mean)
mic7704_KC_anno2 <- mic7704_KC_anno %>% arrange(desc(avg)) %>% 
  distinct(GeneSymbol,.keep_all = TRUE) %>% select(8,2:7) %>% column_to_rownames(var = "GeneSymbol")
Group <- factor(c(rep("NC",3), rep("mic7704_KC",3)), levels = c("NC", "mic7704_KC"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- c("Control","Test")
rownames(design) <- colnames(mic7704_KC_anno2)
design
fit <- lmFit(mic7704_KC_anno2, design)
summary(fit$Amean)
keep <- fit$Amean > 5
fit <- fit[keep,]
contrast.matrix <- makeContrasts(Test-Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix); fit2
fit2 <- eBayes(fit2) 
mic7704_KC_anno2_f <- mic7704_KC_anno2 %>% rownames_to_column(var="GeneSymbol")
mic7704_KC_DEG <- topTable(fit2, n=Inf) %>%
  rownames_to_column(var = "probeID") %>% 
  left_join(.,mic7704_KC_anno2_f,by=c("probeID"= "GeneSymbol")) %>% 
  add_column(Source = "micRNA_7704_KC") %>% arrange(desc(logFC))
fwrite(mic7704_KC_DEG, file = "../DEG_all_7704_mimicsKC.txt", sep="\t")

