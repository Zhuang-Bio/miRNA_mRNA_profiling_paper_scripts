library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(data.table)

rm(list = ls())
setwd("C:/Users/zhuliu/Desktop/miRNA/s5miRNA_TargetPredict/s2_wb_miRNA_targets_DE_modules_enrichment/")

####----step1 read the mRNA modules and DE mRNAs----####
wb_mRNA_wgcna <- read_tsv("C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/s2_wb_mRNA_WGCNA/wb_allmRNA_SignedKME_new.txt")
mRNA_m1 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M1")
mRNA_m2 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M2")
mRNA_m3 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M3")
mRNA_m4 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M4")
mRNA_m5 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M5")
mRNA_m6 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M6")
mRNA_m7 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M7")
mRNA_m8 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M8")
mRNA_m9 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M9")
mRNA_m10 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M10")
mRNA_m11 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M11")
mRNA_m12 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M12")
mRNA_m13 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M13")

##load the DE mRNA results
load("C:/Users/zhuliu/Desktop/miRNA/s3DEanalysis/DEresults_wb_mRNA/wb_allgene_NewDE_mRNAresults.RData")
#cw_skin <- wb_allgene_NewDE_mRNAresults$wb_allmRNA_DEanalysisCWSkin %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) >= 0.58)
#w1_skin <- wb_allgene_NewDE_mRNAresults$wb_allmRNA_DEanalysisWound1Skin %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) >= 0.58)
#w7_w1 <- wb_allgene_NewDE_mRNAresults$wb_allmRNA_DEanalysisWound7Skin %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) >= 0.58)
#cw_DEmRNA_combined <- cw_skin %>% 
#  inner_join(.,w1_skin[,c(1,3,6,7)],by=c("ID"="ID")) %>% 
#  inner_join(.,w7_w1[,c(1,3,6,7)],by=c("ID"="ID"))
#cw_DEmRNA_combined$sigW7sk <- ifelse(cw_DEmRNA_combined$padj.x < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange.x) >= 0.58,ifelse(cw_DEmRNA_combined$log2FoldChange.x >= 0.58,'Down','Up'),'Not')
#cw_DEmRNA_combined$sigW7w1 <- ifelse(cw_DEmRNA_combined$padj.y < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange.y) >= 0.58,ifelse(cw_DEmRNA_combined$log2FoldChange.y >= 0.58,'Down','Up'),'Not')
#cw_DEmRNA_combined$sigcwW7 <- ifelse(cw_DEmRNA_combined$padj < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange) >= 0.58,ifelse(cw_DEmRNA_combined$log2FoldChange >= 0.58,'Down','Up'),'Not')
#cw_DEmRNA_combined_down <- cw_DEmRNA_combined %>% filter(sigW7sk == "Down" & sigW7w1 == "Down" & sigcwW7 == "Down") %>% 
#  left_join(.,wb_allgene_NewDE_mRNAresults$HumanGene_PC_v32[,1:2],by=c("ID"="Geneid"))
##data.table::fwrite(cw_DEmRNA_combined_down,"skin_DEmRNA_combined_down.txt",sep = "\t")
#cw_DEmRNA_combined_up <- cw_DEmRNA_combined %>% filter(sigW7sk == "Up" & sigW7w1 == "Up" & sigcwW7 == "Up") %>% 
#  left_join(.,wb_allgene_NewDE_mRNAresults$HumanGene_PC_v32[,1:2],by=c("ID"="Geneid"))
##data.table::fwrite(cw_DEmRNA_combined_up,"skin_DEmRNA_combined_up.txt",sep = "\t")
cw_up <- fread('wb_cw_DEmRNA_combined_up.txt')
cw_down <- fread('wb_cw_DEmRNA_combined_down.txt')
wound7_down <- fread('w7_DEmRNA_combined_down_noCW.txt')
wound7_up <- fread('w7_DEmRNA_combined_up_noCW.txt')
wound1_down <- fread('w1_DEmRNA_combined_down_noCW.txt')
wound1_up <- fread('w1_DEmRNA_combined_up_noCW.txt')
skin_down <- fread('skin_DEmRNA_combined_down_noCW.txt')
skin_up <- fread('skin_DEmRNA_combined_up_noCW.txt')

##miRNA
load('C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/DEresults_wb_miRNAExpression_new.RData')
#cw_DEmRNA_combined <- allmiRNAdata_wb$DE_wb_miRNA10Exp_CWWound1[,c(1,5,8,9)] %>% 
#  inner_join(.,allmiRNAdata_wb$DE_wb_miRNA10Exp_Wound1Skin[,c(1,5,8,9)],by=c("ID"="ID")) %>% 
#  inner_join(.,allmiRNAdata_wb$DE_wb_miRNA10Exp_Wound7Wound1[,c(1,5,8,9)],by=c("ID"="ID"))
#cw_DEmRNA_combined$sigCWsk <- ifelse(cw_DEmRNA_combined$padj.x < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange.x) >=1,ifelse(cw_DEmRNA_combined$log2FoldChange.x >= 1 ,'Down','Up'),'Not')
#cw_DEmRNA_combined$sigCWw1 <- ifelse(cw_DEmRNA_combined$padj.y < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange.y) >=1,ifelse(cw_DEmRNA_combined$log2FoldChange.y >= 1 ,'Up','Down'),'Not')
#cw_DEmRNA_combined$sigCWw7 <- ifelse(cw_DEmRNA_combined$padj < 0.05 & abs(cw_DEmRNA_combined$log2FoldChange) >=1,ifelse(cw_DEmRNA_combined$log2FoldChange >= 1 ,'Down','Up'),'Not')
#table(cw_DEmRNA_combined$sigCWsk);table(cw_DEmRNA_combined$sigCWw1);table(cw_DEmRNA_combined$sigCWw7)
###Define the CW DE mRNA down or up -regulated (keep the comparison the same trend)
#cw_DEmRNA_combined_down <- cw_DEmRNA_combined %>% filter(sigCWsk == "Down" & sigCWw1 == "Down" & sigCWw7 == "Down") 
#data.table::fwrite(cw_DEmRNA_combined_down,"w1_DEmiRNA_combined_down.txt",sep = "\t")
#cw_DEmRNA_combined_up <- cw_DEmRNA_combined %>% filter(sigCWsk == "Up" & sigCWw1 == "Up" & sigCWw7 == "Up")
#data.table::fwrite(cw_DEmRNA_combined_up,"w1_DEmiRNA_combined_up.txt",sep = "\t")

####---step2 read the module and DE miRNA target genes----####
allmiRNA_targets_ori <- read_tsv("../s1_wb_miRNA_targets_multiMiR/all562miR_target_allTop25_TargetScan_miRtogether.txt")
length(unique(allmiRNA_targets_ori$target_symbol))
allmiRNA_targets <- read_tsv("../s1_wb_miRNA_targets_multiMiR/all562miR_target_allTop25_TargetScan_miRtogether.txt")
length(unique(allmiRNA_targets$target_symbol))

miRNAmodules <-  read_tsv("../s1_wb_miRNA_targets_multiMiR/wb_allmiRNA_TPM_SignedKME.txt")
#(kME > 0.809, 0.833, 0.828, 0.835, 0.812, 0.785, 0.83, 0.837, 0.81 and 0.785 
##  for M2,    M3,    M5,    M6,    M7,    M8,    M9,   M10,   M11 and  M12, respectively)

#miRNAmodules %>% filter(moduleLabel == "M2") %>% summarise_at('MM2', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M3") %>% summarise_at('MM3', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M4") %>% summarise_at('MM4', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M5") %>% summarise_at('MM5', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M6") %>% summarise_at('MM6', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M7") %>% summarise_at('MM7', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M8") %>% summarise_at('MM8', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M9") %>% summarise_at('MM9', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M10") %>% summarise_at('MM10', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M11") %>% summarise_at('MM11', median, na.rm = TRUE)
#miRNAmodules %>% filter(moduleLabel == "M12") %>% summarise_at('MM12', median, na.rm = TRUE)


miRNA_m2 <- miRNAmodules %>% filter(moduleLabel == "M2") %>% filter(MM2 > median(MM2)) %>% arrange(desc(MM2)) #median(miRNA_m2$MM2): 0.83
#miRNA_m2 <- miRNA_m2[miRNA_m2$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m2 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m2$ID,]
tmp <- table(miRNA_m2$target_symbol)
tmp_names <- names(which(tmp >= 4))
miRNA_m2 <- miRNA_m2[miRNA_m2$target_symbol%in%tmp_names,]
length(unique(miRNA_m2$target_symbol))

miRNA_m3 <- miRNAmodules %>% filter(moduleLabel == "M3") %>% filter(MM3 > median(MM3)) %>% arrange(desc(MM3)) #median(miRNA_m3$MM3): 0.83
#miRNA_m3 <- miRNA_m3[miRNA_m3$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m3 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m3$ID,]
tmp <- table(miRNA_m3$target_symbol)
tmp_names <- names(which(tmp >= 4))
miRNA_m3 <- miRNA_m3[miRNA_m3$target_symbol%in%tmp_names,]
length(unique(miRNA_m3$target_symbol))

miRNA_m5 <- miRNAmodules %>% filter(moduleLabel == "M5") %>% filter(MM5 > median(MM5)) %>% arrange(desc(MM5)) 
#miRNA_m5 <- miRNA_m5[miRNA_m5$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m5 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m5$ID,]
tmp <- table(miRNA_m5$target_symbol)
tmp_names <- names(which(tmp >= 3))
miRNA_m5 <- miRNA_m5[miRNA_m5$target_symbol%in%tmp_names, ]
length(unique(miRNA_m5$target_symbol))

miRNA_m6 <- miRNAmodules %>% filter(moduleLabel == "M6") %>% filter(MM6 > median(MM6)) %>% arrange(desc(MM6)) 
#miRNA_m6 <- miRNA_m6[miRNA_m6$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m6 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m6$ID,]
tmp <- table(miRNA_m6$target_symbol)
tmp_names <- names(which(tmp >= 3))
miRNA_m6 <- miRNA_m6[miRNA_m6$target_symbol%in%tmp_names, ]
length(unique(miRNA_m6$target_symbol))

miRNA_m7 <- miRNAmodules %>% filter(moduleLabel == "M7") %>% filter(MM7 > median(MM7)) #median(miRNA_m7$MM7): 0.8123159
#miRNA_m7 <- miRNA_m7[miRNA_m7$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m7 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m7$ID,]
tmp <- table(miRNA_m7$target_symbol)
tmp_names <- names(which(tmp >= 3))
miRNA_m7 <- miRNA_m7[miRNA_m7$target_symbol%in%tmp_names,]
length(unique(miRNA_m7$target_symbol))

miRNA_m8 <- miRNAmodules %>% filter(moduleLabel == "M8") %>% filter(MM8 > median(MM8)) %>% arrange(desc(MM8)) 
#miRNA_m8 <- miRNA_m8[miRNA_m8$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m8 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m8$ID,]
tmp <- table(miRNA_m8$target_symbol)
tmp_names <- names(which(tmp >= 2))
miRNA_m8 <- miRNA_m8[miRNA_m8$target_symbol%in%tmp_names, ]
length(unique(miRNA_m8$target_symbol))

miRNA_m9 <- miRNAmodules %>% filter(moduleLabel == "M9") %>% filter(MM9 > median(MM9))
#miRNA_m9 <- miRNA_m9[miRNA_m9$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m9 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m9$ID,]
tmp <- table(miRNA_m9$target_symbol)
tmp_names <- names(which(tmp >= 2))
miRNA_m9 <- miRNA_m9[miRNA_m9$target_symbol%in%tmp_names,]
length(unique(miRNA_m9$target_symbol))

miRNA_m10 <- miRNAmodules %>% filter(moduleLabel == "M10") %>% filter(MM10 > median(MM10)) %>% arrange(desc(MM10)) 
#miRNA_m10 <- miRNA_m10[miRNA_m10$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m10 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m10$ID,]
tmp <- table(miRNA_m10$target_symbol)
tmp_names <- names(which(tmp >= 2))
miRNA_m10 <- miRNA_m10[miRNA_m10$target_symbol%in%tmp_names, ]
length(unique(miRNA_m10$target_symbol))

miRNA_m11 <- miRNAmodules %>% filter(moduleLabel == "M11") %>% filter(MM11 > median(MM11)) %>% arrange(desc(MM11)) 
#miRNA_m11 <- miRNA_m11[miRNA_m11$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m11 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m11$ID,]
tmp <- table(miRNA_m11$target_symbol)
tmp_names <- names(which(tmp >= 2))
miRNA_m11 <- miRNA_m11[miRNA_m11$target_symbol%in%tmp_names, ]
length(unique(miRNA_m11$target_symbol))

miRNA_m12 <- miRNAmodules %>% filter(moduleLabel == "M12") %>% filter(MM12 > median(MM12))
#miRNA_m12 <- miRNA_m12[miRNA_m12$ID%in%miRNA_CW_DE_onlyfilterFDR$ID,]
miRNA_m12 <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRNA_m12$ID,]
tmp <- table(miRNA_m12$target_symbol)
tmp_names <- names(which(tmp >= 2))
miRNA_m12 <- miRNA_m12[miRNA_m12$target_symbol%in%tmp_names,]
length(unique(miRNA_m12$target_symbol))

#up and down regulated miRNAs
cwmiRNA_up <- fread('wb_CW_DEmiRNA_up.txt')
cw_up_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%cwmiRNA_up$ID,]
tmp <- table(cw_up_targets$target_symbol)
tmp_names <- names(which(tmp >= 3))
cw_up_targets <- cw_up_targets[cw_up_targets$target_symbol%in%tmp_names,]
length(unique(cw_up_targets$target_symbol))

cwmiRNA_down <- fread('wb_CW_DEmiRNA_down.txt')
cw_down_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%cwmiRNA_down$ID,]
tmp <- table(cw_down_targets$target_symbol)
tmp_names <- names(which(tmp >= 2))
cw_down_targets <- cw_down_targets[cw_down_targets$target_symbol%in%tmp_names,]
length(unique(cw_down_targets$target_symbol))

#############################################################################
##------Only Wound7,Wound1 and Skin DEmiRNAs compared with each other------##
#############################################################################
w7miRNA_up <- fread('w7_DEmiRNA_combined_up_noCW.txt')
w7_up_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%w7miRNA_up$ID,]
tmp <- table(w7_up_targets$target_symbol)
tmp_names <- names(which(tmp >= 5))
w7_up_targets <- w7_up_targets[w7_up_targets$target_symbol%in%tmp_names,]
length(unique(w7_up_targets$target_symbol))

w7miRNA_down <- fread('w7_DEmiRNA_combined_down_noCW.txt')
w7_down_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%w7miRNA_down$ID,]
tmp <- table(w7_down_targets$target_symbol)
tmp_names <- names(which(tmp >= 3))
w7_down_targets <- w7_down_targets[w7_down_targets$target_symbol%in%tmp_names,]
length(unique(w7_down_targets$target_symbol))

w1miRNA_up <- fread('w1_DEmiRNA_combined_up_noCW.txt')
w1_up_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%w1miRNA_up$ID,]
tmp <- table(w1_up_targets$target_symbol)
tmp_names <- names(which(tmp >= 3))
w1_up_targets <- w1_up_targets[w1_up_targets$target_symbol%in%tmp_names,]
length(unique(w1_up_targets$target_symbol))

w1miRNA_down <- fread('w1_DEmiRNA_combined_down_noCW.txt')
w1_down_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%w1miRNA_down$ID,]
tmp <- table(w1_down_targets$target_symbol)
tmp_names <- names(which(tmp >= 2))
w1_down_targets <- w1_down_targets[w1_down_targets$target_symbol%in%tmp_names,]
length(unique(w1_down_targets$target_symbol))

skinmiRNA_up <- fread('skin_DEmiRNA_combined_up_noCW.txt')
skin_up_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%skinmiRNA_up$ID,]
tmp <- table(skin_up_targets$target_symbol)
tmp_names <- names(which(tmp >= 4))
skin_up_targets <- skin_up_targets[skin_up_targets$target_symbol%in%tmp_names,]
length(unique(skin_up_targets$target_symbol))

skinmiRNA_down <- fread('skin_DEmiRNA_combined_down_noCW.txt')
skin_down_targets <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%skinmiRNA_down$ID,]
tmp <- table(skin_down_targets$target_symbol)
tmp_names <- names(which(tmp >= 6))
skin_down_targets <- skin_down_targets[skin_down_targets$target_symbol%in%tmp_names,]
length(unique(skin_down_targets$target_symbol))

####----step3 run the fish.test----####
#first calculate the overlapped number of genes
mRNAmodules <- bind_rows(
  data.frame(modules = 'M1', mRNA_m1),
  data.frame(modules = 'M2', mRNA_m2),
  data.frame(modules = 'M3', mRNA_m3),
  data.frame(modules = 'M4', mRNA_m4),
  data.frame(modules = 'M5', mRNA_m5),
  data.frame(modules = 'M6', mRNA_m6),
  data.frame(modules = 'M7', mRNA_m7),
  data.frame(modules = 'M8', mRNA_m8),
  data.frame(modules = 'M9', mRNA_m9),
  data.frame(modules = 'M10', mRNA_m10),
  data.frame(modules = 'M11', mRNA_m11),
  data.frame(modules = 'M12', mRNA_m12),
  data.frame(modules = 'M13', mRNA_m13),
  data.frame(modules = 'cw_up', cw_up[,c(1,17)]),
  data.frame(modules = 'cw_down', cw_down[,c(1,17)]),
  data.frame(modules = 'wound7_down', wound7_down[,c(1,10)]),
  data.frame(modules = 'wound7_up', wound7_up[,c(1,10)]),
  data.frame(modules = 'wound1_down', wound1_down[,c(1,10)]),
  data.frame(modules = 'wound1_up', wound1_up[,c(1,10)]),
  data.frame(modules = 'skin_down', skin_down[,c(1,10)]),
  data.frame(modules = 'skin_up', skin_up[,c(1,10)]))
table(mRNAmodules$modules)
#data.table::fwrite(mRNAmodules, file = 'Top_all_DEmRNAs_.txt', sep = '\t')

miRNAmodules <- bind_rows(
  data.frame(modules = 'M2', miRNA_m2),
  data.frame(modules = 'M3', miRNA_m3),
  data.frame(modules = 'M5', miRNA_m5),
  data.frame(modules = 'M6', miRNA_m6),
  data.frame(modules = 'M7', miRNA_m7),
  data.frame(modules = 'M8', miRNA_m8),
  data.frame(modules = 'M9', miRNA_m9),
  data.frame(modules = 'M10', miRNA_m10),
  data.frame(modules = 'M11', miRNA_m11),
  data.frame(modules = 'M12', miRNA_m12),
  data.frame(modules = 'cw_up', cw_up_targets),
  data.frame(modules = 'cw_down', cw_down_targets),
  data.frame(modules = 'wound7_down', w7_down_targets),
  data.frame(modules = 'wound7_up', w7_up_targets),  
  data.frame(modules = 'wound1_down', w1_down_targets),
  data.frame(modules = 'wound1_up', w1_up_targets),
  data.frame(modules = 'skin_down', skin_down_targets),
  data.frame(modules = 'skin_up', skin_up_targets))
table(miRNAmodules$modules)
#data.table::fwrite(miRNAmodules, file = 'Top_miRNA_targets.txt', sep = '\t')

##
cw_pubGenes <- read.table('clipboard', sep = '\t',header = T)
mRNAmodules <- cw_pubGenes

##########################################################
####------mRNA-miRNA modules enrichment analysis------####
##########################################################
names <- unique(mRNAmodules$modules)
namesmiR <- unique(miRNAmodules$modules)

fisher_re <- list()
fisher_re_overGene <- list()
for (i in 1:length(names)) {
  require(tidyverse)
  module <- names[i]
  moduleData <- mRNAmodules %>% filter(modules == module)
  NmRNA <- length(unique(moduleData$GeneSymbol))
  
  fisher_remi <- list()
  fisher_overlap <- list()
  fisher_overlapGene <- list()  
  for (j in 1:length(namesmiR)) {
    module2 <- namesmiR[j]
    moduleData2 <- miRNAmodules %>% filter(modules == module2)
    NmiRNA <- length(unique(moduleData2$target_symbol))
    fisher_remi[[module2]] <- NmiRNA
    genesOverlapps <- intersect(moduleData2$target_symbol, moduleData$GeneSymbol)
    if(is_empty(genesOverlapps)){
      fisher_overlapGene[[module2]] <- data.frame(miRmodule = module2, OverlapGene = "emptyOverlap")
    }else{
      fisher_overlapGene[[module2]] <- data.frame(miRmodule = module2, OverlapGene = genesOverlapps)
    }
    fisher_overlap[[module2]] <- length(unique(genesOverlapps))
  }
  fisher_remiCom <- do.call("rbind",fisher_remi) %>% as.data.frame()
  fisher_overlapCom <- do.call("rbind",fisher_overlap) %>% as.data.frame()
  fisher_overlapComGene <- do.call("rbind",fisher_overlapGene) %>% as.data.frame()
  fisher_re[[module]] <- data.frame(fisher_overlapCom, mRNA_number = NmRNA, fisher_remiCom) 
  fisher_re_overGene[[module]] <- data.frame(fisher_overlapComGene, mRNAmoudle = module) 
}

fisher_re_com <- do.call("rbind",fisher_re) %>% as.data.frame() %>% rownames_to_column(var = 'mModule_miModule') %>% 
  mutate(Background = 15969)

colnames(fisher_re_com)[c(2,4)] <- c('Overlap','miRNAtargets')
fisher_re_com$mModule_miModule <- gsub('\\.', '-', fisher_re_com$mModule_miModule)

fisher_re_com_f <- fisher_re_com %>% 
  mutate(mRNAnew = mRNA_number - Overlap,
         miRNAnew = miRNAtargets - Overlap,
         BGnew = Background - mRNA_number - miRNAtargets + Overlap) %>% 
  select(2,6:8)
#Extract the overlapped genes
fisher_re_overGene_final <- do.call("rbind",fisher_re_overGene) %>% as.data.frame() 
data.table::fwrite(fisher_re_overGene_final, file = "module_module_VUgenes.txt", sep = '\t')
fisher_test <- list()
for (i in 1:nrow(fisher_re_com_f)) {
  require(tidyverse)
  tmp1 <- fisher_re_com_f[i,] %>% as.character() %>% as.numeric()
  data_fisher <- matrix(tmp1, nrow = 2, byrow = TRUE,
                        dimnames = list(miRNA=c("Yes","No"), mRNA=c("Yes","No")))
  fisher_results <- fisher.test(data_fisher, alternative = "two.sided")
  fisher_results_odds <- fisher_results$estimate
  fisher_results_pvalue <- fisher_results$p.value
  fisher_results_conf1 <- fisher_results$conf.int[[1]]
  fisher_results_conf2 <- fisher_results$conf.int[[2]]
  fisher_test[[i]] <- c(fisher_results_odds, fisher_results_pvalue, fisher_results_conf1, fisher_results_conf2)
}
fisher_test_com <- do.call("rbind",fisher_test)%>% as.data.frame()
padjval <- p.adjust(fisher_test_com$V2, method = 'BH', n=length(fisher_test_com$V2))
fisher_test_comfinal <- data.frame(fisher_test_com, padjVal = padjval) %>% select(1,2,5,3,4)
colnames(fisher_test_comfinal) <- c('Odds_Ratio', 'pvalue', 'padj', 'CI95%L', 'CI95%R')
finalres <- cbind(fisher_re_com, fisher_test_comfinal)
data.table::fwrite(finalres, file = "allfisher_test_sigModules_Enrich_skinW1W7CW.txt",sep = "\t")
#data.table::fwrite(finalres, file = "allfisher_test_sigModules_Enrich_VUgenes.txt",sep = "\t")


########################################################################
####------Individual miRNAs enrichment analysis in each module------####
########################################################################
names <- unique(mRNAmodules$modules)
namesmiR <- unique(miRNAmodules$modules)
##only extract the wound7 modules(M5,M6) and DEmiRNAs of wound7, wound1 and skin
#namesmiR <- namesmiR[c(3,4,13:18)]
#i=1
#j=1
#z=1
fisher_refinal <- list()
for (i in 1:length(namesmiR)) {
  require(tidyverse)
  module <- namesmiR[i]
  moduleData <- miRNAmodules %>% filter(modules == module)
  
  miRls <- unique(moduleData$mature_mirna_id)
  
  fisher_re <- list()
  for (j in 1:length(miRls)) {
    miRtmp <- miRls[j]
    mirData <- allmiRNA_targets %>% filter(mature_mirna_id == miRtmp)
    NmiRNA <- length(unique(mirData$target_symbol))
    
    fisher_remi <- list()
    fisher_overlap <- list()
    for (z in 1:length(names)) {
      module2 <- names[z]
      moduleData2 <- mRNAmodules %>% filter(modules == module2)
      NmRNA <- length(unique(moduleData2$GeneSymbol))
      fisher_remi[[module2]] <- NmRNA
      genesOverlapps <- intersect(moduleData2$GeneSymbol, mirData$target_symbol)
      fisher_overlap[[module2]] <- length(unique(genesOverlapps))
    }
    fisher_remiCom <- do.call("rbind",fisher_remi) %>% as.data.frame()
    fisher_overlapCom <- do.call("rbind",fisher_overlap) %>% as.data.frame()
    fisher_re[[miRtmp]] <- data.frame(ID=miRtmp, Overlap = fisher_overlapCom, mRNA_number=fisher_remiCom, miRNA_number = NmiRNA) %>% rownames_to_column(var = 'mRNAmodules')
  }
  fisher_recom <- do.call("rbind",fisher_re) %>% as.data.frame()
  fisher_refinal[[module]] <- data.frame(miRNAmodule=module, fisher_recom)
}
fisher_re_com <- do.call("rbind",fisher_refinal) %>% as.data.frame() %>% 
  mutate(Background = 15969)
length(unique(fisher_re_com$ID))
colnames(fisher_re_com)[c(3:6)] <- c('miRNA_ID','Overlap','mRNA_number','miRNAtargets')

fisher_re_com_f <- fisher_re_com %>% mutate(mRNAnew = mRNA_number - Overlap, 
                                            miRNAnew = miRNAtargets - Overlap,
                                            BGnew = Background - mRNA_number - miRNAtargets + Overlap) %>% 
  select(4,8:10)

fisher_test <- list()
for (i in 1:nrow(fisher_re_com_f)) {
  require(tidyverse)
  tmp1 <- fisher_re_com_f[i,] %>% as.character() %>% as.numeric()
  data_fisher <- matrix(tmp1, nrow = 2, byrow = TRUE,
                        dimnames = list(miRNA=c("Yes","No"), mRNA=c("Yes","No")))
  fisher_results <- fisher.test(data_fisher, alternative = "two.sided")
  fisher_results_odds <- fisher_results$estimate
  fisher_results_pvalue <- fisher_results$p.value
  fisher_results_conf1 <- fisher_results$conf.int[[1]]
  fisher_results_conf2 <- fisher_results$conf.int[[2]]
  fisher_test[[i]] <- c(fisher_results_odds, fisher_results_pvalue, fisher_results_conf1, fisher_results_conf2)
}
fisher_test_com <- do.call("rbind",fisher_test)%>% as.data.frame()
padjval <- p.adjust(fisher_test_com$V2, method = 'BH', n=length(fisher_test_com$V2))
fisher_test_comfinal <- data.frame(fisher_test_com, padjVal = padjval) %>% select(1,2,5,3,4)
colnames(fisher_test_comfinal) <- c('Odds_Ratio', 'pvalue', 'padj', 'CI95%L', 'CI95%R')
finalres <- cbind(fisher_re_com, fisher_test_comfinal)
data.table::fwrite(finalres, file = "./miRNA_mRNA_RShinyApp/allfisher_test_IndimiRs_Enrich_skinW1W7CW.txt",sep = "\t")


################################################################################
####------Individual miRNAs Targets enrichment analysis in each module------####
################################################################################
names <- unique(mRNAmodules$modules)
namesmiR <- unique(miRNAmodules$modules)
#i=1
#j=1
#z=1
fisher_refinal <- list()
for (i in 1:length(namesmiR)) {
  require(tidyverse)
  module <- namesmiR[i]
  moduleData <- miRNAmodules %>% filter(modules == module)
  
  miRls <- unique(moduleData$mature_mirna_id)
  
  fisher_re <- list()
  for (j in 1:length(miRls)) {
    miRtmp <- miRls[j]
    mirData <- allmiRNA_targets %>% filter(mature_mirna_id == miRtmp)
    
    fisher_overlap <- list()
    for (z in 1:length(names)) {
      module2 <- names[z]
      moduleData2 <- mRNAmodules %>% filter(modules == module2)
      genesOverlapps <- intersect(moduleData2$GeneSymbol, mirData$target_symbol)
      if(is_empty(genesOverlapps)){
        fisher_overlap[[module2]] <- "emptyOverlap"
      }else{
        fisher_overlap[[module2]] <- genesOverlapps
      }
    }
    fisher_overlapCom <- do.call("cbind",fisher_overlap)
    datatmp <- list()
    for (k in seq_along(colnames(fisher_overlapCom))) {
      mrnaName <- colnames(fisher_overlapCom)[k]
      datatmp[[mrnaName]] <- data.frame(mRNAmodule=mrnaName,
                                        OverlapGenes=unique(fisher_overlapCom[,k]))
    }
    datatmp_c <- do.call("bind_rows",datatmp) %>% as.data.frame()
    fisher_re[[miRtmp]] <- data.frame(miRNAid=miRtmp, 
                                      Overlap = datatmp_c)
  }
  fisher_recom <- do.call("bind_rows",fisher_re) %>% as.data.frame()
  fisher_refinal[[module]] <- data.frame(miRNAmodule=module, fisher_recom)
}

fisher_re_com <- do.call("bind_rows",fisher_refinal) %>% as.data.frame()
unique(fisher_re_com$miRNAmodule)
length(unique(fisher_re_com$miRNAid))

colnames(fisher_re_com) <- c('miRNAmodule','miRNA_ID','mRNA_module','target_symbol')
data.table::fwrite(fisher_re_com, file = "miRNA_mRNA_RShinyApp/allfisher_test_IndimiRs_Enrich_skinW1W7CW_OverlapGenes.txt",sep = "\t")

#mRNAmodulestmp <- bind_rows(
#  data.frame(modules = 'cw_up', cw_up[,c(1,17)]),
#  data.frame(modules = 'cw_down', cw_down[,c(1,17)]),
#  data.frame(modules = 'wound7_down', wound7_down[,c(1,10)]),
#  data.frame(modules = 'wound7_up', wound7_up[,c(1,10)]),
#  data.frame(modules = 'wound1_down', wound1_down[,c(1,10)]),
#  data.frame(modules = 'wound1_up', wound1_up[,c(1,10)]),
#  data.frame(modules = 'skin_down', skin_down[,c(1,10)]),
#  data.frame(modules = 'skin_up', skin_up[,c(1,10)]),
#  data.frame(modules = 'M1', mRNA_m1),
#  data.frame(modules = 'M2', mRNA_m2),
#  data.frame(modules = 'M3', mRNA_m3),
#  data.frame(modules = 'M4', mRNA_m4),
#  data.frame(modules = 'M5', mRNA_m5),
#  data.frame(modules = 'M6', mRNA_m6),
#  data.frame(modules = 'M7', mRNA_m7),
#  data.frame(modules = 'M8', mRNA_m8),
#  data.frame(modules = 'M9', mRNA_m9),
#  data.frame(modules = 'M10', mRNA_m10),
#  data.frame(modules = 'M11', mRNA_m11),
#  data.frame(modules = 'M12', mRNA_m12),
#  data.frame(modules = 'M13', mRNA_m13))

#mRNAmodules1 <- mRNAmodulestmp %>% select(1,3) %>% distinct(GeneSymbol, .keep_all = TRUE)
fisher_re_com_f <- fisher_re_com %>% left_join(.,mRNAmodules[,c(1,3)],by=c("target_symbol"="GeneSymbol"))
DEmiRNAsall <- rbind(
  data.frame(miRTrait= "CW", DEmiRNA_specific = "Up", ID=cwmiRNA_up$ID),
  data.frame(miRTrait= "CW", DEmiRNA_specific = "Down", ID=cwmiRNA_down$ID),
  data.frame(miRTrait= "Wound7", DEmiRNA_specific = "Up", ID=w7miRNA_up$ID),
  data.frame(miRTrait= "Wound7", DEmiRNA_specific = "Down", ID=w7miRNA_down$ID),
  data.frame(miRTrait= "Wound1", DEmiRNA_specific = "Up", ID=w1miRNA_up$ID),
  data.frame(miRTrait= "Wound1", DEmiRNA_specific = "Down", ID=w1miRNA_down$ID),
  data.frame(miRTrait= "Skin", DEmiRNA_specific = "Up", ID=skinmiRNA_up$ID),
  data.frame(miRTrait= "Skin", DEmiRNA_specific = "Down", ID=skinmiRNA_down$ID))
data.table::fwrite(DEmiRNAsall, file = "./miRNA_mRNA_RShinyApp/all_DEmiRNAs.txt",sep = "\t")

DEmRNAsall <- rbind(
  data.frame(mRNATrait= "CW", DEmRNA_specific = "Up", ID=cw_up$GeneSymbol),
  data.frame(mRNATrait= "CW", DEmRNA_specific = "Down", ID=cw_down$GeneSymbol),
  data.frame(mRNATrait= "Wound7", DEmRNA_specific = "Up", ID=wound7_up$GeneSymbol),
  data.frame(mRNATrait= "Wound7", DEmRNA_specific = "Down", ID=wound7_down$GeneSymbol),
  data.frame(mRNATrait= "Wound1", DEmRNA_specific = "Up", ID=wound1_up$GeneSymbol),
  data.frame(mRNATrait= "Wound1", DEmRNA_specific = "Down", ID=wound1_down$GeneSymbol),
  data.frame(mRNATrait= "Skin", DEmRNA_specific = "Up", ID=skin_up$GeneSymbol),
  data.frame(mRNATrait= "Skin", DEmRNA_specific = "Down", ID=skin_down$GeneSymbol))
data.table::fwrite(DEmRNAsall, file = "./miRNA_mRNA_RShinyApp/all_DEmRNAs.txt",sep = "\t")

fisher_re_com_f2 <- fisher_re_com_f %>% left_join(.,DEmiRNAsall,by=c('miRNA_ID'='ID')) %>% 
  left_join(.,DEmRNAsall,by=c("target_symbol"="ID"))
data.table::fwrite(fisher_re_com_f2, file = "allfisher_test_IndimiRs_Enrich_skinW1W7CW_OverlapGenes.txt",sep = "\t")


######################################################
####----step4 Draw the plot of DE mRNA results----####
######################################################
##Draw the plot of DE mRNA results of CW
cw_DEmRNA_combined_new <- rbind(cw_DEmRNA_combined_up, cw_DEmRNA_combined_down)
cw_DEmRNA_combined_FPKM <- wb_allgene_NewDE_mRNAresults$FPKM_wb %>% select(-1) %>% select(6:20,1:5)
cw_DEmRNA_combined_FPKM <- cw_DEmRNA_combined_FPKM[rownames(cw_DEmRNA_combined_FPKM)%in%cw_DEmRNA_combined_new$ID,]
DE_zscore <- t(scale(t(cw_DEmRNA_combined_FPKM)))
colnames(DE_zscore) <- gsub("_FPKM","",colnames(DE_zscore))

p1 <- Heatmap(DE_zscore, name = "Z-scoreTPM",
              show_column_names = T, show_row_names = F,
              #left_annotation = ha_row,
              #right_annotation = ha_row_label,
              cluster_columns = FALSE,
              row_km = 2,
              #column_km = 4,
              col= colorRamp2(c(-2,0,2),c("#2171b5","#F7F7F7","#B2182B")), use_raster = TRUE)
print(p1)
pdf("DEmRNA_Skin_Markers_heatmap.pdf",useDingbats = F,width = 9,height = 9)
p1
dev.off()

rm(cw_DEmRNA_combined_FPKM,DE_zscore,p1,cw_DEmRNA_combined,cw_DEmRNA_combined_new,cw_skin,cw_w1,cw_w7)

####----step5 Draw the plot of gene set enrichment results----####
library(WGCNA)
library(RColorBrewer)
gsea_results <- read.table("clipboard", header = T, sep = "\t", row.names = 1)
gsea_oddsratio <- gsea_results[,c(seq(1,13,by=2))]
gsea_pvalue <- gsea_results[,c(seq(2,14,by=2))]
gsea_pvaluelog <- -log10(gsea_pvalue)
colnames(gsea_pvaluelog) <- gsub("\\.", " ", colnames(gsea_oddsratio))

gsea_pvaluelog_f <- t(gsea_pvaluelog)
p1 <- Heatmap(gsea_pvaluelog_f, name = "-log10(Pvalue)",
        show_column_names = T, show_row_names = T,
        row_names_side = "left",
        column_names_rot = 45,
        cluster_columns = FALSE, cluster_rows = FALSE,
        col= colorRamp2(c(0,2,4,6),c("white","#ffd8d8","#ff7676","red")),
        border = TRUE)
p1
pdf("GSEA_alltop25_heatmap.pdf",useDingbats = F,width = 9,height = 9)
p1
dev.off()

##Filter the miRNAS
library(reshape2)

final.mir <- m8_m9miRNAgsea %>% filter(Odds_ratio > 1 & P_value < 0.05)
m89updown <- m8_m9miRNAgsea[m8_m9miRNAgsea$miRNA_ID%in%final.mir$miRNA_ID,]
m89updown.f <- m89updown %>% select(2:3,11:12) %>% distinct()

data.table::fwrite(m89updown.f, file = "fisher_test_results_all_indimiR_M8M9UpDownmiRNAresults_fil.txt",sep = "\t")

m89updown.fodd <- dcast(m89updown.f[,1:3], mRNA_module ~ miRNA_ID)
m89updown.fpvalue <- dcast(m89updown.f[,c(1:2,4)], mRNA_module ~ miRNA_ID)

colnames(m89updown.fpvalue) <- paste0(colnames(m89updown.fpvalue), 'pvalue')

textMatrix = paste(signif(unlist(m89updown.fodd[,-1]), 3), "\n(",
                   signif(unlist(m89updown.fpvalue[,-1]), 4), ")", sep = "")
textMatrix_f <- matrix(textMatrix, nrow=7)
colnames(textMatrix_f) <- colnames(m89updown.fodd)[-1]
textMatrix_f2 <- t(data.frame(ID=m89updown.fodd$mRNA_module, textMatrix_f))
colnames(textMatrix_f2) <- textMatrix_f2[1,]
textMatrix_f2 <- textMatrix_f2[-1,c(1,5,7,6,2:4)] %>% as.data.frame()
rownames(textMatrix_f2) <- gsub('\\.','-',rownames(textMatrix_f2))
textMatrix_f3 <- textMatrix_f2 %>% rownames_to_column(.,var='ID')
textMatri <- modules %>% left_join(.,textMatrix_f3,by=c('miRNA_ID'='ID'))
textMatri_F <- textMatri[,-1:-2]

htdata <- cbind(m89updown.fodd,m89updown.fpvalue[,-1])
htdata_f <- htdata %>% slice(1,5,7,6,2:4)

gsea_oddsratio <- htdata_f %>% select(2:21)
gsea_pvalue <- htdata_f %>% select(22:41)
gsea_pvaluelog <- -log10(gsea_pvalue)
colnames(gsea_pvaluelog) <- colnames(gsea_oddsratio)
rownames(gsea_pvaluelog) <- htdata_f$mRNA_module
gsea_pvaluelog_f <- t(gsea_pvaluelog)
gsea_pvaluelog_f <- data.frame(ID=rownames(gsea_pvaluelog_f), gsea_pvaluelog_f)

modules <- read.table('clipboard',header = T,sep = '\t')
modules_f <- modules %>% left_join(.,gsea_pvaluelog_f,by=c('miRNA_ID'='ID')) %>% column_to_rownames(.,var = 'miRNA_ID') %>% 
  select(-1)
colnames(modules_f) <- gsub('\\.','-',colnames(modules_f))

RowAnn <- data.frame(modules$Module)
colnames(RowAnn) <- c("Module")
col_fun = list("Module"=c("Down"="#1b9e77","M9"="#d95f02","Up"="#7570b3","M8"="#e7298a"))
RowAnn <- HeatmapAnnotation(df=RowAnn, which="row")


p1 <- Heatmap(as.matrix(modules_f), name = "-log10(Pvalue)",
              show_column_names = T, show_row_names = T,
              row_names_side = "right",
              column_names_rot = 45,
              left_annotation = RowAnn,
              cluster_columns = FALSE, cluster_rows = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(textMatri_F[i, j] > 1)
                  grid.text(sprintf("%s", textMatri_F[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              col= colorRamp2(c(0,2,4,6),c("white","#ffd8d8","#ff7676","red")),
              border = TRUE)
p1
pdf("GSEA_m8m9updown_heatmap_ori2.pdf",useDingbats = F,width = 9,height = 9)
p1
dev.off()
