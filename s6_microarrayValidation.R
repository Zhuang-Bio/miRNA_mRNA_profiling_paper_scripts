library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(easyGgplot2)
library(patchwork)
library(data.table)
rm(list = ls())

setwd("C:/Users/zhuliu/Desktop/miRNA/s7miRNA_validation/")
####----Example: mimics 7704_KC microarray data----####
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

####################################################
######------Gene set enrichment analysis------######
####################################################

####----step1 read the mRNA modules and DE mRNAs----####
wb_mRNA_wgcna <- read_tsv("C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/s2_wb_mRNA_WGCNA/wb_allmRNA_SignedKME_new.txt")
mRNA_m5 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M5")
mRNA_m9 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M9")
mRNA_m10 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M10")
mRNA_m11 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M11")
mRNA_m12 <- wb_mRNA_wgcna %>% filter(moduleLabel == "M12")
cw_up <- fread('../s5miRNA_TargetPredict/s2_wb_miRNA_targets_DE_modules_enrichment/wb_cw_DEmRNA_combined_up.txt')
cw_down <- fread('../s5miRNA_TargetPredict/s2_wb_miRNA_targets_DE_modules_enrichment/wb_cw_DEmRNA_combined_down.txt')

mRNAmodules <- bind_rows(
  data.frame(modules = 'M5', mRNA_m5),
  data.frame(modules = 'M9', mRNA_m9),
  data.frame(modules = 'M10', mRNA_m10),
  data.frame(modules = 'M11', mRNA_m11),
  data.frame(modules = 'M12', mRNA_m12),
  data.frame(modules = 'cw_up', cw_up),
  data.frame(modules = 'cw_down', cw_down))
table(mRNAmodules$modules)

####---step2 read the module and DE miRNA target genes----####
allmiRNA_targets <- read_tsv("C:/Users/zhuliu/Desktop/miRNA/s5miRNA_TargetPredict/s1_wb_miRNA_targets_multiMiR/all562miR_target_allTop25_TargetScan_miRtogether.txt")
length(unique(allmiRNA_targets$target_symbol))
min(allmiRNA_targets$score);max(allmiRNA_targets$score)

microarrayVal <- function(inputDEfile, miRname){
  #Only keep the target of interested miRNAs
  miR_targets_valmic <- allmiRNA_targets[allmiRNA_targets$mature_mirna_id%in%miRname,]
  
  #read the DE results of microarray data
  OE_DEG_miRvalmic <- read_tsv(inputDEfile)
  colnames(OE_DEG_miRvalmic)[1] <- "SYMBOL"
  OE_DEG_miRvalmic_overlap <- OE_DEG_miRvalmic[OE_DEG_miRvalmic$SYMBOL%in%miR_targets_valmic$target_symbol,]
  a <- data.frame(Type="The Strongest Targets predicted by TargetScan", OE_DEG_miRvalmic_overlap$logFC)
  colnames(a)[2] <- "FoldChange"
  OE_DEG_miRvalmic_other <- OE_DEG_miRvalmic[!OE_DEG_miRvalmic$SYMBOL%in%OE_DEG_miRvalmic_overlap$SYMBOL,] 
  b <- data.frame(Type="Genes not predicted to be miRNA targets by TargetScan", OE_DEG_miRvalmic_other$logFC)
  colnames(b)[2] <- "FoldChange"
  density_df <- rbind(a,b)
  mean(a$FoldChange); mean(b$FoldChange)
  ##add the wilcox.test p value
  wilcox.test(FoldChange ~ Type, data = density_df, alternative = "two.sided")
  tresult <- t.test(FoldChange ~ Type, data = density_df, alternative = "two.sided")
  p1 <- ggplot2.density(data=density_df, xName='FoldChange', groupName='Type',
                        groupColors=c('black', 'red'), showLegend=FALSE,legendPosition="top",
                        backgroundColor="white",addMeanLine=TRUE,
                        xtitle="mRNA log2(FoldChange)", ytitle="Density", 
                        mainTitle= paste0(miRname, " Overexpression \nt-test ", tresult$p.value),
                        removePanelGrid=TRUE,removePanelBorder=TRUE,
                        axisLine=c(0.3, "solid", "black"))#,xlim=c(-2.3,2), ylim=c(0,2))
  
  p2 <- ggplot(density_df, aes(FoldChange, color = Type)) + 
    stat_ecdf(geom = "line", size = 0.5)+
    scale_y_continuous(labels = scales::percent) +
    #scale_x_continuous(limits = c(-1.5,2.5),breaks = c(-1)) +
    scale_color_manual(values = c('black', 'red')) +
    theme_classic() +
    xlab("mRNA log2(FoldChange)") +
    ylab("Cumulative distribution") +
    theme(legend.position = "none")
  
  pdfout <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position='bottom')
  ggsave(paste0(miRname, "_Overexp.pdf"), plot = pdfout, width = 10, height = 6,units = "in", useDingbats = FALSE)
  
  ##Do the Gene set enrichment analysis
  OE_DEG_miRvalmic_overlap_psig <- OE_DEG_miRvalmic_overlap %>% filter(logFC <= -0.3 & P.Value < 0.05)
  miR_targets_valmic_psig <- miR_targets_valmic[miR_targets_valmic$target_symbol%in%OE_DEG_miRvalmic_overlap_psig$SYMBOL,]
  
  names <- unique(mRNAmodules$modules)
  
  fisher_re <- list()
  fisher_re_overGene <- list()
  for (i in 1:length(names)) {
    require(tidyverse)
    module <- names[i]
    moduleData <- mRNAmodules %>% filter(modules == module)
    NmRNA <- length(unique(moduleData$GeneSymbol))
    
    genesOverlapps <- intersect(miR_targets_valmic_psig$target_symbol, moduleData$GeneSymbol)
    if(is_empty(genesOverlapps)){
      fisher_re_overGene[[module]] <- "emptyOverlap"
    }else{
      fisher_re_overGene[[module]] <- genesOverlapps
    }
    
    overlapN <- length(unique(genesOverlapps))
    NmiRNA <- length(unique(miR_targets_valmic_psig$target_symbol))
    bgGenes = 15969
    fisher_re[[module]] <- c(overlapN, NmRNA, NmiRNA, bgGenes)
  }
  fisher_re_com <- do.call("bind_rows",fisher_re) %>% t() %>% as.data.frame() 
  #add the column names
  colnames(fisher_re_com) <- c('Overlap','mRNA_number','miRNAtargets','Background')
  
  fisher_re_com_f <- fisher_re_com %>% mutate(mRNAnew = mRNA_number - Overlap, 
                                              miRNAnew = miRNAtargets - Overlap,
                                              BGnew = Background - mRNA_number - miRNAtargets + Overlap) %>% 
    select(1,5:7)
  
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
  finalres <- cbind(fisher_re_com, fisher_test_comfinal) %>% rownames_to_column(var = "mRNA_module") %>% 
    mutate(miRNA = miRname) %>% select(11, everything())
  fwrite(finalres, file = paste0(miRname, "_Overexp_GSEAenrichResult.txt"),sep = "\t")
}

microarrayVal(inputDEfile = "DEG_all_218_5p_mimicsKC.txt", miRname = "hsa-miR-218-5p")
microarrayVal(inputDEfile = "DEG_all_149_5p_mimicsKC.txt", miRname = "hsa-miR-149-5p")
#microarrayVal(inputDEfile = "DEG_all_7704_mimicsKC.txt", miRname = "hsa-miR-7704")
microarrayVal(inputDEfile = "DEG_all_424_5p_mimicsKC.txt", miRname = "hsa-miR-424-5p")
microarrayVal(inputDEfile = "DEG_all_450_5p_mimicsKC.txt", miRname = "hsa-miR-450-5p")
microarrayVal(inputDEfile = "DEG_all_517b_3p_mimicsKC.txt", miRname = "hsa-miR-517b-3p")



