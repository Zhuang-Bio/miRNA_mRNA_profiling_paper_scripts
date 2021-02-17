library(tidyverse)
library(WGCNA)
library(data.table)
setwd('C:/Users/zhuliu/Desktop/miRNA/a_Manuscript/AI figures')
##Use the all combined data for GO results
data <- read.table('clipboard', sep = '\t', header = T)

go10bp <- data %>% slice(1:10)
go10bp$description <- factor(go10bp$description, levels = rev(go10bp$description))

plotfun <- function(data1, modules, colorDe){
  plot1 <- ggplot(data = data1) +
    geom_point(aes(x=description, y=enrichmentRatio, 
                   size = overlap), stat='identity', colour = colorDe)  +
    labs(title= modules, x=NULL) + 
    #scale_color_continuous(low = "blue", high = "red") +
    #scale_fill_gradient(low = "white", high = "red")+
    coord_flip() +
    theme(
      plot.background = element_rect(fill = NA), #color="grey50",size=2
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "black",fill = NA,size = 1),
      panel.grid.major = element_line(colour = "grey90",size = 0.5,linetype = 1),
      axis.text.x = element_text(size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      legend.position = "bottom",#top,bottom,left,right
      legend.spacing = unit(0,"cm"))
  return(plot1)
}

names <- unique(data$Module)
#using the mRNA WGCNA module color instruction file
colorlab <- read.table('clipboard', sep = '\t', header = T)
plotsall <- list()
for (i in seq_along(names)) {
  mod <- names[i]
  colorlab_f <- colorlab[colorlab$Label == mod,]$Color
  datause <- data %>% filter(Module == mod)
  plotsall[[mod]] <- plotfun(datause, mod, colorlab_f)
}

plotsall[['M6']]

library(patchwork)
pdf('Combined_enrichment_New.pdf', useDingbats = FALSE, height = 20, width = 22)
plotsall[['M1']] + plotsall[['M2']] + plotsall[['M3']] + plotsall[['M4']] + plotsall[['M5']] + 
  plotsall[['M6']] + plotsall[['M7']] + plotsall[['M8']] + plotsall[['M9']] + plotsall[['M10']] + 
  plotsall[['M11']] + plotsall[['M12']] + plotsall[['M13']] +
  plot_layout(ncol = 3, widths = c(1,1,1)) #, guides = 'collect'
dev.off()

#pdf('Combined_enrichment_New.pdf', useDingbats = FALSE, height = 20, width = 22)
#M1 + M2 + M3 + M4 + M5 + M6 + M7 + M8 + M9 + M10 + M11 + M12 + M13 +
#  plot_layout(ncol = 3, widths = c(1,1,1)) #, guides = 'collect'
#dev.off()


#######################################################################
######------miRNA module to mRNA module targets correlation------######
#######################################################################
##Read the top miRNAs in each module
library(data.table)
miR_tops <- fread('prepares/Top_miRNAs_eachModule.txt')
miR_tops_target <- fread('prepares/Top_miRNA_targets.txt')
##read the expression values of miRNA and mRNA from supplementary Table 
miRNAexp <- readxl::read_xlsx("Table S2 DEresults miRNA.xlsx", sheet = 1)
#miRNAexp_logTransform <- log2(miRNAexp[,-1] + 1)
#miRNAexp_logTransform <- data.frame(ID=miRNAexp$ID, miRNAexp_logTransform)

mRNAexp <- readxl::read_xlsx("Table S3 DEresults mRNA.xlsx", sheet = 1)
#mRNAexp_logTransform <- log2(mRNAexp[,-1:-2] + 1)
#mRNAexp_logTransform <- data.frame(GeneID=mRNAexp$GeneID, GeneSymbol=mRNAexp$GeneSymbol, mRNAexp_logTransform)

module_overGene <- fread('prepares/module_module_overlappedGenes.txt')
miRmod_color <- read.table('C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/s1_wb_miRNA_WGCNA/wb_allmiRNA_TPM_Color_Instruction_for_module.txt', sep = '\t', header = T)
mRNAmod_color <- read.table('C:/Users/zhuliu/Desktop/miRNA/s4WGCNA_analysis/s2_wb_mRNA_WGCNA/wb_allmRNA_Color_Instruction_for_module.txt', sep = '\t', header = T)

#function to calculate the correlation
miR_mod <- 'M3';m_mod <- 'M9'

moduleCor <- function(miR_mod, m_mod){
  tmp_miR <- miR_tops %>% select(1,2) %>% filter(modules == miR_mod)
  tmp_miR_exp <- tmp_miR %>% left_join(., miRmod_color[,1:3], by=c('modules'='Label')) %>% 
    left_join(., miRNAexp, by=c('ID' = 'ID'))
  datExpr_miR <- tmp_miR_exp %>% select(2,5:24) %>% column_to_rownames(var = 'ID') %>% t()
  modcolor <- as.character(tmp_miR_exp$Color)
  MEs_miR <- moduleEigengenes(datExpr_miR, modcolor)$eigengenes
  #tempa <- prcomp(t(tmp_miR_exp[,3:22]), center = TRUE)$x
  #tmpa_scale <- scale(t(tmp_miR_exp[,3:22]))
  #tmpa_svd <- svd(tmpa_scale)
  
  tmp_modGene <- module_overGene %>% filter(miRmodule == miR_mod & mRNAmoudle == m_mod)
  tmp_mRNA_exp <- tmp_modGene %>% left_join(., mRNAmod_color[,1:3], by=c('mRNAmoudle' = 'Label')) %>% 
    left_join(., mRNAexp, by=c('OverlapGene' = 'GeneSymbol')) %>% distinct(OverlapGene, .keep_all = TRUE) 
  datExpr_mRNA <- tmp_mRNA_exp %>% select(2,7:26) %>% column_to_rownames(var = 'OverlapGene') %>% t()
  mRNAmodcolor <- as.character(tmp_mRNA_exp$Color)
  MEs_mRNA <- moduleEigengenes(datExpr_mRNA, mRNAmodcolor)$eigengenes
  #tempb <- prcomp(t(tmp_mRNA_exp[,5:24]), center = TRUE)$x
  #tmpb_scale <- scale(t(tmp_mRNA_exp[,5:24]))
  #tmpb_svd <- svd(tmpb_scale)
  
  pvaluea <- cor.test(MEs_miR[,1],MEs_mRNA[,1])
  #pvaluea <- cor.test(tmpa_svd$u[,1],tmpb_svd$u[,1])
  out <- c(pvaluea$p.value, pvaluea$estimate)
  names(out) <- c('pvalue','cor')
  return(out)

}
moduleCor('M3', 'M9')
moduleCor('M3', 'M10')
moduleCor('M3', 'M11')
moduleCor('M3', 'M12')
moduleCor('M3', 'M5')

moduleCor('M7', 'M9')
moduleCor('M7', 'M10')
moduleCor('M7', 'M11')
moduleCor('M7', 'M12')
moduleCor('M7', 'M5')

moduleCor('M9', 'M9')
moduleCor('M9', 'M10')
moduleCor('M9', 'M11')
moduleCor('M9', 'M12')
moduleCor('M9', 'M5')

moduleCor('M8', 'M5')
moduleCor('M8', 'M9')
moduleCor('M8', 'M10')
moduleCor('M8', 'M11')
moduleCor('M8', 'M12')

moduleCor('M12', 'M5')
moduleCor('M12', 'M9')
moduleCor('M12', 'M10')
moduleCor('M12', 'M11')
moduleCor('M12', 'M12')


########################################################
######------Correlate the DEmiRNAs to DEmRNA------######
########################################################
##using the targets of DEmiRNAs to DEmRNAs
DEmiRNA <- readxl::read_xlsx("Table S6 metacore.xlsx", sheet = 1)
DEmRNA <- readxl::read_xlsx("Table S6 metacore.xlsx", sheet = 2)

tmp_miR <- miR_tops %>% select(1,2) %>% filter(modules == "cw_down") #| modules == "cw_up"
tmp_miR_exp <- tmp_miR %>%
  left_join(., miRNAexp, by=c('ID' = 'ID'))
datExpr_miR <- tmp_miR_exp %>% select(2,3:22) %>% column_to_rownames(var = 'ID') %>% t()
modcolor <- as.character(rep("CWup_down", ncol(datExpr_miR)))
MEs_miR <- moduleEigengenes(datExpr_miR, modcolor)$eigengenes

tmp_modGene <- module_overGene %>% filter(miRmodule == "cw_down") %>% # | miRmodule == "cw_up"
  filter(mRNAmoudle == "cw_up") #mRNAmoudle == "cw_down" | 
tmp_mRNA_exp <- tmp_modGene %>%
  left_join(., mRNAexp, by=c('OverlapGene' = 'GeneSymbol')) %>% distinct(OverlapGene, .keep_all = TRUE) 
datExpr_mRNA <- tmp_mRNA_exp %>% select(2,5:24) %>% column_to_rownames(var = 'OverlapGene') %>% t()
mRNAmodcolor <- as.character(rep("CWup_down", ncol(datExpr_mRNA)))
MEs_mRNA <- moduleEigengenes(datExpr_mRNA, mRNAmodcolor)$eigengenes

pc1 <- cbind(MEs_miR,MEs_mRNA)
colnames(pc1) <- c("miRNA_PC1", "mRNA_PC1")
pvaluea <- cor.test(MEs_miR$MECWup_down, MEs_mRNA$MECWup_down)

plot3 <- ggplot(data = pc1, aes(x=miRNA_PC1, y=mRNA_PC1)) +
  geom_point(size = 2, stat='identity', colour = 'black')  +
  labs(title= paste0(as.character(pvaluea$p.value),' ', as.character(pvaluea$estimate))) + 
  geom_smooth(method="lm", formula = y ~ x, size=1.5, se=TRUE) +
  #scale_x_continuous(limits = c(-5,5),breaks = c(-4,-3,-2,-1,0,1,2,3,4)) +
  #scale_y_continuous(limits = c(-40,20),breaks = c(-40,-20,0,20)) +
  theme_classic()+
  theme(aspect.ratio = 1)

library(patchwork)
pdf('Correlation_miRNA_mRNA.pdf', useDingbats = FALSE, height = 6, width = 12)
plot1 + plot2 + plot3 +
  theme(aspect.ratio = 1)
dev.off()


####################################################
######------GO plot for CW DEmRNAs genes------######
####################################################
data <- read.table('clipboard', sep = '\t', header = T)

go10bp <- data %>% slice(1:10)
go10bp$description <- factor(go10bp$description, levels = rev(go10bp$description))

plot2 <- ggplot(data = go10bp,aes(x=description, y=enrichmentRatio, 
                                  size = overlap, color=-log10(FDR))) +
  geom_point(stat='identity')  +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(0,25)) +
  coord_flip()+
  theme(
    plot.background = element_rect(fill = NA), #color="grey50",size=2
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(colour = "black",fill = NA,size = 1),
    panel.grid.major = element_line(colour = "grey90",size = 0.5,linetype = 1),
    axis.text.x = element_text(size = 12, color = "black"), 
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.position = "bottom",#top,bottom,left,right
    legend.spacing = unit(0,"cm"))
plot2
library(patchwork)
pdf('Combined_enrichment_CWDEmRNA.pdf', useDingbats = FALSE, height = 10, width = 8)
plot1 + plot2 +
  plot_layout(ncol = 1, guides = 'collect')
dev.off()


####################################################
######------GO line plot for WGCNA mRNAs------######
####################################################
goline <- fread('prepares/moduleTraitCor.txt')
library(reshape2)
library(directlabels)
library(patchwork)
goline_t <- melt(goline)
goline_t$variable <- as.factor(goline_t$variable)
goline_t$module <- factor(goline_t$module, levels = c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13"))

goline_t_tmp <- goline_t[goline_t$module %in% c("M2","M4"), ]
goline_t_tmp <- goline_t[goline_t$module %in% c("M1","M3","M5","M8"), ]
goline_t_tmp <- goline_t[goline_t$module %in% c("M7"), ]
goline_t_tmp <- goline_t[goline_t$module %in% c("M9","M11","M12"), ]

p4 <- ggplot(data = goline_t_tmp, aes(x=variable, y=value, group=module, color=module)) +
  geom_line(size=1) + geom_point(size=2) + #,linetype = "dashed"
  scale_y_continuous(limits = c(-1,1),breaks = seq(-1,1,by=0.5)) +
  scale_color_manual(values = c("turquoise","red","black","pink","magenta","purple","greenyellow","tan","salmon","blue","brown","yellow","green")) +
  geom_dl(aes(label = module), method = list(dl.combine("first.points", "last.points")), cex = 0.8) +
  theme_classic() +
  theme(
    legend.position = 'right',
    panel.grid.major = element_line(colour = "grey90",size = 0.5,linetype = 1))
pdf('linePlot_WGCNAmRNA3.pdf', useDingbats = FALSE, height = 8, width = 6)
p1 + p2 + p3 + p4 +
  plot_layout(ncol = 1, guides = 'collect')
dev.off()  


######################################################
######------Violin plot for qRT-PCR miRNAs------######
######################################################
vioplotFun <- function(datadf,aesX, aesY){
  out <-  ggviolin(datadf, x = aesX, y = aesY, fill = aesX,
                   palette = c('#e41a1c','#377eb8','#4daf4a','#984ea3'),
                   add = "boxplot", add.params = list(fill = "white"))
  return(out)
}

miRNAs <- readxl::read_xlsx('../Table S2 DEresults miRNA.xlsx', sheet = 1)
miRNAplot <- read.table('clipboard')
vioplotdata <- miRNAplot %>% left_join(., miRNAs, by=c('V1' = 'ID')) %>% 
  slice(6) %>% 
  select(7:21,2:6) %>% 
  pivot_longer(everything(), names_to = "SampleID" , values_to = "TPM") %>% 
  mutate(Group=factor(rep(c("Skin", "Wound1", "Wound7", "CW"), each=5), levels = c("Skin", "Wound1", "Wound7", "CW")))

vioplotdata$log2TPM <- log2(vioplotdata$TPM + 1)
vioplot <- vioplotFun(
  datadf = vioplotdata, 
  aesX = "Group",
  aesY = "log2TPM")
vioplot_6 <- ggpar(vioplot, ylab="log2(TPM+1)")
vioplot_6

library(patchwork)
pdf('violinPlot_DEmiRs.pdf', useDingbats = FALSE, height = 10, width = 20)
vioplot_1 + vioplot_2 + vioplot_3 + vioplot_4 + vioplot_5 + vioplot_6 +
  plot_layout(ncol = 3, widths = c(1,1,1), guides = 'collect') +
  theme(legend.position = 'right')
dev.off()  

#######################################################
######------Heatmap of GOs in M8 and M11M12------######
#######################################################
library(ComplexHeatmap)
library(circlize)
metacore_results <- read.table("clipboard", header = T, sep = "\t", row.names = 1)
metacore_fdrlog <- -log10(metacore_results)
metacore_fdrlog <- as.matrix(metacore_fdrlog)
p1 <- Heatmap(metacore_fdrlog, name = "-log10(FDR)",
              show_column_names = T, show_row_names = T,
              row_names_side = "left",
              column_names_rot = 45,
              cluster_columns = FALSE, cluster_rows = FALSE,
              col= colorRamp2(c(0,5,10),c("#148aff","#F7F7F7","#ff4000")),
              border = TRUE)
p1
pdf("metacore_heatmap_process.pdf",useDingbats = F,width = 15,height = 5)
p1
dev.off()










