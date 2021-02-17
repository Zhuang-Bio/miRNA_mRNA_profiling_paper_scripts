library(multiMiR)
library(tidyverse)
rm(list = ls())
setwd("C:/Users/zhuliu/Desktop/miRNA/s5miRNA_TargetPredict/s1_wb_miRNA_targets_multiMiR")
allkme_miR <- read_tsv("wb_allmiRNA_TPM_SignedKME.txt")
mir100 <- allkme_miR %>% column_to_rownames(var="ID")
mir100_target <- rownames(mir100)
# Plug miRNA's into multiMiR and getting validated targets
#system.time(
#  multimir_results_all_cons25 <- get_multimir(org     = 'hsa',
#                                                 mirna   = mir100_target,
#                                                 table   = 'predicted',
#                                                 predicted.site = "conserved",
#                                                 predicted.cutoff  = 25,
#                                                 predicted.cutoff.type = "p",
#                                                 summary = FALSE)
#)
#save(multimir_results_all_cons25, file = "multimir_results_mirall_cons25.RData")

system.time(
  multimir_results_all_all25 <- get_multimir(org     = 'hsa',
                                              mirna   = mir100_target,
                                              table   = 'predicted',
                                              predicted.site = "all",
                                              predicted.cutoff  = 25,
                                              predicted.cutoff.type = "p",
                                              summary = FALSE)
)
save(multimir_results_all_all25, file = "multimir_results_mirall_all25.RData")

#load("../s1_wb_miRNA_targets_multiMiR/multimir_results_mirall_cons25.RData")
load("../s1_wb_miRNA_targets_multiMiR/multimir_results_mirall_all25.RData")
all562miR_target_cons25_re <- multimir_results_all_all25@data
length(unique(all562miR_target_cons25_re$target_symbol)) 
all562miR_target_cons25_re_TS <- all562miR_target_cons25_re %>% filter(database == "targetscan") %>%
  distinct(mature_mirna_id, target_symbol,.keep_all = TRUE)
length(unique(all562miR_target_cons25_re_TS$target_symbol)) 
tmp <- table(all562miR_target_cons25_re_TS$target_symbol)
tmp_names <- names(which(tmp >= 2))
all562miR_target_cons25_re_TS <- all562miR_target_cons25_re_TS[all562miR_target_cons25_re_TS$target_symbol%in%tmp_names,]
length(unique(all562miR_target_cons25_re_TS$target_symbol)) #18653
#read the expression file of mRNA (did not filter the expression FPKM greater than 1)
mRNA_expressed10sample <- read_tsv("C:/Users/zhuliu/Desktop/miRNA/s5miRNA_TargetPredict/s1_wb_miRNA_targets_multiMiR/wb_allGenesFPKM_mRNA_nofilterFPKM1_justhalfExpressed.txt")
mRNA_expressed10sample <- mRNA_expressed10sample %>% left_join(., wb_allgene_NewDE_mRNAresults$HumanGene_PC_v32[,1:2],by=c("ID"="Geneid"))
all562miR_target_cons25_re_TS_f <- all562miR_target_cons25_re_TS[all562miR_target_cons25_re_TS$target_symbol%in%mRNA_expressed10sample$GeneSymbol,]
length(unique(all562miR_target_cons25_re_TS_f$target_symbol)) #15969

load("C:/Users/zhuliu/Desktop/miRNA/s3DEanalysis/DEresults_wb_mRNA/wb_allgene_NewDE_mRNAresults.RData")
exp_mRNA <- wb_allgene_NewDE_mRNAresults$readCount_wb %>%
  left_join(.,wb_allgene_NewDE_mRNAresults$HumanGene_PC_v32[,1:2],by=c("ID"="Geneid")) #%>% select(1,22)
##filtering the targets
all562miR_target_cons25_re <- all562miR_target_cons25_re[all562miR_target_cons25_re$target_symbol%in%exp_mRNA$GeneSymbol,]
#check the length of origin backgroud gene list
length(unique(all562miR_target_cons25_re$target_symbol)) #11930

all562miR_target_cons25_re_TS <- all562miR_target_cons25_re %>% filter(database == "targetscan") %>%
  distinct(mature_mirna_id, target_symbol,.keep_all = TRUE)
all562miR_target_cons25_re_TS$score <- parse_double(all562miR_target_cons25_re_TS$score)
min(all562miR_target_cons25_re_TS$score);max(all562miR_target_cons25_re_TS$score)
length(unique(all562miR_target_cons25_re_TS$target_symbol)) #6361
data.table::fwrite(all562miR_target_cons25_re_TS,file = "all562miR_target_allTop25_TargetScan_miRtogether.txt",sep = "\t")
