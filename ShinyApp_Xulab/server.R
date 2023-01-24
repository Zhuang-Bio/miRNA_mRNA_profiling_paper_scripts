server <- function(input, output, session) {
    
    DEmiRNA.df <- readRDS('data/DEmiRNA_df.rds')
    DEmRNA.df <- readRDS('data/DEmRNA_df.rds')
    miRNA.exp <- readRDS('data/miRNA_exp.rds')
    mRNA.exp <- readRDS('data/mRNA_exp.rds')
    miRNAkMEs <- readRDS('data/WGCNA_wb_allmiRNA_TPM_SignedKME.rds')
    mRNAkMEs <- readRDS('data/WGCNA_wb_allmRNA_SignedKME.rds')
    miRNAmodule_trait <- readRDS('data/miRNAmodule_trait.rds')
    miRNAmodule_trait_pvalues <- readRDS('data/miRNAmodule_traitpvalues.rds')
    miRNAmodule_trait_textcor <- readRDS('data/miRNAmodule_trait_textcor.rds')
    mRNAmodule_trait <- readRDS('data/mRNAmodule_trait.rds')
    mRNAmodule_trait_pvalues <- readRDS('data/mRNAmodule_traitpvalues.rds')
    mRNAmodule_trait_textcor <- readRDS('data/mRNAmodule_trait_textcor.rds')

    source("app-functions.R")

    source("server-miRNA-mRNAexpression.R", local = TRUE)
    source("server-DEmiRNA.R", local = TRUE)
    source("server-miRNA-WGCNA.R", local = TRUE)
    source("server-DEmRNA.R", local = TRUE)
    source("server-mRNA-WGCNA.R", local = TRUE)
    
    GotoTab <- function(name){
        shinyjs::show(selector = paste0("a[data-value=\"",name,"\"]"))
       
        shinyjs::runjs("window.scrollTo(0, 0)")
    }
}
