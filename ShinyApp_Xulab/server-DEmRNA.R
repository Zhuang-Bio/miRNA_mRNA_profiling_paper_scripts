observeEvent(input$DEmRNA.Plot,{
  withProgress(message = "Processing , please wait",{

    #Filter the data according to the GroupType of comparisons and different criteria
    if (input$DEmRNA.FCtype == "greatthan") {
      de_list_mRNA <- DEmRNA.df %>% filter(.data$log2FoldChange >= input$DEmRNA.foldChange &
                                       .data$Grouptype == input$DEmRNA.geneData &
                                       .data$padj < input$DEmRNA.pvalue)
    }
    if (input$DEmRNA.FCtype == "lessthan") {
      de_list_mRNA <- DEmRNA.df %>% filter(.data$log2FoldChange <= input$DEmRNA.foldChange &
                                       .data$Grouptype == input$DEmRNA.geneData &
                                       .data$padj < input$DEmRNA.pvalue)
    }
    if (input$DEmRNA.FCtype == "both") {
      de_list_mRNA <- DEmRNA.df %>% filter(abs(.data$log2FoldChange) >= input$DEmRNA.foldChange &
                                       .data$Grouptype == input$DEmRNA.geneData &
                                       .data$padj < input$DEmRNA.pvalue)
    }

    #Further clean the data after filtering and obtain the final expression data
    exptmpmRNA <- mRNA.exp[mRNA.exp$EnsemblID %in% de_list_mRNA$EnsemblID,] %>% 
      distinct(.data$GeneSymbol, .keep_all = TRUE) #remove the duplicates if their gene symbols are same
    #exptmpmRNA <- miRNA.exp %>% left_join(., de_list_mRNA, by=c("miRNA" = "miRNA"))

    #Extract the metadata info (Groups SampleID and gene data)
    groupCommRNA <- str_split(input$DEmRNA.geneData, "vs")[[1]]
    #First group data
    exp1mRNA <- exptmpmRNA[,grep(groupCommRNA[1], colnames(exptmpmRNA))]
    #Second group data
    exp2mRNA <- exptmpmRNA[,grep(groupCommRNA[2], colnames(exptmpmRNA))]
    #Combined the data again
    exp_commRNA <- as.matrix(cbind(exp1mRNA, exp2mRNA))
    rownames(exp_commRNA) <- exptmpmRNA$GeneSymbol
    #Calculate Z scores for rows of expression matrix
    mRNAZscore <- t(scale(t(exp_commRNA), center = TRUE, scale = TRUE))

    if(nrow(mRNAZscore) < 2){
      showNotification("Genes less than 2.", type="error")
      return()
    }
    
    ##########################################complex heatmap#######################################
    # prepare pre-data 
    rowname1 <- input$DEmRNA.rownames
    if (rowname1 == "None") {
      show_rowname1 <- FALSE
    }else {
      show_rowname1 <- TRUE
    }   

    columnname1 <- input$DEmRNA.colnames
    if (columnname1 == "None") {
      show_columnname1 <- FALSE
    }else {
      show_columnname1 <- TRUE
    }   

    #If cluster the gene and cells/samples
    rowcluster1 <- input$DEmRNA.clusterGene
    if (rowcluster1 == "No") {
      cluster_row1 <- FALSE
    }else {
      cluster_row1 <- TRUE
    }   

    columncluster1 <- input$DEmRNA.clusterSample
    if (columncluster1 == "No") {
      cluster_column1 <- FALSE
    }else {
      cluster_column1 <- TRUE
    }   

    #Extract the metadata info (Groups SampleID and gene data)
    mRNAgroup1.name <- as.character(groupCommRNA[1])
    mRNAgroup2.name <- as.character(groupCommRNA[2])
    #Create the annotation groups
    group12name <- c(rep(mRNAgroup1.name, 5), rep(mRNAgroup2.name, 5))
    #Make color for each group
    tmplist <- list(Group = c(mRNAgroup1.name = "#ff7f00", mRNAgroup2.name = "#4daf4a"))
    #Rename the group name
    names(tmplist$Group) <- c(mRNAgroup1.name, mRNAgroup2.name)
    anno_column1=HeatmapAnnotation(Group = group12name,
      col = tmplist,
      annotation_name_side = "left")

    ######################
    ###Draw the heatmap###
    ######################
    DEmRNAheatmap <- DEmRNAheatmap_f(
      matZval = mRNAZscore, 
      columnName = show_columnname1, 
      rowName = show_rowname1, 
      columnClu = cluster_column1, 
      rowClu = cluster_row1, 
      annoColumn = anno_column1
    )
    
    output$demRNAPlot1 <- renderPlot({
      draw(DEmRNAheatmap)
    })
    
    output$demRNATable <- DT::renderDataTable({
      outdata_demRNA <- copy(de_list_mRNA) #data you want to show
      #outdata_demRNA$GeneSymbol <- createLink(outdata_demRNA$GeneSymbol, ensembl_LINK, outdata_demRNA$GeneSymbol)
      #outdata_demRNA <- outdata_demRNA %>% mutate_if(is.numeric, signif, 6) %>% select(1:2, 5:10, 3:4)#set the output decimals
      return(outdata_demRNA)},
      escape = FALSE,server = FALSE,extensions = c("Buttons"), 
      options = list(dom = 'Bfrtip',buttons = c('copy', 'csv')), #'excel': some of the padj values in output are strings. So we suggest to choose the csv format.
      rownames= FALSE)

    output$DEmRNA.UI <- renderUI({
      tabBox(
        title = "", width = 12,
        tabPanel("Heatmap", plotOutput("demRNAPlot1", height = 2000)),
        tabPanel("Data Output", DT::dataTableOutput("demRNATable"), style = "overflow-x: scroll;")
      )
    })

    output$downloadDEmRNAHeatmap <- downloadHandler(
      filename <- function() {
          paste0("Heatmap_DEmRNA", ".pdf")
      }, 
      content <- function(file) {
        pdf(file, width = 8, height = 8)
          print(DEmRNAheatmap)
        dev.off()
      }
    )

  })#withProgress
})#observeEvent



