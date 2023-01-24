mRNAtrait_plot <- Heatmap(
  as.matrix(mRNAmodule_trait), name = "Cor\nr",
  show_column_names = T, show_row_names = T,
  row_names_side = "left",
  column_names_rot = 45,
  cluster_columns = FALSE, cluster_rows = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(mRNAmodule_trait_pvalues[i, j] < 0.05)
      grid.text(sprintf("%s", mRNAmodule_trait_textcor[i, j]), x, y, gp = gpar(fontsize = 14))
  },
  col= colorRamp2(c(-1,0,1),c("#2171b5","#F7F7F7","#B2182B")),
  border = TRUE,
  column_names_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 16, fontface = "bold")
)
output$mRNAmoduleTraitplot <- renderPlot({
  draw(mRNAtrait_plot)
})

output$mRNA_WGCNAimg <- renderUI({
  tabBox(
    title = "", width = 12,
    tabPanel("mRNA WGCNA", plotOutput("mRNAmoduleTraitplot", height = 800))
  )
})

updateSelectizeInput(session, "mRNA.modules", choices = unique(sort(mRNAkMEs$moduleLabel)), server = TRUE)

observeEvent(input$mRNAWGCNAPlot1,{
  withProgress(message = "Processing , please wait",{ 

    #select the top N kME mRNAs to draw the heatmap
    if (input$mRNA.KMEs == 10) {
        moduleDatamRNA <- mRNAkMEs %>% filter(.data$moduleLabel == input$mRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$mRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:10) %>% #select the top N rows data
            select(1:4, grep(pattern = paste0("\\b", input$mRNA.modules, "\\b") , x =  colnames(mRNAkMEs))) %>% ##grep the index of selected module
            left_join(., mRNA.exp[,-1], by=c("ID"="EnsemblID")) #add the mRNA expression data
        colnames(moduleDatamRNA)[5] <- paste0(colnames(moduleDatamRNA)[5], "_kMEs")
    }
    if (input$mRNA.KMEs == 20) {
        moduleDatamRNA <- mRNAkMEs %>% filter(.data$moduleLabel == input$mRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$mRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:20) %>% #select the top N rows data
            select(1:4, grep(pattern = paste0("\\b", input$mRNA.modules, "\\b") , x =  colnames(mRNAkMEs))) %>% ##grep the index of selected module
            left_join(., mRNA.exp[,-1], by=c("ID"="EnsemblID")) #add the mRNA expression data
            colnames(moduleDatamRNA)[5] <- paste0(colnames(moduleDatamRNA)[5], "_kMEs")
    }
    if (input$mRNA.KMEs == 30) {
        moduleDatamRNA <- mRNAkMEs %>% filter(.data$moduleLabel == input$mRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$mRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:30) %>% #select the top N rows data
            select(1:4, grep(pattern = paste0("\\b", input$mRNA.modules, "\\b") , x =  colnames(mRNAkMEs))) %>% ##grep the index of selected module
            left_join(., mRNA.exp[,-1], by=c("ID"="EnsemblID")) #add the mRNA expression data
            colnames(moduleDatamRNA)[5] <- paste0(colnames(moduleDatamRNA)[5], "_kMEs")
    }
    if (input$mRNA.KMEs == 40) {
        moduleDatamRNA <- mRNAkMEs %>% filter(.data$moduleLabel == input$mRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$mRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:40) %>% #select the top N rows data
            select(1:4, grep(pattern = paste0("\\b", input$mRNA.modules, "\\b") , x =  colnames(mRNAkMEs))) %>% ##grep the index of selected module
            left_join(., mRNA.exp[,-1], by=c("ID"="EnsemblID")) #add the mRNA expression data
            colnames(moduleDatamRNA)[5] <- paste0(colnames(moduleDatamRNA)[5], "_kMEs")
    }    
    
    mRNAhtplot <- as.matrix(moduleDatamRNA[,7:26], rownames = moduleDatamRNA$GeneSymbol)
    mRNAhtplotZscore <- t(scale(t(mRNAhtplot), center = TRUE, scale = TRUE))

    anno_columnmRNA <- data.frame(Group = factor(rep(c("Skin", "Wound1", "Wound7", "VU"), each = 5)), 
        row.names = colnames(mRNAhtplot))

    mRNAWGCNAheatmap <- DEcomplexheatmap(
      matZval = mRNAhtplotZscore, 
      columnName = TRUE, 
      rowName = TRUE, 
      columnClu = FALSE, 
      rowClu = TRUE, 
      annoColumn = anno_columnmRNA#, 
      #annoCol = anno_color
    )
    
    output$mRNAmoduleplot <- renderPlot({
      mRNAWGCNAheatmap
    })

    #Set plot height using gtable package functions
    gtable_h3 <- convertHeight(gtable_height(mRNAWGCNAheatmap$gtable),"inches",valueOnly = TRUE)
    legend_row3 <- length(unique(anno_columnmRNA$Group))
    h3 <- max(gtable_h3, legend_row3/2 + 1.5)

    plot_height3 <- h3 * 70 + 150
    plotDE_pdf_height3 <<- h3 + 0.5 #the output heatmap

    output$mRNAWGCNATable <- DT::renderDataTable({
      outdatamRNA <- copy(moduleDatamRNA) #data you want to show
      outdatamRNA$GeneSymbol <- createLink(outdatamRNA$GeneSymbol, ensembl_LINK, outdatamRNA$GeneSymbol)
      outdatamRNA <- outdatamRNA %>% mutate_if(is.numeric, signif, 8) #set the output decimals
      colnames(outdatamRNA)[7:26] <- paste0(colnames(outdatamRNA)[7:26], "_FPKM")
      return(outdatamRNA)},
      escape = FALSE,server = FALSE,extensions = c("Buttons"), 
      options = list(dom = 'Bfrtip',buttons = c('copy', 'csv')), #'excel': some of the padj values in output are strings. So we suggest to choose the csv format.
      rownames= FALSE)

    output$mRNAWGCNA.UI <- renderUI({
      tabBox(
        title = "", width = 12,
        tabPanel("Module Heatmap", plotOutput("mRNAmoduleplot", height = plot_height3)),
        tabPanel("Data Output", DT::dataTableOutput("mRNAWGCNATable"), style = "overflow-x: scroll;") #height:800px; overflow-y: scroll;
      )
    })

    output$downloadmRNAwgcnaPlot <- downloadHandler(
      filename <- function() {
          paste0("Heatmap_mRNAmodules", ".pdf")
      }, 
      content <- function(file) {
        pdf(file, width = 10, height = plotDE_pdf_height3)
          grid::grid.newpage()
          grid::grid.draw(mRNAWGCNAheatmap$gtable)
        dev.off()
      }
    )

  })#withProgress
})#observeEvent
