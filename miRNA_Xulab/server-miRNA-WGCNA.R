miRNAtrait_plot <- Heatmap(
  as.matrix(miRNAmodule_trait), name = "Cor\nr",
  show_column_names = T, show_row_names = T,
  row_names_side = "left",
  column_names_rot = 45,
  cluster_columns = FALSE, cluster_rows = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(miRNAmodule_trait_pvalues[i, j] < 0.05)
      grid.text(sprintf("%s", miRNAmodule_trait_textcor[i, j]), x, y, gp = gpar(fontsize = 14))
  },
  col= colorRamp2(c(-1,0,1),c("#2171b5","#F7F7F7","#B2182B")),
  border = TRUE,
  column_names_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 16, fontface = "bold")
)
output$miRmoduleTraitplot <- renderPlot({
  draw(miRNAtrait_plot)
})

output$miRNA_WGCNAimg <- renderUI({
  tabBox(
    title = "", width = 12,
    tabPanel("miRNA WGCNA", plotOutput("miRmoduleTraitplot", height = 800))
  )
})


updateSelectizeInput(session, "miRNA.modules", choices = unique(sort(miRNAkMEs$moduleLabel)), server = TRUE)

observeEvent(input$miRNAWGCNAPlot1,{
  withProgress(message = "Processing , please wait",{ 

    #select the top N kME miRNAs to draw the heatmap
    if (input$miRNA.KMEs == 10) {
        moduleData <- miRNAkMEs %>% filter(.data$moduleLabel == input$miRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$miRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:10) %>% #select the top N rows data
            select(1:3, grep(pattern = paste0("\\b", input$miRNA.modules, "\\b") , x =  colnames(miRNAkMEs))) %>% ##grep the index of selected module
            left_join(., miRNA.exp, by=c("ID"="miRNA")) #add the miRNA expression data
        colnames(moduleData)[4] <- paste0(colnames(moduleData)[4], "_kMEs")
    }
    if (input$miRNA.KMEs == 20) {
        moduleData <- miRNAkMEs %>% filter(.data$moduleLabel == input$miRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$miRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:20) %>% #select the top N rows data
            select(1:3, grep(pattern = paste0("\\b", input$miRNA.modules, "\\b") , x =  colnames(miRNAkMEs))) %>% ##grep the index of selected module
            left_join(., miRNA.exp, by=c("ID"="miRNA")) #add the miRNA expression data
            colnames(moduleData)[4] <- paste0(colnames(moduleData)[4], "_kMEs")
    }
    if (input$miRNA.KMEs == 30) {
        moduleData <- miRNAkMEs %>% filter(.data$moduleLabel == input$miRNA.modules) %>% #extract the module data
            arrange(desc(eval(parse(text = input$miRNA.modules)))) %>% ##order the selected module according to kMEs
            slice(1:30) %>% #select the top N rows data
            select(1:3, grep(pattern = paste0("\\b", input$miRNA.modules, "\\b") , x =  colnames(miRNAkMEs))) %>% ##grep the index of selected module
            left_join(., miRNA.exp, by=c("ID"="miRNA")) #add the miRNA expression data
            colnames(moduleData)[4] <- paste0(colnames(moduleData)[4], "_kMEs")
    }
    
    miRhtplot <- as.matrix(moduleData[,5:24], rownames = moduleData$ID)
    miRhtplotZscore <- t(scale(t(miRhtplot), center = TRUE, scale = TRUE))

    anno_columnmiR <- data.frame(Group = factor(rep(c("Skin", "Wound1", "Wound7", "VU"), each = 5)), 
        row.names = colnames(miRhtplot))

    miRWGCNAheatmap <- DEcomplexheatmap(
      matZval = miRhtplotZscore, 
      columnName = TRUE, 
      rowName = TRUE, 
      columnClu = FALSE, 
      rowClu = TRUE, 
      annoColumn = anno_columnmiR#, 
      #annoCol = anno_color
    )
    
    output$miRmoduleplot <- renderPlot({
      miRWGCNAheatmap
    })

    #Set plot height using gtable package functions
    gtable_h2 <- convertHeight(gtable_height(miRWGCNAheatmap$gtable),"inches",valueOnly = TRUE)
    legend_row2 <- length(unique(anno_columnmiR$Group))
    h2 <- max(gtable_h2, legend_row2/2 + 1.5)

    plot_height2 <- h2 * 70 + 150
    plotDE_pdf_height2 <<- h2 + 0.5 #the output heatmap

    output$miRWGCNATable <- DT::renderDataTable({
      outdata <- copy(moduleData) #data you want to show
      outdata$ID <- createLink(outdata$ID, miRBase_LINK, outdata$ID)
      outdata <- outdata %>% mutate_if(is.numeric, signif, 8) #set the output decimals
      colnames(outdata)[5:24] <- paste0(colnames(outdata)[5:24], "_TPM")
      return(outdata)},
      escape = FALSE,server = FALSE,extensions = c("Buttons"), 
      options = list(dom = 'Bfrtip',buttons = c('copy', 'csv')), #'excel': some of the padj values in output are strings. So we suggest to choose the csv format.
      rownames= FALSE)

    output$miRNAWGCNA.UI <- renderUI({
      tabBox(
        title = "", width = 12,
        tabPanel("Module Heatmap", plotOutput("miRmoduleplot", height = plot_height2)),
        tabPanel("Data Output", DT::dataTableOutput("miRWGCNATable"), style = "overflow-x: scroll;") #height:800px; overflow-y: scroll;
      )
    })

    output$downloadmiRwgcnaPlot <- downloadHandler(
      filename <- function() {
          paste0("Heatmap_miRNAmodules", ".pdf")
      }, 
      content <- function(file) {
        pdf(file, width = 10, height = plotDE_pdf_height2)
          grid::grid.newpage()
          grid::grid.draw(miRWGCNAheatmap$gtable)
        dev.off()
      }
    )

  })#withProgress
})#observeEvent
