
updateSelectizeInput(session, "miRNA.names", choices = unique(sort(miRNA.exp$miRNA)), server = TRUE)
updateSelectizeInput(session, "mRNA.geneName", choices = unique(sort(mRNA.exp$GeneSymbol)), server = TRUE)

observeEvent(input$miRNAExpPlot1,{
  withProgress(message = "Processing , please wait",{ 

    upmiRNAname <- toupper(input$miRNA.names) #transform all input data into uppercase
    miRNAexp.cp <- copy(miRNA.exp) %>% as.data.frame()
    rownames(miRNAexp.cp) <- miRNAexp.cp$miRNA #add the rownames to the data
    miRNAexp.cp$miRNA <- toupper(miRNAexp.cp$miRNA) #transform all data into uppercase
    miRNAexp.cp <- miRNAexp.cp %>% distinct(miRNA, .keep_all = TRUE) #Check if it has duplicates

    vioplotdata <- miRNAexp.cp[miRNAexp.cp$miRNA %in% upmiRNAname, -1] #this data is for the rownames
    viodata <- vioplotdata %>% 
      pivot_longer(everything(), names_to = "SampleID" , values_to = "TPM") %>% #prepare the data for plot
      mutate(Group=factor(rep(c("Skin", "Wound1", "Wound7", "VU"), each=5), levels = c("Skin", "Wound1", "Wound7", "VU")))

    if(nrow(vioplotdata)==0){
      return()
    }  

    my_comparisons <- list(c("VU", "Skin"), c("VU", "Wound1"), c("VU", "Wound7"),
      c("Wound1", "Skin"), c("Wound7", "Skin"), c("Wound7", "Wound1"))

    if (input$miRNA.TPM == "TPM") {
      vioplot <- vioplotFun(
        datadf = viodata, 
        aesX = "Group",
        aesY = "TPM",
        compares = my_comparisons)
    }
    if (input$miRNA.TPM == "log2(TPM+1)") {
      viodata$log2TPM <- log2(viodata$TPM + 1)
      vioplot <- vioplotFun(
        datadf = viodata, 
        aesX = "Group",
        aesY = "log2TPM",
        compares = my_comparisons)
      vioplot <- ggpar(vioplot, ylab="log2(TPM+1)")
    }    

    output$vioplotout <- renderPlot({
      vioplot
    })

    output$miRNAExp.UI <- renderUI({
      tabBox(
        title = "", width = 12,
        tabPanel(rownames(vioplotdata), plotOutput("vioplotout"))
      )
    })

    output$downloadmiRNAPlot <- downloadHandler(
      filename <- function() {
          paste0("miRNAexpression_violin", ".pdf")
      }, 
      content <- function(file) {
        pdf(file, width = 8, height = 8)
          print(vioplot)
        dev.off()
      }
    )

  })#withProgress
})#observeEvent


observeEvent(input$mRNAExpPlot1,{
  withProgress(message = "Processing , please wait",{ 

    upmRNAname <- toupper(input$mRNA.geneName) #transform all input data into uppercase
    mRNAexp.cp <- copy(mRNA.exp) %>% as.data.frame() %>% distinct(GeneSymbol, .keep_all = TRUE) #Check if it has duplicates
    rownames(mRNAexp.cp) <- mRNAexp.cp$GeneSymbol #add the rownames to the data
    mRNAexp.cp$GeneSymbol <- toupper(mRNAexp.cp$GeneSymbol) #transform all data into uppercase

    vioplotmRNA <- mRNAexp.cp[mRNAexp.cp$GeneSymbol %in% upmRNAname, -c(1,2,3)] #this data is for the rownames
    viomRNA <- vioplotmRNA %>% 
      pivot_longer(everything(), names_to = "SampleID" , values_to = "FPKM") %>% #prepare the data for plot
      mutate(Group=factor(rep(c("Skin", "Wound1", "Wound7", "VU"), each=5), levels = c("Skin", "Wound1", "Wound7", "VU")))

    if(nrow(vioplotmRNA)==0){
      return()
    }  

    my_comparisons <- list(c("VU", "Skin"), c("VU", "Wound1"), c("VU", "Wound7"),
      c("Wound1", "Skin"), c("Wound7", "Skin"), c("Wound7", "Wound1"))
    
    if (input$mRNA.FPKM == "FPKM") {
      vioplot2 <- vioplotFun(
        datadf = viomRNA, 
        aesX = "Group",
        aesY = "FPKM",
        compares = my_comparisons)
    }
    if (input$mRNA.FPKM == "log2(FPKM+1)") {
      viomRNA$log2FPKM <- log2(viomRNA$FPKM + 1)
      vioplot2 <- vioplotFun(
        datadf = viomRNA, 
        aesX = "Group",
        aesY = "log2FPKM",
        compares = my_comparisons)
      vioplot2 <- ggpar(vioplot2, ylab="log2(FPKM+1)")
    }

    output$vioplotout2 <- renderPlot({
      vioplot2
    })

    output$mRNAExp.UI <- renderUI({
      tabBox(
        title = "", width = 12,
        tabPanel(rownames(vioplotmRNA), plotOutput("vioplotout2"))
      )
    })

    output$downloadmRNAPlot <- downloadHandler(
      filename <- function() {
          paste0("mRNAexpression_violin", ".pdf")
      }, 
      content <- function(file) {
        pdf(file, width = 8, height = 8)
          print(vioplot2)
        dev.off()
      }
    )

  })#withProgress
})#observeEvent
