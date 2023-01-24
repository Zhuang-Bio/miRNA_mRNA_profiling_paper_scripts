tabItem(tabName = "TabDEmRNA",
   fluidRow(
    box(
      width = 12, solidHeader = TRUE, collapsible = TRUE,
      fluidRow(
        width = 12,
        column(6, 
          selectInput("DEmRNA.geneData", "Comparisons of Differential Expression Analysis",
            list("VU vs. Skin" = "VUvsSkin",
                 "VU vs. Wound1" = "VUvsWound1",
                 "VU vs. Wound7" = "VUvsWound7",
                 "Wound1 vs. Skin" = "Wound1vsSkin",
                 "Wound7 vs. Skin" = "Wound7vsSkin",
                 "Wound7 vs. Wound1" = "Wound7vsWound1"))),
        column(2, selectInput("DEmRNA.clusterGene","Row Clustering", c("Yes","No"))),
        column(2, selectInput("DEmRNA.clusterSample","Column Clustering", c("Yes","No")))
      ),      
     fluidRow(
        width = 12,
        column(2,numericInput("DEmRNA.pvalue", "P-adjust value threshold", value = 0.05, min = 0, max = 1, step = 0.01)), 
        column(2,numericInput("DEmRNA.foldChange", "Log2(FC) threshold", value = 0.58, min = -Inf, max = Inf, step=0.1)),
        column(2,selectInput("DEmRNA.FCtype", "Fold Change Type", list("abs(FC) >= " = "both", ">=" = "greatthan", "<=" =  "lessthan"))),
        column(2,selectInput("DEmRNA.rownames", "Show Gene Names", c("None","Right"))),
        column(2,selectInput("DEmRNA.colnames", "Show Sample Names", c("Bottom","None"))),
        column(2, style = "margin-top: 25px;", 
          actionButton("DEmRNA.Plot", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
      ),
    ) #box
  ), #fluidRow
  #hr(),
  fluidRow(
    uiOutput("DEmRNA.UI")
  ),
  fixedPanel(
    downloadButton("downloadDEmRNAHeatmap", label = "Download Heatmap", class = "button-primary", icon = icon("download")),
    right = 15,
    bottom = 15
  )
)

