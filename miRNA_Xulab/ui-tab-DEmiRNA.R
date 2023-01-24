tabItem(tabName = "TabDEmiRNA",
   fluidRow(
    box(
      width = 12, solidHeader = TRUE, collapsible = TRUE,
      fluidRow(
        width = 12,
        column(6, 
          selectInput("DEmiRNA.geneData", "Comparisons of Differential Expression Analysis",
            list("VU vs. Skin" = "VUvsSkin",
                 "VU vs. Wound1" = "VUvsWound1",
                 "VU vs. Wound7" = "VUvsWound7",
                 "Wound1 vs. Skin" = "Wound1vsSkin",
                 "Wound7 vs. Skin" = "Wound7vsSkin",
                 "Wound7 vs. Wound1" = "Wound7vsWound1"))),
        column(2, selectInput("DEmiRNA.clusterGene","Row Clustering", c("Yes","No"))),
        column(2, selectInput("DEmiRNA.clusterSample","Column Clustering", c("Yes","No")))
      ),      
     fluidRow(
        width = 12,
        column(2,numericInput("DEmiRNA.pvalue", "P-adjust value threshold", value = 0.05, min = 0, max = 1, step = 0.01)), 
        column(2,numericInput("DEmiRNA.foldChange", "Log2(FC) threshold", value = 1, min = -Inf, max = Inf, step=0.1)),
        column(2,selectInput("DEmiRNA.FCtype", "Fold Change Type", list("abs(FC) >= " = "both", ">=" = "greatthan", "<=" =  "lessthan"))),
        column(2,selectInput("DEmiRNA.rownames", "Show Gene Names", c("Right","None"))),
        column(2,selectInput("DEmiRNA.colnames", "Show Sample Names", c("Bottom","None"))),
        column(2, style = "margin-top: 25px;", 
          actionButton("DemiRNA.Plot", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
      ),
    ) #box
  ), #fluidRow
  #hr(),
  fluidRow(
    uiOutput("DEmiRNA.UI")
  ),
  fixedPanel(
    downloadButton("downloadDEHeatmap", label = "Download Heatmap", class = "button-primary", icon = icon("download")),
    right = 15,
    bottom = 15
  )
)

