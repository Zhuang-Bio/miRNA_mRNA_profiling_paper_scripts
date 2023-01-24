tabItem(tabName = "miRNAWGCNA",
  fluidRow(      
    box(
      width = 6, height = 1000, title="miRNA WGCNA trait-module relationships", solidHeader = TRUE, collapsible = TRUE,
      column(12, uiOutput("miRNA_WGCNAimg"))
    ),
    box(
      width = 6, title="miRNA WGCNA", solidHeader = TRUE, collapsible = TRUE,
      column(4, selectizeInput("miRNA.modules", "miRNA modules", choices = NULL, selected=NULL, options = list(openOnFocus = FALSE))),
      column(4, selectInput("miRNA.KMEs", "Top(N) kMEs:", c(20, 10, 30))),
      column(2, style = "margin-top: 25px;", actionButton("miRNAWGCNAPlot1", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
    ),
    box(
      title="", width = 6,
      uiOutput("miRNAWGCNA.UI")
    ),
  ), 
  fixedPanel(
    downloadButton("downloadmiRwgcnaPlot", label = "Download Heatmap", class = "button-primary", icon = icon("download")),
    right = 5,
    bottom = 5
  )
)

