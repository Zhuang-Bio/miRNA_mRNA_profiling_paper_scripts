tabItem(tabName = "mRNAWGCNA",
  fluidRow(      
    box(
      width = 6, height = 1000, title="mRNA WGCNA trait-module relationships", solidHeader = TRUE, collapsible = TRUE,
      column(12, uiOutput("mRNA_WGCNAimg"))
    ),
    box(
      width = 6, title="mRNA WGCNA", solidHeader = TRUE, collapsible = TRUE,
      column(4, selectizeInput("mRNA.modules", "mRNA modules", choices = NULL, selected=NULL, options = list(openOnFocus = FALSE))),
      column(4, selectInput("mRNA.KMEs", "Top(N) kMEs:", c(20, 10, 30, 40))),
      column(2, style = "margin-top: 25px;", actionButton("mRNAWGCNAPlot1", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
    ),
    box(
      title="", width = 6,
      uiOutput("mRNAWGCNA.UI")
    ),
  ), 
  fixedPanel(
    downloadButton("downloadmRNAwgcnaPlot", label = "Download Heatmap", class = "button-primary", icon = icon("download")),
    right = 5,
    bottom = 5
  )
)

