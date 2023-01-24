tabItem(tabName = "RNAexpressions",
  fluidRow(
    box(
      width = 6, title="miRNA Expression", solidHeader = TRUE, collapsible = TRUE,
      column(4, selectInput("miRNA.TPM", "Type:", c("log2(TPM+1)", "TPM"))),
      column(4, selectizeInput("miRNA.names", "miRNA Name", choices = NULL, selected=NULL, options = list(openOnFocus = FALSE))),
      column(2, style = "margin-top: 25px;", actionButton("miRNAExpPlot1", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
    ),
    box(
      width = 6, title="mRNA Expression", solidHeader = TRUE, collapsible = TRUE,
      column(4, selectInput("mRNA.FPKM", "Type:", c("log2(FPKM+1)", "FPKM"))),
      column(4, selectizeInput("mRNA.geneName", "mRNA Name", choices = NULL, selected=NULL, options = list(openOnFocus = FALSE))),
      column(2, style = "margin-top: 25px;", actionButton("mRNAExpPlot1", "VISUALIZATION", class = "button-3d button-block button-pill button-primary"))
    )
  ), 
  #hr(),
  fluidRow(
    box(
      title="miRNA Plot", width = 6,
      uiOutput("miRNAExp.UI")
    ),
    box(
      title="mRNA Plot", width = 6,
      uiOutput("mRNAExp.UI")
    )
  ),
  fixedPanel(
    downloadButton("downloadmiRNAPlot", label = "Download miRPlot", class = "button-primary", icon = icon("download")),
    downloadButton("downloadmRNAPlot", label = "Download mRNAPlot", class = "button-primary", icon = icon("download")),
    right = 5,
    bottom = 5
  )
)

