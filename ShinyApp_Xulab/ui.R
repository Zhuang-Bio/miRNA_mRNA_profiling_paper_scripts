library(shiny)
library(shinydashboard)
library(shinyjs)
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(V8) #required by shinyjs
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(zip)
library(foreach)
library(DT)
library(dplyr)
library(gtable)
library(grid)
library(ggpubr)



ui <- tagList(
    dashboardPage(
        dashboardHeader(
            title = "Integrative small and long RNA-omics analysis of human healing and non-healing wounds discovers cooperating microRNAs as therapeutic targets",
            titleWidth = "98%"),
        dashboardSidebar(
            sidebarMenu(
                id = "tabs", #if id="tabs", then input$tabs will be the tabName of the currently-selected tab. If you want to be able to bookmark and restore the selected tab, an id is required.
                menuItem("Quick start tutorial", tabName = "startTutorial"),
                menuItem("miRNA-mRNA expression", tabName = "RNAexpressions", icon = icon("dice-one","fa-1.5x")), # 1.5x normal size                
                menuItem("DE miRNA analysis", tabName = "TabDEmiRNA", icon = icon("dice-two","fa-1.5x")), 
                menuItem("miRNA WGCNA", tabName = "miRNAWGCNA", icon = icon("dice-three","fa-1.5x")),
                menuItem("DE mRNA analysis", tabName = "TabDEmRNA", icon = icon("dice-four","fa-1.5x")),
                menuItem("mRNA WGCNA", tabName = "mRNAWGCNA", icon = icon("dice-five","fa-1.5x"))
            )#,
            #tags$style(HTML(".main-sidebar { font-size: 15px; }")) #set the sidebar names font sizes
        ),
        dashboardBody(
            shinyjs::useShinyjs(),
            extendShinyjs(script = "www/custom.js", functions=c()),
            tags$head(
                tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"), 
                tags$link(rel = "stylesheet", type = "text/css", href = "buttons.css")
            ),
            #type: Specifies the media type of the linked document
            #rel:  Required. Specifies the relationship between the current document and the linked document. rel = "stylesheet" specifies a persistent or preferred style while rel="alternate" defines an alternate style. 
            #href: Specifies the location of the linked document (URL).
            tabItems(
                source("ui-tab-tutorial.R", local = TRUE)$value,
                source("ui-tab-miRNA-mRNAexpression.R", local = TRUE)$value,
                source("ui-tab-DEmiRNA.R", local = TRUE)$value,
                source("ui-tab-miRNA-WGCNA.R", local = TRUE)$value,
                source("ui-tab-DEmRNA.R", local = TRUE)$value,
                source("ui-tab-mRNA-WGCNA.R", local = TRUE)$value
            )
        )
    ),
    tags$footer(
        HTML(
          '
           <p align="center" style="margin:10px;font-size:15px">© 2021-2023 Zhuang Liu<br>
                <a href="https://www.xulandenlab.com/data" target="_blank">Xu Land&#233;n Laboratory | Karolinska Institutet</a> 
           </p>
          '
          )
    )
)
