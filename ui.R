library(shiny)
library(DT)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("SB SCC mciroarray data - 21 samples"),

  sidebarPanel(
    selectInput("colors", "Color scheme:",
                list("Green/Red" = "greenred",
                     "Blue/Red" = "bluered",
                     "Blue/Yellow" = "blueyellow",
                     "Matlab" = "matlab")),
    
    selectInput("list", "Gene list:",
                list("Proliferation" = "gps",
                     "EMT" = "emt",
                     "Immune Response" = "th1",
                     "Gene correlation" = "gcor",
                     "Top n DE genes" = "topg",
                     "User Defined" = "user")),
    
    conditionalPanel(
      condition = "input.list == 'user'",
      textInput("genes","Gene list", 
                value = "Zmiz1 Ereg Ptgs2")
      ),

    conditionalPanel(
      condition = "input.list == 'topg'",
      numericInput("ngenes","Number of genes", 
                value = 50)
    ),
    
    conditionalPanel(
      condition = "input.list == 'gcor'",
      textInput("corGene","Gene", 
                value = "Zmiz1"),
      numericInput("corThres","Correlation threshold", 
                value = 0.7, min=0, max=1, step=0.05)
      ),
    
    selectInput("order", "Sample order:",
                list("Metagene" = "mg",
                     "Clustering" = "hc",
                     "Zmiz1 status" = "mt")
                ),
      
    conditionalPanel(
      condition = "input.order == 'hc'",
      radioButtons("radioDist","Distance metric",
                   c("Euclidean","Correlation"),
                   selected="Euclidean")
      )
      
#    checkboxInput("slider", "Slider?", value = TRUE),
    
#    conditionalPanel(
#        condition = "input.slider",
#        sliderInput("survSlider","Survival cut",0,1,0.5)
#    ),
    
#    selectInput("samp", "Samples:",
#                list("All samples" = "all",
#                     "BRAF mutants"="braf",
#                     "NRAS mutants"="nras",
#                     "Primary tumours" = "pri",
#                     "All metastases" = "met",
#                     "Distant metastases" = "dmt",
#                     "Regional Lymph Nodes" = "rln",
#                     "Regional tissue (metastasis)"="rtm")
#    )
    
         
    ,width=2
    ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(h4("Expression Heatmap"), plotOutput("expMap", height = 700, width = 1500)),
#      tabPanel(h4("Survival Plot"), h4(textOutput("cph")), plotOutput("survPlot", height = 600, width = 1100)),
      tabPanel(h4("Gene info (DiffExp: Zmiz1 status)"), dataTableOutput("giDisplay", width=600)),
      tabPanel(h4("Gene list"), tableOutput("glDisplay"))
    )
    ,width=10
  )
))
