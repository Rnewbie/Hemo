#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("function.R")


library(shiny)
library(shinythemes)
library(protr)
library(markdown)
library(shiny)
library(protr)
library(RWeka)
library(seqinr)
library(caret)
library(randomForest)

df <- read.csv("alpha_data_final.csv")
fit <- J48(Affinity~., data = df)

# Define UI for application that draws a histogram
ui <- fluidPage(title="HOBBIT: A webserver for prediction the oxygen affinity of human hemoglobin", theme=shinytheme("cerulean"),
                navbarPage(strong("FPOP"),
                           tabPanel("Submit Job", titlePanel("HOBBIT: A webserver for prediction the oxygen affinity of human hemoglobin"),
                                    sidebarLayout(
                                      wellPanel(
                                        radioButtons(inputId = "Selection",label = "Select the Sequence", choices = c("Alpha", "Beta")),
                                        tags$label("Enter your input sequence(s) in FASTA format",style="float: none; width: 100%;"),
                                        actionLink("addlink", "Insert example data"),
                                        tags$textarea(id="Sequence", rows=5, cols=100, style="float: none; width:100%;", ""),
                                        #actionLink("addlink", "Insert example data"),
                                        #tags$label("or",style="float: none; width: 100%;"),
                                        fileInput('file1', 'or upload file',accept=c('text/FASTA','FASTA','.fasta','.txt')),
                                        # tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
                                        actionButton("submitbutton", "Submit", class = "btn btn-primary"),
                                        actionButton("clearbutton", "Clear", class = "btn btn-danger")
                                      ), #wellPanel
                                      
                                      mainPanel(
                                        verbatimTextOutput('contents'),
                                        downloadButton('downloadData', 'Download CSV')
                                        
                                      )  
                                    ) #sidebarLayout
                           )#, #tabPanel Submit Job
                           
                          # tabPanel("About", titlePanel("Fluorescent protein oligomerization"), div(includeMarkdown("about.md"), align="justify")),
                          # tabPanel("Citing Us", titlePanel("Citing Us"), includeMarkdown("citingus.md")),
                          # tabPanel("Contact", titlePanel("Contact"), includeMarkdown("contact.md"))	
                           
                ) #navbarPage
) #fluidPage	


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  observe({
    FASTADATA <- ''
#    fastaexample <- '>wild_type(Normal)
#    VLDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHPDASTVSKYR
#    >Hb_Lyon_Bron(Decrease)
#    ALDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHPDASTVSKYR
#    >Hb_Diamant(Increase)
#    VLDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHLDASTVSKYR'
    fastaexample_alpha <- '>HBA_Human Hemoglobin subunit alpha
    MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'
    fastaexample_beta <- '>HBA_Human Hemoglobin subunit beta
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
    '
    if(input$addlink>0) {
      isolate({
        if (input$Selection == "Alpha") {
          FASTADATA <- fastaexample_alpha
          updateTextInput(session, inputId = "Sequence", value = FASTADATA)
        } else if (input$Selection == "Beta") {
        FASTADATA <- fastaexample_beta
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
        }
      })
    }
  })
  
  observe({
    emptyDATA <- ""
    if(input$clearbutton>0) {
      isolate({
        updateTextInput(session, inputId = "Sequence", value = emptyDATA) 
        is.null(datasetInput())
        
      })
    }
  })
  
  
  datasetInput <- reactive({
    
    inFile <- input$file1 
    inTextbox <- input$Sequence
    
    
    if (is.null(inTextbox)) {
      return("Please insert/upload sequence in FASTA format")
    } else {
      if (is.null(inFile)) {
        x <- inTextbox
        write.fasta(sequence = x, names = names(x),
                    nbchar = 80, file.out = "text.fasta")
       # x <- readFASTA("text.fasta")
      #  x <- x[(sapply(x, protcheck))]
       # DPC <- t(sapply(x, bossrequestedDes))
        DPC <- alpha_z3("text.fasta")
        test <- data.frame(DPC)
        Prediction <- predict(fit, test)
        Prediction <- as.data.frame(Prediction)
        Protein <- cbind(Name = rownames(DPC, Prediction))
        results <- cbind(Protein, Prediction)
        results <- data.frame(results, row.names=NULL)
        print(results)
      } 
      else {     
        #x <- readFASTA(inFile$datapath)
        #x <- x[(sapply(x, protcheck))]
        #DPC <- t(sapply(x, bossrequestedDes))
        DPC <- alpha_z3(inFile$datapath)
        test <- data.frame(DPC)
        Prediction <- predict(fit, test)
        Prediction <- as.data.frame(Prediction)
        Protein <- cbind(Protein = rownames(DPC, Prediction))
        results <- cbind(Protein, Prediction)
        results <- data.frame(results, row.names=NULL)
        print(results)
        
      }
    }
    
    
  })
  
  
  
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } else {
      if (input$clearbutton>0) {
        isolate(is.null(datasetInput()))
      } else {
        return("Please insert/upload sequence in FASTA format")
      }
    }
  })
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('Predicted_Results', '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
    })

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

