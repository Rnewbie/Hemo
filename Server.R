library(shiny)
library(protr)
library(RWeka)
library(seqinr)
library(caret)
library(randomForest)


extractCTDD = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
                normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
                polarity        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
                polarizability  = c('G', 'A', 'S', 'D', 'T'),
                charge          = c('K', 'R'),
                secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
                solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
  
  group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
                normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
                polarity        = c('P', 'A', 'T', 'G', 'S'),
                polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
                charge          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 
                                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
                secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
                solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
  
  group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
                normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
                polarity        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
                polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
                charge          = c('D', 'E'),
                secondarystruct = c('G', 'N', 'P', 'S', 'D'),
                solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))
  
  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)
  
  G = vector('list', 7)
  for (i in 1:7) G[[i]] = rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:7) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- 'G1')
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- 'G2')
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- 'G3')
  }
  
  # Compute Distribution
  
  D = vector('list', 7)
  for (i in 1:7) D[[i]] = matrix(ncol = 5, nrow = 3)
  
  for (i in 1:7) {
    inds = which(G[[i]] == 'G1')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][1, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G2')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][2, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G3')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][3, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
  }
  
  D = do.call(rbind, D)
  D = as.vector(t(D))
  
  names(D) = paste(rep(paste('prop', 1:7, sep = ''), each = 15),
                   rep(rep(c('.G1', '.G2', '.G3'), each = 5), times = 7),
                   rep(paste('.residue', c('0', '25', '50', '75', '100'), 
                             sep = ''), times = 21), sep = '')
  
  return(D)
  
}

x <- readFASTA("alpha.fasta")
data <- read.csv("data.csv", header = TRUE)
Affinity <- data$Affinity
bossrequestedDes <- function(x) {
  c(extractCTDC(x), ### composition
    extractCTDT(x), ### transition
    extractCTDD(x)) ### distribution
}
alpha <- t(sapply(x, bossrequestedDes))
rownames(alpha) <- NULL
rownames(Affinity) <- NULL
alpha <- data.frame(alpha)
data <- cbind(Affinity, alpha)
fit <- J48(Affinity~., data = data)



shinyServer(function(input, output, session) {
  
  
  observe({
    FASTADATA <- ''
    fastaexample <- '>wild_type(Normal)
VLDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHPDASTVSKYR
>Hb_Lyon_Bron(Decrease)
ALDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHPDASTVSKYR
>Hb_Diamant(Increase)
VLDAWKYARLKTPHFDSHAVKGHGKADANDDMPLSSDLHAHKRDPVNAHLDASTVSKYR'
    if(input$addlink>0) {
      isolate({
        FASTADATA <- fastaexample
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
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
        x <- readFASTA("text.fasta")
        x <- x[(sapply(x, protcheck))]
        DPC <- t(sapply(x, bossrequestedDes))
        test <- data.frame(DPC)
        Prediction <- predict(fit, test)
        Prediction <- as.data.frame(Prediction)
        Protein <- cbind(Name = rownames(DPC, Prediction))
        results <- cbind(Protein, Prediction)
        results <- data.frame(results, row.names=NULL)
        print(results)
      } 
      else {     
        x <- readFASTA(inFile$datapath)
        x <- x[(sapply(x, protcheck))]
        DPC <- t(sapply(x, bossrequestedDes))
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
  
  
  
  })






