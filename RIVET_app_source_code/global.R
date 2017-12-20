# ui-------------------------------------------------------------------

# allows user to choose samples
sampleUploadUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,
      #makes checkbox of control choices
      uiOutput(ns("ctlSamples"))
      ),
      column(6,
      #makes checkbox of experiment choices
      uiOutput(ns("exptSamples"))
      )
    )
  )
}

# download matrix
DownloadButtonUI <- function(id){
  ns <- NS(id)
  downloadButton(ns('download'), 'Download')
}

# set fold-change and statistic threshold
sliderInputFcUI <- function(id){
  ns <- NS(id)
  tagList(
    sliderInput(ns('logFC'), 'Log Fold-Change', 
                min=0, max=5, value=min(0), 
                step=0.2, round=FALSE),
    radioButtons(ns('prob'), 'P-val or Adj. p-val',
                 choices = c(
                   'p-value' = 'P.Value',
                   'adj p-val' = 'adj.P.Val'
                 ),'P.Value'),
    sliderInput(ns('pVal'), 'Statistic Threshold', 
                min=0, max=.99, value=min(0.05), 
                step=0.01, round=FALSE)
  )
}

# Download non-ggplot graphs
DownloadPicUI <- function(id){
  ns <- NS(id)
  downloadButton(ns('dndPlot'), 'Download Pic')
}

# Download ggplot graphs
DownloadPicGGUI <- function(id){
  ns <- NS(id)
  downloadButton(ns('dndPlotGG'), 'Download Pic')
}

# Download heatmap graphs
heatmapUI <-function(id){
  ns <- NS(id)
  downloadButton(ns('dwnHeatmap'), 'Download Pic')
}

# server------------------------------------------------------------

# create normalized filtered object, used as input for voom and edgeR
dge <- function(input, output, session, expMatrix){
normalize <- reactive({
  req(expMatrix())
  dge <- DGEList(counts=expMatrix())
  filtered <- rowSums(cpm(dge) > 1) >= 2
  dge <- dge[filtered,]
  dge <- calcNormFactors(dge)
  return(dge)
})
return(normalize)
}

# creates input boxes for control and experiment samples
sampleUpload <- function(input,output,session, expMatrix){
  
  # user clicks control samples
  output$ctlSamples <- renderUI({
    req(expMatrix())  
    ns <- session$ns
    selectInput(ns("ctlallSamples"), 'Control', colnames(expMatrix()), multiple=TRUE, selectize=TRUE)
  })

  # user clicks experiment samples
  output$exptSamples <- renderUI({
    req(expMatrix())
    ns <- session$ns
    selectInput(ns("exptallSamples"), 'Experiment', colnames(expMatrix()), multiple=TRUE, selectize=TRUE)
  })
  
  # makes a list of control samples, experiment samples that can be used to reorder
  # original matrix
  listcond <- reactive({
    req(input$ctlallSamples)
    req(input$exptallSamples)
    cond1_t = input$ctlallSamples
    cond2_t = input$exptallSamples
    transcription_conditions = list(cond1_t, cond2_t)
    return(transcription_conditions)
    })
  return(listcond)
  
}

# limma statistical test
limma <- function(input,output,session, expMatrix, contrast_matrix, design_matrix,compare){
  calc_limma <- reactive({
    req(expMatrix())
    req(contrast_matrix())
    req(design_matrix())
    req(compare())
    exprs = expMatrix()
    contrast.matrix = contrast_matrix()
    designMat = design_matrix()
    comparison = compare()
    logset <- new("ExpressionSet", exprs = exprs)
    fit <- lmFit(logset, designMat)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    all = topTable(fit2, coef=comparison, adjust="BH", number=Inf)
    all_tofile = all[,c('logFC', 'P.Value', 'adj.P.Val')]
    colnames(all_tofile) <- lapply(colnames(all_tofile), paste, comparison, sep='_')
    all_tofile = all_tofile[order(rownames(all_tofile)),]
    return(all_tofile)
  })
  return(calc_limma)
}

# edgeR statistical test
edgeR <- function(input,output,session, expMatrix, contrast_matrix, design_matrix, compare){
  calcEdgeR <- reactive({
    req(expMatrix())
    req(contrast_matrix())
    req(design_matrix())
    req(compare())
    expMatrix = expMatrix()
    design_matrix = design_matrix()
    contrast_matrix = contrast_matrix()
    comparison = compare()
    expMatrix <- estimateDisp(expMatrix, design_matrix)
    fit <- glmQLFit(expMatrix, design_matrix)
    qlf_test <- glmQLFTest(fit, contrast=contrast_matrix[,comparison])
    stats2 <- topTags(qlf_test, n=Inf)
    stats_tofile <- stats2[,c('logFC', 'PValue', 'FDR')]
    stats_tofile <- stats_tofile[order(rownames(stats_tofile)),]
    colnames(stats_tofile$table) <- lapply(c('logFC', 'P.Value', 'adj.P.Val'), paste, comparison, sep='_')
    return(stats_tofile)
  })
  return(calcEdgeR)
}

# volcano plot for transcription, translation, and translational efficiency tab
volcanoPlot <- function(input, output, session, whole_dataset, logFC, pval, option){
  
  volcanoInput <- reactive({
    
    req(whole_dataset())
    req(logFC())
    req(pval())
    req(option())
    
    whole_data <- whole_dataset()
    pval <- pval()
    logFC <- logFC()
    option <- option()
    
    # isolate logFC, pval, adjusted pval names
    splitname <- function(name){
      unlist(strsplit(name, '_'))[1]
    }
    colnames(whole_data) <- lapply(colnames(whole_data), splitname)
    
    # create threshold column based on user selected fc and statistical thresholds
    whole_data$threshold = factor(whole_data[,option] <= pval & abs(whole_data[,'logFC'])>=logFC)
    
    # Set log threshold for graph
    minFC = abs(min(whole_data[,'logFC']))
    maxFC = abs(max(whole_data[,'logFC']))
    if (minFC>maxFC) {
      lim = minFC
    } else {
      lim = maxFC}
    
    # Construct the plot object
    ggplot(data=whole_data, aes(x=whole_data[,'logFC'], y=-log10(whole_data[,option]), colour=threshold)) +
      geom_point(alpha=0.4, size=1.75) +
      xlim(c(-lim, lim)) + ylim(c(0, max(-log10(whole_data[,option])))) +
      xlab("log2 fold change") + ylab("-log10 p-value") + theme(legend.position='none')
    
    })
  
  return(volcanoInput)
}

# volcano plot for download, some sort of bug with downloading input directly
volcanoDwnInput <- function(input, output, session, whole_dataset, logFC, pval, option){
  
  # this is a function instead of a reactive object
  volcanoInput <- function(){
    
    req(whole_dataset())
    req(logFC())
    req(pval())
    req(option())
    
    whole_data <- whole_dataset()
    pval <- pval()
    logFC <- logFC()
    option <- option()
    splitname <- function(name){
      unlist(strsplit(name, '_'))[1]
    }
    colnames(whole_data) <- lapply(colnames(whole_data), splitname)
    return(whole_data[,option])
    whole_data$threshold = factor(whole_data[,option] <= pval & abs(whole_data[,'logFC'])>=logFC)
    
    ##Set log threshold for graph
    minFC = abs(min(whole_data[,'logFC']))
    maxFC = abs(max(whole_data[,'logFC']))
    if (minFC>maxFC) {
      lim = minFC
    } else {
      lim = maxFC}
    
    ##Construct the plot object
    ggplot(data=whole_data, aes(x=whole_data[,'logFC'], y=-log10(whole_data[,option]), colour=threshold)) +
      geom_point(alpha=0.4, size=1.75) +
      xlim(c(-lim, lim)) + ylim(c(0, max(-log10(whole_data[,option])))) +
      xlab("log2 fold change") + ylab("-log10 p-value") + theme(legend.position='none')
    
  }
  return(volcanoInput)
}

# heatmap for translation and translational efficiency
heatmap <- function(input, output, session, translation_foldchange){
  
  heatmapInput <- reactive({
    
    req(translation_foldchange())
    whole_data <- translation_foldchange()
    
    ##Heatmap object
    myColors = brewer.pal(n=11, name="RdBu")
    myColors = colorRampPalette(myColors)(50)
    heatmap.2(whole_data, dendrogram = 'row', Colv=FALSE, tracecol='black', labRow=FALSE, col=myColors, cexCol=2)
    
  })
  
  return(heatmapInput)
}

# Heatmap download plot
heatmapDwnInput <- function(input, output, session, translation_foldchange){
  
  heatmapInput <- function(){
    
    req(translation_foldchange())
    
    whole_data <- translation_foldchange()
    
    ##Heatmap object
    myColors = brewer.pal(n=11, name="RdBu")
    myColors = colorRampPalette(myColors)(50)
    heatmap.2(whole_data, dendrogram = 'row', Colv=FALSE, tracecol='black', labRow=FALSE, col=myColors, cexCol=2)
    
  }
  
  return(heatmapInput)
}

# slider for determining fc, pvalue, adj pval thresholds
sliderInputFC <- function(input,output, session){
  logFC <- reactive(input$logFC)
  pval <- reactive(input$pVal)
  input_choice <- reactive(input$prob)
  params <- reactive(list(logFC(), pval(), input_choice()))
  return(params)
}

# threshold stat matrix based on user input
top_genes <- function(input, output, session, whole_dataset, logFC, pval, option){
  
  top <- reactive({
    req(whole_dataset())
    req(logFC())
    req(pval())
    req(option())
    logFC <- logFC()
    pval <- pval()
    option <- option()
    whole_data <- whole_dataset()
    splitname <- function(name){
      unlist(strsplit(name, '_'))[1]
    }
    colnames(whole_data) <- lapply(colnames(whole_data), splitname)
    top_part = whole_data[whole_data[,option] <= pval & abs(whole_data[,'logFC'])>=logFC,]
    top_tofile = top_part[,c('logFC', 'P.Value', 'adj.P.Val')]
    return(top_tofile)
  })
  
  return(top)
}

# Download ggplot pic
DownloadPicGG <- function(input, output, session, plotInput)
  output$dndPlotGG <- downloadHandler(
    filename = function() {'pic.pdf'},
    content = function(file) {
      ggsave(file,plotInput())
    }
  )

# Download non-ggplot pic
DownloadPic <- function(input, output, session, plotInput)
  output$dndPlot <- downloadHandler(
    filename = function(){'pic.pdf'},
    content = function(file){
      pdf(file)
      print(plotInput())
      dev.off()
    }
  )

# Download files
DownloadButton <- function(input, output, session, finalfile){    
  output$download <- downloadHandler(
    filename = 'test.txt',
    content = function(file) {
      write.table(finalfile(), file, sep='\t', col.names=NA)
    }
  )
}






