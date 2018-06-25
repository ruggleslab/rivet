library(shiny)
library("limma")
library('ggplot2')
library("edgeR")
library('gplots')
library('Biobase')
library(RColorBrewer)
library(shinythemes)
library(shinyjs)

# file size
options(shiny.maxRequestSize=50*1024^2)

server <- function(input , output, session){
  
  observeEvent(input$testme,{
    if (input$testme) {
      hide('user_input')
    } else {
      show('user_input')
    }
  })
  
  # input matrix
  expMatrix <- reactive({
    if (input$testme) {
      # test file
      inputfile <- read.table("www/Geter_et_al.txt", header=TRUE, sep='\t', row.names=1)
    } else {
      # user inputs count matrix
      req(input$file1)
      inFile <- input$file1
      inputfile <- read.table(inFile$datapath, sep=input$sep, header=TRUE, row.names=1)
     }
    exprs = as.matrix(inputfile)
    return(exprs)
  })
  
  #test whether contrasts object is reset when reset button pushed
  observeEvent(input$reset, {
    reset("norm_tab")
  })
  
  # update polysome input and label input if the test me data is chosen
  observeEvent(input$testme, {
    if (input$testme){
      output$tutorialGroup <- renderText({"Groups have been named and assigned!\n Please continue to the next Tab"})
      updateNumericInput(session, "numInputs", value = 2)
      updateTextInput(session, "control_name", value= 'NS')
      updateTextInput(session, "exp_name", value= 'E4')
    } else {
      output$tutorialGroup <- renderText({""})
      updateNumericInput(session, "numInputs", value = 1) 
      updateTextInput(session, "control_name", value = "")
      updateTextInput(session, "exp_name", value = "")
    }
  })
  
  # input transcription scamples, ui
  observeEvent(input$testme,{
    
    output$txgroup <- renderUI({
      # test samples
      if (input$testme) {
        ctrltx = c('R1_TamR_NS_T','R2_TamR_NS_T')
        exptx = c('R1_TamR_4E_T', 'R2_TamR_4E_T')
        fluidRow(
          fluidRow(
            column(6, 
                   selectInput("tx1", label = "",
                   choices = ctrltx, multiple = TRUE, selected = ctrltx)
            ),
            column(6, 
                   selectInput("tx2", label = "", choices = exptx, 
                   multiple = TRUE, selected = exptx)
            )
          )
        )
       } else {
         # samples from file
         sampleUploadUI("tx")
       }
    })
  })
  
  # collect transcription samples
  transcription <- c()
  makeReactiveBinding('transcription')
  observeEvent(input$testme,{
    if (input$testme) {
      ctrltx = c('R1_TamR_NS_T','R2_TamR_NS_T')
      exptx = c('R1_TamR_4E_T', 'R2_TamR_4E_T')
      transcription <<- reactive({list(ctrltx,exptx)})
    } else {
      transcription <<- callModule(sampleUpload, "tx", expMatrix)
    }
  })
  
  # input translation samples based on number of polysome fractions, output reordered matrix  
  results <- c()
  makeReactiveBinding('results')
  
  # observe changes in "numInputs" (polysome input), and create corresponding number of inputs
  observeEvent(input$numInputs, {
      
    output$inputGroup = renderUI({
      if (input$testme){
        ctrl1 = c('R1_TamR_NS_HD','R2_TamR_NS_HD')
        exp1 = c('R1_TamR_4E_HD', 'R2_TamR_4E_HD')
        ctrl2 = c('R1_TamR_NS_LD','R2_TamR_NS_LD')
        exp2 = c('R1_TamR_4E_LD','R2_TamR_4E_LD')
        fluidRow(
          fluidRow(column(6, 
                        selectInput("input1", label = "",
                        choices = ctrl1, multiple = TRUE, selected = ctrl1)
                 ),
                 column(6, 
                        selectInput("input2", label = "", choices = exp1, 
                        multiple = TRUE, selected = exp1)
                 )
          ),
          fluidRow(column(6, 
                        selectInput("input3", label = "", choices = ctrl2, 
                        multiple = TRUE, selected = ctrl2)
                  ),
                   column(6, 
                        selectInput("input4", label = "", choices = exp2, 
                        multiple = TRUE, selected = exp2)
                   )
          ) 
        )
      } else {
        # warning to user to upload file if samples aren't present
        validate(
          need(input$file1, 'Please provide a file in the Upload Tab')
          ) 
        input_list <- lapply(1:input$numInputs, function(i) {
          inputName <- paste0("input", i)
          sampleUploadUI(inputName)
        })
      
      }
    })
    
    # combine all polysome fractions into a list per fraction of sublists of ctrl and exp
    results <<- lapply(1:input$numInputs, function(i) {
      inputName <- paste0("input", i)
      callModule(sampleUpload, inputName, expMatrix)
    })
  })
  
  # translation results
  translation <- reactive({
    if (input$testme) {
      # example data
      ctrl1 = c('R1_TamR_NS_HD','R2_TamR_NS_HD')
      exp1 = c('R1_TamR_4E_HD', 'R2_TamR_4E_HD')
      ctrl2 = c('R1_TamR_NS_LD','R2_TamR_NS_LD')
      exp2 = c('R1_TamR_4E_LD','R2_TamR_4E_LD')
      g1 = list(ctrl1,exp1)
      g2 = list(ctrl2,exp2)
      glist = list(g1,g2)
      glist
    } else {
    # user input  
    lapply(1:input$numInputs, function(i) {
      results[[i]]()
    })
    }
  })
  
  # warning the user of duplicate sample assignment
  output$group_warning <- renderText({
    req(transcription())
    req(translation())
    tx_list = c(unlist(transcription()), unlist(translation()))
    tx_dups = duplicated(tx_list)
    if ('TRUE' %in% tx_dups) {
      message = c("Warning: A sample has been assigned to more than one group!")
    } else {
      message = c(" ")
    }
  })
  
  # reorder the input matrix based on input transcription and translation samples
  reorder_mat <- reactive({
    req(transcription())
    req(translation())
    exprs_reorder = expMatrix()[,c(unlist(transcription()), unlist(translation()))]
    return(exprs_reorder)
  })  
  
  # normalize for MDS plot input
  normalize <- callModule(dge, "mds_plot", expMatrix)
  
  # normalize for statistics
  normalize_reorder <- callModule(dge, "edgeRinput", reorder_mat)
  
  #MDS plot of normalized data
  MDSInput <- reactive({
    req(normalize())
    plotMDS(normalize(), cex=0.6)
  })
  
  observeEvent(input$reset, {
    output$MDSPlot <- renderPlot({
      print(MDSInput())
    })
  })
  
  # displays MDS plot
  output$MDSPlot <- renderPlot({
    print(MDSInput())
  })
  
  # creates input for pic download
  MDSInput2 <- function(){
    req(normalize())
    plotMDS(normalize(), cex=0.6)
  }
  
  # download MDS plot
  callModule(DownloadPic, 'MDS', MDSInput2)
  
  # voom normalize for limma
  voomnorm <- reactive({
    req(transcription())
    req(translation())
    req(normalize_reorder())
    norm = voom(normalize_reorder(), plot=FALSE)
    return(norm$E)
  })

  # treatments - ordered list of labels
  treatments <- reactive({
    req(input$numInputs)
    req(input$control_name)
    req(input$exp_name)
    polynames <- sapply(1:input$numInputs, function(i){
      c(paste(input$control_name, 'translation', i, sep='.'),
      paste(input$exp_name, 'translation', i, sep='.'))
    })
    treat <- c(paste(input$control_name, 'transcription', sep='.'), paste(input$exp_name, 'transcription', sep='.'),
               polynames)
    return(treat)
  })
  
  # targets frame
  designMatrix <- reactive({
    req(transcription())
    req(translation())
    req(input$numInputs)
    req(reorder_mat())
    req(treatments())
    tx <- sapply(transcription(), length)
    tl <- sapply(1:input$numInputs, function(i){
      sapply(translation()[[i]], length)
    })
    sample_num <- c(tx, tl)
    #targets are created by setting the levels to be equal to 1+number of poly inputs
    #and then replicate number is the times argument
    targets = data.frame(samples = colnames(reorder_mat()), levels = rep(1:((1+input$numInputs)*2), times=sample_num))
    design <- model.matrix(~0+factor(targets$levels))
    colnames(design) = treatments()
    rownames(design) <- targets$samples
    return(design)
  })

  # contrasts 
  contrasts <- reactive({
    req(input$exp_name)
    req(input$control_name)
    req(input$numInputs)
    req(treatments())
    tx_label <- paste(paste(input$exp_name, 'transcription',sep='.'), paste(input$control_name, 'transcription',sep='.'), sep='-')
    tl_label <- sapply(1:input$numInputs, function(i){
      paste(paste(input$exp_name, 'translation',i,sep='.'), paste(input$control_name, 'translation',i,sep='.'), sep='-')
    })
    diff <- sprintf("(%s)-(%s)", tl_label, tx_label)
    contrasts_label <- c(tx_label, tl_label, diff)
    contrast.matrix <- makeContrasts(contrasts = contrasts_label,  levels=treatments())
    return(contrast.matrix)
  })

  # limma stat
  limmastat <- list()
  makeReactiveBinding('limmastat')
    
  observeEvent(contrasts(), {
    limmastat <<- lapply(1:length(colnames(contrasts())), function(i){
      compare <- reactive(colnames(contrasts())[i])
      tagname <- paste0("stat", i)
      switch(input$platform,
        ma = callModule(limma, tagname, reorder_mat, contrasts, designMatrix, compare),
        rs = callModule(limma, tagname, voomnorm, contrasts, designMatrix, compare)
      )
    })
  })

  # dataframe of limma results for all genes
  limmatest <- reactive({
    result <- lapply(1:(length(colnames(contrasts()))-input$numInputs), function(i) {
      limmastat[[i]]()
    })
    as.data.frame(result)
  })

  # edgeRstat
  edgeRstat <- list()
  makeReactiveBinding('edgeRstat')

  observeEvent(contrasts(), {
  
    edgeRstat <<- lapply(1:length(colnames(contrasts())), function(i){
      compare <- reactive(colnames(contrasts())[i])
      tagname <- paste0("stat_edge", i)
      callModule(edgeR, tagname, normalize_reorder, contrasts, designMatrix, compare)
    })
  })

  # dataframe of edgeR results for all genes
  edgetest <- reactive({
    result <-lapply(1:length(colnames(contrasts())), function(i) {
      edgeRstat[[i]]()
    })
    as.data.frame(result)
  })

  # choose between limma and edgeR
  norm_options <- reactive({
    req(edgetest())
    req(limmatest())
    switch(input$stat,
      eR = edgetest(),
      limma = limmatest()
    )
  })  

  # to download all results for transcription and translation
  callModule(DownloadButton, 'stat_test', norm_options)
  
# transcription-------------------------------------------------------------------------------
  
  transcription_results <- reactive({
    req(norm_options())
    norm_options()[,1:3]
  })
  
  # params for transcription
  param_tx <- callModule(sliderInputFC, 'tx_parameters')
  
  # fold-change transcription
  fc_tx <- reactive({
    req(param_tx())
    param_tx()[[1]]
  })
  
  # statistic transcription, number
  pval_tx <- reactive({
    req(param_tx())
    param_tx()[[2]]
  })
  
  # either will be 'pvalue' or 'adjusted pval'
  option_tx <- reactive({
    param_tx()[[3]]
  })
  
  # volcano plot transcription
  tx_plot <- callModule(volcanoPlot, 'tx_vol', transcription_results, fc_tx, pval_tx, option_tx)
  output$transcriptionPlot <- renderPlot({
    req(tx_plot())
    print(tx_plot())
  })
  callModule(DownloadPicGG, 'reg_tx', tx_plot)
  
  # top genes for download for transcription
  top_tx <- callModule(top_genes, 'tx_log_FC', transcription_results, fc_tx, pval_tx, option_tx)
  callModule(DownloadButton, 'top_tx_dl', top_tx)
   
# translation--------------------------------------------------------------------------------------
  
  # translation parameters from slider
  param_tl <- callModule(sliderInputFC, 'tl_parameters')
  
  # fold-change translation
  fc_tl <- reactive({
    req(param_tl())
    param_tl()[[1]]
    })
  
  # statistic translation
  pval_tl <- reactive({
    req(param_tl())
    param_tl()[[2]]
  })
  
  # 'p-value' or 'adjusted p-value' translation, bug fix********************
  option_tl <- reactive({
    param_tl()[[3]]
  })
  
  # results matrix for translation
  translation_results <- reactive({
    req(norm_options())
    norm_options()[,4:length(colnames(norm_options()))]
  })

  # creates swtich for user to choose between polysome fractions
  tl_tags <- reactive({
    req(input$numInputs)
    values <- c()
    choices <- c()
    for (i in seq(1, input$numInputs)){
      choices[[i]] <- paste0('Poly',i) 
      values[[i]] <- paste0('poly',i)
    }
    buttons <- setNames(values, choices)
    return(buttons)
  })
 
  # saves the translation subframe as a named list of poly fractions dataframes  
  translation_sub <- reactive({
    req(translation_results())
    subframes <-lapply(seq(1, length(colnames(translation_results())),3), function(i){
      translation_results()[,seq(i,i+2)]
    }
    )
    names(subframes) <- as.vector(tl_tags())
    return(subframes)
  })  
  
  # isolates top genes for translation
  tl_top <- reactive({
    req(translation_results())
    req(fc_tl())
    req(pval_tl())
    j = 0  
    results <- c()
    for (i in seq(1, length(colnames(translation_results())),3)){
      j = j + 1
      tl_compare <- translation_results()[,seq(i,i+2)]
      tagname <- paste0("tl_top", i)
      top <- callModule(top_genes, tagname, reactive(tl_compare), fc_tl, pval_tl,option_tl)
      results[[j]] <- rownames(top())
    }
    topnames = unique(unlist(results))
    return(topnames)
  })
  
  # top translation across all fractions
  regulated <- reactive({
    req(translation_results())
    as.matrix(translation_results()[tl_top(),])
  })

  # top translation across all fractions for heatmap
  regulated_foldchange <- reactive({
    req(translation_results())
    fc <- as.matrix(translation_results()[tl_top(),seq(1, length(colnames(translation_results())),3)])
    colnames(fc) <- sapply(colnames(fc), function (x) substr(x=x,start=nchar(x),stop=nchar(x)) )
    return(fc)
  })

  # number_of_inputs
  num_inputs <- reactive({input$numInputs})

  # bugfix for tl volcano plot *****************************************************************
  # either volcano or heatmap depending on input number of polysome fractions, generates plot
  tlPlot <- reactive({
    req(num_inputs())
    if (num_inputs()>1){
      tl_plot <- callModule(heatmap, 'tl_heat', regulated_foldchange)
    } else {
      tl_plot <- callModule(volcanoPlot, 'tl_vol', translation_results, fc_tl, pval_tl, option_tl)
    }
    return(tl_plot())
  })
  
  #generates figure for download
  observeEvent(input$numInputs,{
    
    output$reg_tl <- renderUI({
      if (input$numInputs > 1){
        DownloadPicUI('reg_tl_heatmap')
      } else {
        DownloadPicGGUI('reg_tl_volcano')
      }
    })
    
    if (input$numInputs > 1){
      tl_dwnld_htmap <- callModule(heatmapDwnInput, 'tl_ht', regulated_foldchange)
       callModule(DownloadPic, 'reg_tl_heatmap', tl_dwnld_htmap)
    } else {
      callModule(DownloadPicGG, 'reg_tl_volcano', tlPlot)
    }
    
  })

  # display translation plot
  output$translationPlot <- renderPlot({print(tlPlot())})
  
  # download translation top genes
  callModule(DownloadButton, 'top_tl_dl', regulated)

# TE log2 ratio method--------------------------------------------------------------------------------------------------

  # unlist translation samples so can be used in ratio
  te_tl <- reactive({
    req(translation())
    tl<-voomnorm()[,unlist(translation())]
    return(tl)
  })
  
 # unlist transcription samples to be used in ratio
  te_tx <- reactive({
    req(transcription())
    req(input$numInputs)
    tx_te <- rep(unlist(transcription()), input$numInputs)
    tx <- voomnorm()[,tx_te]
    return(tx)
  })
  
  # log ratio
  log2ratio <- reactive({
    req(te_tl())
    req(te_tx())
    if (length(colnames(te_tl())) == length(colnames(te_tx()))){
      ratio <- te_tl() - te_tx()
      return(ratio)
    } 
  })
  
  # targets frame for translational efficiency
  designMatrix_te <- reactive({
    req(input$numInputs)
    req(log2ratio())
    req(input$control_name)
    req(input$exp_name)
    tl <- sapply(1:input$numInputs, function(i){
      sapply(translation()[[i]], length)
    })
    tx <- sapply(transcription(), length)
    targets = data.frame(samples = colnames(log2ratio()), levels = rep(1:(input$numInputs*2), times=tl))
    design <- model.matrix(~0+factor(targets$levels))
    colnames(design) <- sapply(1:input$numInputs, function(i){
      c(paste(input$control_name, 'te', i, sep='.'),
        paste(input$exp_name, 'te', i, sep='.'))
    })
    rownames(design) <- targets$samples
    return(design)
  })
  
  # contrasts translational efficiency 
  contrasts_te <- reactive({
    req(input$exp_name)
    req(input$control_name)
    req(designMatrix_te())
    contrasts_label <- sapply(1:input$numInputs, function(i){
      paste(paste(input$exp_name, 'te',i,sep='.'), paste(input$control_name, 'te',i,sep='.'), sep='-')
    })
    contrast.matrix <- makeContrasts(contrasts = contrasts_label,  levels=colnames(designMatrix_te()))
    return(contrast.matrix)
  }) 
  
  # limma_te
  limma_te <- list()
  makeReactiveBinding('limma_te')
  
  observeEvent(contrasts_te(), {
    limma_te <<- lapply(1:length(colnames(contrasts_te())), function(i){
      compare_te <- reactive(colnames(contrasts_te())[i])
      tagname_te <- paste0("stat_te", i)
      callModule(limma, tagname_te, log2ratio, contrasts_te, designMatrix_te, compare_te)
    }
    )
  })
  
  # translational efficiency interaction, this dataframe is calculated using contrast matrix defined above
  limmatest_int <- reactive({
    result <- lapply((length(colnames(contrasts()))-input$numInputs):length(colnames(contrasts())), function(i) {
      limmastat[[i]]()
    })
    as.data.frame(result)
  })
  
  # limma results for log ratio
  limmaTE_frame <- reactive({
    req(limma_te)
    req(contrasts_te())
    
    result_te <- lapply(1:length(colnames(contrasts_te())), function(i) {
      limma_te[[i]]()
    })
    as.data.frame(result_te)
  })
  
  # interaction or log ratio option choices
  limmaTEtest <- reactive({
    switch(input$teChoice,
           teRatio = limmaTE_frame(),
           int = limmatest_int())
  })
  
  # slider bar for log2 ratio, fixed bug option_te**********************************************
  param_te <- callModule(sliderInputFC, 'te_parameters')
  fc_te <- reactive({
    req(param_te())
    param_te()[[1]]
    })
  pval_te <- reactive({
    req(param_te())
    param_te()[[2]]
    })
  option_te <- reactive({
    param_te()[[3]]
  })
  
  # top te genes
  te_top <- reactive({
    req(limmaTEtest())
    req(fc_te())
    req(pval_te())
    j = 0  
    results <- c()
    for (i in seq(1, length(colnames(limmaTEtest())),3)){
      j = j + 1
      tl_compare <- limmaTEtest()[,seq(i,i+2)]
      tagname <- paste0("tl_top", i)
      top <- callModule(top_genes, tagname, reactive(tl_compare), fc_te, pval_te, option_te)
      results[[j]] <- rownames(top())
    }
    topnames = unique(unlist(results))
    return(topnames)
  })
  
  # top te across all fractions
  regulated_te <- reactive({
    req(translation_results())
    as.matrix(translation_results()[te_top(),])
  }) 
  
  callModule(DownloadButton, 'te_reg', regulated_te)
  
  # top translational effeciency genes across all polysome fractions for heatmap
  regulated_foldchange_te <- reactive({
    req(translation_results())
    fc <- as.matrix(translation_results()[te_top(),seq(1, length(colnames(translation_results())),3)])
    colnames(fc) <- sapply(colnames(fc), function (x) substr(x=x,start=nchar(x),stop=nchar(x)) )
    return(fc)
  })
  
  # te for polysome tab, if more than 1 poly fraction will make heatmap otherwise volcano plot
  tePlot <- reactive({
    req(num_inputs())
    if (num_inputs()>1){
      te_plot <- callModule(heatmap, 'te_heat', regulated_foldchange_te)
    } else {
      te_plot <- callModule(volcanoPlot, 'te_vol', limmaTEtest, fc_te, pval_te, option_te)
    }
    return(te_plot())
  })
  
  #generates figure for download, te
  observeEvent(input$numInputs,{
    
    output$reg_te <- renderUI({
      if (input$numInputs > 1){
        DownloadPicUI('reg_te_heatmap')
      } else {
        DownloadPicGGUI('reg_te_volcano')
      }
    })
    
    # bugfix regulated_foldchange_te*********************************************************************
    if (input$numInputs > 1){
      te_dwnld_htmap <- callModule(heatmapDwnInput, 'te_ht', regulated_foldchange_te)
      callModule(DownloadPic, 'reg_te_heatmap', te_dwnld_htmap)
    } else {
      callModule(DownloadPicGG, 'reg_te_volcano', tePlot)
    }
    
  })
  
  output$tePlot <- renderPlot({print(tePlot())})

  #bugfix this is duplicated, commented out********************************************************************
  #callModule(DownloadButton, 'te_reg', regulated)


# translational regulation--------------------------------------------------------------------

  # display different graphs for different poly fractions    
  output$polyGroup = renderUI({
    req(tl_tags())
    radioButtons('polynum', 'Fraction',
                 choices = tl_tags(),
                 'poly1')
  })

  # grabs data for chosen poly fraction
  poly <- reactive({
    req(input$polynum)
    req(translation_sub())
    poly_frac <- {input$polynum}
    translation_sub()[[poly_frac]]
  })

  # combines transcription results and chosen poly results into a matrix
  bigmat <- reactive({
    req(transcription_results())
    req(poly())
    transcrip <- transcription_results()
    transl <- poly()
    cbind(transcrip, transl)
  })


  #translational efficiency scatter matrix, for 1 poly fraction alone
  #creates regulation scatter matrix based on user selection
  scatter_mat <- reactive({
    req(bigmat())
    req(input$col)
    req(te_top())
  
    bigmat <- bigmat()
  
    # whether to use p-value or adj pvalue transcription
    if (option_tx() == 'P.Value'){
      tx_stat = bigmat[,2]
    } else {
      tx_stat = bigmat[,3]
    }
    
    # whether to use p-value or adj pvalue translation  
    if (option_tl() == 'P.Value'){
      tl_stat = bigmat[,5]
    } else {
      tl_stat = bigmat[,6]
    }
  
    # add a categorical column to transcription/translation mat to
    # describe the type of translational regulation
    bigmat$Regulation = rep('none',nrow(bigmat))
  
    # Collecting types of regulation
    txUp = bigmat[,1]>= fc_tx() & tx_stat < pval_tx()
    txDown = bigmat[,1]<= -(fc_tx()) & tx_stat < pval_tx()
    txNone = tx_stat > pval_tx()
  
    polyUp = bigmat[,4]>= fc_tl() & tl_stat < pval_tl()
    polyDown = bigmat[,4]<= -(fc_tl()) & tl_stat < pval_tl()
    polyNone = tl_stat > pval_tl()
  
    # transcription/translation
    downUp =  txDown & polyUp
    upDown = txUp & polyDown
    upUp = txUp & polyUp
    downDown = txDown & polyDown
    noneUp = txNone & polyUp
    noneDown = txNone & polyDown
    upNone = txUp & polyNone
    downNone = txDown & polyNone
  
    # categories
    noReg = 1:length(bigmat$Regulation)
    opps = downUp | upDown
    all_change = upUp | downDown
    tl_alone = noneUp | noneDown
    tx_alone = upNone | downNone
    te_alone = rownames(bigmat()) %in% te_top()
  
    # all_reg
    bigmat_all = bigmat
    bigmat_all$Regulation[opps]<-"opposite"
    bigmat_all$Regulation[all_change] <- "transcription and translation"
    bigmat_all$Regulation[tl_alone] <- "translation alone"
    bigmat_all$Regulation[tx_alone] <- "transcription alone"
  
    # user chooses type of regulation
    bigmat_options <- switch(input$col,
                           none= list(noReg,'none'),
                           opp = list(opps, 'opposite'),
                           tx = list(tx_alone,'transcription alone'),
                           tl = list(tl_alone,'translation alone'),
                           txtl = list(all_change, 'transcription and translation'),
                           all = list(1:length(bigmat_all$Regulation), bigmat_all$Regulation),
                           te = list(te_alone, 'translational efficiency'),
                           none)
  
   bigmat$Regulation[bigmat_options[[1]]]<-bigmat_options[[2]]
   return(bigmat)
  })

  # Download only the genes that correspond to translational regulation
  sub_scatter_mat <- reactive({
    req(scatter_mat())
    sub <- scatter_mat()
    genes_reg <- if(length(sub$Regulation[sub$Regulation!='no change'])<length(sub$Regulation)){
      sub[sub$Regulation!='no change',]
    } else {
      sub
    }
    return(genes_reg)
  })

  callModule(DownloadButton, 'regulation', sub_scatter_mat)

  # Input for scatter graph
  ScatterInput <- reactive({ 
    req(scatter_mat())
    mat <- scatter_mat()
  
    # Set log Regulation for graph, creates symmetrical axes
    minFC = abs(min(c(mat[,1], mat[,4])))
    maxFC = abs(max(c(mat[,1], mat[,4])))
    if (minFC>maxFC) {
      lim = minFC
    } else {
      lim = maxFC}
  
    # Set the color pallete for the graph
    cols <- c('opposite' = '#E69F00', 'transcription alone' = '#56B4E9',
              'transcription and translation' = '#009E73', 
              'translational efficiency' = '#F0E442', 'none' = "#000000",
              'translation alone'='#CC79A7')
  
    # Construct the plot object
    ggplot(data=mat, aes(x=mat[,1], y=mat[,4],colour=Regulation)) +
      geom_point(alpha=0.4, size=1.75) +
      scale_colour_manual(values = cols) +
      xlim(c(-lim, lim)) + ylim(c(-lim,lim)) +
      xlab("Transcription FC (Log2)") + 
      ylab("Translation FC (Log2)") 
      
  })

  output$ScatterPlot <- renderPlot({
    req(ScatterInput())
    print(ScatterInput())
  })

  # Download scatter plot
  callModule(DownloadPicGG, 'reg_scatter', ScatterInput)

  # barplot of regulation
  barInput <- reactive({
  
    # Set the color pallete for the graph
    cols <- c('none' = '#000000', 'opposite' = '#E69F00', 'transcription alone' = '#56B4E9',
              'transcription and translation' = '#009E73', 
              'translational efficiency' = '#F0E442', 'none' = "#000000",
              'translation alone'='#CC79A7')
  
    ggplot(scatter_mat(), aes(x=factor(Regulation), fill=Regulation)) +
      geom_bar(stat="count") +
      scale_fill_manual(values = cols) +
      theme(axis.title.x=element_blank(), 
            axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
            legend.position="none") +
      geom_text(stat='count',aes(label=..count.., vjust=-0.5)) 
  })

  output$barPlot <- renderPlot({
    req(barInput())
    print(barInput())
  })

  # download bar plot
  callModule(DownloadPicGG, 'reg_bar', barInput)

}