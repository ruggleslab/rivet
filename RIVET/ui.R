library(shinythemes)
library(shinyjs)

ui <- fluidPage(
  shinyUI(
    tagList(
      useShinyjs(),
      div(
        id = "norm_tab",
    navbarPage(
      'RIVET',
      theme = shinytheme("sandstone"),
      tabPanel('Home',
        fluidPage(
          titlePanel("RIVET"),
          fluidRow(
            column(7,
              'Ribosomal Investigation and Visualization to Evaluate Translation is a simple to use',
              'tool to automate the statistical analysis of RNA seq data acquired from polysome profile',
              'or ribosome footprinting experiments.',
              h3("Methods Workflow"),
              wellPanel(
                img(src='for_website3.png', height=350)
              )
            ),
            column(5,
              tabsetPanel(
                tabPanel("About",
                  h3('Resources'),
                  wellPanel(
                    'Example files:',
                            tags$li(
                              a('Geter et al.', href = 'https://github.com/ruggleslab/rivet/blob/master/example_files/Geter_et_al.txt', target = "_blank")
                              ),
                            tags$li(
                              a('Hsieh et al.', href = 'https://github.com/ruggleslab/rivet/blob/master/example_files/Hseih_et_al.txt', target = '_blank')
                              ),
                            tags$li(
                              a('Silvera et al.', href = 'https://github.com/ruggleslab/rivet/blob/master/example_files/Silvera_et_al.txt', target = '_blank')
                            ),
                    br(),
                    'Citations for example files:',
                            tags$li(
                              'Geter et al. (2017).  Hyperactive mTOR and MNK1 phosphorylation of eIF4E confer tamoxifen resistance and estrogen independence through selective mRNA translation reprogramming. Genes Dev. 31:2235–49.',
                              a('Geter et al. paper', href = 'http://genesdev.cshlp.org/content/31/22/2235.long', target='blank')
                            ),
                            tags$li(
                              'Hsieh et al. (2012), The translational landscape of mTOR signalling steers cancer initiation and metastasis. Nature. 485:55–61.',
                              a('Hsieh et al. paper', href = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3663483/', target='blank')
                            ),
                            tags$li(
                              'Silvera et al. (2017).  mTORC1 and -2 Coordinate Transcriptional and Translational Reprogramming in Resistance to DNA Damage and Replicative Stress in Breast Cancer Cells.  Mol Cell Biol. 37:e00577-16.',
                              a('Silvera et al. paper', href = 'http://mcb.asm.org/content/37/5/e00577-16.long)', target='blank')
                            )
                  ),
                  h3('External Information'),
                  wellPanel(
                    p('RIVET is an open-source project developed in the Ruggles Lab',
                      'at NYU Langone Medical Center.',
                      br(),
                      a('Ruggles lab homepage', href='http://www.ruggleslab.org',target="_blank")
                    ),
                    p('Source code for RIVET can be found here.',
                      br(),
                      a('Ruggles lab github', href='http://www.github.com/ruggleslab/rivet', target='_blank')
                    ),
                    p('For additional information on using RIVET or to cite RIVET in your work,',
                      'please refer to the following paper:')
                  )
                ),
                tabPanel("Instructions",
                  tags$ol(h3('To Use RIVET:'),
                    tags$li('Upload data into Normalization tab'), 
                    tags$li('Select statistic and input experiment and control names.'),
                    tags$li('Input number of polysome fractions to be analyzed.'),
                    tags$li('Using sample tab, select samples that are transcription and samples that belong to each polysome fraction.'),
                    tags$li('Data will be calculated.',
                            'To download results move to Transcription, Polysome,',
                            'and / or Translational Efficiency Tabs.')
                  )
                )
              )
            )      
          )
        )
      ),
      tabPanel('Normalization',
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              
              # test dataset
              checkboxInput("testme", "Try a Sample Dataset!", 
                            value = FALSE),
              
              # upload input file
              tagList(
                div(id='user_input',
                  fileInput('file1', 'Choose CSV File',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv')
                  ),
                  # file separator
                  radioButtons('sep', 'Separator',
                           c(Tab='\t', Comma=','),
                           '\t')
                )
              ),
              
              #reset analysis button
              actionButton("reset", "Reset form"),
              
              # platform choice
              radioButtons('platform', 'RNAseq or microarray',
                           choices = c(
                             'Microarray' = 'ma',
                             'RNA seq' = 'rs'
                           ),'rs'),
              
              # statistical choice for RNA seq
              conditionalPanel(
                condition = "input.platform == 'rs'", radioButtons('stat', 'Statistical Method',
                                                                   choices = c('EdgeR' = 'eR',
                                                                               'Limma' = 'limma'),
                                                                   'limma')
              ),
              
              # defines number polysome fractions
              numericInput("numInputs", "How many polysome fractions?", 1),
              
              # downloads results stats matrix of all genes
              DownloadButtonUI('stat_test')
              
            ),
            mainPanel(
              textOutput("tutorialGroup"),
              tabsetPanel(
                
                # MDS plot
                tabPanel("Plot", 
                         plotOutput("MDSPlot"),
                         DownloadPicUI('MDS')
                         ),
                
                # Transcription and translation sample input
                tabPanel("Samples",
                  
                  fluidRow(
                    h4("Sample Labels"),
                      fluidRow(
                        column(6,
                              # defines name of control
                               textInput("control_name", label="Control label")
                          ),
                             column(6,
                                    # defines experiment name
                                    textInput("exp_name", label="Experiment label")
                             )
                      )
                  ),
                  hr(),
                  fluidRow(
                    h4("Transcription"),
                    uiOutput("txgroup")
                  ),
                  hr(),
                  fluidRow(
                    h4("Translation"),
                    # place to hold dynamic inputs
                    uiOutput("inputGroup")
                  ),
                  textOutput("group_warning")
                )
              )
            )
          )
        )
      ),
      tabPanel('Transcription',
        fluidPage(
          titlePanel("Top Gene Selection"),
          sidebarLayout(
            sidebarPanel(
              
              # slider threshold
              sliderInputFcUI('tx_parameters'),
              
              # download top genes matrix
              DownloadButtonUI('top_tx_dl')
              
            ),
            mainPanel(
              
              # busy calculating gif
              tagList(
                tags$head(
                  tags$link(rel="stylesheet", type="text/css",href="style.css"),
                  tags$script(type="text/javascript", src = "busy.js")
                )
              ),
              div(class = "busy",  
                  p("Calculation in progress.."), 
                  img(src="giphy.gif")
              ),
              div(class = "span8",
                  
                  # volcano plot
                  plotOutput('transcriptionPlot'),
                  DownloadPicGGUI('reg_tx')
                  
              )
            )
          )
        )
      ),
      tabPanel('Translation',
        fluidPage(
          titlePanel("Top Gene Selection"),
          sidebarLayout(
            sidebarPanel(
              
              # slider for top genes translation
              sliderInputFcUI('tl_parameters'),
              
              # download top genes
              DownloadButtonUI('top_tl_dl')
              
            ),
            mainPanel(
              tagList(
                tags$head(
                  tags$link(rel="stylesheet", type="text/css",href="style.css"),
                  tags$script(type="text/javascript", src = "busy.js")
                )
              ),
              div(class = "busy",  
                  p("Calculation in progress.."), 
                  img(src="giphy.gif")
              ),
              div(class = "span8", 
                  plotOutput('translationPlot'),
                  uiOutput('reg_tl')
              )
            )
          )
        )
      ),
      tabPanel('Translational Efficiency',
        fluidPage(
          titlePanel("Top Gene Selection"),
          sidebarLayout(
            sidebarPanel(
              
              # choice of translation efficiency method
              radioButtons('teChoice', 'Translational Efficiency Method',
                           choices = c('Log2 ratio' = 'teRatio',
                                       'Interaction' = 'int'),
                           'int'),
              
              # parameter selection for translational efficiency
              sliderInputFcUI('te_parameters'),
              
              # download top genes translational efficiency
              DownloadButtonUI("te_reg")
              
            ),
            mainPanel(
              tagList(
                tags$head(
                  tags$link(rel="stylesheet", type="text/css",href="style.css"),
                  tags$script(type="text/javascript", src = "busy.js")
                )
              ),
              div(class = "busy",  
                  p("Calculation in progress.."), 
                  img(src="giphy.gif")
              ),
              div(class = "span8", 
                  plotOutput('tePlot'),
                  uiOutput('reg_te')
              )
            )
          )
        )
      ),
      tabPanel('Translational Regulation',
        sidebarLayout(
          sidebarPanel(
            
            # output which poly group is displayed as checkbox
            uiOutput("polyGroup"),
            
            # choose type of translational regulation for scatter and download
            radioButtons('col', 'Regulation',
                         choices = c('None' = 'none',
                                     'Transcription Alone'= 'tx',
                                     'Translation Alone' = 'tl',
                                     'Opposite' = 'opp',
                                     'Transcription and Translation' = 'txtl',
                                     'All Regulation' = 'all',
                                     'Translational efficiency' = 'te'),
                         'none'),
            
            # only download selected type of translational regulation
            DownloadButtonUI("regulation")
            
          ),
          mainPanel(
            tagList(
              tags$head(
                tags$link(rel="stylesheet", type="text/css",href="style.css"),
                tags$script(type="text/javascript", src = "busy.js")
              )
            ),
            div(class = "busy",  
                p("Calculation in progress.."), 
                img(src="giphy.gif")
            ),
            div(class = "span8", 
                splitLayout(cellWidths = c("70%", "30%"),
                            plotOutput("ScatterPlot"),
                            conditionalPanel(
                              condition = "input.col != 'none'", plotOutput('barPlot')
                            )
                ),
                splitLayout(cellWidths = c("70%", "30%"),
                            DownloadPicGGUI("reg_scatter"),
                            conditionalPanel(
                              condition = "input.col != 'none'", DownloadPicGGUI('reg_bar')
                            )
                )
            )
          )
        )
      )
    )
      )
    )
  )
)
          