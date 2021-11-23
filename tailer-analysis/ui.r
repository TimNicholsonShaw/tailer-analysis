# Modular Pieces

options_ui <- function(id, dots=F, position=T, analysis_window=F, pdisplay=F, mature_end=F, PT_tail=F){
  ns <- NS(id)
  tags<- tagList(
    textInput(ns("gene_name"), "Gene Name"),
    checkboxInput(ns("multi_loc"), "Multi-Locus Gene"))
  
  if (position==T){
    tags <- tagAppendChild(tags, numericInput(ns("x_min"), "Minimum Position", value=-10))
    tags <- tagAppendChild(tags, numericInput(ns("x_max"), "Maximum Position", value=10))
  }
  if (PT_tail){
    tags <- tagAppendChild(tags, numericInput(ns("y_min"), "Minimum Y-Value", value=0))
    tags <- tagAppendChild(tags, numericInput(ns("y_max"), "Maximum Y-Value", value=5))
  }
  else{
    tags <- tagAppendChild(tags, numericInput(ns("y_min"), "Minimum Y-Value", value=0))
    tags <- tagAppendChild(tags, numericInput(ns("y_max"), "Maximum Y-Value", value=1))

  }

  if(analysis_window==T) {
    tags <- tagAppendChild(tags, numericInput(ns("analysis_min"), "Start of Analysis Window", value=-100))
    tags <- tagAppendChild(tags, numericInput(ns("analysis_max"), "End of Analysis Window", value=100))
  }
  if(mature_end){
    tags <- tagAppendChild(tags, numericInput(ns("mature_end"), "Mature 3' End", value=0))
  }
  if (pdisplay==T){
    tags <- tagAppendChild(tags, checkboxInput(ns("pdisplay"), "p values"))
  }
  if (dots==T) {
    tags <- tagAppendChild(tags, checkboxInput(ns("dots"), "Dots?"))
    }
  tags <- tagAppendChild(tags, actionButton(ns("make_plot"), "Make Plot"))
  tags
}
plot_ui <- function(id){
  ns<-(NS(id))
  tagList(
    plotOutput(ns("plot"))%>% withSpinner(color="#0dc5c1"),
    downloadButton(ns("download_plot"), "Download Plot"),
    downloadButton(ns("download_data"), "Download Data"),
    fluidRow(
      column(6,numericInput(ns("height"), "Plot Height", value=2500)),
      column(6, numericInput(ns("width"), "Plot Width", value=3500))
    ),
    fluidRow(
      tableOutput(ns("n_table"))
    )
  )
}

whole_plot_page_ui <- function(id, dots=F, pdisplay=F, analysis_window=F, mature_end=F, position=T, PT_tail=F){
  fluidPage(
    sidebarLayout(
      sidebarPanel(options_ui(id, dots=dots, pdisplay=pdisplay, analysis_window=analysis_window, mature_end=mature_end, position=position, PT_tail=PT_tail)),
      mainPanel(plot_ui(id))
    )
  )
}

#Pieces
ui_files <- fluidPage(
  #Input files and give them experimental groupings
  fileInput("tail_files", "Upload Tail Files", multiple=TRUE, accept=c(".csv")),
  DTOutput("sample_table"),
  actionButton("make_df_button", "Format Data"),
  verbatimTextOutput("df_preview") %>% withSpinner(color="#0dc5c1"),
  uiOutput("sample_order")
)
ui_candidate_finder <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            uiOutput("can_con1"),
            uiOutput("can_con2"),
            actionButton("find_cans_button", "Find Candidates"),
            downloadButton("download_cans", "Download Data"),
            p(" "),
            p("Dynamically updated options:"),
            numericInput("min_cans", "Minimum Reads For a Candidate", 10, min=1),
            numericInput("delN_cutoff", "Average Δ End PositionCutoff", value=0),
            numericInput("pval_cutoff", "p-value End Position Cutoff", value=1)
        ),
        mainPanel(
            DTOutput("candidates") %>% withSpinner(color="#0dc5c1"),
            DTOutput("n_candidates")
            )
    )
)

ui_statistics_page <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      textInput("stats_gene_name", "Gene Name"),
      uiOutput("con1"),
      uiOutput("con2"),
      actionButton("get_stats", "Get Stats"),
      checkboxInput("stat_multi_loc", "Multi Locus")
    ),
    mainPanel(
      uiOutput("stat_gene_name"),
      p(strong("End Position KS-test Matrix")),
      tableOutput("end_position_tab"),
      p(strong("Pooled End Position KS-test")),
      textOutput("end_position_text"),
      p(" "),
      p(strong("Tail Length KS-test Matrix")),
      tableOutput("tail_len_tab"),
      p(strong("Pooled Tail Length KS-test")),
      textOutput("tail_len_text"),
      p(" "),
      p(strong("Number of Observations")),
      tableOutput("n_df_table")
      )
    )
)




# Actual UI declaration
ui <- fluidPage(
    titlePanel(
        "Tailer Analysis"
    ),
    tabsetPanel(
        tabPanel("File Upload", ui_files),
        tabPanel("Candidate Finder", ui_candidate_finder),
        tabPanel("Tail Bar Graph", whole_plot_page_ui("tail_bar", analysis_window=T, mature_end=T)),
        tabPanel("Cumulative Plot", whole_plot_page_ui("cum_plot", dots=T, analysis_window=T, mature_end=T)),
        tabPanel("Tail Logo Plot", whole_plot_page_ui("tail_logo", analysis_window=T, mature_end=T)),
        tabPanel("Post-Transcriptional Tailing", whole_plot_page_ui("pt_tail", pdisplay=T, analysis_window=T, mature_end=T, position=F, PT_tail=T)),
        tabPanel("Statistics", ui_statistics_page)
    )
)

