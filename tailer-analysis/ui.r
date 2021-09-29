# Modular Pieces

options_ui <- function(id, dots=F, position=T, analysis_window=F){
  ns <- NS(id)
  tags<- tagList(
    textInput(ns("gene_name"), "Gene Name"),
    checkboxInput(ns("multi_loc"), "Multi-Locus Gene"))
  
  if (position==T){
    tags <- tagAppendChild(tags, numericInput(ns("x_min"), "Minimum Position", value=-10))
    tags <- tagAppendChild(tags, numericInput(ns("x_max"), "Maximum Position", value=10))
  }
  
  tags <- tagAppendChild(tags, numericInput(ns("y_min"), "Minimum Fraction", value=0))
  tags <- tagAppendChild(tags, numericInput(ns("y_max"), "Maximum Fraction", value=1))
  
  if(analysis_window==T) {
    tags <- tagAppendChild(tags, numericInput(ns("analysis_min"), "Start of Analysis Window", value=-100))
    tags <- tagAppendChild(tags, numericInput(ns("analysis_max"), "End of Analysis Window", value=100))
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
    plotOutput(ns("plot")),
    downloadButton(ns("download_plot"), "Download Plot"),
    fluidRow(
      column(6,numericInput(ns("height"), "Plot Height", value=400)),
      column(6, numericInput(ns("width"), "Plot Width", value=400))
    )
  )
}

whole_plot_page_ui <- function(id, dots=F){
  fluidPage(
    sidebarLayout(
      sidebarPanel(options_ui(id, dots=dots)),
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
  verbatimTextOutput("df_preview") %>% withSpinner(color="#0dc5c1")
)
ui_candidate_finder <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            numericInput("min_cans", "Minimum Reads For a Candidate", 10, min=1),
            actionButton("find_cans_button", "Find Candidates")
        ),
        mainPanel(
            DTOutput("candidates") %>% withSpinner(color="#0dc5c1")
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
        tabPanel("Cumulative Plot", whole_plot_page_ui("cum_plot", dots=T)),
        tabPanel("Tail Bar Graph", whole_plot_page_ui("tail_bar")),
        tabPanel("Tail Logo Plot", whole_plot_page_ui("tail_logo")),
        tabPanel("Post-Transcriptional Tailing", whole_plot_page_ui("pt_tail")),
        tabPanel("Test", plotOutput("test"))
    )
)

