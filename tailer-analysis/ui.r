# Modular Pieces
download_button_ui <- function(id){
  ns <- NS(id)
  tagList(
    downloadButton(ns("plot_download"), "Download Plot"),
    numericInput(ns("height"), "height", value=1000),
    numericInput(ns("width"), "width", value=1000)
  )
}

test_options_ui <- function(id, dots=F){
  ns <- NS(id)
  tags<- tagList(
    textInput(ns("gene_name"), "Gene Name"),
    checkboxInput(ns("multi_loc"), "Multi-Locus Gene"),
    numericInput(ns("x_min"), "Minimum Position", value=-10),
    numericInput(ns("x_max"), "Maximum Position", value=10),
    numericInput(ns("y_min"), "Minimum Fraction", value=0),
    numericInput(ns("y_max"), "Maximum Fraction", value=1),
    )
  if (dots==T) {tags <- tagAppendChild(tags, checkboxInput("dots", "Dots?"))}
  tags
}

test_plot_ui <- function(id){
  ns<-(NS(id))
  tagList(
    actionButton(ns("make_plot"), "Make Plot"),
    plotOutput(ns("plot")),
    downloadButton(ns("download_plot"), "Download Plot"),
    fluidRow(
      column(6,numericInput(ns("height"), "Plot Height", value=400)),
      column(6, numericInput(ns("width"), "Plot Width", value=400))
    )
  )
}

plot_page_ui <- function(id){
  ns <-NS(id)
  tagList(
    sidebarLayout(
    sidebarPanel(
      textInput(ns("gene_name"), "Gene Name"),
      numericInput(ns("x_min"), "Minimum Position", value=-10),
      numericInput(ns("x_max"), "Maximum Position", value=10),
      numericInput(ns("y_min"), "Minimum Fraction", value=0),
      numericInput(ns("y_max"), "Maximum Fraction", value=1),
      actionButton(ns("make_plot"), "Make Plot")
    ),
    mainPanel(
      plotOutput(ns("plot"))%>% withSpinner(color="#0dc5c1"),
      downloadButton(ns("download_plot"), "Download Plot"),
      fluidRow(
        column(6,numericInput(ns("height"), "Plot Height", value=400)),
        column(6, numericInput(ns("width"), "Plot Width", value=400))
        )
      )
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

ui_cumulative_tail_plot <- fluidPage(
    sidebarLayout(
        sidebarPanel( "Options",
            textInput("cum_plot_gene_name",
                label="Gene Name"),
            checkboxInput("cum_plot_multiloc",
                label="Multi Locus Gene",
                value=F),
            numericInput("cum_plot_x_min", 
                label="Minimum Position",
                value=-10
            ),
            numericInput("cum_plot_x_max",
                label="Maximum Position",
                value=10
            ),
            numericInput("cum_plot_y_min",
                label="Cumulative Fraction Minimum",
                value=0,
                step=0.1
            ),
            numericInput("cum_plot_y_max",
                label="Cumulative Fraction Maximum",
                value=1,
                step=0.1
            ), 
            checkboxInput("cum_plot_dots",
                label="Individual Experiment Dots",
                value=F
            ),
            actionButton("make_cum_plot",
                label="Make Plot"
            )
        ),
        mainPanel(
            plotOutput("cum_plot") %>% withSpinner(color="#0dc5c1"),
            download_button_ui("cum_plot_download")
        )
    )
)

ui_tail_bar_graph <- fluidPage(
    sidebarLayout(
        sidebarPanel( "Options",
            textInput("tail_bar_gene_name",
                label="Gene Name"),
            checkboxInput("tail_bar_multiloc",
                label="Multi Locus Gene",
                value=F),
            numericInput("tail_bar_x_min", 
                label="Minimum Position",
                value=-10
            ),
            numericInput("tail_bar_x_max",
                label="Maximum Position",
                value=10
            ),
            numericInput("tail_bar_y_min",
                label="Fraction Minimum",
                value=0,
                step=0.1
            ),
            numericInput("tail_bar_y_max",
                label="Fraction Maximum",
                value=1,
                step=0.1
            ), 

            actionButton("make_tail_bar",
                label="Make Plot"
            )
        ),
        mainPanel(
            plotOutput("tail_bar") %>% withSpinner(color="#0dc5c1")
        )
    )
)

ui_tail_logo_plot <- fluidPage(
    sidebarLayout(
        sidebarPanel("Options",
          textInput("tail_logo_gene_name",
            label="Gene Name"
          ),
        checkboxInput("tail_logo_multiloc",
                label="Multi Locus Gene",
                value=F),
          numericInput("tail_logo_xmin",
            label="Tail Position Start",
            min=1,
            value=1
          ),
          numericInput("tail_logo_xmax",
            label="Tail Position End",
            min=1,
            value=10
          ),
          numericInput("tail_logo_ymin",
            label="Fraction Minimum",
            min=0,
            max=1,
            step=0.1,
            value=0
          ),
          numericInput("tail_logo_ymax",
            label="Fraction Maximum",
            min=0,
            max=1,
            step=0.1,
            value=1
          ),
          actionButton("make_tail_logo",
            label="Make Plot"
          )

        ),
        mainPanel(
            plotOutput("tail_logo") %>% withSpinner(color="#0dc5c1")
        )
    )
)

ui_pt_graph <- fluidPage(
    sidebarLayout(
        sidebarPanel("Options",
            textInput("pt_gene_name",
                label="Gene Name"
            ),
            checkboxInput("pt_multiloc",
                label="Multi Locus Gene",
                value=F),
            numericInput("pt_ymin",
                label="Fraction Minimum",
                min=0,
                max=1,
                value=0,
                step=0.1
            ),
            numericInput("pt_ymax",
                label="Fraction Maximum",
                min=0,
                max=1,
                value=1,
                step=0.1
            ),
            checkboxInput("pt_pdisplay",
                label="Display P Values",
                value=F
            ),
            actionButton("make_pt_graph",
                label="Make Plot"
            )
        ),
        mainPanel(
            plotOutput("pt_graph") %>% withSpinner(color="#0dc5c1")
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
        tabPanel("Cumulative Plot", plot_page_ui("cum_plot")),
        tabPanel("Tail Bar Graph", ui_tail_bar_graph),
        tabPanel("Tail Logo Plot", ui_tail_logo_plot),
        tabPanel("Post-Transcriptional Tailing", ui_pt_graph),
        tabPanel("test", fluidPage(sidebarLayout(
          sidebarPanel(
            test_options_ui("test", dots=T)
          ),
          mainPanel(test_plot_ui("test_plot"),
                    verbatimTextOutput("debug"))
          
          #mainPanel(plotOutput("test_plot"))
        )))
    )
)

