
#Pieces
ui_files <- fluidPage(
  #Input files and give them experimental groupings
  fileInput("tail_files", "Upload Tail Files", multiple=TRUE, accept=c(".csv")),
  DTOutput("sample_table"),
  actionButton("make_df_button", "Format Data"),
  DTOutput("df_preview") %>% withSpinner(color="#0dc5c1")
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
        mainPanel("Cumulative Plot",
            plotOutput("cum_plot") %>% withSpinner(color="#0dc5c1")
        )
    )
)

ui_tail_bar_graph <- fluidPage(
    sidebarLayout(
        sidebarPanel( "Options",
            textInput("tail_bar_gene_name",
                label="Gene Name"),
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
        mainPanel("Tail Bar Graph",
            plotOutput("tail_bar") %>% withSpinner(color="#0dc5c1")
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
        tabPanel("Cumulative Plot", ui_cumulative_tail_plot),
        tabPanel("Tail Bar Graph", ui_tail_bar_graph),
        tabPanel("Graph3"),
        tabPanel("Graph4")
    )
)