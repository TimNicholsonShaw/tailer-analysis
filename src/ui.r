
#Pieces
ui_files <- fluidPage(
  #Input files and give them experimental groupings
  fileInput("tail_files", "Tail Files", multiple=TRUE),
  DTOutput("sample_table"),
  actionButton("make_df_button", "Format Data")
  #could add custom colors here
)

ui_candidate_finder <- fluidPage(
  DTOutput("test") %>% withSpinner(color="#0dc5c1"),
  numericInput("min_cans", "Minimum Reads", 10, min=10),
  actionButton("find_cans_button", "Find Candidates"),
  DTOutput("candidates") %>% withSpinner(color="#0dc5c1"),
)

# Actual UI declaration
ui <- fluidPage(
    titlePanel(
        "3' End Analysis"
    ),
    ui_files,
    ui_candidate_finder
)