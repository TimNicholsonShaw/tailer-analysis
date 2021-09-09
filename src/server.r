server <- function(input, output, session){

  #Put file dataframe into a reactive value
  values <- reactiveValues(df=NULL) 

  #populates and reset reactive dataframe when uploading files
  observeEvent(input$tail_files, {
    values$df <- input$tail_files
    values$df$group <- "" 
  })

  #Captures any edits made to the dataframe
  observeEvent(input$sample_table_cell_edit, {
    info <- input$sample_table_cell_edit
    values$df[info$row, 5] <- info$value
  })
  
  #Outputs the tail_File table with just the file name and grouping
  output$sample_table <- renderDT({
    req(input$tail_files)
    select(values$df, 1:5)
  }, editable='cell', rownames=FALSE)

  df <- reactive({
    req(input$make_df_button)
    dfBuilder(values$df$datapath, values$df$group)
  })

  output$test <- renderDT({
    req(input$make_df_button)
    head(isolate(df()))
  })

  ########## Candidate finder backend ##############
  cans <- reactive({
    req(input$find_cans_button)
    discover_candidates(isolate(df()), min=isolate(input$min_cans))
  })
  output$candidates <- renderDT({
      req(input$find_cans_button)
      isolate(cans())
  })

}