server <- function(input, output, session){

############### File input backend #######################
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
  #Allows metadata aenrty
  output$sample_table <- renderDT({
    req(input$tail_files)
    select(values$df, 1,5)
    }, rownames=F, editable=list(target = "cell", disable = list(columns = c(0)))
  )


  # Reactive dataframe to hold the combined data formatted by dfBuilder
  df <- reactive({
    req(input$make_df_button)
    dfBuilder(values$df$datapath, values$df$group)
  })

  # Preview the data frame
  output$df_preview <- renderDT({
    req(input$make_df_button)
    head(isolate(df()))
    }, options = list(scrollX = TRUE)
  )

########## Candidate finder backend ##############

  #  reactive to hold candidates in
  cans <- reactive({
    req(input$find_cans_button)
    discover_candidates(isolate(df()), min=isolate(input$min_cans))
  })
  # Output to dataframe 
  output$candidates <- renderDT({
      req(input$find_cans_button)
      isolate(cans())
  })

################### Cumulative Plotter #####################

  observeEvent(input$make_cum_plot, {
    output$cum_plot <- renderPlot({
      cumulativeTailPlotter(
        df(),
        isolate(input$cum_plot_gene_name),
        start=isolate(input$cum_plot_x_min),
        stop=isolate(input$cum_plot_x_max),
        ymin=isolate(input$cum_plot_y_min),
        ymax=isolate(input$cum_plot_y_max),
        dots=isolate(input$cum_plot_dots)
      )
    })
  })
################## Tail Bar Graph ########################

  observeEvent(input$make_tail_bar, {
    output$tail_bar <- renderPlot({
      tail_bar_grapher(
        df(),
        isolate(input$tail_bar_gene_name),
        start=isolate(input$tail_bar_x_min),
        stop=isolate(input$tail_bar_x_max),
        ymin=isolate(input$tail_bar_y_min),
        ymax=isolate(input$tail_bar_y_max),
      )
    })
  })

################### Tail Logo Graph #####################

  observeEvent(input$make_tail_logo)

}