# Modular Pieces


plot_page_server <- function(id, dataframe, plot_fun, ...){
  moduleServer(
    id,
    function(input, output, session){
      output$plot <- renderPlot({NULL})
      plt <- reactive({
        plot_fun(
          dataframe,
          isolate(input$gene_name),
          x_min=isolate(input$x_min),
          x_max=isolate(input$x_max),
          ymin=isolate(input$y_min),
          ymax=isolate(input$y_max)
        )
      })
      observeEvent(input$make_plot,{
        output$plot <- renderPlot(plt())
      })
    }
  )
}
  

options_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      return(list(
        gene_name=input$gene_name,
        multi_loc=input$multi_loc,
        x_min=input$x_min,
        x_max=input$x_max,
        y_min=input$y_min,
        y_max=input$y_max,
        dots=input$dots,
        analysis_min=input$analysis_min,
        analysis_max=input$analysis_max
      ))
    }
  )
}

plot_server <- function(id, plt){
  moduleServer(
    id,
    function(input, output, session){
      observeEvent(input$make_plot, {
        output$plot <- renderPlot({isolate(plt())})
        
      })

    }
  )
}


server <- function(input, output, session){
  
################### Cumulative Plotter #####################
  #plot_page_server("cum_plot", df(), cumulativeTailPlotter, ymax=0.5)
  cum_plot_options <- reactive({options_server("cum_plot")})
  
  
  cum_plot <- reactive({
    print(cumulativeTailPlotter(df(), 
                          cum_plot_options()$gene_name,
                          start=cum_plot_options()$x_min,
                          stop=cum_plot_options()$x_max,
                          ymin=cum_plot_options()$y_min,
                          ymax=cum_plot_options()$y_max,
                          dots=cum_plot_options()$dots,
                          multi_locus=cum_plot_options()$multi_loc
                        ))
  })
  plot_server("cum_plot", cum_plot)
  
  output$test <- renderPlot(cum_plot())

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
  output$df_preview <- renderText({
    req(input$make_df_button)
    out=""
    for (group in unique(df()$Grouping)){
      out <- paste0(out,group,
                    " Condition:\n",
                    length(unique(filter(df(), Grouping==group)$Sample)), 
                    " Sample(s)\n")
      }
    out
  })

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


  
################## Tail Bar Graph ########################
  output$tail_bar <- renderPlot({NULL})
  observeEvent(input$make_tail_bar, {
    output$tail_bar <- renderPlot({
      tail_bar_grapher(
        df(),
        isolate(input$tail_bar_gene_name),
        start=isolate(input$tail_bar_x_min),
        stop=isolate(input$tail_bar_x_max),
        ymin=isolate(input$tail_bar_y_min),
        ymax=isolate(input$tail_bar_y_max),
        multi_locus=isolate(input$tail_bar_multiloc)
      )
    })
  })

################### Tail Logo Graph #####################
  output$tail_logo <- renderPlot({NULL})
  observeEvent(input$make_tail_logo, {
    output$tail_logo <- renderPlot({
      tail_logo_grapher(
        df(),
        isolate(input$tail_logo_gene_name),
        xmin=isolate(input$tail_logo_xmin),
        xmax=isolate(input$tail_logo_xmax),
        ymin=isolate(input$tail_logo_ymin),
        ymax=isolate(input$tail_logo_ymax),
        multi_locus=isolate(input$tail_logo_multiloc)
      )

    })
  })

####################### PT Nuc Graph ########################
  output$pt_graph <- renderPlot({NULL})
  observeEvent(input$make_pt_graph, {
    output$pt_graph <- renderPlot({
      tail_pt_nuc_grapher(
        df(),
        isolate(input$pt_gene_name),
        ymin=isolate(input$pt_ymin),
        ymax=isolate(input$pt_ymax),
        pdisplay=isolate(input$pt_pdisplay),
        multi_locus=isolate(input$pt_multiloc)
      )
    })
  })
  opt <- reactive({test_options_server("test")})
  output$test_text <- renderText(opt())
  
}  


