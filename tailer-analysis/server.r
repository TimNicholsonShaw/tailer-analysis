# Modular Pieces
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
        analysis_max=input$analysis_max,
        pdisplay=input$pdisplay,
        mature_end=input$mature_end
      ))
    }
  )
}

plot_server <- function(id, plt, plottype=""){
  moduleServer(
    id,
    function(input, output, session){
      flag <- F
      if (flag==F){
        output$plot <- renderPlot({NULL})
        flag <-T
      }
      observeEvent(input$make_plot, {
        output$plot <- renderPlot({isolate(plt()$plot)})
      })
      
      output$download_plot <- downloadHandler(
        filename=paste(input$gene_name, plottype, ".tiff", sep=""),
        content=function(file){
          tiff(file, height=input$height, width=input$width, res=600)
          print(plt())
          dev.off()
        }
      )
      output$download_data <- downloadHandler(
        filename=paste(input$gene_name, "_data", "csv", sep="."),
        content=function(file){
          write.csv(plt()$data, file, row.names=F)
        }
      )
    }
  )
}


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
    isolate(dfBuilder(values$df$datapath, values$df$group))
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

  # Set sample order
  output$sample_order <- renderUI({
    orderInput("sample_order", "Drag to reorder samples", unique(df()$Grouping))
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
  }, rownames=FALSE)
  
  output$download_cans<- downloadHandler(
    filename="candidates.csv",
    content=function(file){
      write.csv(cans(), file, row.names=F)
    }
  )
################### Cumulative Plotter #####################
  cum_plot_options <- reactive({options_server("cum_plot")})
  
  
  cum_plot <- reactive({
    cumulativeTailPlotter(df(), 
                          cum_plot_options()$gene_name,
                          start=cum_plot_options()$x_min,
                          stop=cum_plot_options()$x_max,
                          ymin=cum_plot_options()$y_min,
                          ymax=cum_plot_options()$y_max,
                          dots=cum_plot_options()$dots,
                          multi_locus=cum_plot_options()$multi_loc,
                          analysis_min=cum_plot_options()$analysis_min,
                          analysis_max=cum_plot_options()$analysis_max,
                          mature_end=cum_plot_options()$mature_end,
                          order=input$sample_order
                            
    )
  })
  plot_server("cum_plot", cum_plot, plottype="cumulative_plot")
  
################## Tail Bar Graph ########################
  tail_bar_options <- reactive({options_server("tail_bar")})
    
  tail_bar_plot <- reactive({
    print(
      tail_bar_grapher(df(),
                       tail_bar_options()$gene_name,
                       start=tail_bar_options()$x_min,
                       stop=tail_bar_options()$x_max,
                       ymin=tail_bar_options()$y_min,
                       ymax=tail_bar_options()$y_max,
                       #dots=tail_bar_options()$dots,
                       multi_locus=tail_bar_options()$multi_loc,
                       analysis_min=tail_bar_options()$analysis_min,
                       analysis_max=tail_bar_options()$analysis_max,
                       mature_end=tail_bar_options()$mature_end,
                       order=input$sample_order
                       )
    )
  })
  
  plot_server("tail_bar", tail_bar_plot, plottype="tail_bar")

################### Tail Logo Graph #####################
  tail_logo_options <- reactive({options_server("tail_logo")})
  
  tail_logo_plot <- reactive({
    print(
      tail_logo_grapher(df(),
                        tail_logo_options()$gene_name,
                        xmin=1,
                        xmax=tail_logo_options()$x_max,
                        ymin=tail_logo_options()$y_min,
                        ymax=tail_logo_options()$y_max,
                        multi_locus = tail_logo_options()$multi_loc,
                        order=input$sample_order
                        )
    )
  })
  plot_server("tail_logo", tail_logo_plot, plottype="tail_logo")

####################### PT Nuc Graph ########################
  pt_tail_options <- reactive({options_server("pt_tail")})
  
  pt_tail_plot <- reactive({
    print(
      tail_pt_nuc_grapher(df(),
                          pt_tail_options()$gene_name,
                          ymin=pt_tail_options()$y_min,
                          ymax=pt_tail_options()$y_max,
                          multi_locus=pt_tail_options()$multi_loc,
                          pdisplay=pt_tail_options()$pdisplay,
                          order=input$sample_order
                          )
    )
  })
  
  plot_server("pt_tail", pt_tail_plot, plottype="pt_tail")
  
#################### Statistics Page #######################
  output$con1 <- renderUI({
    selectInput("con1", "Condition 1", unique(df()$Grouping))
  })
  output$con2 <- renderUI({
    selectInput("con2", "Condition 2", unique(df()$Grouping))
  })
  
  stats <- reactive({stat_matrix_maker(df(), input$stats_gene_name, input$con1, input$con2)})
  output$end_position_tab <- renderTable(stats()$p.mat_end_position, rownames=T, digits=-6)
  output$end_position_text <- renderText(stats()$p.end_position_total)
  output$tail_len_tab <- renderTable(stats()$p.mat_tail_len, rownames=T, digits=-6)
  output$tail_len_text <- renderText(stats()$p.tail_len_total)
  

}


