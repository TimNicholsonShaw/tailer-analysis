library(shiny)
library(DT)
library(tidyverse)
source("./tailer-analysis_functions.r")

options(shiny.maxRequestSize = 1000*1024^2) # raises file limit 

source("./ui.r")
source("./server.r")

shinyApp(ui, server)