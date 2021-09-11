library(shiny)
library(DT)
library(tidyverse)
source("./scripts/tailer-analysis_functions.r")

options(shiny.maxRequestSize = 1000*1024^2) # raises file limit 
