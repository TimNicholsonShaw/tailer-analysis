library(tidyverse)

dfBuilder <- function(files, grouping){
  # There's gotta be a better way to do this, please tell me at timnicholsonshaw@gmail.com
  # Prep dataframe
  out <- data.frame(Sequence="",
                    Count="",
                    EnsID="",
                    Gene_Name="",
                    Three_End="",
                    Tail_Length="",
                    Tail_Sequence="",
                    Sample="",
                    Grouping="")
  
  for (i in c(1:length(files))){
    temp_df <- read_csv(files[i]) 
    temp_df$Sample <- files[i] # Add file name as sample name, maybe do prefix?
    temp_df$Grouping <- grouping[i] # Add sample metadata
    out <- rbind(out, temp_df)
  }
  out$Count <- as.numeric(out$Count)
  out$Three_End <- as.numeric(out$Three_End)
  out$Tail_Length <- as.numeric(out$Tail_Length)
  return (out[-1,])
}

cumulativeTailPlotter <- function(df, gene, start, stop){
  # Only interested in a single gene
  # Maybe add some flexibility for multi-mappers
  df <- filter(df, Gene_Name==gene) 


  # Pre-populate dataframe
  out = data.frame(Pos="", Percentage="",  Sample="", Condition="")

  for (sample in unique(df$Sample)){

    for (x in c(start:stop)){
      total <- sum(filter(df, Sample==sample)$Count)
      cum_sum <- sum(filter(df, Sample==sample, Three_End <= x)$Count)/total
      out = rbind(out, c(x, cum_sum, sample, filter(df, Sample==sample)$Grouping[1]))
    }
  }
  out$Percentage <- as.numeric(out$Percentage)
  out$Pos <- as.numeric(out$Pos)
  out=out[-1,]
  out %>%
  dplyr::group_by(Condition, Pos) %>%
  dplyr::summarise(
    average = mean(Percentage)
  ) %>%
  ggplot(aes(x=Pos, y=average, color=Condition)) + 
    geom_step() + 
    geom_point(data=out, aes(x=Pos, y=Percentage))
}