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

cumulativeTailPlotter <- function(df, gene, start, stop, gimme=FALSE, show_legend=TRUE, ymin=0, ymax=1, dots=FALSE, linecolors=""){
  # Only interested in a single gene
  # Maybe add some flexibility for multi-mappers
  df <- filter(df, Gene_Name==gene) 

  # Pre-populate dataframe
  out = data.frame(Pos="", Total_Percentage="", Three_End_Percentage="",  Sample="", Condition="")

  for (sample in unique(df$Sample)){ 
    total <- sum(filter(df, Sample==sample)$Count) # sum all reads for this gene

    # Loop through positions, find cumulative percentage
    for (x in c(start:stop)){ 
      # Three_End + Tail_Length gives total tail
      total_cum_sum <- sum(filter(df, Sample==sample, Three_End + Tail_Length <= x)$Count)/total
      three_end_cum_sum <- sum(filter(df, Sample==sample, Three_End <= x)$Count)/total
      #Add to dataframe
      out = rbind(out, c(x, total_cum_sum, three_end_cum_sum, sample, filter(df, Sample==sample)$Grouping[1]))
    }
  }
  #Convert columns to numeric type
  out$Total_Percentage <- as.numeric(out$Total_Percentage)
  out$Three_End_Percentage <- as.numeric(out$Three_End_Percentage)
  out$Pos <- as.numeric(out$Pos)
  out=out[-1,] # drop that weird first row

  if (gimme) {
    return(out)
  }

  ######################Plot Output#######################
  plt <- out %>%
    dplyr::group_by(Condition, Pos) %>%
    dplyr::summarise(
      total_avg = mean(Total_Percentage),
      three_end_avg = mean(Three_End_Percentage)
    ) %>% 
  ggplot(aes(x=Pos, y=total_avg, color=Condition)) + 

  # Add gray rectangle indicating mature length
  geom_rect(aes(xmin=start, 
    xmax=0,
		ymin=0, 
    ymax=1), 
		linetype=0, 
		fill="grey90", 
		alpha=0.05) +
  
  # Step Plots for full tail and 3' end mapping only
  geom_step() + 
  geom_step(aes(x=Pos, y=three_end_avg), linetype=3) +

  # Rectangle to denote post-transcriptional tailing
  geom_rect(aes(xmin=Pos, xmax=Pos+0.99,
		ymin=total_avg, ymax=three_end_avg, 
		color=Condition, fill=Condition),
		linetype=0, 
		alpha=0.2,
    show.legend=FALSE) +

  # Axes
  scale_x_continuous(name="3' end position", 
		breaks=seq(start, stop, 2),
		expand=c(0,0),
		limits=c(start,stop)) +
  scale_y_continuous(name="Cumulative Percent", 
		breaks=seq(0.1, 1, 0.2), 
		expand=c(0,0), 
		limits=c(ymin, ymax),
		position = "right"
		) +  

  # Theming
  theme_classic() +
  theme(axis.text=element_text(family="Helvetica", size=12),
		axis.line=element_line(size=0.5),
		axis.ticks.length=unit(0.1, "cm"),
		axis.ticks=element_line(size=0.5),
		axis.title=element_text(family="Helvetica", face="bold",
		size=16),
		legend.text=element_text(family="Helvetica", size=8),
		legend.title=element_blank(),
		legend.key.height=unit(0.1, "cm"),
		legend.key.width=unit(0.4, "cm")) +

    theme(strip.background = element_rect(fill="grey", linetype=1,
		size=0.8), 
		strip.text = element_text(size=16, family="Helvetica",
		color="black"), 
		strip.text.x = element_text(margin = margin(0.05, 0,
			0.05, 0, "cm")),
		strip.text.y = element_text(margin = 
			margin(0, 0.08, 0, 0.08, "cm"))
		) +
  
  #Colors
  #scale_fill_manual(values=linecolors) +

  
  ######## Plot Options #########

  if (show_legend == FALSE){
    plt <- plt + theme(legend.position = "none")
  }

  if (dots) {
    plt <- plt +   geom_point(data=out, aes(x=Pos, y=Total_Percentage))
  }

  return(plt)

     
  }