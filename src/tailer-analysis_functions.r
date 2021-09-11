library(tidyverse)
library(ggseqlogo)
library(cowplot)
library(progress)
library(shiny)
library(shinycssloaders)



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

cumulativeTailPlotter <- function(df, gene, start=-10, stop=10, gimme=FALSE, show_legend=TRUE, ymin=0, ymax=1, dots=FALSE, linecolors=""){
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

  if (gimme) { # for debugging
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

  # Pretty segments
  	geom_segment(aes(x=stop,xend=stop,
      y=ymin, 
      yend=ymax), 
      size=0.5) +

  
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

tail_bar_grapher <- function(df, gene, start=-10, stop=10, gimme=F, ymin=0, ymax=1, AUCGcolors="", show_legend=TRUE, dots=FALSE) {
  df <- filter(df, Gene_Name==gene) # Take one gene

  #Pre-populate dataframe
  out <- data.frame(Pos="", reads_at_end="", reads_start_A="", reads_start_C="", reads_start_G="", reads_start_T="", Sample="", Condition="")

  for (sample in unique(df$Sample)){
    #save total number of reads for that gene in that sample
    total <- sum(filter(df, Sample==sample)$Count)

    for (i in c(start:stop)){ #loop through region of interest
      #save total number of reads that mapped to that end
      reads_at_end <- sum(filter(df, Sample==sample, Three_End==i)$Count)
      
      #save how many reads have a tail that starts with each nucleotide
      reads_start_A <- sum(filter(df, Sample==sample, Three_End==i, startsWith(Tail_Sequence, "A"))$Count)
      reads_start_C <- sum(filter(df, Sample==sample, Three_End==i, startsWith(Tail_Sequence, "C"))$Count)
      reads_start_G <- sum(filter(df, Sample==sample, Three_End==i, startsWith(Tail_Sequence, "G"))$Count)
      reads_start_T <- sum(filter(df, Sample==sample, Three_End==i, startsWith(Tail_Sequence, "T"))$Count)
      
      #add to out df
      out <- rbind(out, c(i, reads_at_end/total, reads_start_A/total, reads_start_C/total, 
                          reads_start_G/total, reads_start_T/total, 
                          sample, filter(df, Sample==sample)$Grouping[1]))

      # Convert columns to numeric types
      out$reads_at_end = as.numeric(out$reads_at_end)
      out$reads_start_A = as.numeric(out$reads_start_A)
      out$reads_start_C = as.numeric(out$reads_start_C)
      out$reads_start_G = as.numeric(out$reads_start_G)
      out$reads_start_T = as.numeric(out$reads_start_T)
      out$Pos = as.numeric(out$Pos)
    }
  
  }
  if (gimme){ # for debugging
    return(out)
  }

    #################### Data Pre-processing #######################
  plt <- out[-1,] %>%
    dplyr::group_by(Condition, Pos) %>%
    dplyr::summarise(reads_at_end=mean(reads_at_end),
                     reads_start_A=mean(reads_start_A),
                     reads_start_C=mean(reads_start_C),
                     reads_start_G=mean(reads_start_G),
                     reads_start_T=mean(reads_start_T)) %>%

    pivot_longer(cols=c(reads_start_A, reads_start_C, reads_start_G, reads_start_T, reads_at_end))

    plt$name <- factor(plt$name, levels=c('reads_start_A', 'reads_start_C', 'reads_start_G', 'reads_start_T', 'reads_at_end'))
    plt$name <- recode_factor(plt$name, reads_start_A="A",
                              reads_start_C="C",
                              reads_start_G="G",
                              reads_start_T="T",
                              reads_at_end="No Tail")
    
  ######################### Plotting ############################
  plt <- plt %>%

    # main bar graph
    ggplot(aes(x=Pos, y=value, color=name, fill=name)) +
    geom_bar(stat='identity') +
    facet_grid(rows=vars(Condition), switch="y") +

    # themeing
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
                                          size=0.3), 
          strip.text = element_text(size=16, family="Helvetica",
                                    color="black"), 
          strip.text.x = element_text(margin = margin(0.05, 0,
                                                      0.05, 0, "cm")),
          strip.text.y = element_text(margin = 
                                        margin(0, 0.08, 0, 0.08, "cm"))
    ) +
  
  # Axes
  scale_x_continuous(name="Position", 
		breaks=seq(start, stop, 2),
		expand=c(0,0),
		limits=c(start,stop)) +
  scale_y_continuous(name="Fraction", 
		breaks=seq(0.1, 1, 0.2), 
		expand=c(0,0), 
		limits=c(ymin, ymax),
		position = "right"
		) +
  
  # Some segments
  geom_segment(aes(x=stop,xend=stop,
                     y=ymin, 
                     yend=ymax), 
                 size=0.3) +
    geom_segment(aes(x=start,xend=start,
                     y=ymin, 
                     yend=ymax), 
                 size=0.3) +
    geom_segment(aes(x=start,xend=stop,
                     y=ymin, 
                     yend=ymin), 
                 size=0.3) +
    geom_segment(aes(x=start,xend=stop,
                     y=ymax, 
                     yend=ymax), 
                 size=0.3) +

  # Custom colors
  # scale_fill_manual(values=AUCGcolors)


  # options
  if (show_legend == FALSE){
    plt <- plt + theme(legend.position = "none")
  }

  if (dots) { #doesn't work yet
    #plt <- plt + geom_point(data=out[-1,], aes(x=Pos, y=reads_at_end))
  }

  return(plt) 
}

tail_logo_grapher <- function(df, gene, xmin=1, xmax=10, ymin=0, ymax=1, gimme=F) {

  # So percentages can be entered as parameters as well as fractions
  if (ymax>1){
    ymax = ymax/100
  }
  if (ymin>1){
    ymin = ymin/100
  }

  ############ data pre-processing #################
  df <- filter(df, Gene_Name==gene) # Take one gene

  out <- data.frame(Pos="", Sample="", Condition="", Nuc="", Frequency="")

  for (sample in unique(df$Sample)){
    #save total number of reads for that gene in that sample
    total <- sum(filter(df, Sample==sample)$Count)

    for (i in c(xmin:xmax)) {#loop through tail lengths of interest
      for (nuc in c("A", "C", "G", "T")) {
        # Find frequency of each nucleotide at each position
        frequency <- sum(filter(df, Sample==sample, substr(Tail_Sequence, i, i) == nuc)$Count)/total

        # Add to out df
        out <- rbind(out, c(i, sample, filter(df, Sample==sample)$Grouping[1], nuc, frequency))
      }
    }
  }
  # futz with columns
  out=out[-1,]
  out$Pos <- as.numeric(out$Pos)
  out$Frequency <- as.numeric(out$Frequency)
  
  if(gimme){ # debugging
    return(out)
  }


  # calculate out nucleotide frequencies
  out <-out %>%
    group_by(Condition, Pos, Nuc) %>% 
    summarise(freq_avg = mean(Frequency)) %>%
    ungroup()
  plot_list <- list() # to be added to for cowplot

  conditions <- c(unique(out$Condition))

  for (i in 1:length(conditions)){
    plot_df <- out %>% 
                filter(Condition==conditions[i]) %>% # Look at singular condition
                pivot_wider(id_cols=Nuc, names_from=Pos, values_from=freq_avg) %>% 
                {. ->> nucs} %>% # Save intermediate to get nucleotide names
                select(-Nuc) %>% # Remove Nuc column
                as.matrix(nrow=4, ncol=ncol(nucs)-1) # ggseqlogo needs a probability matrix
    rownames(plot_df) <- nucs$Nuc # Names from Nuc intermediate makes sure we don't get them out of order

    # Add plot to plotlist
    plot_list[[i]] <- ggplot() + 
      geom_logo(data=plot_df, method='custom', font='roboto_bold') + 
      coord_cartesian() +

      # Themeing
      theme_classic() +
      theme(
        axis.text=element_text(family="Helvetica", size=12),
        axis.text.x=element_text(hjust=0.5),
        axis.line=element_line(size=0.3),
        axis.ticks.length=unit(0.1, "cm"),
        axis.ticks=element_line(size=0.3),
        axis.title=element_text(family="Helvetica", 
          face="bold", size=16),
        legend.text=element_text(family="Helvetica", size=8),
        legend.title=element_blank(), 
        legend.key.height=unit(0.3, "cm"),
        legend.key.width=unit(0.3, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.background=element_rect(fill="white", linetype=1,
          size=0.3, color="black")
	    ) + 

      # Axes
      scale_x_continuous(
        name="Tail length",
        breaks = seq(1, xmax-xmin+1, 1),
        labels = seq(xmin, xmax, 1),
        expand=c(0,0),
        limits=c(0.5, xmax-xmin+1.5)
	    ) +
      scale_y_continuous(
        name="% Nucleotide",
        breaks=seq(ymin, ymax, 0.1),
        labels=seq(ymin*100, ymax*100, 10),
        expand=c(0,0),
        limits=c(0, ymax), 
        position = "right"
	    ) +
      ggtitle(conditions[i])
  }

  return(plot_grid(plotlist = plot_list, ncol=1)) # Use cowplot to place them in a grid
}

tail_pt_nuc_grapher <- function(df, gene, gimme=F, ymin=0, ymax=1, pdisplay=F){
  df <- filter(df, Gene_Name==gene) # Take one gene
  out <- data.frame(Sample="", Condition="", Nuc="", Frequency="")

  for (sample in unique(df$Sample)){
    #save total number of reads for that gene in that sample
    total <- sum(filter(df, Sample==sample)$Count)

    for (nuc in c("A", "C", "G", "T")) {
      # Find frequency of each nucleotide at each position
      frequency <- sum(filter(df, Sample==sample, substr(Tail_Sequence, 1, 1) == nuc)$Count)/total
      # Add to out df
      out <- rbind(out, c(sample, filter(df, Sample==sample)$Grouping[1], nuc, frequency))
      }
    
  }
  # Change to numeric types
  out=out[-1,]
  out$Frequency <- as.numeric(out$Frequency)

  out_summed <-out %>%
    group_by(Condition, Nuc) %>% 
    summarise(freq_avg = mean(Frequency), sd= sd(Frequency), n=n(), se=sd/sqrt(n)) %>%
    ungroup()

  if(gimme){ # debugging
    return(out)
  }

  # t-testing
  pvals = c()
  conditions = unique(out$Condition)
  if (pdisplay){
    if (length(conditions) != 2){
      return("Can't do the stats on a different number than 2 conditions because statistics is hard") # Better error handling here
    }
    for (nuc in unique(out$Nuc)){
      temp_val <- t.test(filter(out, Nuc==nuc, Condition==conditions[1])$Frequency, filter(out, Nuc==nuc, Condition==conditions[2])$Frequency)$p.val
      pvals <- append(pvals, temp_val)
    }
  }

  ################ Plotting ##########################
  plt <- ggplot(out) + 
    coord_cartesian() +
    theme_classic() +

    # Individual dots
    geom_jitter(aes(y=Frequency, color=Condition, x=Nuc)) +

    # Lines for mean
    geom_linerange(data=out_summed, aes(y=freq_avg, 
      xmax=c(1.4, 2.4, 3.4, 4.4, 1.4, 2.4, 3.4, 4.4), 
      xmin=c(0.6, 1.6, 2.6, 3.6, 0.6, 1.6, 2.6, 3.6), 
      color=Condition)) +


    geom_errorbar(data=out_summed, 
      aes(x=Nuc, ymin=freq_avg-se, ymax=freq_avg+se, color=Condition, width=0.2)) +

    # Themeing
    theme(strip.background.x = element_rect(
      fill="grey",
      linetype=1, size=0.3),
      strip.background.y = element_rect(
      fill="grey",
      linetype=1, size=0.3), 
      strip.text = element_text(size=16, family="Helvetica",
        color="black"), 
      strip.text.x = element_text(margin = 
        margin(0.05, 0, 0.05, 0, "cm")),
      strip.text.y = element_text(margin = 
        margin(0, 0.08, 0, 0.08, "cm"))
    ) +
		
	theme(axis.text=element_text(family="Helvetica", size = 12),
		axis.text.x=element_text(angle = 45, 
			vjust = 0.8, hjust = 0.8),
		axis.line=element_line(size=0.3),
		axis.ticks.length=unit(0.1, "cm"),
		axis.ticks=element_line(size=0.3),
		axis.title=element_text(family="Helvetica", face="bold",
		size=16),
		legend.text=element_text(family="Helvetica", size=8),
		legend.title=element_blank(),
		legend.key.height=unit(0.1, "cm"),
		legend.key.width=unit(0.4, "cm")
  ) +

  #Axes
  scale_y_continuous( 
		limits=c(ymin-0.01, ymax), 
		position = "right")


  ######## options ############

  if (pdisplay){
    pvals <- lapply(pvals, round, digits=4)
    pvals <- paste0("p=", pvals)
    plt <- plt + geom_text(data=out_summed, aes(x=Nuc, y=ymax, label=append(pvals, c("", "", "", "")))) # this is dumb, but it works
  }


  return(plt)
}

discover_candidates <- function(df, min=1){

  conditions <- unique(df$Grouping)
  if (length(conditions) != 2){ # Statistics is hard
    return("Candidate discovery doesn't support this number of conditions :(, try with just 2")
  }
  genes <- unique(df$Gene_Name)

  pb <- progress_bar$new(total=length(genes),
          format = "Finding candidates [:bar] :percent | :elapsed | Estimated eta: :eta",
          clear=FALSE)

  out <- data.frame(Gene=NA, pval_three_end=NA, pval_tail_length=NA)
  
 
  withProgress( message="Finding Candidates", value=0, { # Progress bar wrapper for shiny
    for (i in 1:length(genes)) {
      pb$tick()
      incProgress(1/length(genes))
      # toss if it doesn't meet the minimum number of reads
      if (genes[i] == "None") next
      if (nrow(filter(df, Grouping==conditions[1],Gene_Name==genes[i])) < min) next
      if (nrow(filter(df, Grouping==conditions[2],Gene_Name==genes[i])) < min) next


      statistic_three_end <- tryCatch({
        suppressWarnings(ks.test(filter(df, Grouping==conditions[1],Gene_Name==genes[i])$Three_End,
                          filter(df, Grouping==conditions[2],Gene_Name==genes[i])$Three_End))
                          },
        error = function(e) e
      )
      
      #Erros to NA
      if(inherits(statistic_three_end, "error")) statistic_three_end$p.value <- NA


      statistic_tail_length <- tryCatch({
        suppressWarnings(ks.test(filter(df, Grouping==conditions[1],Gene_Name==genes[i])$Tail_Length,
                          filter(df, Grouping==conditions[2],Gene_Name==genes[i])$Tail_Length))
                          },
        error = function(e) e
      )
      
      # Errors to NA
      if(inherits(statistic_tail_length, "error")) statistic_tail_length$p.value <- NA
      

      out <- rbind(out, c(genes[i], statistic_three_end$p.value, statistic_tail_length$p.value))
    }
  }) 

  out<-out[-1,] # Remove crap line
  # Change to numeric types
  out$pval_three_end <- as.numeric(out$pval_three_end)
  out$pval_tail_length <- as.numeric(out$pval_tail_length)

  out <- out[order(out$pval_three_end),]
  return (out)
  

}

 