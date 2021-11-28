library(tidyverse)
library(ggseqlogo)
library(cowplot)
library(progress)
library(shiny)
library(shinycssloaders)
library(ggthemes)
library(lemon)
library(ggprism)
library(grid)
library(stringr)


#################### Graphing Functions ##########################
cumulativeTailPlotter <- function(df, gene, start=-10, stop=10, gimme=FALSE, show_legend=TRUE, ymin=0, ymax=1, dots=FALSE, multi_locus=F, analysis_min=-100, analysis_max=100, mature_end=0, order=""){
  # Only interested in a single gene
  # Maybe add some flexibility for multi-mappers
  if (grepl("|", gene, fixed=T)) multi_locus=T #if entering multiple genes, automatically does multilocus search
  if (multi_locus) df <- multimap_gene_subsetter(df, gene) #optional slower more exhaustivesubsetter
  else if(startsWith(gene, "ENS")) df <- filter(df, EnsID==gene) # handling for ENSIDs
  else df <- filter(df, Gene_Name==gene) 
  df <- filter(df, End_Position>=analysis_min, End_Position<=analysis_max) # throws out samples out of analysis range
  df$End_Position <- df$End_Position - mature_end # mature end correction

  if (nrow(df)==0) return(NULL) # Return a null if you didn't find any of that gene

  # Pre-populate dataframe
  out = data.frame(Pos="", Total_Percentage="", Three_End_Percentage="",  Sample="", Condition="")

  for (sample in unique(df$Sample)){ 
    total <- sum(filter(df, Sample==sample)$Count) # sum all reads for this gene

    # Loop through positions, find cumulative percentage
    for (x in c(start:stop)){ 
      # Three_End + Tail_Length gives total tail
      total_cum_sum <- sum(filter(df, Sample==sample, End_Position <= x)$Count)/total
      three_end_cum_sum <- sum(filter(df, Sample==sample, End_Position - Tail_Length <= x)$Count)/total
      #Add to dataframe
      out = rbind(out, c(x, total_cum_sum, three_end_cum_sum, sample, filter(df, Sample==sample)$Grouping[1]))
    }
  }
  #Convert columns to numeric type
  out$Total_Percentage <- as.numeric(out$Total_Percentage)
  out$Three_End_Percentage <- as.numeric(out$Three_End_Percentage)
  out$Pos <- as.numeric(out$Pos)
  if (order!="") out$Condition <- factor(out$Condition, levels=order)
  out=out[-1,] # drop that weird first row




  ######################Plot Output#######################
  plt <- out %>%
    dplyr::group_by(Condition, Pos) %>%
    dplyr::summarise(
      total_avg = mean(Total_Percentage),
      three_end_avg = mean(Three_End_Percentage)
    ) 

  plt$top_label=gene
  plt$left_label=""

  # save the data for output
  data_out<-plt

  plt<-plt %>%
  ggplot(aes(x=Pos, y=total_avg, color=Condition)) + 

  # Add gray rectangle indicating mature length
  geom_rect(aes(xmin=start, 
    xmax=0,
		ymin=0, 
    ymax=1), 
		linetype=0, 
		fill="grey90", 
		alpha=0.02) +
  
  # Step Plots for full tail and 3' end mapping only
  geom_step(size=.7) + 
  geom_step(aes(x=Pos, y=three_end_avg), linetype=3, size=.7) +

  # Rectangle to denote post-transcriptional tailing
  geom_rect(aes(xmin=Pos, xmax=Pos+0.99,
		ymin=total_avg, ymax=three_end_avg, 
		color=Condition, fill=Condition),
		linetype=0, 
		alpha=0.2,
    show.legend=FALSE) +

  # Axes
  scale_x_continuous(name="Position", 
    expand=c(0,0),
		limits=c(start,stop),
    guide = guide_prism_minor()) + # add minor ticks
  scale_y_continuous(name="\nCumulative Fraction\n",  
		expand=c(0,0), 
		limits=c(ymin, ymax),
		position = "right",
    guide = guide_prism_minor() # add minor ticks
		) +  

  # Theming
    common_theme() +
    scale_color_manual(values=c(jens_colors), aesthetics=c("color", "fill")) + 
    facet_grid(left_label~top_label, switch='y')

  ######## Plot Options #########

  if (show_legend == FALSE){
    plt <- plt + theme(legend.position = "none")
  }

  if (dots) { # individual replicate dots
    plt <- plt +   geom_point(data=out, aes(x=Pos+0.5, y=Total_Percentage))
  }

  return(list(plot=plt, data=data_out, n_table=n_condition_reporter(df, gene, multimap=multi_locus)))
}

tail_bar_grapher <- function(df, gene, start=-10, stop=10, gimme=F, ymin=0, ymax=1, AUCGcolors=AUCGcolors, show_legend=TRUE, dots=FALSE, multi_locus=F, analysis_min=-100, analysis_max=100, mature_end=0, order="") {
  if (grepl("|", gene, fixed=T)) multi_locus=T #if entering multiple genes, automatically does multilocus search
  if (multi_locus) df <- multimap_gene_subsetter(df, gene) # slower but more exhaustive subsetter
  else if(startsWith(gene, "ENS")) df <- filter(df, EnsID==gene) # handle ENSIDs
  else df <- filter(df, Gene_Name==gene) 
  df <- filter(df, End_Position>=analysis_min, End_Position<=analysis_max) # throws out anything outside of analysis window
  df$End_Position <- df$End_Position - mature_end # user correction of mature-end location


  #Pre-populate dataframe
  out <- data.frame(Pos="", Nuc="", Freq="", Sample="", Condition="")
  
  
  ######## Make matrices and dataframes ##############

  for (sample in unique(df$Sample)){
    temp_matrix <- tail_end_matrix_maker(filter(df, Sample==sample), xmin=start, xmax=stop)
    out <- rbind(out, tail_end_matrix_to_df(temp_matrix, sample, filter(df, Sample==sample)$Grouping[1]))
  }
  
  out <- out[-1,]
  if (order!="") out$Condition <- factor(out$Condition, levels=order)
  if (gimme){ # for debugging
    return(out)
  }

    #################### Data Pre-processing #######################
  
  # Fix T -> U and put in the correct order
  out$Nuc <- factor(out$Nuc, levels=c("A", "T", "C", "G", "X"))
  out$Nuc <- recode_factor(out$Nuc, T="U", X="Terminal Encoded")
  out$Nuc <- factor(out$Nuc, levels=c("A", "U", "C", "G", "Terminal Encoded"))

  # Fix numeric type
  out$Pos <- as.numeric(out$Pos)
  out$Freq <- as.numeric(out$Freq)

  plt <- out %>%
    group_by(Condition, Pos, Nuc) %>% 
    summarise(freq_avg = mean(Freq)) %>%
    ungroup()

  data <- plt # save data for output
  ######################### Plotting ############################
  plt$Nucleotide <- plt$Nuc # This is the easiest way to fix the legend unfortunately
  plt <- plt %>%

    # main bar graph
    ggplot(aes(x=Pos, y=freq_avg, color=Nucleotide, fill=Nucleotide)) +
    geom_bar(stat='identity') +
    facet_rep_grid(rows=vars(Condition), cols=vars(toString(gene)), switch="y") +
    

    # themeing
    common_theme() +
    theme(panel.spacing = unit(2, "lines")) + # Increase spacing between facets
    scale_color_manual(values=c("#7EC0EE", "#1f78b4", "#33a02c", "#A2CD5A", "grey"), aesthetics=c("color", "fill")) +
  # Axes
  scale_x_continuous(name="Position", 
		expand=c(0,0),
		limits=c(start,stop),
    guide = guide_prism_minor()) + # minor ticks
  scale_y_continuous(name="\nFraction\n", 
		expand=c(0,0), 
		limits=c(ymin, ymax),
		position = "right",
    guide = guide_prism_minor() # minor ticks
		) + 
    theme(panel.spacing = unit(0.1, "lines")) 

  # options
  if (show_legend == FALSE){
    plt <- plt + theme(legend.position = "none")
  }

  if (dots) { #doesn't work yet, but probably too messy
    #plt <- plt + geom_point(data=out[-1,], aes(x=Pos, y=reads_at_end))
  }

  return(list(data=data, plot=plt, n_table=n_condition_reporter(df, gene, multimap=multi_locus)))
}

tail_logo_grapher <- function(df, gene, xmin=1, xmax=10, ymin=0, ymax=1, gimme=F, multi_locus=F, analysis_min=-100, analysis_max=100, mature_end=0, order="") {
  if (grepl("|", gene, fixed=T)) multi_locus=T #if entering multiple genes, automatically does multilocus search
  ############ data pre-processing #################
  if (multi_locus) df <- multimap_gene_subsetter(df, gene) # slower but more exhaustive subsetter
  else if(startsWith(gene, "ENS")) df <- filter(df, EnsID==gene) # handling for ensembl IDs
  else df <- filter(df, Gene_Name==gene) 
  df <- filter(df, End_Position>=analysis_min, End_Position<=analysis_max) # throws out things that are outside of analysis range
  df$End_Position <- df$End_Position - mature_end # user correction for mature end


  df <- mutate(df, Grouping=replace(Grouping, Grouping=="", " ")) # needed blank grouping to be a character or else matrices blow up

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
  plot_list <- list() # to be added to for faceting

  out$Nuc <- recode_factor(as.factor(out$Nuc), "T"="U") # Fix T -> U

  conditions <- c(unique(out$Condition))
  if (order!="") conditions<-order

  for (i in 1:length(conditions)){
    plot_df <- out %>% 
                filter(Condition==conditions[i]) %>% # Look at singular condition
                pivot_wider(id_cols=Nuc, names_from=Pos, values_from=freq_avg) %>% 
                {. ->> nucs} %>% # Save intermediate to get nucleotide names
                select(-Nuc) %>% # Remove Nuc column
                as.matrix(nrow=4, ncol=ncol(nucs)-1) # ggseqlogo needs a probability matrix

    rownames(plot_df) <- nucs$Nuc # Names from Nuc intermediate makes sure we don't get them out of order
    plot_list[[conditions[i]]] <- plot_df
  }
  data <- plot_list # save data for output
  cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'U'), 
                      cols=c("#7EC0EE","#33a02c", "#A2CD5A", "#1f78b4")) 
 
  	

    plt <- ggplot() + 
      geom_logo(data=plot_list, method='custom', font='roboto_bold', col_scheme=cs1) +
      facet_rep_grid(rows=vars(seq_group), cols=vars(toString(gene)), switch='y') +


      # Themeing
      common_theme() +   
      theme(panel.spacing = unit(0.1, "lines")) + # Increase spacing between facets  

      # Axes
      scale_x_continuous(
        name="Tail length",
        breaks = seq(1, xmax-xmin+1, 1),
        labels = seq(xmin, xmax, 1),
        expand=c(0,0),
        limits=c(0.5, xmax-xmin+1.5)
	    ) +
      scale_y_continuous(
        name="\nFraction Nucleotide\n",
        breaks=seq(ymin, ymax, 0.2),
        labels=seq(ymin, ymax, .2),
        expand=c(0,0),
        limits=c(0, ymax), 
        position = "right",
        guide = guide_prism_minor() # add minor ticks
	    ) 

      return(list(data=data, plot=plt, n_table=n_condition_reporter(df, gene, multimap=multi_locus)))
}

tail_pt_nuc_grapher <- function(df, gene, gimme=F, ymin=0, ymax=6, pdisplay=F, multi_locus=F, analysis_min=-100, analysis_max=100, mature_end=0, order=""){
  if (grepl("|", gene, fixed=T)) multi_locus=T #if entering multiple genes, automatically does multilocus search
  if (multi_locus) df <- multimap_gene_subsetter(df, gene) # slower but more exhaustive subsetter
  else if(startsWith(gene, "ENS")) df <- filter(df, EnsID==gene) # handling for ensembl IDs
  else df <- filter(df, Gene_Name==gene) 
  df <- filter(df, End_Position>=analysis_min, End_Position<=analysis_max) # throws out observations outside of the analysis window
  df$End_Position <- df$End_Position - mature_end # user correction for mature end

  out <- data.frame(Sample="", Condition="", Nuc="", Frequency="")

  for (sample in unique(df$Sample)){
    #save total number of reads for that gene in that sample
    total <- sum(filter(df, Sample==sample)$Count)

    for (nuc in c("A", "C", "G", "T")) {
      # Find frequency of each nucleotide at each position
      #frequency <- sum(filter(df, Sample==sample, substr(Tail_Sequence, 1, 1) == nuc)$Count)/total
      frequency <- sum(
        unlist(
          lapply(
            filter(df, Sample==sample)$Tail_Sequence, 
            fraction_nuc_tail, pattern=nuc)
              )
            ) / sum(filter(df, Sample==sample)$Count)
      # Add to out df
      out <- rbind(out, c(sample, filter(df, Sample==sample)$Grouping[1], nuc, frequency))
      }
    
  }
  # Change to numeric types
  out=out[-1,]
  out$Frequency <- as.numeric(out$Frequency)
  out$Nuc <- recode_factor(out$Nuc, T="U")
  out$Nuc <- factor(out$Nuc, levels=c("A", "C", "G", "U"))
  if (order!="") out$Condition <- factor(out$Condition, levels=order)

  out_summed <-out %>%
    group_by(Condition, Nuc) %>% 
    summarise(freq_avg = mean(Frequency), sd= sd(Frequency), n=n(), se=sd/sqrt(n)) %>%
    ungroup()

  if(gimme){ # debugging
    return(out)
  }
  data <- out
  # t-testing
  pvals = c()
  conditions = unique(out$Condition)
  if (pdisplay){
    if (length(conditions) != 2){
      return("Can't do the stats on a different number than 2 conditions because statistics is hard") # Better error handling here
    }
    for (nuc in unique(out$Nuc)){
      temp_val <- t.test(filter(out, Nuc==nuc, Condition==conditions[1])$Frequency, filter(out, Nuc==nuc, Condition==conditions[2])$Frequency)$p.val
      pvals <- append(pvals, round(temp_val,3))
    }
  }
  ################ Plotting ##########################
  plt <- ggplot(out) + 
    coord_cartesian() +
    theme_classic() +

    
    # Individual dots
    geom_point(aes(y=Frequency, color=Condition, x=Nuc),size=1.5, position=position_jitter(width=.3)) +

    # Lines for mean
    geom_linerange(data=out_summed, size=1.5, aes(y=freq_avg,  
      xmin=rep(c(.6, 1.6, 2.6, 3.6), length(conditions)),
      xmax=rep(c(1.4,2.4,3.4,4.4), length(conditions)),
      color=Condition))+

    facet_grid(cols=vars(toString("PT-Tail")), rows=vars(toString(gene)), switch='y') +

    geom_errorbar(data=out_summed, size=1,
      aes(x=Nuc, ymin=freq_avg-se, ymax=freq_avg+se, width=0)) +

    # Themeing
    common_theme() +
    scale_color_manual(values=c(jens_colors), aesthetics=c("color", "fill")) + 
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
		

    #Axes
    scale_y_continuous(name="\nAvg #Nucs Per RNA\n",
      limits=c(ymin-0.01, ymax), 
      position = "right",
      guide = guide_prism_minor()) + # minor ticks

    scale_x_discrete(name="PT-Tail")



  ######## options ############

  if (pdisplay){
    pvals <- lapply(pvals, round, digits=4)
    pvals <- paste0("p=", pvals)
    plt <- plt + geom_text(data=out_summed, size=4.5, aes(x=Nuc, y=ymax, label=append(pvals, c("", "", "", "")))) # this is dumb, but it works
  }


  return(list(data=data, plot=plt, n_table=n_condition_reporter(df, gene, multimap=multi_locus)))
}

############## Helper Functions #####################
 multimap_gene_subsetter <- function(df, query){
  if (startsWith(query,"ENSG")){ # handler for ensids
      query <- str_split(query, "\\|")[[1]] #splits query by | character
      
      return (df[unlist(lapply(df$EnsID, function(x) length(
          intersect(
            str_split(x, "\\|")[[1]], 
            query)
          )>0)),])
    
  }
  else{
      query <- str_split(query, "\\|")[[1]]
      #if (length(query)==1) return(filter(df, Gene_Name==query))
    
      df[unlist(lapply(df$Gene_Name, function(x) length(
          intersect(
            str_split(x, "\\|")[[1]], 
            query)
          )>0)),]
  }
}
dfBuilder <- function(files, grouping){
  # There's gotta be a better way to do this, please tell me at timnicholsonshaw@gmail.com
  # Prep dataframe
  out <- data.frame(Count="",
                    EnsID="",
                    Gene_Name="",
                    End_Position="",
                    Tail_Length="",
                    Tail_Sequence="",
                    Sample="",
                    Grouping="")
  for (i in c(1:length(files))){
    temp_df <- read_csv(files[i]) 
    if ("Sequence" %in% names(temp_df)) temp_df <- select(temp_df, -Sequence) # prevent problems with sequence containing tail files
    temp_df$Sample <- files[i] # Add file name as sample name, maybe do prefix?
    temp_df$Grouping <- grouping[i] # Add sample metadata

    out <- rbind(out, temp_df)
  }
  out$Count <- as.numeric(out$Count)
  out$Tail_Length <- as.numeric(out$Tail_Length)
  out$End_Position <- as.numeric(out$End_Position)

  return (out[-1,])
}
discover_candidates <- function(df, min=1, conditions=NA, isShiny=T){

  if (is.na(conditions)){
    conditions <- unique(df$Grouping)
    }
  if (length(conditions) != 2){ # Statistics is hard
    return("Candidate discovery doesn't support this number of conditions :(, try with just 2")
  }
  genes <- unique(df$Gene_Name)

  pb <- progress_bar$new(total=length(genes),
          format = "Finding candidates [:bar] :percent | :elapsed | Estimated eta: :eta",
          clear=FALSE)

  out <- data.frame(Gene=NA, pval_end_position=NA, pval_PT_tail=NA,delta_end_pos=NA, delta_PT_tail=NA, n_con1=NA, n_con2=NA)
  
  if (isShiny){
    withProgress( message="Finding Candidates", value=0, { # Progress bar wrapper for shiny
      for (i in 1:length(genes)) {
        pb$tick()
        incProgress(1/length(genes))
        # toss if it doesn't meet the minimum number of reads
        if (genes[i] == "None") next
        con1 = filter(df, Grouping==conditions[1],Gene_Name==genes[i])
        con2 = filter(df, Grouping==conditions[2],Gene_Name==genes[i])
        if (nrow(con1) < min) next
        if (nrow(con1) < min) next


        statistic_three_end <- tryCatch({
          suppressWarnings(ks.test(con1$End_Position,
                            con2$End_Position))
                            },
          error = function(e) e
        )
        
        #Erros to NA
        if(inherits(statistic_three_end, "error")) statistic_three_end$p.value <- NA


        statistic_tail_length <- tryCatch({
          suppressWarnings(ks.test(con1$Tail_Length,
                            con2$Tail_Length))
                            },
          error = function(e) e
        )
        
        # Errors to NA
        if(inherits(statistic_tail_length, "error")) statistic_tail_length$p.value <- NA
        # report how many observations in each condition
        n <- c(sum(con1$Count), sum(con2$Count))

        # calculate mean differences 

        del_end <-  mean(con2$End_Position) - mean(con1$End_Position)
        del_tail <- mean(con2$Tail_Len) - mean(con1$Tail_Len)

        out <- rbind(out, c(genes[i], statistic_three_end$p.value, statistic_tail_length$p.value, del_end,del_tail, n[1], n[2]))
      }
    }) 
  }
  else{
    for (i in 1:length(genes)) {
        # toss if it doesn't meet the minimum number of reads
        if (genes[i] == "None") next

        con1 = filter(df, Grouping==conditions[1],Gene_Name==genes[i])
        con2 = filter(df, Grouping==conditions[2],Gene_Name==genes[i])
        if (nrow(con1) < min) next
        if (nrow(con1) < min) next


        statistic_three_end <- tryCatch({
          suppressWarnings(ks.test(con1$End_Position,
                            con2$End_Position))
                            },
          error = function(e) e
        )
        
        #Erros to NA
        if(inherits(statistic_three_end, "error")) statistic_three_end$p.value <- NA


        statistic_tail_length <- tryCatch({
          suppressWarnings(ks.test(con1$Tail_Length,
                            con2$Tail_Length))
                            },
          error = function(e) e
        )
        
        # Errors to NA
        if(inherits(statistic_tail_length, "error")) statistic_tail_length$p.value <- NA

        n <- c(sum(con1$Count), sum(con2$Count))


        del_end <-  mean(con2$End_Position) - mean(con1$End_Position)
        del_tail <- mean(con2$Tail_Len) - mean(con1$Tail_Len)

        out <- rbind(out, c(genes[i], statistic_three_end$p.value, statistic_tail_length$p.value,del_end,del_tail, n[1], n[2]))
     }
  }

  out<-out[-1,] # Remove crap line
  # Change to numeric types
  out$pval_end_position <- as.numeric(out$pval_end_position)
  out$pval_PT_tail <- as.numeric(out$pval_PT_tail)
  out$delta_end_pos <- as.numeric(out$delta_end_pos)
  out$delta_PT_tail <- as.numeric(out$delta_PT_tail)
  out$n_con1 <- as.numeric(out$n_con1)
  out$n_con2 <- as.numeric(out$n_con2)

  out <- out[order(out$pval_end_position),]

  # reorder columns and make more user readable
  out <- select(out, Gene, delta_end_pos, pval_end_position, delta_PT_tail, pval_PT_tail, n_con1, n_con2) %>%
    rename(
      "Δ End Position"=delta_end_pos, 
      "p-value End Position"=pval_end_position,
      "Δ PT Tail"=delta_PT_tail,
      "p-value PT Tail"=pval_PT_tail,
      "# Reads Condition 1"=n_con1,
      "# Reads Condition 2"=n_con2
      )
  return (out)
  

}
tail_end_matrix_maker <- function(df, xmin, xmax){
  positions <- c(xmin:xmax)
  tail_end_matrix <- matrix(0,5, length(positions))
  rownames(tail_end_matrix) <- c("A", "C", "G", "T", "X")
  colnames(tail_end_matrix) <- positions
  total <- sum(df$Count)
  
  for (row in 1:nrow(df)) {
    End_Position <- df[row, "End_Position"]
    Tail_Sequence  <- df[row, "Tail_Sequence"]
    count <- df[row, "Count"]
    Tail_Length <- df[row, "Tail_Length"]
    
    if (is.na(Tail_Sequence)){
      if(End_Position> xmax | End_Position < xmin) next
      tail_end_matrix["X", toString(End_Position)] <- tail_end_matrix["X", toString(End_Position)] + count 
      next
    }
    
    start <- End_Position - Tail_Length
    if(start<=xmax & start>=xmin){
      tail_end_matrix["X", toString(start)] <- tail_end_matrix["X", toString(start)] + count
    }
    

    for (i in 1:Tail_Length){
      nuc <- substr(Tail_Sequence, i, i)
      if (nuc=='N') next
      pos <- i+start
      if (pos > xmax | pos < xmin) next
      tail_end_matrix[nuc, toString(pos)] <- tail_end_matrix[nuc, toString(pos)] + count

    }
  
  }
  tail_end_matrix/total
}
tail_end_matrix_to_df <- function(matrix, sample, condition){
  out <- data.frame(Pos="", Nuc="", Freq="", Sample="", Condition="")
  for (i in rownames(matrix)){
    for (j in colnames(matrix)){
      out <- rbind(out, c(j, i, matrix[i,j], sample, condition))
    }
  }
  out <- out[-1,]
  out$Pos <- as.numeric(out$Pos)
  out$Freq <- as.numeric(out$Freq)
  out
}

stat_matrix_maker <- function(df, gene, con1, con2, multi_locus=F){
  if (multi_locus) df <- multimap_gene_subsetter(df, gene) # slower but more exhaustive subsetter
  else if(startsWith(gene, "ENS")) df <- filter(df, EnsID==gene) # handle ENSIDs
  else df <- filter(df, Gene_Name==gene) 

  # Pulling out samples for the listed conditions
  samples1 = unique(filter(df, Grouping==con1)$Sample) 
  samples2 = unique(filter(df, Grouping==con2)$Sample)

  # Matrix to hold results
  end_position_matrix <- matrix(ncol=length(samples1), nrow=length(samples2))
  #rownames(end_position_matrix) <- samples2
  #colnames(end_position_matrix) <- samples1
  custom_col_names <- c()
  custom_row_names <- c()
  for (i in 1:length(samples1)){
    custom_col_names <- append(custom_col_names, paste0(con1,"_", i))
  }
  for (i in 1:length(samples2)){
    custom_row_names <- append(custom_row_names, paste0(con2,"_", i))
  }
  colnames(end_position_matrix) <- custom_col_names
  rownames(end_position_matrix) <- custom_row_names


  tail_len_matrix <- matrix(ncol=length(samples1), nrow=length(samples2))
  colnames(tail_len_matrix) <- custom_col_names
  rownames(tail_len_matrix) <- custom_row_names

  print(length(samples1))
  print(length(samples2))
  
  for (i in c(1:length(samples1))){
    for (j in c(1:length(samples2))){
      a <- filter(df, Sample==samples1[i])$End_Position
      b <- filter(df, Sample==samples2[j])$End_Position

      p <- suppressWarnings(ks.test(a,b)$p.value)

      if (p<1e-50) p="<1e-50"

      end_position_matrix[i,j] <- p
    }
  }

  for (i in c(1:length(samples1))){
    for (j in c(1:length(samples2))){
      a <- filter(df, Sample==samples1[i])$Tail_Len
      b <- filter(df, Sample==samples2[j])$Tail_Len

      p <- suppressWarnings(ks.test(a,b)$p.value)
      if (p<1e-50) p="<1e-50"
      tail_len_matrix[i,j] <- p
    }
  }

  p.end_position_total <- suppressWarnings(ks.test(filter(df, Sample %in% samples1)$End_Position,filter(df, Sample %in% samples2)$End_Position)$p.value)
  if (p.end_position_total < 1e-50) p.end_position_total <- "<1e-50"
  p.tail_len_total <- suppressWarnings(ks.test(filter(df, Sample %in% samples1)$Tail_Len,filter(df, Sample %in% samples2)$Tail_Len)$p.value)
  if (p.tail_len_total < 1e-50) p.tail_len_total <- "<1e-50"

  # Also return df of number of observations
  con1_counts <- df %>% filter(Grouping==con1) %>%
                    group_by(Sample) %>% 
                    summarise(n = sum(Count))
  con1_counts <- con1_counts$n

  con2_counts <- df %>% filter(Grouping==con2) %>%
                    group_by(Sample) %>% 
                    summarise(n = sum(Count))
  con2_counts <- con2_counts$n


  n.df <- data.frame(Condition=c(custom_col_names, custom_row_names), n=c(con1_counts, con2_counts))
  return (list(
    p.mat_end_position=end_position_matrix, 
    p.end_position_total=p.end_position_total,
    p.mat_tail_len=tail_len_matrix,
    p.tail_len_total=p.tail_len_total,
    n.df=n.df
                ))
}

n_sample_reporter <- function(df, gene, multimap=F){
  # Utility to report number of observed reads for a particular gene by sample
  if (multimap) temp_df <- multimap_gene_subsetter(df, gene)
  else temp_df <- filter(df, Gene_Name==gene)
  groupings <- unique(df$Grouping)
  sample_designations <- c()
  totals <- c()
  for (i in 1:length(groupings)){
    samples <- unique(filter(df, Grouping==groupings[i])$Sample)
    for (j in 1:length(samples)){
      sample_designations <- c(sample_designations, paste0(groupings[i], "_", j))
      totals <- c(totals, sum(filter(temp_df, Grouping==groupings[i], Sample==samples[j])$Count))
    }
  }

  return (data.frame(Sample=sample_designations, n=totals))
}

n_condition_reporter <- function(df, gene, multimap=F){
  if (multimap) temp_df <- multimap_gene_subsetter(df, gene)
  else temp_df <- filter(df, Gene_Name==gene)
  groupings <- unique(df$Grouping)
  sample_designations <- c()
  totals <- c()
  for (i in 1:length(groupings)){
    samples <- unique(filter(df, Grouping==groupings[i])$Sample)
    for (j in 1:length(samples)){
      sample_designations <- c(sample_designations, paste0(groupings[i], "_", j))
      totals <- c(totals, sum(filter(temp_df, Grouping==groupings[i], Sample==samples[j])$Count))
    }
  }
  return (data.frame(Sample=sample_designations, n=totals))
}
fraction_nuc_tail <- function(tail_sequence, pattern=""){
  if (is.na(tail_sequence)) return(0)
  str_count(tail_sequence, pattern=pattern)
}
################### Pretty Making Functions ################
common_theme <- function() { 
  theme_base() +
    theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
     panel.background = element_rect(colour = NA),
     plot.background = element_rect(colour = NA),
     panel.border = element_rect(colour = "black", size=1),
     axis.title = element_text(face = "bold",size = rel(1)),
     axis.title.y = element_text(angle=90,vjust =2),
     axis.title.x = element_text(vjust = -0.2),
     axis.text = element_text(size=rel(1.1)), 
     axis.line = element_line(colour="black"),
     axis.ticks = element_line(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     legend.key = element_rect(colour = NA),
     legend.position = "right",
     legend.key.size= unit(0.6, "cm"),
     legend.margin = unit(0, "cm"),
     legend.title = element_text(face="italic", size=rel(1)),
     plot.margin=unit(c(10,5,5,5),"mm"),
     strip.background=element_rect(colour="#000000",fill="#b3b3b3"),
     strip.text = element_text(face="bold"),
    )
}
jens_colors <- c(
  "#0073c2", 
  "#efc000", 
  "#868686",
  "#cd534c", 
  "#7aa6dc", 
  "#003c67", 
  "#8f7700", 
  "#3b3b3b",
	"#a73030", 
  "#4a6990"
)