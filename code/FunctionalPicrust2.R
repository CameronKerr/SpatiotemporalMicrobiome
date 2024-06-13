# Load KEGG orthology data and format to numeric instead of integer
ko_df <- read.table(file = "data/feature-table-ko.tsv", sep="\t", header = TRUE)
ko_df <- as.data.frame(cbind(ko_df[1], lapply(ko_df[2:ncol(ko_df)], as.numeric)))
ko_df <- cbind(ko_df[1], ko_df[3:ncol(ko_df)])

# Load pathways data and format
path_df <- read.table(file = "data/feature-table-path.tsv", sep="\t", header = TRUE)
path_df <- cbind(path_df[1], path_df[3:ncol(ko_df)])

###########
# Barplot #
###########

barplot_ko <- function(sample_list, ko_df, ko_number){
  # Subset data according to only samples in sample_list and to ko_number
  plot_df <- ko_df[which(ko_df$X == ko_number),]
  plot_df <- plot_df[2:length(plot_df)]
  plot_df <- plot_df[c(sample_list)]
  
  # Get columns sums and standardize plot_df to relative abundance
  ko_sums <- colSums(ko_df[c(sample_list)])
  for(samp in 1:length(sample_list)){
      plot_df[samp] <- plot_df[samp]/ko_sums[i]
  }
  
  # Plot abundance of ko_number for all samples
  p_df <- data.frame('samples' = colnames(plot_df), 'abundance' = unlist(plot_df[1,]))
  p <- ggplot(p_df, aes(x=samples, y=abundance)) + geom_bar(stat = "identity", width=0.7) + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.margin = unit(c(0.2,0,0,0),"cm")) +
    labs(fill = "Class")
  print(p)
}