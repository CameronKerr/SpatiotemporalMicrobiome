####################
# Overlap Analysis #
####################

# Morisita index calculation function
morisitaindex <- function(sampleid1, sampleid2, data){
  samples_data <- data[c(sampleid1, sampleid2)]
  cleaned_comb_data <- samples_data[!(samples_data[1]==0 & samples_data[2]==0),]
  frac_data <- cbind(cleaned_comb_data[1]/sum(cleaned_comb_data[1]), cleaned_comb_data[2]/sum(cleaned_comb_data[2]))
  numerator <- 2 * sum(frac_data[1]*frac_data[2])
  simpson1 <- simpson_index(sampleid1, data)
  simpson2 <- simpson_index(sampleid2, data)
  
  numerator/(simpson1 + simpson2)
}
# Generating morisita overlap matrix
overlapmatrix_gen <- function(data){
  all_samples <- colnames(data)
  n_samp <- length(all_samples)
  overlapmatrix <- as.data.frame(matrix(data = NA, nrow = n_samp, ncol = n_samp))
  colnames(overlapmatrix) <- all_samples
  rownames(overlapmatrix) <- all_samples
  pb <- txtProgressBar(min = 0,max = n_samp^2,style = 3,   width = 100, char = "=")
  i <- 0
  for(n in 1:n_samp){
    for(m in 1:n_samp){
      i <- i + 1
      setTxtProgressBar(pb, i)
      if(n < m){ 
        overlapmatrix[n,m] <- morisitaindex(all_samples[n], all_samples[m], data)
        overlapmatrix[m,n] <- overlapmatrix[n,m]
      }
    }
  }
  close(pb)
  overlapmatrix <- as.matrix(overlapmatrix)
  diag(overlapmatrix) <- 1
  overlapmatrix
} 
# Generating heatmap function
overlapheatmap_gen <- function(overlapmatrix){
  annotation <- as.data.frame(matrix(NA, nrow=nrow(overlapmatrix)))
  rownames(annotation) <- rownames(overlapmatrix)
  
  Urban <- c("Cedarvale", "Humber", "Grenadier")
  Rural <- c("Preston", "KSR", "Terracotta")
  Suburban <- c("Rouge", "Aviemore", "Bruce")
  annotation$type <- NA
  annotation$location <- NA
  annotation$timepoint <- NA
  annotation$area <- NA
  for(i in 1:length(rownames(annotation))){
    split_r <- strsplit(rownames(annotation)[i], "_")[[1]]
    annotation$type[i] <- split_r[1]
    annotation$location[i] <- split_r[2]
    if(annotation$type[i] == "C"){
      annotation$timepoint[i] <- split_r[4]
    } else {
      annotation$timepoint[i] <- split_r[3]
    }
    if(annotation$location[i] %in% Urban){
      annotation$area[i] <- "Urban"
    } else if(annotation$location[i] %in% Rural){
      annotation$area[i] <- "Rural"
    } else {annotation$area[i] <- "Suburban"}
  }
  annotation <- annotation[2:length(annotation)]
  ann_colors <- list(
    type = c(C = "#FFFF00", D = "#1CE6FF"),
    location = c(Aviemore = "#FF34FF", BruceH = "#FF4A46", Cedarvale = "#008941", 
                 Grenadier = "#006FA6",Humber = "#A30059", KSR = "#FFDBE5", 
                 Preston = "#7A4900", Rouge = "#FDC086", Terracotta = "#BEAED4"),
    timepoint = c(`1` = "darkred", `2` = "#FC8D59", `3` = "#FEF0D9"),
    area = c(Urban = "#ff7e1c", Suburban = "#4b4ebf", Rural = "#4bbf77")
  )
  overlapheatmap <<- pheatmap(overlapmatrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, fontsize_col = 6, annotation_col = annotation, annotation_row = annotation, 
                              annotation_colors = ann_colors, cellheight=6, cellwidth=6)
}

# Generating heatmaps
otu <- as.data.frame(otu_table(physeq_rarefy))
moroverlap_df <- overlapmatrix_gen(otu)
unifracdist <- as.matrix(distance(physeq_rarefy, method="unifrac"))

overlapheatmap_gen(unifracdist)
overlapheatmap_gen(moroverlap_df)
