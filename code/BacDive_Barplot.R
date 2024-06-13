###############################
# Functional Microbe Analysis #
###############################

# Generate list of unique genus
g_unique_genus = c(unique(tax_table(physeq_rarefy)[,6]))[-c(1)]
unique_genus = c()
for(i in 1:length(g_unique_genus)){unique_genus[i] <- c(strsplit(g_unique_genus[i], "__")[[1]][2])}

# Generate list of unique species
s_unique_species = c(unique(tax_table(physeq_rarefy)[,7]))[-c(1)]
unique_species = c()
for(i in 1:length(s_unique_species)){unique_species[i] <- c(strsplit(s_unique_species[i], "__")[[1]][2])}

# Load anaerobe data
anaerobedata = unique(read.table("data/anaerobe_data.csv", sep = ",", header = TRUE)[,5:6])
anaerobedata$genus = NA
anaerobedata$Species = NA
for(i in 1:nrow(anaerobedata)){
  anaerobedata$genus[i] <- c(strsplit(anaerobedata$species[i], " ")[[1]][1])
  anaerobedata$Species[i] <- c(strsplit(anaerobedata$species[i], " ")[[1]][2])
  }

# Adding in anaerobe data to unique_species
species_df <- data.frame('species' = unique_species, 'ox_req' = NA)
for(i in 1:length(unique_species)){
  if(species_df$species[i] %in% anaerobedata$species){
    species_df$ox_req[i] <- anaerobedata$Oxygen.tolerance[which(anaerobedata$species == species_df$species[i])][1]
  }
}




relevant_anaerobedata <- anaerobedata[which(anaerobedata$species %in% unique_species),]

barplot_anaerobic <- function(sample_list, phyloseq_object){
  otu <- as.data.frame(otu_table(phyloseq_object))
  tax <- as.data.frame(tax_table(phyloseq_object))
  otu_df <- otu[, sample_list]
  otu_df["Species"] <- NA
  for(i in 1:nrow(otu)){
    otu_df[i, "Species"] <- tax[which(rownames(tax) == rownames(otu)[i]), "Species"]
    otu_df[i, "Order"] <- tax[which(rownames(tax) == rownames(otu)[i]), "Order"]
  }
  otu_df <- otu_df[which(rowSums(otu_df[1:(ncol(otu_df) - 2)]) != 0),]
  for(c in sample_list){
    otu_df[c] <- otu_df[c]/sum(otu_df[c])
  }
  # Create column depicting oxygen dependency
  otu_df$Oxygen.Tolerance <- NA
  for(i in 1:nrow(otu_df)){
    specie = strsplit(otu_df$Species[i], "__")[[1]][2]
    if(specie %in% relevant_anaerobedata$species){
      relevant = relevant_anaerobedata[which(relevant_anaerobedata$species == specie),]
      otu_df$Oxygen.Tolerance[i] = relevant$Oxygen.tolerance[1]
    }
    else {
      otu_df$Oxygen.Tolerance[i] = "Unknown"
    }
  }
  # Converting from wide to long for bar plot
  otu_bar_df <- pivot_longer(otu_df, cols = sample_list[1]:sample_list[length(sample_list)],
                             names_to = "SampleId", values_to = "RelativeAbundance")
  otu_bar_df$SampleId <- factor(otu_bar_df$SampleId, levels = c(sample_list))
  levels(otu_bar_df$SampleId) <- c(sample_list)
  
  p <- ggplot(data = otu_bar_df, aes(x = SampleId, y = RelativeAbundance, fill = Oxygen.Tolerance)) + 
    geom_bar(stat = "identity", width=0.7) + #scale_fill_manual(values = unique(otu_df$Oxygen.Tolerance)) +
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
  
barplot_anaerobic(c("D_Humber_1_A", "D_Humber_1_B"), physeq_rarefy)