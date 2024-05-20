#############
# Bar plots #
#############

# Color list
color = colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# Function for creating barplot for a specific subset of samples and taxonomic rank
barplot <- function(sample_list, phyloseq_object, taxonomic_rank){
  otu <- as.data.frame(otu_table(phyloseq_object))
  tax <- as.data.frame(tax_table(phyloseq_object))
  otu_df <- otu[, sample_list]
  otu_df[taxonomic_rank] <- NA
  for(i in 1:nrow(otu)){
    otu_df[i, taxonomic_rank] <- tax[which(rownames(tax) == rownames(otu)[i]), taxonomic_rank]
  }
  otu_df <- otu_df[which(rowSums(otu_df[1:(ncol(otu_df) - 1)]) != 0),]
  for(c in sample_list){
    otu_df[c] <- otu_df[c]/sum(otu_df[c])
  }
  # Creating fill for bar plot
  n_taxa <- nrow(unique(otu_df[taxonomic_rank]))
  color_list <- color[floor(seq(1, 433, length.out=n_taxa))]
  # Converting from wide to long for bar plot
  otu_bar_df <- pivot_longer(otu_df, cols = sample_list[1]:sample_list[length(sample_list)],
                             names_to = "SampleId", values_to = "RelativeAbundance")
  otu_bar_df$SampleId <- factor(otu_bar_df$SampleId, levels = c(sample_list))
  levels(otu_bar_df$SampleId) <- c(sample_list)
  
  p <- ggplot(data = otu_bar_df, aes(x = SampleId, y = RelativeAbundance, fill = eval(as.name(taxonomic_rank)))) + 
    geom_bar(stat = "identity", width=0.7) + scale_fill_manual(values = color_list) +
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

# Lists of sample by location and sample typ
D_BruceH <- c("D_BruceH_1_A", "D_BruceH_1_B", 
             "D_BruceH_2_A", "D_BruceH_2_B")
C_BruceH <- c("C_BruceH_YMA_1", "C_BruceH_YMA_2",
              "C_BruceH_TY_2")

D_Cedarvale <- c("D_Cedarvale_1_A", "D_Cedarvale_1_B",
                 "D_Cedarvale_2_A", "D_Cedarvale_2_B",
                 "D_Cedarvale_3_A", "D_Cedarvale_3_B")
C_Cedarvale <- c("C_Cedarvale_LB_1", "C_Cedarvale_TY_1", "C_Cedarvale_YMA_1",
                 "C_Cedarvale_LB_2", "C_Cedarvale_TY_2", "C_Cedarvale_YMA_2")

D_Grenadier <- c("D_Grenadier_1_A", "D_Grenadier_1_B",
                 "D_Grenadier_2_A", "D_Grenadier_2_B",
                 "D_Grenadier_3_A", "D_Grenadier_3_B")
C_Grenadier <- c("C_Grenadier_LB_1", "C_Grenadier_TY_1", "C_Grenadier_YMA_1",
                 "C_Grenadier_LB_2", "C_Grenadier_TY_2", "C_Grenadier_YMA_2")

D_Humber <- c("D_Humber_1_A", "D_Humber_1_B",
              "D_Humber_2_A", "D_Humber_2_B",
              "D_Humber_3_A", "D_Humber_3_B")
C_Humber <- c("C_Humber_LB_1", "C_Humber_TY_1", "C_Humber_YMA_1",
              "C_Humber_LB_2", "C_Humber_TY_2", "C_Humber_YMA_2")

D_KSR <- c("D_KSR_1_A", "D_KSR_1_B", 
           "D_KSR_2_A", "D_KSR_2_B")
C_KSR <- c("C_KSR_LB_1", "C_KSR_TY_1", "C_KSR_YMA_1",
           "C_KSR_LB_2", "C_KSR_TY_2", "C_KSR_YMA_2")

D_Preston <- c("D_Preston_1_A", "D_Preston_1_B", 
               "D_Preston_2_A", "D_Preston_2_B")
C_Preston <- c("C_Preston_LB_1", "C_Preston_TY_1", "C_Preston_YMA_1",
               "C_Preston_LB_2", "C_Preston_TY_2", "C_Preston_YMA_2")

D_Rouge <- c("D_Rouge_1_A", "D_Rouge_1_B",
             "D_Rouge_2_A", "D_Rouge_2_B")
C_Rouge <- c("C_Rouge_LB_1", "C_Rouge_TY_1", "C_Rouge_YMA_1",
             "C_Rouge_LB_2", "C_Rouge_TY_2", "C_Rouge_YMA_2")

D_Terracotta <- c("D_Terracotta_1_A", "D_Terracotta_1_B",
                  "D_Terracotta_2_A", "D_Terracotta_2_B")
C_Terracotta <- c("C_Terracotta_LB_1", "C_Terracotta_TY_1", "C_Terracotta_YMA_1",
                  "C_Terracotta_TY_2", "C_Terracotta_YMA_2")

D_Aviemore <- c("D_Aviemore_1_A", "D_Aviemore_1_B", 
                "D_Aviemore_2_A", "D_Aviemore_2_B")


