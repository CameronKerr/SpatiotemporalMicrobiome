################################
# Principal Component Analysis #
################################

# Functions for testing
gen_permanovadf <- function(dist_mat, sampdata, runtype = T, runculture = T){
  # Generates a dataframe which summarizes the adonis tests for all categorical
  # variables
  all_adonis <- data.frame(matrix(data=NA, nrow=5, ncol=4))
  rownames(all_adonis) <- c("Sample Type", "Location", "Timepoint", "Location Type", 
                            "Culture")
  colnames(all_adonis) <- c("SSq", "R2", "F", "p")
  
  # If 'runtype' is True an adonis test is completed grouping samples by on-site vs cultured-samples
  if(runtype){
    adonis_type <- adonis2(as.dist(dist_mat)~type, data=as(sample_data(sampdata), "data.frame"), 
                           permutations = 1000)
    all_adonis[1,] <- adonis_type[1, c(2,3,4,5)]
  }
  adonis_location <- adonis2(as.dist(dist_mat)~location, data=as(sample_data(sampdata), "data.frame"), 
                             permutations = 1000)
  all_adonis[2,] <- adonis_location[1, c(2,3,4,5)]
  adonis_timepoint <- adonis2(as.dist(dist_mat)~timepoint, data=as(sample_data(sampdata), "data.frame"), 
                              permutations = 1000)
  all_adonis[3,] <- adonis_timepoint[1, c(2,3,4,5)]
  adonis_location_type <- adonis2(as.dist(dist_mat)~location_type, data=as(sample_data(sampdata), "data.frame"), 
                                  permutations = 1000)
  all_adonis[4,] <- adonis_location_type[1, c(2,3,4,5)]
  # If 'runtype' is True an adonis test is completed grouping samples by the type of media
  if(runculture){
    adonis_culture <- adonis2(as.dist(dist_mat)~culture, data=as(sample_data(sampdata), "data.frame"), 
                              permutations = 1000)
    all_adonis[5,] <- adonis_culture[1, c(2,3,4,5)]
  }
  all_adonis
}

gen_betadf <- function(dist_mat, sampdata, runtype = T, runculture = T){
  # Generates a datframe which summarizes the beta dispersal tests for all 
  # categorical variables
  all_beta <- data.frame(matrix(data=NA, nrow=5, ncol=4))
  rownames(all_beta) <- c("Sample Type", "Location", "Timepoint", "Location Type", 
                            "Culture")
  colnames(all_beta) <- c("SSq", "MeanSSq", "F", "p")
  if(runtype){
    all_beta[1,] <- anova(betadisper(dist_mat, sampdata$type))[1,c(2,3,4,5)]
  }
  all_beta[2,] <- anova(betadisper(dist_mat, sampdata$location))[1,c(2,3,4,5)]
  all_beta[3,] <- anova(betadisper(dist_mat, sampdata$timepoint))[1,c(2,3,4,5)]
  all_beta[4,] <- anova(betadisper(dist_mat, sampdata$location_type))[1,c(2,3,4,5)]
  if(runculture){
    all_beta[5,] <- anova(betadisper(dist_mat, sampdata$culture))[1,c(2,3,4,5)]
  }
  all_beta
}

# All samples

## PCA Ordination
ordu = ordinate(physeq_rarefy, method="PCoA", distance = "unifrac")
myp <- plot_ordination(physeq_rarefy, ordu, color = "type") + 
  geom_point(size = 2.5) + ggtitle("All samples - Unifrac PCoA") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(legend.key=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA)) +
  labs(color = "Sample Type")
## Permanova
sampdata <- sample_data(physeq_rarefy)
dist_mat <- distance(physeq_rarefy, method="unifrac")

gen_permanovadf(dist_mat, sampdata)
gen_betadf(dist_mat, sampdata)

# Only on-site
d_physeq <- subset_samples(physeq_rarefy, type == "D")
d_physeq <- subset_samples(d_physeq, timepoint != 3)

## PCA Ordination
d_ordu = ordinate(d_physeq, method="PCoA", distance = "unifrac")
d_myp <- plot_ordination(d_physeq, d_ordu, color = "timepoint", shape = "location_type") + 
  geom_point(size = 2.5) + ggtitle("On-site samples - Unifrac PCoA") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(legend.key=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA)) +
  labs(color = "Timepoint", shape = "Location Type")
## Permanova
d_sampdata <- sample_data(d_physeq)
d_dist_mat <- distance(d_physeq, method="unifrac")

gen_permanovadf(d_dist_mat, d_sampdata, F, F)
gen_betadf(d_dist_mat, d_sampdata, F, F)

# Only cultured
c_physeq <- subset_samples(physeq_rarefy, type == "C")
## PCA Ordination
c_ordu = ordinate(c_physeq, method="PCoA", distance = "unifrac")
c_myp <- plot_ordination(c_physeq, c_ordu, color = "timepoint", shape = "culture") + 
  geom_point(size = 2.5) + ggtitle("Cultured samples - Unifrac PCoA") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(legend.key=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA)) +
  labs(color = "Timepoint", shape = "Culture Type")

## Permanova
c_sampdata <- sample_data(c_physeq)
c_dist_mat <- distance(c_physeq, method="unifrac")

gen_permanovadf(c_dist_mat, c_sampdata, F, T)
gen_betadf(c_dist_mat, c_sampdata, F, T)
