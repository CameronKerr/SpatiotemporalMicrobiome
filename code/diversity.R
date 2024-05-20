######################
# Diversity Analysis #
######################

# Computation functions
shannon_index <- function(sampleid, data) {
  # Compute the Shannon diversity of a sample
  sample_data <- data[sampleid][data[sampleid] != 0]
  sample_frac <- sample_data/sum(sample_data)
  -sum(sample_frac*log(sample_frac))
}
richness <- function(sampleid, data) {
  # Compute the richness of a sample
  sample_data <- data[sampleid][data[sampleid] != 0]
  length(sample_data)
}
eveness <- function(sampleid, data) {
  # Compute the eveness of a sample
  # Note that eveness = 1 - clonality
  sample_data <- data[sampleid][data[sampleid] != 0]
  shannon <- shannon_index(sampleid, data)
  richness <- richness(sampleid, data)
  shannon/log(richness)
}
simpson_index <- function(sampleid, data) {
  # Compute the Simpson diversity of a sample
  sample_data <- data[sampleid][data[sampleid] != 0]
  sample_frac <- sample_data/sum(sample_data)
  sum(sample_frac^2)
}
div_table_gen <- function(data){
  all_samples <- colnames(data)
  shannon_list <- c()
  richness_list <- c()
  eveness_list <- c()
  simpson_list <- c()
  for(samp in all_samples){
    shannon_list <- c(shannon_list, shannon_index(samp, data))
    richness_list <- c(richness_list, richness(samp, data))
    eveness_list <- c(eveness_list, eveness(samp,data))
    simpson_list <- c(simpson_list, 1/simpson_index(samp,data))
  }
  div_table <- data.frame ( Shannon = shannon_list,
                            Richness = richness_list,
                            Eveness = eveness_list,
                            InverseSimpson = simpson_list)
  rownames(div_table) <- all_samples
  div_table
}



# Compute diversity statistics for all samples
otu <- as.data.frame(otu_table(physeq_rarefy))
output_table <- div_table_gen(otu)
output_table <- cbind(output_table, sample_data(physeq_rarefy))

# Diversity boxplots
par(mfrow=c(1,4), mar = c(5.1, 4.1, 4.7, 2.1))
boxplot(output_table$Shannon~output_table$type, ylab = "Shannon Diversity Index",
        names = c("Culture", "On-Site"), xlab = "Sample Type", 
        main = "Shannon Diversity \n by Sample Type")
boxplot(output_table$Eveness~output_table$type, ylab = "Eveness Index", 
        names = c("Culture", "On-Site"), xlab = "Sample Type", 
        main = "Eveness by \n Sample Type")
boxplot(output_table$Richness~output_table$type, ylab = "Richness", 
        names = c("Culture", "On-Site"), xlab = "Sample Type", 
        main = "Richness by \n Sample Type")
boxplot(output_table$InverseSimpson~output_table$type, 
        ylab = "Inverse Simpson Diversity Index", names = c("Culture", "On-Site"), 
        xlab = "Sample Type", main = "Inverse Simpson \n by Sample Type")

# Compute diversity statistics for rural, urban, suburban
output_table_D <- output_table[which(output_table$type == "D"),]

# Diversity boxplots
par(mfrow=c(1,3), mar = c(5.1, 4.1, 4.7, 2.1))
boxplot(output_table_D_t$Shannon~output_table_D_t$timepoint, ylab = "Shannon Diversity Index", xlab = "Location Type", 
        main = "Shannon Diversity \n by Location Type")
boxplot(output_table_D_t$Eveness~output_table_D_t$timepoint, ylab = "Eveness Index", xlab = "Location Type", 
        main = "Eveness by \n Location Type")
boxplot(output_table_D_t$Richness~output_table_D_t$timepoint, ylab = "Richness", xlab = "Location Type", 
        main = "Richness by \n Location Type")
boxplot(output_table_D$InverseSimpson~output_table_D$location_type, 
        ylab = "Inverse Simpson Diversity Index", 
        xlab = "Location Type", main = "Inverse Simpson \n by Location Type")

# Compute diversity statistics for culture type
output_table_C <- output_table[which(output_table$type == "C"),]


# Diversity boxplots
par(mfrow=c(1,3), mar = c(5.1, 4.1, 4.7, 2.1))
boxplot(output_table_C$Shannon~output_table_C$culture, ylab = "Shannon Diversity Index", xlab = "Culture Type", 
        main = "Shannon Diversity \n by Culture Type")
boxplot(output_table_C$Eveness~output_table_C$culture, ylab = "Eveness Index", xlab = "Culture Type", 
        main = "Eveness by \n Culture Type")
boxplot(output_table_C$Richness~output_table_C$culture, ylab = "Richness", xlab = "Culture Type", 
        main = "Richness by \n Culture Type")
boxplot(output_table_C$InverseSimpson~output_table_C$culture, 
        ylab = "Inverse Simpson Diversity Index", 
        xlab = "Culture Type", main = "Inverse Simpson \n by Culture Type")

# Blocked ANOVA (t-test because 2 variable) for on-site vs culture
fit <- lm(Richness~type + timepoint + location,data=output_table)
anova(fit)

# Blocked ANOVA for culture type diversity
shapiro.test(output_table_C$Shannon)
leveneTest(Shannon ~ culture, output_table_C)
fit <- lm(Shannon~culture + location + timepoint,data=output_table_C)
anova(fit)

shapiro.test(output_table_C$Eveness)
leveneTest(Eveness ~ culture, output_table_C)
fit <- lm(Eveness~culture + location + timepoint,data=output_table_C)
anova(fit)

shapiro.test(output_table_C$Richness)
leveneTest(Richness ~ culture, output_table_C)
fit <- lm(Richness~location_type + timepoint,data=output_table_C)
anova(fit)

# Blocked ANOVA location type
leveneTest(Shannon ~ location_type, output_table_D)
fit <- lm(Shannon~location_type + timepoint,data=output_table_D)
anova(fit)
TukeyHSD( aov(Shannon~location_type + timepoint,data=output_table_D), conf.leve=0.95)

leveneTest(Shannon ~ location_type, output_table_C)
fit <- lm(Shannon~location_type + timepoint,data=output_table_C)
anova(fit)

leveneTest(Eveness ~ location_type, output_table_D)
fit <- lm(Eveness~location_type + timepoint,data=output_table_D)
anova(fit)
TukeyHSD( aov(Eveness~location_type + timepoint,data=output_table_D), conf.leve=0.95)

leveneTest(Eveness ~ location_type, output_table_C)
fit <- lm(Eveness~location_type + timepoint,data=output_table_C)
anova(fit)

leveneTest(Richness ~ location_type, output_table_D)
fit <- lm(Richness~location_type + timepoint,data=output_table_D)
anova(fit)
TukeyHSD( aov(Richness~location_type + timepoint,data=output_table_D), conf.leve=0.95)

leveneTest(Richness ~ location_type, output_table_C)
fit <- lm(Richness~location_type + timepoint,data=output_table_C)
anova(fit)

# Blocked ANOVA timepoint
output_table_D_t <- output_table_D[which(output_table_D$timepoint %in% c("1", "2")),]
output_table_C_t <- output_table_C[which(output_table_C$timepoint %in% c("1", "2")),]

leveneTest(Shannon ~ timepoint, output_table_D_t)
fit <- lm(Shannon~timepoint + location ,data=output_table_D_t)
anova(fit)

leveneTest(Shannon ~ timepoint, output_table_C_t)
fit <- lm(Shannon~timepoint + location + culture,data=output_table_C_t)
anova(fit)

leveneTest(Eveness ~ timepoint, output_table_D_t)
fit <- lm(Eveness~timepoint + location ,data=output_table_D_t)
anova(fit)

leveneTest(Eveness ~ timepoint, output_table_C_t)
fit <- lm(Eveness~timepoint + location + culture,data=output_table_C_t)
anova(fit)

leveneTest(Richness ~ timepoint, output_table_D_t)
fit <- lm(Richness~timepoint + location ,data=output_table_D_t)
anova(fit)

leveneTest(Richness ~ timepoint, output_table_C_t)
fit <- lm(Richness~timepoint + location + culture,data=output_table_C_t)
anova(fit)
