###############
# Rarefaction #
###############

# Plotting rarefaction curve
set.seed(1)
subsamples <- seq(0, 77600, by=500)[-1]
p <- plot_alpha_rcurve(physeq, index="observed",
                       subsamples <- subsamples,
                       lower.conf = 0.025,
                       upper.conf = 0.975,
                       group="samples",
                       label.color="brown2",
                       label.size = 3,
                       label.min=TRUE)
mycols <- c(rep("brown3", 52), rep("steelblue", 42))

myp <- p + scale_color_manual(values = mycols) + 
  xlab("Number of Sequencing Reads") + ylab("Number of Observed OTUs") +
  scale_fill_manual(values = mycols) + ggtitle("Rarefaction Curve")+
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "blank")

# Sequencing read distribution
seq_dist <- colSums(otu_table(physeq))
type <- substring(names(seq_dist), first = 1, last = 1)
seq_df <- data.frame("n_seqs" = seq_dist, "type" = type)
D <- seq_df

ggplot(aes(x=n_seqs), data=D) +
  geom_histogram(binwidth = 1500) + 
  xlab("Number of Sequencing Reads") + ylab("Number of Samples") + 
  ggtitle("Sequncing Reads Histogram") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "blank")

# Sequencing read scatterplot
ggplot(aes(x=1, y=n_seqs), data=D) +
  geom_jitter() + scale_y_log10() + 
  xlab("") + ylab("Number of Sequencing Reads") + 
  ggtitle("Sequencing Reads Scatterplot") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "blank")

# Sequencing reads cumulative plot
D %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs)) + 
  geom_line() + 
  xlab("Number of Samples") + ylab("Number of Sequencing Reads") + 
  ggtitle("Sequencing Reads Cumulative Plot") +
  theme(axis.title.y = element_text(size = 13, vjust=0.5), 
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "blank")

# Rarefying (One sample is to be removed)
physeq_rarefy <- rarefy_even_depth(physeq, sample.size = 48014, rngseed = 1,
                                   replace = FALSE)
