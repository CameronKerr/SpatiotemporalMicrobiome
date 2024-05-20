#############################
# Data Loading & Processing #
#############################

# Loading packages
require(stringr)
require(phyloseq)
require(microbiomeutilities)
require(readr)
require(qiime2R)
require(vegan)
require(pheatmap)
require(grDevices)
require(tidyr)

# Loading data
tax_metadata <- read.table(file = "metadata.tsv", sep="\t", header = TRUE)
table_output <- read.delim(file = "table-medb30.tsv", skip = 1)
tax <- merge(tax_metadata, table_output, by.x = "Feature.ID", by.y = "X.OTU.ID")

new_columns <- str_split_fixed(as.character(tax$Taxon), ";", 7)
colnames(new_columns) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat <- cbind(new_columns, tax)

qiimegen_tree <- qza_to_phyloseq(tree = "rooted-tree-medb30.qza")

# Creating phyloseq packages

## Sample data
type <- substring(colnames(taxmat[12:ncol(taxmat)]), first = 1, last = 1)
location <- c()
for(i in 1:94){location <- c(location, strsplit(colnames(taxmat[12:ncol(taxmat)]), "_")[[i]][2])}
timepoint <- parse_number(colnames(taxmat[12:ncol(taxmat)]))
samples <- colnames(taxmat[12:ncol(taxmat)])
sampdata <- cbind(type, location, timepoint, samples)
rownames(sampdata) <- colnames(taxmat[12:ncol(taxmat)])
SAM <- sample_data(data.frame(sampdata, stringsAsFactors=FALSE))

rural <- c("KSR", "Preston", "Terracotta")
suburban <- c("Aviemore", "Rouge", "BruceH")
urban <- c("Cedarvale", "Grenadier", "Humber")
location_type <- c()
for(i in 1:length(SAM$location)){
  if(SAM$location[i] %in% rural){
    location_type[i] <- "rural"
  }
  else if(SAM$location[i] %in% suburban){
    location_type[i] <- "suburban"
  }
  else if(SAM$location[i] %in% urban){
    location_type[i] <- "urban"
  }
}
SAM$location_type <- location_type

culture <- c()
for(i in 1:length(SAM$samples)){
  if(SAM$type[i] == "C"){
    culture[i] <- strsplit(SAM$samples[i], "_")[[1]][3]
  }
  else {
    culture[i] <- "D"
  }
}
SAM$culture <- culture

## Taxonomy table
rownames(taxmat) <- tax$Feature.ID
taxmat <- as.matrix(taxmat)
colnames(taxmat) <- c("Domain", colnames(taxmat)[2:ncol(taxmat)])
TAX = tax_table(taxmat)

## OTU table
otumat <- cbind(tax[12:length(tax)])
rownames(otumat) = tax$Feature.ID
OTU = otu_table(otumat, taxa_are_rows = TRUE)

physeq = phyloseq(OTU, TAX, SAM, qiimegen_tree)

