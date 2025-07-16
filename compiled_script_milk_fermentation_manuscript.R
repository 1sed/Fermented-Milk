######################### PHYLOSEQ OBJECT CREATION #################

library(phyloseq)
library(readxl)
library(ape)

# Set working directory
setwd

# 1. Import OTU table
# Load OTU table
otu <- read.table("feature-table.tsv",
                  header = TRUE,
                  sep = "\t",
                  row.names = 1,
                  skip = 1,               # Skip the "# Constructed from biom file" line
                  comment.char = "",     # Prevent "#" from removing header
                  check.names = FALSE)   # Keep sample names as-is (e.g. "sample1")

# Optional: Remove taxonomy column if it exists
if ("taxonomy" %in% colnames(otu)) {
  otu <- otu[, -which(colnames(otu) == "taxonomy")]
}

# Convert to matrix
otu_mat <- as.matrix(otu)

# Create phyloseq OTU table object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)


tax <- read.table("exported-taxonomy/taxonomy.tsv", 
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1, 
                  comment.char = "", 
                  quote = "", 
                  stringsAsFactors = FALSE)

# Split taxonomy into levels
tax_split <- strsplit(tax$Taxon, ";\\s*")
max_ranks <- max(sapply(tax_split, length))
tax_matrix <- t(sapply(tax_split, function(x) c(x, rep(NA, max_ranks - length(x)))))

# Assign biological rank names
colnames(tax_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:max_ranks]

# ❗ Assign original ASV hashes as row names
rownames(tax_matrix) <- rownames(tax)

# Create taxonomy table
TAX <- tax_table(tax_matrix)



# 3. Import tree
tree <- read_tree("exported-tree/tree.nwk")



# Step 1: Read metadata
metadata <- read_excel("metadata_milk.xlsx")

# Step 2: Convert to data frame
metadata <- as.data.frame(metadata)

# Step 3: Set row names to "SAMPLE ID"
rownames(metadata) <- metadata$`SAMPLE ID`

# Step 4: Remove the "SAMPLE ID" column now that it's in the rownames
metadata$`SAMPLE ID` <- NULL

# Step 5: Convert to phyloseq sample_data object
META <- sample_data(metadata)


# 7. Create phyloseq object
ps <- phyloseq(OTU, TAX, META, tree)



############ CORE MICROBIAL VISUALIZATION ##################################
setwd

#####load phyloseq object


ps <- readRDS("phyloseq_object.rds")



#############load relevant libraries 


library(microbiome)  # for transform()
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(stringr)


#####################PHYLUM#################################
####################PHYLUM##################################

# Normalize to relative abundance
physeq_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Agglomerate to Phylum level
physeq_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")

# Convert to long format for ggplot
df <- psmelt(physeq_phylum)

# Convert to percentage
df$Abundance <- df$Abundance * 100

# Optionally: round for clean labels
df$label <- ifelse(df$Abundance > 1, paste0(round(df$Abundance, 1), "%"), "")

# Plot
ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  facet_wrap(~ Product, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  labs(
    title = "Phylum-Level Composition by Product (with Percentages)",
    x = "Samples (Grouped by Product)",
    y = "Relative Abundance (%)"
  )


####################GENUS#######################
####################GENUS#######################

# Transform to relative abundance
ps_rel <- transform(ps, "compositional")

# Agglomerate at Genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# Melt to long format
df <- psmelt(ps_genus)

# Add Product info from sample_data
df$Product <- sample_data(ps_genus)$Product[match(df$Sample, rownames(sample_data(ps_genus)))]

# Identify top 10 genera by overall abundance
top10_genera <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Identify top 5 for labeling
top5_genera <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

# Filter and group
df_top10 <- df %>%
  mutate(Genus = as.character(Genus),
         Genus = ifelse(Genus %in% top10_genera, Genus, "Other")) %>%
  group_by(Product, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Percentage = Abundance * 100,
         Label = ifelse(Genus %in% top5_genera & Percentage > 1,
                        paste0(round(Percentage, 1), "%"), ""))

# Optional: wrap long product names
df_top10$Product <- str_wrap(df_top10$Product, width = 10)

# Plot with percentage labels for top 5 genera
ggplot(df_top10, aes(x = Product, y = Percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  theme_minimal() +
  labs(title = "Top 10 Genera by Product (with % for Top 5)",
       x = "Product", y = "Mean Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold")) +
  scale_fill_brewer(palette = "Set3")




#################SPECIES PLOT WITHOUT USING PHYLOSEQ OBJECT#############
#########################################################################


# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

# Step 1: Load Excel file
setwd("C:/Users/DAVID/Desktop/phyloseq/collapsed_tables_freq/")
data <- read_excel("taxa_freq_table7.xlsx")

# Step 2: Convert to long format
df_long <- data %>%
  pivot_longer(cols = -Species, names_to = "Product", values_to = "Abundance") %>%
  filter(!is.na(Abundance))

# Step 3: Calculate relative abundance (percentage) per product
df_long <- df_long %>%
  group_by(Product) %>%
  mutate(TotalProduct = sum(Abundance),
         Percentage = (Abundance / TotalProduct) * 100) %>%
  ungroup()

# Step 4: Identify top 10 species per product
top10_per_product <- df_long %>%
  group_by(Product) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE) %>%
  mutate(IsTop10 = TRUE)

# Step 5: Identify top 3 species per product for labeling
top3_labels <- df_long %>%
  group_by(Product) %>%
  slice_max(order_by = Abundance, n = 3, with_ties = FALSE) %>%
  select(Species, Product)

# Step 6: Keep only top 10 per product (plus “Other”) and label top 3
df_labeled <- df_long %>%
  left_join(top10_per_product %>% select(Species, Product, IsTop10),
            by = c("Species", "Product")) %>%
  mutate(Species = ifelse(is.na(IsTop10), "Other", Species)) %>%
  select(-IsTop10) %>%
  group_by(Product, Species) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  mutate(Label = ifelse(paste(Species, Product) %in%
                          paste(top3_labels$Species, top3_labels$Product),
                        paste0(round(Percentage, 1), "%"), "")) %>%
  filter(Species %in% top10_per_product$Species | Species == "Other")

# Step 7: Custom color palette (ensure enough colors for unique species)
custom_palette <- c(
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#fdbf6f", "#cab2d6", "#ff7f00", "#6a3d9a", "#1b9e77",
  "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#a6cee3", "#fb9a99", "#b15928",
  "#ffffb3", "#bebada"
)

# Step 8: Plot pie charts with facets and percentage labels
p <- ggplot(df_labeled, aes(x = "", y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") +
  facet_wrap(~ Product, scales = "free") +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Step 9: Save high-resolution image
ggsave("top_species_pie_facets.png", plot = p, width = 12, height = 8, dpi = 800, bg = "white")






################### ALPHA DIVERSITY PHYLUM LEVEL#############################

# Load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

# Set working directory and load phyloseq object
setwd("C:/Users/......")
ps <- readRDS("phyloseq_object.rds")

# ----------------------------------------------------
# 1. Agglomerate to Phylum level
# ----------------------------------------------------
ps_phylum <- tax_glom(ps, taxrank = "Phylum")

# ----------------------------------------------------
# 2. Estimate alpha diversity metrics
# ----------------------------------------------------
alpha_div <- estimate_richness(ps_phylum, measures = c("Observed", "Shannon", "Chao1"))
alpha_div$SampleID <- rownames(alpha_div)

# ----------------------------------------------------
# 3. Merge with Product info
# ----------------------------------------------------
meta <- data.frame(sample_data(ps_phylum))
meta$SampleID <- rownames(meta)
alpha_merged <- left_join(alpha_div, meta, by = "SampleID")

# ----------------------------------------------------
# 4. Save alpha diversity values to CSV
# ----------------------------------------------------
write.csv(alpha_merged, "phylum_level_alpha_diversity_values.csv", row.names = FALSE)
write.csv(alpha_merged, "phylum_level_alpha_diversity_values.csv", row.names = FALSE)
# ----------------------------------------------------
# 5. Perform Kruskal-Wallis tests
# ----------------------------------------------------
kw_observed <- kruskal.test(Observed ~ Product, data = alpha_merged)
kw_shannon  <- kruskal.test(Shannon ~ Product, data = alpha_merged)
kw_chao1    <- kruskal.test(Chao1 ~ Product, data = alpha_merged)

# Save test results
kw_results <- data.frame(
  Metric = c("Observed", "Shannon", "Chao1"),
  p_value = c(kw_observed$p.value, kw_shannon$p.value, kw_chao1$p.value)
)

write.csv(kw_results, "kruskal_test_results_phylum_level.csv", row.names = FALSE)

# ----------------------------------------------------
# 6. Plot with p-values
# ----------------------------------------------------
alpha_long <- melt(alpha_merged, id.vars = c("SampleID", "Product"),
                   measure.vars = c("Observed", "Shannon", "Chao1"))

colors <- RColorBrewer::brewer.pal(n = 5, name = "Set2")
names(colors) <- unique(alpha_long$Product)

ggplot(alpha_long, aes(x = Product, y = value, fill = Product)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Product), width = 0.2, alpha = 0.8, size = 4) +  # Larger dots
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Product (Phylum Level) with Kruskal-Wallis Tests",
       x = "Product", y = "Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)



################ALPHA DIVERSITY GENUS LEVEL#################################
############################################################################



# ----------------------------------------------------
# 1. Agglomerate to Genus level
# ----------------------------------------------------
ps_genus <- tax_glom(ps, taxrank = "Genus")

# ----------------------------------------------------
# 2. Estimate alpha diversity metrics
# ----------------------------------------------------
alpha_div <- estimate_richness(ps_genus, measures = c("Observed", "Shannon", "Chao1"))
alpha_div$SampleID <- rownames(alpha_div)

# ----------------------------------------------------
# 3. Merge with Product info
# ----------------------------------------------------
meta <- data.frame(sample_data(ps_genus))
meta$SampleID <- rownames(meta)
alpha_merged <- left_join(alpha_div, meta, by = "SampleID")

# ----------------------------------------------------
# 4. Save alpha diversity values to CSV
# ----------------------------------------------------
write.csv(alpha_merged, "genus_level_alpha_diversity_values.csv", row.names = FALSE)

# ----------------------------------------------------
# 5. Perform Kruskal-Wallis tests
# ----------------------------------------------------
kw_observed <- kruskal.test(Observed ~ Product, data = alpha_merged)
kw_shannon  <- kruskal.test(Shannon ~ Product, data = alpha_merged)
kw_chao1    <- kruskal.test(Chao1 ~ Product, data = alpha_merged)

# Save test results
kw_results <- data.frame(
  Metric = c("Observed", "Shannon", "Chao1"),
  p_value = c(kw_observed$p.value, kw_shannon$p.value, kw_chao1$p.value)
)

write.csv(kw_results, "kruskal_test_results_genus_level.csv", row.names = FALSE)

# ----------------------------------------------------
# 6. Plot with p-values
# ----------------------------------------------------
alpha_long <- melt(alpha_merged, id.vars = c("SampleID", "Product"),
                   measure.vars = c("Observed", "Shannon", "Chao1"))

colors <- RColorBrewer::brewer.pal(n = 5, name = "Set2")  # or any palette you prefer
names(colors) <- unique(alpha_long$Product)

ggplot(alpha_long, aes(x = Product, y = value, fill = Product)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Product), width = 0.2, alpha = 0.8, size = 4) +  # Increased dot size
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Product (Genus Level) with Kruskal-Wallis Tests",
       x = "Product", y = "Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)







#################DENDROGRAM COLOUR##############################################
#######script to extract unifract distnace######################################

library(phyloseq)
library(vegan)
library(dendextend)

# Set working directory
setwd("C:/Users/DAVID/Desktop/phyloseq/")

# Load phyloseq object
ps <- readRDS("phyloseq_object.rds")

# Clean object
ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

# Compute Weighted UniFrac distance
unifrac_dist <- distance(ps_clean, method = "unifrac", weighted = TRUE)

# ✅ Save UniFrac distance matrix as CSV
unifrac_matrix <- as.matrix(unifrac_dist)
write.csv(unifrac_matrix, file = "weighted_unifrac_distance_matrix.csv")
cat("✅ UniFrac distance matrix saved as 'weighted_unifrac_distance_matrix.csv'\n")

# Perform hierarchical clustering
hc <- hclust(as.dist(unifrac_dist), method = "average")

# Extract metadata
meta <- data.frame(sample_data(ps_clean))
product_labels <- as.character(meta$Product)
names(product_labels) <- rownames(meta)

# Convert hc to dendrogram
dend <- as.dendrogram(hc)

# Match product labels to dendrogram tip order
labels_in_dend <- labels(dend)
tip_labels <- product_labels[labels_in_dend]

# Create consistent color map
unique_products <- unique(tip_labels)
product_colors <- setNames(RColorBrewer::brewer.pal(length(unique_products), "Set2"), unique_products)

# Apply product labels and colors to dendrogram
dend <- dend %>%
  set("labels", tip_labels) %>%
  set("labels_colors", product_colors[tip_labels])  # exact color mapping

# Save plot
png("products_dendrogram_unifrac_named.png", width = 1000, height = 800)
par(cex = 1.8)

plot(dend,
     main = "Dendrogram (Weighted UniFrac) Labeled by Product",
     ylab = "Weighted UniFrac Distance")

legend("topright",
       legend = names(product_colors),
       col = product_colors,
       pch = 15, cex = 1.2)

dev.off()

cat("✅ Plot saved as 'products_dendrogram_unifrac_named.png'\n")


############################################################################
################################bray-curtis analysis########################

# Clean phyloseq object: remove empty samples and taxa
ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

# Compute Bray-Curtis distance matrix
dist_all <- distance(ps_clean, method = "bray")

bray_matrix <- as.matrix(dist_all)
write.csv(bray_matrix, file = "bray_curtis_distance_matrix.csv")
cat("✅ Bray distance matrix saved as 'bray_distance_matrix.csv'\n")


# Perform hierarchical clustering
hc <- hclust(as.dist(dist_all), method = "average")

# Extract sample metadata
sample_info <- data.frame(sample_data(ps_clean))
product_labels <- as.character(sample_info$Product)
names(product_labels) <- rownames(sample_info)

# Get labels in dendrogram order
labels_in_dend <- labels(as.dendrogram(hc))
tip_labels <- product_labels[labels_in_dend]

# Create consistent color mapping
unique_products <- unique(tip_labels)
product_colors <- setNames(RColorBrewer::brewer.pal(length(unique_products), "Set2"), unique_products)

# Apply labels and colors
dend <- as.dendrogram(hc)
dend <- dend %>%
  set("labels", tip_labels) %>%
  set("labels_colors", product_colors[tip_labels])

# Save dendrogram to PNG
png("products_dendrogram_bray_named.png", width = 1000, height = 800)
par(cex = 1.8)  # Increase font size of tip labels

plot(dend,
     main = "Dendrogram of Samples (Bray-Curtis) Labeled by Product",
     ylab = "Bray-Curtis Dissimilarity")

legend("topright",
       legend = names(product_colors),
       col = product_colors,
       pch = 15, cex = 1.2)  # Increase legend font size

dev.off()

cat("✅ Dendrogram saved as 'products_dendrogram_bray_named.png'\n")




##################FUNCTIONAL ANALYSIS###############################
#################set working directory##############################


setwd("C:/Users/..........")
library(readxl)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)

ps <- readRDS("phyloseq_object.rds")

# Read PICRUSt2 output files
ec_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/picrust_result_old/EC_metagenome_out/pred_metagenome_unstrat.xlsx")
ko_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/picrust_result_old/KO_metagenome_out/pred_metagenome_unstrat.xlsx")
pathways_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/picrust_result_old/pathways_out/path_abun_unstrat.xlsx")

# Convert from wide to long format for ggplot2

# Explicitly rename first column to "function"
colnames(ec_unstra)[1] <- "function"
colnames(ko_unstra)[1] <- "function"
colnames(pathways_unstra)[1] <- "function"

ec_unstra_long <- ec_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")

ko_unstra_long <- ko_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")

pathways_unstra_long <- pathways_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")


metadata <- data.frame(sample_data(ps)) %>%
  tibble::rownames_to_column("SampleID")
ko_unstra_joined <- left_join(ko_unstra_long, metadata, by = "SampleID")
ec_unstra_joined <- left_join(ec_unstra_long, metadata, by = "SampleID")
pathways_unstra_joined <- left_join(pathways_unstra_long, metadata, by = "SampleID")



################################################################################
###########PATHWAY COMPARATIVE ANALYSIS#########################################

# Get top 10 pathways
top_pathways_unstra <- pathways_unstra_long %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance)) %>%
  top_n(10, Total) %>%
  pull(`function`)

# Define custom colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Prepare data
plot_data <- pathways_unstra_long %>%
  filter(`function` %in% top_pathways_unstra) %>%
  left_join(metadata, by = "SampleID") %>%
  group_by(Product, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Create the plot
p <- ggplot(plot_data, aes(x = mean_abundance, y = fct_reorder(`function`, mean_abundance), fill = Product)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 15) +  # Slightly larger base size
  labs(
    title = "Top 10 Predicted Pathways by Product",
    x = "Mean Abundance",
    y = "Pathway",
    fill = "Product"
  ) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save as high-resolution PNG
ggsave("top10_pathways_barplot.png", plot = p, width = 10, height = 6, dpi = 400)

readr::write_csv(plot_data, paste0("pathway_plot_data.csv"))

# Optional: also display in R session
print(p)





#############heatmap_ko Comparative Analysis ####################
#############heatmap-ko Comparative Analysis ###################
library(dplyr)
library(tidyr)
library(forcats)
library(pheatmap)
library(RColorBrewer)

# Define custom product colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Ensure Product is a factor with desired order
ko_unstra_joined <- ko_unstra_joined %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Get top 10 most abundant KOs
top_kos <- ko_unstra_joined %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  top_n(10, Total) %>%
  pull(`function`)

# Filter for top 10 KOs
top_ko_data <- ko_unstra_joined %>%
  filter(`function` %in% top_kos)

# Pivot to KO × Sample matrix (mean across replicates)
heatmap_data <- top_ko_data %>%
  group_by(SampleID, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames(var = "function")

# Prepare annotation for columns (samples)
sample_metadata <- ko_unstra_joined %>%
  select(SampleID, Product) %>%
  distinct() %>%
  column_to_rownames("SampleID")

# Plot heatmap with clustering
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_col = sample_metadata,
  annotation_colors = list(Product = custom_colors),
  color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 11,
  
  main = "Heatmap of Top 10 Predicted KOs Across Samples"
)




# Step 1: Calculate total KO abundance per sample to compute relative abundance
relative_ko_data <- top_ko_data %>%
  group_by(SampleID) %>%
  mutate(TotalSampleAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(RelativeAbundance = Abundance / TotalSampleAbundance * 100)  # in percent

# Step 2: Calculate mean abundance and mean relative abundance per Product and KO
ko_summary <- relative_ko_data %>%
  group_by(Product, `function`) %>%
  summarise(
    mean_abundance = mean(Abundance),
    mean_relative_abundance = mean(RelativeAbundance),
    .groups = "drop"
  )

# Step 3: Save to CSV
write.csv(ko_summary, "mean_abundance_relative_abundance_top_kos.csv", row.names = FALSE)

# Optional: View part of the summary
head(ko_summary)



#############################################################################
############ EC EC HEATMAP Comparative Analysis##############################

# Load required packages
library(dplyr)
library(tidyr)
library(pheatmap)
library(readr)
library(RColorBrewer)
library(tibble)

# Define custom product colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Step 1: Join EC predictions with metadata
ec_unstra_joined <- left_join(ec_unstra_long, metadata, by = "SampleID") %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Step 2: Identify top 10 most abundant ECs
top_ecs <- ec_unstra_joined %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  top_n(10, Total) %>%
  pull(`function`)

# Step 3: Filter for top 10 ECs
top_ec_data <- ec_unstra_joined %>%
  filter(`function` %in% top_ecs)

# Step 4: Prepare wide matrix (EC x SampleID)
heatmap_data <- top_ec_data %>%
  group_by(`function`, SampleID) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames("function")

# Step 5: Prepare annotation data for samples
sample_metadata <- ec_unstra_joined %>%
  select(SampleID, Product) %>%
  distinct() %>%
  column_to_rownames("SampleID")

# Step 6: Save top ECs summary (optional)
top_ec_table <- ec_unstra_joined %>%
  filter(`function` %in% top_ecs) %>%
  group_by(Product, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

write_csv(top_ec_table, "top10_predicted_ECs_by_product.csv")

# Step 7: Plot heatmap with clustering and dendrogram
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_col = sample_metadata,
  annotation_colors = list(Product = custom_colors),
  scale = "row",  # normalize across rows (ECs)
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("white", "red"))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Heatmap of Top 10 Predicted EC Numbers Across Samples"
)


################################################################################
######################### SAMPLE SPECIFIC FUNCTIONAL ANALYSIS###################
################ GRID GRID KEGG GENES###########################################

# Load required libraries
library(pheatmap)
library(gridExtra)
library(grid)

# === Grab each heatmap as grob ===

# Sample 1 – Nono
heatmap_plot_1 <- grid.grabExpr(pheatmap(heatmap_matrix,
                                         main = "Top KO Functions in Nono",
                                         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 2 – Wara
heatmap_plot_2 <- grid.grabExpr(pheatmap(heatmap_matrix_2,
                                         main = "Top KO Functions in Wara",
                                         color = colorRampPalette(c("black", "yellow", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 3 – Kwerionik
heatmap_plot_3 <- grid.grabExpr(pheatmap(heatmap_matrix_3,
                                         main = "Top KO Functions in Kwerionik",
                                         color = colorRampPalette(c("pink", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 4 – Ghee
heatmap_plot_4 <- grid.grabExpr(pheatmap(heatmap_matrix_4,
                                         main = "Top KO Functions in Ghee",
                                         color = colorRampPalette(c("black", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 5 – Nunu
heatmap_plot_5 <- grid.grabExpr(pheatmap(heatmap_matrix_5,
                                         main = "Top KO Functions in Nunu",
                                         color = colorRampPalette(c("blue", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Create an empty grob for slot 6 (optional)
empty_plot <- grid.rect(gp = gpar(col = NA))  # Transparent placeholder


# Combine plots 1, 2, and 5 into a single frame
frame_1 <- grid.arrange(
  heatmap_plot_1,  # Nono
  heatmap_plot_2,  # Wara
  heatmap_plot_5,  # Nunu
  ncol = 3
)


# Combine plots 3 and 4 into another frame
frame_2 <- grid.arrange(
  heatmap_plot_3,  # Kwerionik
  heatmap_plot_4,  # Ghee
  ncol = 2
)



################################################################################
################################## grid grid Enzyme Commission##################


heatmap_ec_plot_1 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec1,
           main = "Top EC Functions in Nono",
           color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_2 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec2,
           main = "Top EC Functions in Wara",
           color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_3 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec3,
           main = "Top EC Functions in Kwerionik",
           color = colorRampPalette(c("navy", "white", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_4 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec4,
           main = "Top EC Functions in Ghee",
           color = colorRampPalette(c("black", "yellow", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)

heatmap_ec_plot_5 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec5,
           main = "Top EC Functions in Nunu",
           color = colorRampPalette(c("orange", "white", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)

# Frame 1
grid.arrange(heatmap_ec_plot_1, heatmap_ec_plot_2, heatmap_ec_plot_5, ncol = 3)

# Frame 2
grid.arrange(heatmap_ec_plot_3, heatmap_ec_plot_4, ncol = 2)



################################################################################
########################### grid pathway #######################################
# Load required libraries
library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)

# === Create all 5 pathway heatmaps and capture gtables ===

# Sample 1 - Nono
p1 <- pheatmap(heatmap_matrix_pw1,
               main = "Top Pathways in Nono",
               color = colorRampPalette(brewer.pal(9, "BuPu"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g1 <- p1$gtable

# Sample 2 - Wara
p2 <- pheatmap(heatmap_matrix_pw2,
               main = "Top Pathways in Wara",
               color = colorRampPalette(brewer.pal(9, "Reds"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g2 <- p2$gtable

# Sample 3 - Kwerionik
p3 <- pheatmap(heatmap_matrix_pw3,
               main = "Top Pathways in Kwerionik",
               color = colorRampPalette(brewer.pal(9, "Purples"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g3 <- p3$gtable

# Sample 4 - Ghee
p4 <- pheatmap(heatmap_matrix_pw4,
               main = "Top Pathways in Ghee",
               color = colorRampPalette(brewer.pal(9, "Oranges"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g4 <- p4$gtable

# Sample 5 - Nunu
p5 <- pheatmap(heatmap_matrix_pw5,
               main = "Top Pathways in Nunu",
               color = colorRampPalette(brewer.pal(9, "RdPu"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g5 <- p5$gtable

# === Arrange into frames ===

# Frame 1: Nono, Wara, Nunu
frame1 <- grid.arrange(g1, g2, g5, ncol = 3)

# Frame 2: Kwerionik, Ghee
frame2 <- grid.arrange(g3, g4, ncol = 2)

# === Optional: Save to file ===
 ggsave("frame1_pathways.png", plot = frame1, width = 18, height = 6, dpi = 300)
 ggsave("frame2_pathways.png", plot = frame2, width = 12, height = 6, dpi = 300)



























