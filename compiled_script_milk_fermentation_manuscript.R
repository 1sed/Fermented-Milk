######################### PHYLOSEQ OBJECT CREATION #################

library(phyloseq)
library(readxl)
library(ape)

# Set working directory where all your files are located
setwd()

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


                       
#---------------------------------------------------------------------------
############ CORE MICROBIAL VISUALIZATION ##################################
#---------------------------------------------------------------------------
#####load phyloseq object

ps <- readRDS("phyloseq_object.rds")

#############load relevant libraries 
#----------------------------------------------------------------------------
#----------------------script genus plot-------------------------------------
#----------------------------------------------------------------------------
library(phyloseq)
library(microbiome)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

# 1. Transform to relative abundance
ps_rel <- microbiome::transform(ps, "compositional")

# 2. Agglomerate at Genus level for relative abundance
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# 3. Melt relative abundance table
df <- psmelt(ps_genus)

# 4. Clean genus names and add Product
df <- df %>%
  mutate(
    Genus = as.character(Genus),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
    Product = sample_data(ps_genus)$Product[
      match(Sample, rownames(sample_data(ps_genus)))
    ]
  )

# 5. Identify overall top 10 genera across all samples
top10_genera <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# 6. Build sample-level relative abundance table
df_sample_abundance <- df %>%
  mutate(Genus = ifelse(Genus %in% top10_genera, Genus, "Other")) %>%
  group_by(Sample, Product, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 7. Agglomerate ORIGINAL phyloseq object at Genus level for raw counts
ps_genus_counts <- tax_glom(ps, taxrank = "Genus")

# 8. Melt raw count table
df_counts <- psmelt(ps_genus_counts)

# 9. Clean genus names and add Product
df_counts <- df_counts %>%
  mutate(
    Genus = as.character(Genus),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
    Product = sample_data(ps_genus_counts)$Product[
      match(Sample, rownames(sample_data(ps_genus_counts)))
    ]
  )

# 10. Build sample-level raw count table using the same top10 grouping
df_sample_counts <- df_counts %>%
  mutate(Genus = ifelse(Genus %in% top10_genera, Genus, "Other")) %>%
  group_by(Sample, Product, Genus) %>%
  summarise(Count = sum(Abundance), .groups = "drop")

# 11. Merge sample-level relative abundance and raw counts
df_sample_merged <- df_sample_abundance %>%
  left_join(df_sample_counts, by = c("Sample", "Product", "Genus"))

# 12. Final mean table including counts
df_top10 <- df_sample_merged %>%
  group_by(Product, Genus) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    MeanCount = mean(Count),
    TotalCount = sum(Count),
    N_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(Percentage = MeanAbundance * 100)

# 13. Create labels for top 5 genera WITHIN each product
df_top10 <- df_top10 %>%
  group_by(Product) %>%
  arrange(desc(Percentage), .by_group = TRUE) %>%
  mutate(
    RankWithinProduct = row_number(),
    Label = ifelse(
      RankWithinProduct <= 5 & Percentage > 1,
      paste0(round(Percentage, 1), "%"),
      ""
    )
  ) %>%
  ungroup()

# 14. Optional: order products
product_order <- df_top10 %>%
  group_by(Product) %>%
  summarise(Total = sum(Percentage), .groups = "drop") %>%
  pull(Product)

df_top10 <- df_top10 %>%
  mutate(
    Product = factor(Product, levels = unique(product_order)),
    Product_wrapped = str_wrap(as.character(Product), width = 10)
  )

# 15. Optional: order fill legend by overall abundance
genus_order <- df_top10 %>%
  group_by(Genus) %>%
  summarise(Total = sum(Percentage), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

df_top10 <- df_top10 %>%
  mutate(Genus = factor(Genus, levels = genus_order))

# 16. Plot
gn <- ggplot(df_top10, aes(x = Product_wrapped, y = Percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "Top 10 Genera by Product",
    subtitle = "Labels show top 5 genera within each product",
    x = "Product",
    y = "Mean Relative Abundance (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "Set3")

print(gn)

ggsave(
  "Genus_Level_Composition_Top10_Top5Labels.svg",
  plot = gn,
  width = 10,
  height = 8,
  dpi = 1200,
  bg = "white"
)

# 17. Save final mean table with counts
write.csv(df_top10, "top_10_taxa_by_genus_with_counts.csv", row.names = FALSE)

# 18. Optional: save sample-level table actually used for calculating the means
write.csv(df_sample_merged, "top_10_taxa_sample_level_counts_and_abundance.csv", row.names = FALSE)



#-------------------------------------------------------------------------------
############################### ALPHA DIVERSITY#################################
#-------------------------------------------------------------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# ----------------------------------------------------
# 1. Agglomerate at Genus level
# ----------------------------------------------------
ps_genus <- tax_glom(ps, taxrank = "Genus")

# ----------------------------------------------------
# 2. Estimate alpha diversity metrics
# ----------------------------------------------------
alpha_div <- estimate_richness(
  ps_genus,
  measures = c("Observed", "Shannon", "Simpson", "Chao1")
)

alpha_div$SampleID <- rownames(alpha_div)

# ----------------------------------------------------
# 3. Merge with Product metadata
# ----------------------------------------------------
meta <- data.frame(sample_data(ps_genus))
meta$SampleID <- rownames(meta)

alpha_merged <- left_join(alpha_div, meta, by = "SampleID")

# ----------------------------------------------------
# 4. Save alpha diversity values to CSV
# ----------------------------------------------------
write.csv(alpha_merged, "genus_level_alpha_diversity_values.csv", row.names = FALSE)

# ----------------------------------------------------
# 5. Reshape for plotting
#    Main indices for manuscript emphasis: Shannon, Simpson
# ----------------------------------------------------
alpha_main <- alpha_merged %>%
  select(SampleID, Product, Shannon, Simpson) %>%
  pivot_longer(
    cols = c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  )

# Optional supporting richness metrics
alpha_richness <- alpha_merged %>%
  select(SampleID, Product, Observed, Chao1) %>%
  pivot_longer(
    cols = c(Observed, Chao1),
    names_to = "Metric",
    values_to = "Value"
  )

# ----------------------------------------------------
# 6. Define colors
# ----------------------------------------------------
products <- unique(alpha_merged$Product)
colors <- brewer.pal(n = max(3, min(length(products), 8)), name = "Set2")
colors <- colors[seq_along(products)]
names(colors) <- products

# ----------------------------------------------------
# 7. Plot Shannon and Simpson (recommended main figure)
#    Use points because there is only one sample per product
# ----------------------------------------------------
p_main <- ggplot(alpha_main, aes(x = Product, y = Value, color = Product)) +
  geom_point(size = 5) +
  geom_text(
    aes(label = round(Value, 2)),
    color = "black",
    size = 5,
    fontface = "bold",
    vjust = -1.2
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Genus-Level Alpha Diversity Across Fermented Dairy Products",
    subtitle = "Descriptive diversity indices (single sample per product)",
    x = "Product",
    y = "Diversity index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) +
  #scale_color_manual(values = colors)
  scale_fill_brewer(palette = "Set3")

print(p_main)

ggsave(
  filename = "alpha_diversity_shannon_simpson_genus_level2.png",
  plot = p_main,
  width = 10,
  height = 8,
  dpi = 1200
)


#-----------------------------------------------------------------------                       
##################FUNCTIONAL ANALYSIS###################################
#################set working directory##################################
#-----------------------------------------------------------------------

setwd("C:/Users/..........")
library(readxl)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)

ps <- readRDS("phyloseq_object.rds")

# Read PICRUSt2 output files
ec_unstra <- read_excel("C:/Users/provide the path to the file on your computer")
ko_unstra <- read_excel("C:/Users/provide the path to the file on your computer")
pathways_unstra <- read_excel("C:/Users/provide the path to the file on your computer")

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



























