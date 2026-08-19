# Inferred Functional Potential of Microbial Communities in Traditional African Fermented Dairy Products Using 16S rRNA

## This repository contains the scripts, workflows, and analysis pipelines used for the manuscript investigating the microbial diversity and predicted functional potential of traditionally fermented dairy products using 16S rRNA amplicon sequencing, QIIME 2, PICRUSt2, and phyloseq.

Repositroy structure
├── qiime_script.txt
├── picrust_script.txt
├── compiled_script_milk_fermentation_manuscript.R
├── metadata/
├── exported_table/
├── exported_taxonomy/
├── exported-tree/
├── picrust2_out/
├── figures/
├── results/
└── README.md


## Overview of Workflow
The analysis pipeline consists of:
* Raw sequence import and preprocessing
* Denoising and ASV generation
* Taxonomic classification
* Phylogenetic tree construction
* Construction of phyloseq object
* Diversity analyses
* Taxonomic visualization
* Functional prediction using PICRUSt2
* KEGG Orthology (KO), Enzyme Commission (EC), and pathway analyses
* Heatmap and clustering visualization



## Software and Dependencies
### QIIME2 Environment
### Recommended version:
* QIIME 2 2024.10
* Python ≥ 3.9
* Required QIIME2 plugins
* dada2
* feature-classifier
* phylogeny
* diversity


### PICRUSt2 Environment
### Recommended version:
* PICRUSt2
* Python ≥ 3.8

### R Environment
### Recommended R version:
* R ≥ 4.2
```
# Required R packages:
install.packages(c(
  "phyloseq",
  "microbiome",
  "ggplot2",
  "dplyr",
  "tidyr",
  "reshape2",
  "ggpubr",
  "readxl",
  "ape",
  "vegan",
  "dendextend",
  "RColorBrewer",
  "pheatmap",
  "gridExtra",
  "grid",
  "forcats",
  "scales",
  "tibble",
  "readr"
))
```
## Step-by-Step Workflow

```
# Import Raw Reads into QIIME2
# Import paired-end FASTQ files using a manifest file
qiime tools import \
--type "SampleData[PairedEndSequencesWithQuality]" \
--input-format PairedEndFastqManifestPhred33V2 \
--input-path manifest-paired-end.tsv \
--output-path demuxed-dss.qza
# Reference workflow:qiime_script
```
```
# Visualize Demultiplexed Reads
qiime demux summarize \
--i-data demuxed-dss.qza \
--o-visualization demuxed-dss.qza.qzv
```
```
# Denoising Using DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# Reference workflow:qiime_script
# Expected Outputs:
# 1. ASV table
# 2. Representative sequences
# 3. Denoising statistics
```
```
# Taxonomic Classification
# Extract SILVA reference reads
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 250 \
  --p-n-jobs 8 \
  --o-reads silva-138-99-515F-806R.qza
```

```
# Train Classifier 
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-515F-806R.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-515F-806R-classifier.qza
```
```
#Assign Taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515F-806R-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
# Reference workflow: qiime_script
```
```
# Construct Phylogenetic Tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```
```
# Export Files for Downstream R Analysis
# Export taxonomy
qiime tools export \
    --input-path taxonomy.qza \
    --output-path exported_taxonomy
```

```
# Export feature table
qiime tools export \
    --input-path table.qza \
    --output-path exported_table
```
```
# Convert BIOM to TSV
biom convert \
    -i exported_table/feature-table.biom \
    -o exported_table/feature-table.tsv \
    --to-tsv

# Reference workflow:qiime_script
# Expected outputs:
# 1. taxonomy.tsv
# 2. feature-table.tsv
# 3. tree.nwk
# 4. metadata file
```
## Functional Prediction Using PICRUSt2
### Activate PICRUSt2 environment:
* source ~/miniconda3/etc/profile.d/conda.sh
* conda activate picrust2

```
# Run pipeline:
picrust2_pipeline.py \
-s dna-sequences.fasta \
-i feature-table.tsv \
-o picrust2_out \
-p 4 \
--stratified

# Reference workflow:picrust_script.txt
# Outputs include:
# 1. KO predictions
# 2. EC predictions
# 3. Pathway abundance tables
```

## Create Phyloseq Object in R
The phyloseq object was generated using:
* feature table
* taxonomy table
* rooted phylogenetic tree
* sample metadata

Reference Core script:compiled_script_milk_fermentation

Main output:
* phyloseq_object.rds

Generate Taxonomic Composition Figures
Phylum-level bar plots
Script section:
Relative abundance normalization
Taxonomic agglomeration
Percentage labeling

Reference:compiled_script_milk_fermentation

Output:
* Phylum-Level Composition by Product

Genus-level stacked bar plots
Reference:compiled_script_milk_fermentation

Output:
* Top 10 Genera by Product

Species-level pie charts
Reference:compiled_script_milk_fermentation

Output:
* Top_species_pie_facets.png

Alpha Diversity Analysis
Genus-level diversity

Metrics:
* Shannon
* Simpson


Reference:compiled_script_milk_fermentation
Outputs:
* alpha diversity CSV files
* diversity plots

Functional Analysis and Visualization
Pathway analysis
Output:
* Top10_pathways_barplot.png
* KO heatmaps
* EC heatmaps
* Reference:compiled_script_milk_fermentation

Sample-Specific Functional Heatmaps
The repository includes:
* KO heatmaps
* EC heatmaps
* pathway heatmaps
Generated using:
* pheatmap
* gridExtra
* grid graphics
* KO heatmaps
* EC heatmaps
* pathway heatmaps


Notes and Limitations
Functional predictions generated by PICRUSt2 represent inferred functional potential and not direct measurements of metabolic activity.
Taxonomic resolution is limited by 16S rRNA amplicon sequencing.
Analyses are descriptive due to the absence of biological replicates.

Citation

If using this repository, please cite:

Contact
For questions regarding the analysis workflow or scripts, please contact the corresponding author.
