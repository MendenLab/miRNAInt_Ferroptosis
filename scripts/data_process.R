## Data download and pre-process

# RNA-seq data
library(readr)
raw_counts <- read_csv("RNA-seq data/raw_counts.csv") # raw RNA-Seq data
rownames <- raw_counts$...1
raw_counts <- raw_counts[,-1]
rownames(raw_counts) <- rownames

# miRNA data
library(readxl)
miRNA_stats <- read_excel("Data/DE_miRNA.xlsx")

# proteomics data
library(readxl)
proteins_stats <- read_excel("Data/DE_proteomics.xlsx", sheet = "Volcano LFQ intensity")

# Integration datasets
# -> mirTarBase
mirTarBase_PPI_top10 <- read_excel("Data/mirTarBase_PPI_top10.xlsx")

# -> TargetScan
TargetScan_PPI_top10 <- read_excel("Data/TargetScan_PPI_top10.xlsx")
