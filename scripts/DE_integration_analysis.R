###### Script to run differential expression and integration analysis ####

# packages needed
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
library(DESeq2)
library('biomaRt')
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Column data with the condition corresponding to each ID in count matrix
condition <- as.data.frame(cbind(c("",colnames(raw_counts)),
                                 c("condition",rep("CTRL",length(grep("control",colnames(raw_counts)))),
                                   rep("miRNA_KO",length(grep("miRNA",colnames(raw_counts)))))))
col.condition <- condition
colnames(col.condition) <- as.character(unlist(col.condition[1,]))
col.condition <- col.condition[-1,]
rownames(col.condition) <- col.condition[,1]
col.condition[,1] <- NULL

# -> Creation of DESeq2 Data Matrix
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = col.condition,
                              design = ~ condition)

# Eliminate the zero counts from mRNA data (pre-filter)
mRNA.count.tokeep <- rowSums(counts(dds)) >= 10
dds.new <- dds[mRNA.count.tokeep,]

# Factor levels
dds.new$condition <- factor(dds.new$condition,
                                        levels = c("CTRL","miRNA_KO"))

# DE analysis (miRNA KO vs CTRL)
dds.new <- DESeq(dds.new)
# default: p-value 0.1
res <- results(dds.new)

# convert ENSEMBL IDs to Gene Symbols
results_matrix <- res %>% as.data.frame() %>% mutate(ensembl_gene_id = row.names(.))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(results_matrix)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=genes,mart= mart)
results_matrix <- merge(results_matrix,G_list,by = "ensembl_gene_id")
results_matrix$hgnc_symbol <- ifelse(results_matrix$hgnc_symbol == "SLC11A2","DMT1",results_matrix$hgnc_symbol)

# genes and miRNA of interest 
genes_interest <- c("CHAC1", "LPCAT3", "NCOA4","DMT1","ACSL4","GPX4")
miRNA_interest <- "hsa-miR-940"

writexl::write_xlsx(results_matrix,"RNA-seq data/DE_RNA-seq.xlsx")

results_matrix <- results_matrix %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.01 ~ "up",
                               log2FoldChange <= (-1) & padj <= 0.01 ~ "down",
                               TRUE ~ "ns")) 
results_gene_interes <- results_matrix %>%
  filter(hgnc_symbol %in% genes_interest)

# volcano plot of transcriptomics
cols <- c("up" = "#e10f98", "down" = "#2378ba", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
ggplot(data = results_matrix,
       aes(x = log2FoldChange,
           y = -log10(padj))) + 
  geom_point(aes(colour = gene_type), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = results_gene_interes,
             shape = 21,
             size = 2,aes(fill = gene_type)) +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed") + 
  annotate(geom="text", x=-7.2, y=-log10(0.01) + 8, label="FDR = 1%") +
  geom_vline(xintercept = c(- 1, 1),
             linetype = "dashed") +
  geom_label_repel(data = results_gene_interes,
                   aes(label = hgnc_symbol),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),     
                     limits = c(-8, 8)) +
  labs(title = "Transcriptomics differential expression miRNA-940 vs control",
       x = "log2(fold change)",
       y = "-log10(adjusted p-value)",
       colour = "Differential \nExpression") +
  theme_classic() + # Select theme with a white background  
  guides(fill=FALSE)

# volcano plot of proteomics
proteins_stats <- proteins_stats %>% 
  rename(log2FoldChange = `log2 ratio miRNA940 / Control`,
         padj = `LFQ Intensity ttest miRNA/Control`) %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.01 ~ "up",
                               log2FoldChange <= (-1) & padj <= 0.01 ~ "down",
                               TRUE ~ "ns")) 
results_protein_interest <- proteins_stats %>%
  filter(`Gene Name` %in% genes_interest)

ggplot(data = proteins_stats,
       aes(x = log2FoldChange,
           y = -log10(padj))) + 
  geom_point(aes(colour = gene_type), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = results_protein_interest,
             shape = 21,
             size = 2,aes(fill = gene_type)) +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed") + 
  annotate(geom="text", x=-7.2, y=-log10(0.01) + 0.5, label="FDR = 1%") +
  geom_vline(xintercept = c(- 1, 1),
             linetype = "dashed") +
  geom_label_repel(data = results_protein_interest,
                   aes(label = `Gene Name`),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),     
                     limits = c(-8, 8)) +
  labs(title = "Proteomics differential expression miRNA-940 vs control",
       x = "log2(fold change)",
       y = "-log10(adjusted p-value)",
       colour = "Differential \nExpression") +
  theme_classic() + # Select theme with a white background  
  guides(fill=FALSE)


# boxplots of genes of interest
ensembl.interest <- results_matrix %>%
  filter(hgnc_symbol %in% genes_interest) %>%
  pull(ensembl_gene_id)
counts <- data.frame(assay(dds.new)[rownames(assay(dds.new)) %in% ensembl.interest,]) 
counts <- counts %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  left_join(results_matrix %>%
              select(hgnc_symbol,ensembl_gene_id))
m2 <- data.frame(ACSL4 = t(counts[1,1:8]),
                 NCOA4 = t(counts[2,1:8]),
                 LPCAT3 = t(counts[3,1:8]),
                 DMT1 = t(counts[4,1:8]),
                 CHAC1 = t(counts[5,1:8]),
                 GPX4 = t(counts[6,1:8]),
                 group = as.factor(col.condition$condition))
colnames(m2) <- c("ACSL4","NCOA4","LPCAT3","DMT1","CHAC1","GPX4","group")

# boxplot function
ACSL4 <- ggplot(m2,aes(group,ACSL4,fill = group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0))+
  labs(x="",y = "Normalized Counts", title = "ACSL4") +
  theme_classic() +
  scale_x_discrete(labels = c("CTRL","miRNA")) +
  scale_fill_manual(values = c("#2378ba","#e10f98"), breaks = c("CTRL","miRNA_KO"),labels = c("CTRL","miRNA"))

ggarrange(ACSL4,
          NCOA4,
          LPCAT3,
          DMT1,
          CHAC1,
          GPX4,
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2)

# create final dataset with miRNA stats, gene stats and protein stats datasets
miRNA_final <- miRNA_stats %>%
  rename(miRNA =  `Gene symbol`,
         FDR_miRNA = `FDR(BH)`,
         log2FC_miRNA = `log2(Fold change)`) %>%
  select(miRNA,FDR_miRNA,log2FC_miRNA)

gene_final <- results_matrix %>%
  rename(Gene = hgnc_symbol,
         FDR_gene = padj,
         log2FC_gene = log2FoldChange) %>%
  select(Gene,FDR_gene,log2FC_gene)

protein_final <- proteins_stats %>%
  rename(Gene = `Gene Name`,
         FDR_protein = padj,
         log2FC_protein = log2FoldChange) %>%
  select(Gene,FDR_protein,log2FC_protein)

# data with wppi
mirTarBase_PPI_top10 <- mirTarBase_PPI_top10 %>%
  select(source,target) %>%
  mutate(source = toupper(source)) %>%
  mutate(dataset_int = rep("miRTarBase",nrow(.)))

TargetScan_PPI_top10 <- TargetScan_PPI_top10 %>%
  select(source,target) %>%
  mutate(source = toupper(source)) %>%
  mutate(dataset_int = rep("TargetScan",nrow(.)))

miRNA_final$miRNA <- toupper(miRNA_final$miRNA)

mirTarBase_PPI_top10 <- merge(as.data.frame(mirTarBase_PPI_top10),as.data.frame(miRNA_final),
      by.x = "source", by.y = "miRNA",
      all.x = TRUE)

TargetScan_PPI_top10 <- merge(as.data.frame(TargetScan_PPI_top10),as.data.frame(miRNA_final),
                              by.x = "source", by.y = "miRNA",
                              all.x = TRUE)

# all together
Gene_miRNA_integration_PPI <- rbind(mirTarBase_PPI_top10,
                                    TargetScan_PPI_top10) %>%
  mutate(source = ifelse(source == "SLC11A2","DMT1",source), # different protein name 
         target = ifelse(target == "SLC11A2","DMT1",target))

Gene_miRNA_integration_PPI <- merge(Gene_miRNA_integration_PPI,gene_final,
      by.x = "source", by.y = "Gene",
      all.x = TRUE) %>%
  rename(FDR_gene_source = FDR_gene,
         log2FC_gene_source = log2FC_gene) 

Gene_miRNA_integration_PPI <- merge(Gene_miRNA_integration_PPI,protein_final,
                                    by.x = "source", by.y = "Gene",
                                    all.x = TRUE) %>%
  rename(FDR_protein_source = FDR_protein,
         log2FC_protein_source = log2FC_protein) 

Gene_miRNA_integration_PPI <- merge(Gene_miRNA_integration_PPI,gene_final,
                                    by.x = "target", by.y = "Gene",
                                    all.x = TRUE) %>%
  rename(FDR_gene_target = FDR_gene,
         log2FC_gene_target = log2FC_gene) 

Gene_miRNA_integration_PPI <- merge(Gene_miRNA_integration_PPI,protein_final,
                                    by.x = "target", by.y = "Gene",
                                    all.x = TRUE) %>%
  rename(FDR_protein_target = FDR_protein,
         log2FC_protein_target = log2FC_protein) 

writexl::write_xlsx(Gene_miRNA_integration_PPI,"Gene_miRNA_integration_PPI_all.xlsx")

# remove genes or proteins which are not significant, i.e. FDR < 1% and |log2FC| >= 1, and focus only on miRNA direct
Gene_miRNA_integration_PPI_miRNA_rows <- Gene_miRNA_integration_PPI %>%
  filter(grepl("HSA-MIR",source)) %>%
  filter((FDR_gene_target < 0.01 & abs(log2FC_gene_target) >= 1) | (FDR_protein_target < 0.01 & abs(log2FC_protein_target) >= 1)) %>%
  mutate(gene_interest = ifelse(target %in% genes_interest,"interest","no interest"),
         miRNA_interest = ifelse(source %in% toupper(miRNA_interest),"interest","no interest")) %>%
  mutate(gene_DE = ifelse(log2FC_gene_target >= 1, "up", "down")) %>% 
  distinct()

writexl::write_xlsx(Gene_miRNA_integration_PPI_miRNA_rows,"integration data/Gene_miRNA_PPI_direct.xlsx")

# only significant and with miRNA of interest 
Gene_miRNA_interest_integration_PPI_miRNA_significant <- Gene_miRNA_integration_PPI_miRNA_rows %>%
  dplyr::filter(miRNA_interest == "interest") %>%
  mutate(gene_DE = ifelse(log2FC_gene_target >= 1, "up", "down"))

writexl::write_xlsx(Gene_miRNA_interest_integration_PPI_miRNA_significant,"Gene_miRNA_interest_PPI_integration.xlsx")
