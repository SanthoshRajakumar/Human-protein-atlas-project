library(readxl)
library(dplyr)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(tibble)
library(writexl)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)


file_path <- "C:/Users/adith/Desktop/hpa/rt/HPA fetal tissues.xlsx" 
data <- read_excel(file_path)
View(data)
colnames(data)

df <- data

df <- df %>% dplyr::select(ensg_id, tissue, age, sample_comment, ntpm)

View(df)
# Creating a column based on tissue and sample_comment
df <- df %>%
  mutate(dynamic_column = paste0("fetal_", tissue, "_", gsub(",|\\s+", "_", sample_comment), "_ntpm"))

# Reshaping to wide format
df_wide <- df %>%
  dplyr::select(ensg_id, dynamic_column, ntpm) %>%
  pivot_wider(names_from = dynamic_column, values_from = ntpm)
View(df_wide)
colnames(df_wide)

# numeric
df_wide <- df_wide %>%
  mutate(across(starts_with("fetal_"), as.numeric))

# Calculate the averages for different tissues
df_avg <- df_wide %>%
  mutate(
    # Calculate the average for fetal kidney columns
    fetal_kidney_avg = rowMeans(dplyr::select(df_wide, starts_with("fetal_fetal kidney_")), na.rm = TRUE),
    
    # Calculate the average for fetal lung columns
    fetal_lung_avg = rowMeans(dplyr::select(df_wide, starts_with("fetal_fetal lung_")), na.rm = TRUE),
    
    # Calculate the average for fetal CNS columns
    fetal_CNS_avg = rowMeans(dplyr::select(df_wide, starts_with("fetal_fetal CNS_")), na.rm = TRUE)
  )


View(df_avg)

#file path of adult data
file_path_1 <- "C:/Users/adith/Desktop/adult_data.xlsx" 
data_1 <- read_excel(file_path_1)
View(data_1)
colnames(data_1)

#changing the column name
colnames(data_1)[colnames(data_1) == "adult_id"] <- "ensg_id"
View(data_1)


# Merging the two data frames using left join based on ensg_id
merged_data <- df_avg %>%
  left_join(data_1, by = "ensg_id")

View(merged_data)
colnames(merged_data)

#adding a column based on the ratio 
comparison_data <- merged_data %>%
  mutate(
    ratio_kidney = as.numeric(fetal_kidney_avg) / as.numeric(`Consensus Kidney (nTPM)`),
    ratio_lung = as.numeric(fetal_lung_avg) / as.numeric(`Consensus Lung (nTPM)`),
    ratio_CNS = as.numeric(fetal_CNS_avg) / as.numeric(`Consensus Brain (nTPM)`)  # Assuming Consensus Brain for CNS comparison
  )
View(comparison_data)
colnames(comparison_data)





write.csv(comparison_data, "restructured data.csv", row.names = FALSE)




#threshold
threshold <- 2  

# Filter significant genes for all tissue types
significant_genes_all <- comparison_data %>%
  filter(ratio_kidney > threshold | ratio_lung > threshold | ratio_CNS > threshold) %>%
  dplyr::select(ensg_id, `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
         `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`, 
         `fetal_fetal lung_Fetal__Lung_12028_14.6_ntpm`, 
         `fetal_fetal lung_Fetal__Lung_8836_17.4_ntpm`,
         `fetal_fetal CNS_Fetal__CNS_8597_10_ntpm`, 
         `fetal_fetal CNS_Fetal__CNS_8350_13_ntpm`,
         `fetal_fetal CNS_Fetal__CNS_8838_17.4_ntpm`, 
         `fetal_kidney_avg`, `fetal_lung_avg`, `fetal_CNS_avg`,
         `Consensus Kidney (nTPM)`, `Consensus Brain (nTPM)`, `Consensus Lung (nTPM)`,
         `Consensus Distribution`, `Consensus Max`,
         ratio_kidney, ratio_lung, ratio_CNS)
View(significant_genes_all)
colnames(significant_genes_all)



# Filter significant genes for kidney tissue only
significant_genes_kidney <- comparison_data %>%
  filter(ratio_kidney > threshold) %>%
  dplyr::select(ensg_id, 
         `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
         `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`,  
         `fetal_kidney_avg`, 
         `Consensus Kidney (nTPM)`, 
         `Consensus Distribution`, 
         `Consensus Max`, 
         `ratio_kidney`)
View(significant_genes_kidney)

significant_genes_kidney_sorted <- significant_genes_kidney %>%
  arrange(desc(ratio_kidney))

View(significant_genes_kidney_sorted)

# Filter rows where at least one of the fetal kidney nTPM values is >= 1
significant_genes_kidney_1 <- significant_genes_kidney_sorted %>%
  filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` >= 1 | 
           `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` >= 1)



View(significant_genes_kidney_1)
summary(significant_genes_kidney_1)


# Filter significant genes for lung tissue only
significant_genes_lung <- comparison_data %>%
  filter(ratio_lung > threshold) %>%
  dplyr::select(ensg_id, 
         `fetal_fetal lung_Fetal__Lung_12028_14.6_ntpm`, 
         `fetal_fetal lung_Fetal__Lung_8836_17.4_ntpm`,  
         `fetal_lung_avg`, 
         `Consensus Lung (nTPM)`, 
         `Consensus Distribution`, 
         `Consensus Max`, 
         `ratio_lung`)
View(significant_genes_lung)

significant_genes_lung_sorted <- significant_genes_lung %>%
  arrange(desc(ratio_lung))


View(significant_genes_lung_sorted)


# Filter rows where at least one of the fetal kidney nTPM values is >= 1
significant_genes_lung_1 <- significant_genes_lung_sorted %>%
  filter(`fetal_fetal lung_Fetal__Lung_12028_14.6_ntpm` >= 1 | 
           `fetal_fetal lung_Fetal__Lung_8836_17.4_ntpm` >= 1)


View(significant_genes_lung_1)


# Filter significant genes for cns tissue only
significant_genes_cns <- comparison_data %>%
  filter(ratio_CNS > threshold) %>%
  dplyr::select(ensg_id, 
         `fetal_fetal CNS_Fetal__CNS_8597_10_ntpm`, 
         `fetal_fetal CNS_Fetal__CNS_8350_13_ntpm`,
         `fetal_fetal CNS_Fetal__CNS_8838_17.4_ntpm`,  
         `fetal_CNS_avg`, 
         `Consensus Brain (nTPM)`, 
         `Consensus Distribution`, 
         `Consensus Max`, 
         `ratio_CNS`)
View(significant_genes_cns)

significant_genes_cns_sorted <- significant_genes_cns %>%
  arrange(desc(ratio_CNS))

View(significant_genes_cns_sorted)


# Filter rows where at least one of the fetal cns nTPM values is >= 1
significant_genes_cns_1 <- significant_genes_cns_sorted %>%
  filter(`fetal_fetal CNS_Fetal__CNS_8597_10_ntpm` >= 1 | 
           `fetal_fetal CNS_Fetal__CNS_8350_13_ntpm` >= 1 | 
         `fetal_fetal CNS_Fetal__CNS_8838_17.4_ntpm` >= 1)


View(significant_genes_cns_1)




write.csv(significant_genes_kidney_1, "significant_genes_kidney_1.csv", row.names = FALSE)
write.csv(significant_genes_cns_1, "significant_genes_cns_1.csv", row.names = FALSE)
write.csv(significant_genes_lung_1, "significant_genes_lung_1.csv", row.names = FALSE)


# numeric
significant_genes_kidney_1 <- significant_genes_kidney_1 %>%
  mutate(`Consensus Kidney (nTPM)` = as.numeric(`Consensus Kidney (nTPM)`))
significant_genes_cns_1 <- significant_genes_cns_1 %>%
  mutate(`Consensus Brain (nTPM)` = as.numeric(`Consensus Brain (nTPM)`))
significant_genes_lung_1 <- significant_genes_lung_1 %>%
  mutate(`Consensus Lung (nTPM)` = as.numeric(`Consensus Lung (nTPM)`))
# Load ggplot2 for plotting
library(ggplot2)

# Scatter plot for Consensus Kidney (nTPM) vs fetal kidney nTPM values
plot <- ggplot(significant_genes_kidney_1, aes(x = `Consensus Kidney (nTPM)`)) +
  geom_point(aes(y = `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`, color = "Fetal Kidney week 15"), alpha = 0.6) +
  geom_point(aes(y = `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`, color = "Fetal Kidney week 18"), alpha = 0.6) +
  scale_color_manual(values = c("Fetal Kidney week 15" = "blue", "Fetal Kidney week 18" = "red")) +
  labs(
    x = "Consensus Kidney (nTPM)",
    y = "Fetal Kidney nTPM",
    color = "Sample Type",
    title = "Scatter Plot of Adult kidney nTPM vs Fetal Kidney nTPM"
  ) +
  theme_minimal()

# Print the plot
print(plot)

# Load ggplot2 for plotting
library(ggplot2)

# Scatter plot for Consensus Kidney (nTPM) vs fetal kidney nTPM values
plot <- ggplot(significant_genes_kidney_1, aes(x = `Consensus Kidney (nTPM)`)) +
  geom_point(aes(y = `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`, color = "Fetal Kidney week 15"), alpha = 0.6) +
  geom_point(aes(y = `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`, color = "Fetal Kidney week 18"), alpha = 0.6) +
  geom_point(aes(y = `Consensus Kidney (nTPM)`, color = "Adult Kidney"), alpha = 0.6) +
  scale_color_manual(values = c(
    "Fetal Kidney week 15" = "blue", 
    "Fetal Kidney week 18" = "red", 
    "Adult Kidney" = "green"
  )) +
  labs(
    x = "Adult Kidney (nTPM)",
    y = "Fetal Kidney (nTPM)",
    color = "Sample Type",
    title = "Scatter Plot of Adult Kidney (nTPM) vs Fetal Kidney (nTPM) "
  ) +
  theme_minimal()

# Print the plot
print(plot)


#scatter plot 
ggplot(significant_genes_kidney_1, aes(x = fetal_kidney_avg, y = `Consensus Kidney (nTPM)`)) +
  geom_point(alpha = 0.6) +
  scale_fill_viridis_c() +
  geom_density_2d() +
  theme_minimal() +
  labs(title = "Scatter Plot of Fetal vs Consensus Kidney (nTPM)",
       x = "Fetal Kidney Avg (nTPM)",
       y = "Consensus Kidney (nTPM)")

ggplot(significant_genes_cns_1, aes(x = fetal_CNS_avg, y = `Consensus Brain (nTPM)`)) +
  geom_point(alpha = 0.6) +
  scale_fill_viridis_c() +
  geom_density_2d() +
  theme_minimal() +
  labs(title = "Scatter Plot of Fetal vs Consensus CNS (nTPM)",
       x = "Fetal cns Avg (nTPM)",
       y = "Consensus Brain (nTPM)")

ggplot(significant_genes_lung_1, aes(x = fetal_lung_avg, y = `Consensus Lung (nTPM)`)) +
  geom_point(alpha = 0.6) +
  scale_fill_viridis_c() +
  geom_density_2d() +
  theme_minimal() +
  labs(title = "Scatter Plot of Fetal vs Consensus Lung (nTPM)",
       x = "Fetal lung Avg (nTPM)",
       y = "Consensus lung (nTPM)")




# top 200 genes with the highest ratio_kidney
top_200_genes_kidney <- significant_genes_kidney_1 %>%
  top_n(200, wt = ratio_kidney) %>%
  arrange(desc(ratio_kidney)) 
top_200_genes_lung <- significant_genes_lung_1 %>%
  top_n(200, wt = ratio_lung) %>%
  arrange(desc(ratio_lung)) 
top_200_genes_cns <- significant_genes_cns_1 %>%
  top_n(200, wt = ratio_CNS) %>%
  arrange(desc(ratio_CNS)) 

View(top_200_genes_kidney)
View(top_200_genes_lung)
View(top_200_genes_cns)

write.csv(top_200_genes_kidney, "top_200_genes_kidney.csv", row.names = FALSE)
write.csv(top_200_genes_lung, "top_200_genes_lung.csv", row.names = FALSE)
write.csv(top_200_genes_cns, "top_200_genes_cns.csv", row.names = FALSE)

#Gene ontology


#KIDNEY
list_gene_kidney <- top_200_genes_kidney$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "ALL",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)

View(top_200_genes_kidney)
kidney200_go <- as.data.frame(go)
View(kidney200_go)
colnames(kidney200_go)

# dot plot
ggplot(kidney200_go, aes(x = FoldEnrichment, 
                      y = reorder(Description, FoldEnrichment),  # Reorder terms by enrichment
                      size = Count,                             # Dot size by gene count
                      color = p.adjust)) +                      # Dot color by adjusted p-value
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +             # Blue (significant) to red (less significant)
  labs(title = "GO Enrichment Analysis (Top 200 kidney Genes-ALL)", 
       x = "Fold Enrichment", 
       y = "GO Term", 
       size = "Gene Count", 
       color = "Adj. p-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# top 20 GO terms with the highest number of expressed genes
top_20_high_expression_kidney <- kidney200_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_20_high_expression_kidney, aes(x = Count, 
                                   y = reorder(Description, Count),   # Reorder by gene count
                                   size = FoldEnrichment,             # Dot size by Fold Enrichment
                                   color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY- ALL= 200", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




#BP 200
list_gene_kidney <- top_200_genes_kidney$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


View(top_200_genes_kidney)
kidney200_go <- as.data.frame(go)
View(kidney200_go)
colnames(kidney200_go)

# dot plot
ggplot(kidney200_go, aes(x = FoldEnrichment, 
                         y = reorder(Description, FoldEnrichment),  # Reorder terms by enrichment
                         size = Count,                             # Dot size by gene count
                         color = p.adjust)) +                      # Dot color by adjusted p-value
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +             # Blue (significant) to red (less significant)
  labs(title = "GO Enrichment Analysis (Top 200 kidney Genes- BP 200)", 
       x = "Fold Enrichment", 
       y = "GO Term", 
       size = "Gene Count", 
       color = "Adj. p-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# top 20 GO terms with the highest number of expressed genes
top_20_high_expression_kidney <- kidney200_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_20_high_expression_kidney, aes(x = Count, 
                                          y = reorder(Description, Count),   # Reorder by gene count
                                          size = FoldEnrichment,             # Dot size by Fold Enrichment
                                          color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY- BP 200", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))





#MF 200
list_gene_kidney <- top_200_genes_kidney$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "MF",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)

View(top_200_genes_kidney)
kidney200_go <- as.data.frame(go)
View(kidney200_go)
colnames(kidney200_go)

# dot plot
ggplot(kidney200_go, aes(x = FoldEnrichment, 
                         y = reorder(Description, FoldEnrichment),  # Reorder terms by enrichment
                         size = Count,                             # Dot size by gene count
                         color = p.adjust)) +                      # Dot color by adjusted p-value
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +             # Blue (significant) to red (less significant)
  labs(title = "GO Enrichment Analysis (Top 200 kidney Genes- MF 200)", 
       x = "Fold Enrichment", 
       y = "GO Term", 
       size = "Gene Count", 
       color = "Adj. p-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# top 20 GO terms with the highest number of expressed genes
top_20_high_expression_kidney <- kidney200_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_20_high_expression_kidney, aes(x = Count, 
                                          y = reorder(Description, Count),   # Reorder by gene count
                                          size = FoldEnrichment,             # Dot size by Fold Enrichment
                                          color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY- MF 200", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



#for 2000 genes
list_gene_kidney_2000 <- significant_genes_kidney_1$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney_2000,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "ALL",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


kidney2000_go <- as.data.frame(go)
View(kidney2000_go)
colnames(kidney2000_go)
write.csv(kidney2000_go, "kidney2000_go.csv", row.names = FALSE)


# top 20 GO terms with the highest number of expressed genes
top_2000_high_expression_kidney <- kidney2000_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_2000_high_expression_kidney, aes(x = Count, 
                                          y = reorder(Description, Count),   # Reorder by gene count
                                          size = FoldEnrichment,             # Dot size by Fold Enrichment
                                          color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY (total= ALL)", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#BP- 2000
list_gene_kidney_2000 <- significant_genes_kidney_1$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney_2000,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


kidney2000_go <- as.data.frame(go)
View(kidney2000_go)
colnames(kidney2000_go)


# top 20 GO terms with the highest number of expressed genes
top_2000_high_expression_kidney <- kidney2000_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_2000_high_expression_kidney, aes(x = Count, 
                                            y = reorder(Description, Count),   # Reorder by gene count
                                            size = FoldEnrichment,             # Dot size by Fold Enrichment
                                            color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY (total= BP)", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#MF- 2000
list_gene_kidney_2000 <- significant_genes_kidney_1$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_kidney_2000,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "MF",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


kidney2000_go <- as.data.frame(go)
View(kidney2000_go)
colnames(kidney2000_go)


# top 20 GO terms with the highest number of expressed genes
top_2000_high_expression_kidney <- kidney2000_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_2000_high_expression_kidney, aes(x = Count, 
                                            y = reorder(Description, Count),   # Reorder by gene count
                                            size = FoldEnrichment,             # Dot size by Fold Enrichment
                                            color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms KIDNEY (total= MF)", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




#CNS
list_gene_cns <- top_200_genes_cns$ensg_id 
# Perform GO enrichment analysis (Biological Process - BP)
go <- enrichGO(
  gene          = list_gene_cns,   # Your ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "ALL",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)

View(top_200_genes_cns)
cns200_go <- as.data.frame(go)
View(cns200_go)
colnames(cns200_go)

# dot plot
ggplot(cns200_go, aes(x = FoldEnrichment, 
                      y = reorder(Description, FoldEnrichment),  # Reorder terms by enrichment
                      size = Count,                             # Dot size by gene count
                      color = p.adjust)) +                      # Dot color by adjusted p-value
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +             # Blue (significant) to red (less significant)
  labs(title = "GO Enrichment Analysis (Top 200 CNS Genes)", 
       x = "Fold Enrichment", 
       y = "GO Term", 
       size = "Gene Count", 
       color = "Adj. p-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# top 20 GO terms with the highest number of expressed genes
top_20_high_expression <- cns200_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_20_high_expression, aes(x = Count, 
                                   y = reorder(Description, Count),   # Reorder by gene count
                                   size = FoldEnrichment,             # Dot size by Fold Enrichment
                                   color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms CNS", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# LUNGS
list_gene_lung <- top_200_genes_lung$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = list_gene_lung,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "ALL",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)

View(top_200_genes_lung)
lung200_go <- as.data.frame(go)
View(lung200_go)
colnames(lung200_go)

# dot plot
ggplot(lung200_go, aes(x = FoldEnrichment, 
                         y = reorder(Description, FoldEnrichment),  # Reorder terms by enrichment
                         size = Count,                             # Dot size by gene count
                         color = p.adjust)) +                      # Dot color by adjusted p-value
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") +             # Blue (significant) to red (less significant)
  labs(title = "GO Enrichment Analysis (Top 200 lung Genes)", 
       x = "Fold Enrichment", 
       y = "GO Term", 
       size = "Gene Count", 
       color = "Adj. p-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# top 20 GO terms with the highest number of expressed genes
top_20_high_expression_lung <- lung200_go %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)

# Create a dot plot for terms with higher expressed genes
ggplot(top_20_high_expression_lung, aes(x = Count, 
                                          y = reorder(Description, Count),   # Reorder by gene count
                                          size = FoldEnrichment,             # Dot size by Fold Enrichment
                                          color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "Top 20 GO Terms LUNG", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




#DESEQ
#library(DESeq2)
#View(merged_data)
#colnames(merged_data)

# Convert only valid numeric Consensus columns
#merged_data <- merged_data %>%
#  mutate(across(c(`Consensus Kidney (nTPM)`, 
                  #`Consensus Brain (nTPM)`, 
                  #`Consensus Lung (nTPM)`), as.numeric))

# library(dplyr)
# colnames(merged_data)
# head(merged_data$`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`)
# 
# 
# count_matrix <- merged_data %>%
#   select(ensg_id,
#          `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
#          `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`,
#          `Consensus Kidney (nTPM)`) %>%
#   column_to_rownames(var = "ensg_id")
# 
# colnames(count_matrix) <- c("fetal_1", "fetal_2", "adult")
# 
# # Ensure row names in coldata match column names of count_matrix
# rownames(coldata) <- coldata$sample
# 
# 
# str(count_matrix)
# 
# # Create DESeqDataSet
# dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_matrix),
#                               colData = coldata,
#                               design = ~ condition)
# 
# dds <- DESeq(dds)
# results <- results(dds)
# plotMA(results, main = "MA Plot", ylim = c(-2, 2))
# sig_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)
# 
# 
# 
# 
# # Create sample information based on the updated count_data
# num_fetal <- length(grep("^fetal_", colnames(merged_data)))
# num_adult <- length(grep("^Consensus", colnames(merged_data)))
# 
# sample_info <- data.frame(
#   condition = c(rep("fetal", num_fetal), rep("adult", num_adult))
# )
# rownames(sample_info) <- colnames(merged_data)
# 
# # Check the sample info
# print(sample_info)
# 
# dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_data)),
#                               colData = sample_info,
#                               design = ~ condition)
# View(dds)
# 
# # Run DESeq2 analysis
# dds <- DESeq(dds)
# 
# # Extract and view results
# results_dds <- results(dds)
# summary(results_dds)
# print(results_dds)
# 
# # View top differentially expressed genes
# head(results_dds[order(results_dds$padj), ])
# # MA plot
# plotMA(results_dds, ylim = c(-5, 5))
# library(ggplot2)
# volcano_data <- as.data.frame(results_dds)
# ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(color = padj < 0.1), alpha = 0.6) +
#   scale_color_manual(values = c("gray", "red")) +
#   theme_minimal() +
#   labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
# 
# library(ggplot2)
# 
# # Remove rows with NA in log2FoldChange or padj
# results_dds_clean <- results_dds[!is.na(results_dds$log2FoldChange) & !is.na(results_dds$padj), ]
# print(results_dds_clean)
# 
# # Calculate fold change from log2FoldChange
# results_dds_clean$foldChange <- 2^results_dds_clean$log2FoldChange
# 
# # Filter for upregulated genes with Fold Change > 2 and padj < 0.1
# upregulated_genes_fetal_vs_adult <- results_dds_clean[results_dds_clean$foldChange > 2 & results_dds_clean$padj < 0.1, ]
# 
# # View the top upregulated genes
# print(upregulated_genes_fetal_vs_adult)
# 
# 
# ggplot(results_dds_clean, aes(x = foldChange, y = -log10(padj))) +
#   geom_point(aes(color = foldChange > 2 & padj < 0.1), alpha = 0.6) +
#   scale_color_manual(values = c("gray", "red")) +
#   theme_minimal() +
#   xlim(-5, 5) +  # Adjust based on your data distribution
#   ylim(0, 50) +  # Adjust if needed
#   labs(title = "Volcano Plot (Fold Change > 2 and padj < 0.1)", 
#        x = "Log2 Fold Change", 
#        y = "-Log10 Adjusted P-value")
# 
# 
# # Write the upregulated genes to a CSV file
# write.csv(upregulated_genes_fetal_vs_adult, "upregulated_genes_fetal_vs_adult.csv")
# 
# colnames(upregulated_genes_fetal_vs_adult)
# print(upregulated_genes_fetal_vs_adult)


file_path_3 <- "C:/Users/adith/Desktop/ape.xlsx" 
data3 <- read_excel(file_path_3)
View(data3)
colnames(data3)

colnames(data3)[colnames(data3) == "ENSG"] <- "ensg_id" 
ape <- comparison_data %>%
  left_join(data3, by = "ensg_id")

# View the merged data
View(ape)
colnames(ape)


View(significant_genes_kidney_1)
ape_kidney <- significant_genes_kidney_1 %>%
  left_join(data3, by = "ensg_id")
ape_cns <- significant_genes_cns_1 %>%
  left_join(data3, by = "ensg_id")
ape_lung <- significant_genes_lung_1 %>%
  left_join(data3, by = "ensg_id")

View(significant_genes_kidney_1)


View(ape_kidney)

write.csv(ape_kidney, "ape_kidney.csv", row.names = FALSE)
write.csv(ape_cns, "ape_cns.csv", row.names = FALSE)
write.csv(ape_lung, "ape_lung.csv", row.names = FALSE)

library(dplyr)

# Select only the specified columns and store them in kidney_data
kidney_data <- comparison_data %>%
  dplyr::select(
    ensg_id,
    `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
    `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`,
    fetal_kidney_avg,
    `Consensus Kidney (nTPM)`,
    `Consensus Distribution`,
    `Consensus Max`,
    ratio_kidney
  )

# View the new data frame
View(kidney_data)

kidneysdata_ape <- kidney_data %>%
  left_join(data3, by = "ensg_id")


View(kidneysdata_ape)

F15_F18 <- kidneysdata_ape  %>%
  mutate(F15_F18 = `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` / `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`)
View(F15_F18)
summary(F15_F18)



master_kidney_data <-  F15_F18 %>%
  mutate(F18_F15 = `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` / `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`)
View(master_kidney_data)
colnames(master_kidney_data)


master_kidney_data <- master_kidney_data %>%
  mutate(
    F15_high_expression = ifelse(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` > 
                                   `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`, 
                                 "Yes", "No"),
    F18_high_expression = ifelse(`fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` > 
                                   `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`, 
                                 "Yes", "No")
  )

# View the updated dataframe
View(master_kidney_data)
colnames(master_kidney_data)
head(master_kidney_data)

write.csv(master_kidney_data, "master_kidney_data.csv", row.names = FALSE)


file_path_4 <- "C:/Users/adith/Desktop/hpa/rt/literaturelist.xlsx" 
data4 <- read_excel(file_path_4)
View(data4)
summary(data4)
colnames(data4)

colnames(data4)[colnames(data4) == "ENSG"] <- "ensg_id" 

# Merge data4 with master_kidney_data based on the "ensg_id" column
literaturelistmasterdata <- merge(
  x = data4, 
  y = master_kidney_data, 
  by = "ensg_id", 
  all.x = TRUE
)



View(literaturelistmasterdata )

colnames(literaturelistmasterdata)
head(literaturelistmasterdata)



literaturelistmasterdata <- literaturelistmasterdata %>%
  mutate(
    `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` = as.numeric(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`),
    `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` = as.numeric(`fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`),
    `Consensus Kidney (nTPM)` = as.numeric(`Consensus Kidney (nTPM)`)
  )


# Filter rows where both F15 and F18 are greater than the adult value
fetal_higher_than_adult <- literaturelistmasterdata %>%
  filter(
    `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` > `Consensus Kidney (nTPM)`  |
      `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` > `Consensus Kidney (nTPM)`
  )

# View the filtered data
View(fetal_higher_than_adult)
colnames(fetal_higher_than_adult)
write.csv(fetal_higher_than_adult, "literaature_masterdata.csv", row.names = FALSE)


#young
threshold <- 2 
kidney_young <- master_kidney_data %>%
  filter(F15_F18 > threshold) %>%
  dplyr::select(ensg_id, 
         `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
         `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`,  
         `fetal_kidney_avg`, 
         `Consensus Kidney (nTPM)`, 
         `Consensus Distribution`, 
         `Consensus Max`, 
         `ratio_kidney`,
         `APE reliability score`,
         `F15_F18`,
         `F18_F15`)
View(kidney_young)

# kidney_young_sorted <- significant_genes_kidney %>%
#   arrange(desc(ratio_kidney))
# 
# View(significant_genes_kidney_sorted)

# Filter rows where at least one of the fetal kidney nTPM values is >= 1
kidney_young_1 <- kidney_young %>%
  filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` >= 1 | 
           `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` >= 1)



View(kidney_young_1)






threshold <- 2 
kidney_old <- master_kidney_data %>%
  filter(F18_F15 > threshold) %>%
  dplyr::select(ensg_id, 
         `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`,
         `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`,  
         `fetal_kidney_avg`, 
         `Consensus Kidney (nTPM)`, 
         `Consensus Distribution`, 
         `Consensus Max`, 
         `ratio_kidney`,
         `APE reliability score`,
         `F15_F18`,
         `F18_F15`)
View(kidney_old)

# kidney_young_sorted <- significant_genes_kidney %>%
#   arrange(desc(ratio_kidney))
# 
# View(significant_genes_kidney_sorted)

# Filter rows where at least one of the fetal kidney nTPM values is >= 1
kidney_old_1 <- kidney_old %>%
  filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` >= 1 | 
           `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` >= 1)



View(kidney_old_1)







# overexpressed_genes <- ape_not_between %>%
#   filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` > `Consensus Kidney (nTPM)`)
# 
# 
# View(overexpressed_genes)
# 
# 
# overexpressed_genes_1 <- overexpressed_genes %>%
#   filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` >= 1)
# 
# 
# View(overexpressed_genes_1)
# 
# 
# write.csv(overexpressed_genes_1, "ape_overexpressed_genes_1.csv", row.names = FALSE)
# 


kidney_young_list<- kidney_young$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = kidney_young_list,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)



kidney_young_go <- as.data.frame(go)
View(kidney_young_go)
head(kidney_young_go)

# top 20 GO terms with the highest number of expressed genes
kidney_young_go_TOP <- kidney_young_go  %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:50)


# Create a dot plot for terms with higher expressed genes
ggplot(kidney_young_go_TOP, aes(x = Count, 
                                          y = reorder(Description, Count),   # Reorder by gene count
                                          size = FoldEnrichment,             # Dot size by Fold Enrichment
                                          color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "GO Enrichment Analysis (young)", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



ggplot(kidney_young_go_TOP, aes(x = RichFactor, 
                                y = reorder(Description, RichFactor),   # Reorder by RichFactor
                                color = FoldEnrichment)) +             # Color by Fold Enrichment
  geom_segment(aes(x = 0, xend = RichFactor, 
                   y = reorder(Description, RichFactor), 
                   yend = reorder(Description, RichFactor)), 
               color = "gray") +                                      # Add line segment for lollipop
  geom_point(size = 4) +                                              # Add points at the end of each line
  scale_color_gradient(low = "blue", high = "red") +                  # Blue (low) to red (high)
  labs(title = "GO Enrichment Analysis (Young) - Lollipop Plot", 
       x = "RichFactor", 
       y = "GO Term", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



#old

kidney_old_list<- kidney_old$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          = kidney_old_list,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


kidney_old_go <- as.data.frame(go)
View(kidney_old_go)

# top 20 GO terms with the highest number of expressed genes
kidney_old_go_TOP <- kidney_old_go  %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:50)


# Create a dot plot for terms with higher expressed genes
ggplot(kidney_old_go_TOP, aes(x = Count, 
                                y = reorder(Description, Count),   # Reorder by gene count
                                size = FoldEnrichment,             # Dot size by Fold Enrichment
                                color = FoldEnrichment)) +         # Dot color by Fold Enrichment
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +                # Green (low) to red (high)
  labs(title = "GO Enrichment Analysis (old)", 
       x = "Gene Count", 
       y = "GO Term", 
       size = "Fold Enrichment", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



#YOUNG
# Filtering rows for the GO terms and extract geneID info
gene_ids <- kidney_young_go %>%
  filter(ID %in% c("GO:0072001", "GO:0001822")) %>%
  dplyr::select(ID, geneID)

# View the resulting gene IDs for the specified GO terms
print(gene_ids)


# Extract gene lists for the GO terms
genes_GO_0072001 <- unlist(strsplit(kidney_young_go["GO:0072001", "geneID"], "/"))
genes_GO_0001822 <- unlist(strsplit(kidney_young_go["GO:0001822", "geneID"], "/"))

# Find common genes
common_genes <- intersect(genes_GO_0072001, genes_GO_0001822)

# Find non-common genes
non_common_genes_GO_0072001 <- setdiff(genes_GO_0072001, genes_GO_0001822)
non_common_genes_GO_0001822 <- setdiff(genes_GO_0001822, genes_GO_0072001)

# Print results
cat("Common Genes:\n", common_genes, "\n\n")
cat("Non-Common Genes in GO:0072001:\n", non_common_genes_GO_0072001, "\n\n")
cat("Non-Common Genes in GO:0001822:\n", non_common_genes_GO_0001822, "\n\n")




common_genes_df <- data.frame(
  ENSG = common_genes
)

colnames(common_genes_df)[colnames(common_genes_df ) == "ENSG"] <- "ensg_id" 

View(common_genes_df)
colnames(common_genes_df)
# Check for matches with master data
matched_genes <- common_genes_df %>% 
  inner_join(master_kidney_data, by = "ensg_id")


view(matched_genes)
# Display results
print("Matched Genes:")
print(matched_genes)


write.csv(matched_genes, "matched_genes.csv", row.names = FALSE)


# Find common genes between matched_genes and literaturelistmasterdata
common_genes_in_both <- matched_genes %>%
  inner_join(fetal_higher_than_adult, by = "ensg_id")


View(common_genes_in_both)
colnames(common_genes_in_both)

write.csv(common_genes_in_both, "common_genes_in_both.csv", row.names = FALSE)

cleaned_common_genes <- common_genes_in_both %>%
  dplyr::select(
    ensg_id, 
    Gene, 
    `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm.x`,
    `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm.x`,
    fetal_kidney_avg.x, 
    `Consensus Kidney (nTPM).x`, 
    `Consensus Distribution.x`,
    `Consensus Max.x`, 
    ratio_kidney.x,
    `APE reliability score.x`, 
    F15_F18.x, 
    F18_F15.x,
    F15_high_expression.x, 
    F18_high_expression.x
  )

# Rename columns for clarity (if desired)
colnames(cleaned_common_genes) <- sub("\\.x$", "", colnames(cleaned_common_genes))

# View the cleaned data
View(cleaned_common_genes)

write.csv(cleaned_common_genes, "cleaned_common_genes.csv", row.names = FALSE)



common_genes_in_both_list<- common_genes_in_both$ensg_id 
# GO enrichment analysis 
go <- enrichGO(
  gene          =common_genes_in_both_list,   # ENSEMBL gene list
  OrgDb         = org.Hs.eg.db,        # Human annotation database
  keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
  ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
  pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
  qvalueCutoff = 0.2 ,               # q-value cutoff
  readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
)


common_genes_in_both_go <- as.data.frame(go)
View(common_genes_in_both_go)

# top 20 GO terms with the highest number of expressed genes
common_genes_in_both_go_TOP <- common_genes_in_both_go  %>%
  arrange(desc(Count)) %>%   # Sorting by gene count higher count = higher expression
  slice(1:20)
View(common_genes_in_both_go_TOP)
head(common_genes_in_both_go_TOP)


ggplot(common_genes_in_both_go_TOP, aes(x = RichFactor, 
                                y = reorder(Description, RichFactor),   # Reorder by RichFactor
                                color = FoldEnrichment)) +             # Color by Fold Enrichment
  geom_segment(aes(x = 0, xend = RichFactor, 
                   y = reorder(Description, RichFactor), 
                   yend = reorder(Description, RichFactor)), 
               color = "gray") +                                      # Add line segment for lollipop
  geom_point(size = 4) +                                              # Add points at the end of each line
  scale_color_gradient(low = "blue", high = "red") +                  # Blue (low) to red (high)
  labs(title = "GO Enrichment Analysis (Young) -common_genes_in_both_go", 
       x = "RichFactor", 
       y = "GO Term", 
       color = "Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))






# Load necessary library
library(dplyr)
library(readxl)

# Load the Excel data
data5 <- read_excel("C:/Users/adith/Downloads/HPA DATA.xlsx", sheet = "raw data from literature")
View(data5)
# Group by ENSG and Gene, combine Cell types
result <- data5 %>%
  group_by(ENSG, Gene) %>%
  summarise(Cell_type = paste(unique(`Cell type`), collapse = ", "), .groups = "drop")

# View the result
View(result)
colnames(result)

# Optionally write to a new Excel file
write.csv(result, "processed_ENSG_Gene_CellType.csv", row.names = FALSE)


merged_data_egc <- master_kidney_data %>%
  left_join(result %>% dplyr::select(ENSG, Gene, Cell_type), by = c("ensg_id" = "ENSG"))

# View the resulting dataframe
  View(merged_data_egc)
  
  
  merged_list<- merged_data_egc$ensg_id 
  # GO enrichment analysis 
  go <- enrichGO(
    gene          = merged_list,   # ENSEMBL gene list
    OrgDb         = org.Hs.eg.db,        # Human annotation database
    keyType       = "ENSEMBL",           # Specify ENSEMBL ID type
    ont           = "BP",                # Ontology: BP (Biological Process), CC (Cellular Component), or MF (Molecular Function)
    pAdjustMethod = "BH",               # Adjust p-values using Benjamini-Hochberg (BH) method
    qvalueCutoff = 0.2 ,               # q-value cutoff
    readable      = FALSE                 # Convert ENSEMBL IDs to gene symbols in output
  )
  
  
  merged_list_go <- as.data.frame(go)
  View( merged_list_go)
  summary( merged_list_go)
  
  #Select the top 20 most significant GO terms based on adjusted p-value
  top_20_go <- merged_list_go %>%
    arrange(p.adjust) %>%
    slice(1:20)
  
  # Create a column for the rich factor (ratio of gene count to term size)
  top_20_go <- top_20_go %>%
    mutate(rich_factor = Count / as.numeric(gsub("/.*", "", BgRatio)))  # Extract numerator of BgRatio
  
  # Generate the lollipop plot
  ggplot(top_20_go, aes(x = reorder(Description, rich_factor), y = rich_factor)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    geom_segment(aes(xend = Description, yend = 0), linetype = "dashed") +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
    scale_size_continuous(name = "Gene Count") +
    coord_flip() +
    labs(
      title = "Top 20 GO Terms by Enrichment",
      x = "GO Term",
      y = "Rich Factor"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  
  library(ggplot2)
  
  
  
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
 
  
  # Step 1: Prepare merged_list_go
  # Split 'geneID' column into separate rows to associate individual genes with GO terms
  go_terms_expanded <- merged_list_go %>%
    separate_rows(geneID, sep = "/") %>% # Split 'geneID' where multiple genes are present
    dplyr::select(geneID, ID, Description)     # Keep relevant columns: ENSG ID, GO Term (ID), and Description
  
  # Step 2: Merge GO terms with merged_data
  # Perform a left join to add GO terms and descriptions based on ENSG IDs
  final_merged_data <- merged_data_egc %>%
    left_join(go_terms_expanded, by = c("ensg_id" = "geneID"))
  
  # Step 3: Optional - Combine multiple GO terms into single rows
  # Aggregate GO Terms and Descriptions for genes that map to multiple GO terms
  final_aggregated_data <- final_merged_data %>%
    group_by(ensg_id) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = "; ")))
  
  # View the resulting dataframe
  print(head(final_aggregated_data))
  View(final_aggregated_data)
  
  # Step 4: Export the data if needed
  write.csv(final_aggregated_data, "merged_data_with_go_terms.csv", row.names = FALSE)
  

  
dataPE <- read_excel("C:/Users/adith/Downloads/HPA DATA.xlsx", sheet = "PE")
View(dataPE)  
colnames(dataPE)
dataPE <- dataPE %>%
  rename(ensg_id = ENSG)

# Step 2: Perform a left join
# Add 'Gene' and 'PE' columns from dataPE to final_aggregated_data
final_merged_with_PE <- final_aggregated_data %>%
  left_join(dataPE, by = "ensg_id")

# Step 3: View the result
print(head(final_merged_with_PE))
View(final_merged_with_PE)

# Step 4: Export the final dataframe if needed
write.csv(final_merged_with_PE, "final_data_with_PE.csv", row.names = FALSE)


# final_merged_with_PE <- final_merged_with_PE %>%
#   mutate(
#     `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` = as.numeric(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm`),
#     `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` = as.numeric(`fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm`),
#     `Consensus Kidney (nTPM)` = as.numeric(`Consensus Kidney (nTPM)`)
#   )
# 
# 
# # Filter rows where both F15 and F18 are greater than the adult value
# final_merged_with_PE_g <- final_merged_with_PE %>%
#   filter(
#     `fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` > `Consensus Kidney (nTPM)`  |
#       `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` > `Consensus Kidney (nTPM)`
#   )
# 
# # Filter rows where at least one of the fetal kidney nTPM values is >= 1
# final_merged_with_PE_final <- final_merged_with_PE %>%
#   filter(`fetal_fetal kidney_Fetal__Kidney_12276_15_ntpm` >= 1 | 
#            `fetal_fetal kidney_Fetal__Kidney_2643_18.4_ntpm` >= 1)
# View(final_merged_with_PE_final)

 

final_merged_with_PE <- final_merged_with_PE %>%
  mutate(ratio_kidney_greaterthan_2 = ifelse(ratio_kidney > 2, TRUE, FALSE))

final_merged_with_PE <- final_merged_with_PE %>%
  mutate(F15_F18_greaterthan_2 = ifelse(F15_F18 > 2, TRUE, FALSE))

final_merged_with_PE <- final_merged_with_PE %>%
  mutate(F18_F15_greaterthan_2 = ifelse(F18_F15 > 2, TRUE, FALSE))

# View the updated dataframe
View(final_merged_with_PE)
colnames(final_merged_with_PE)
write.csv(final_merged_with_PE, "Final_Data.csv", row.names = FALSE)
