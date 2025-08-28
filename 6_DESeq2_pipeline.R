##### Code by Leon Hofmann


# Load required libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(apeglm)
library(ggplot2)
library(gplots)
library(RColorBrewer)

####################################################################################################

# Set working path
setwd('/media/wagner/Data_Storage3/Leon_H/Thesis/RNAseq_Zpa796_Hmurinum/RNAseq/counts/DESeq2')

### Data preparation ###
# Import a count table from featureCounts output and extract gene lengths and renames the columns of the count table to match sample names.
count_table <- read.table("../Quantified_Zpa796_project6448_all_invitro_inplanta.txt", header = TRUE, row.names = "Geneid")
gene_length <- count_table[,5]
names(gene_length) <- rownames(count_table)
count_table <- count_table[, -(1:5)]
rownames(count_table) <- sub("^file_1_file_1_", "", rownames(count_table))

print(colnames(count_table))
print(ncol(count_table))

# Concatenate
count_table <- count_table %>% # select columns with unique part in name that should be concatenated to new column
  mutate(invitro_1 = paste(rowSums(select(., contains('.01')))),
         invitro_2 = paste(rowSums(select(., contains('.02')))),
         invitro_3 = paste(rowSums(select(., contains('.04')))),
         dpi_10_A = paste(rowSums(select(., contains('10dpi.pool2')))),
         dpi_10_B = paste(rowSums(select(., contains('10dpi.pool3')))),
         dpi_10_C = paste(rowSums(select(., contains('10dpi.pool4')))),
         dpi_4_A = paste(rowSums(select(., contains('4dpi.pool1')))),
         dpi_4_B = paste(rowSums(select(., contains('4dpi.pool3')))),
         dpi_4_C = paste(rowSums(select(., contains('4dpi.pool4')))),
         dpi_7_A = paste(rowSums(select(., contains('7dpi.pool1')))),
         dpi_7_B = paste(rowSums(select(., contains('7dpi.pool2')))),
         dpi_7_C = paste(rowSums(select(., contains('7dpi.pool4'))))) %>%
  select(-ends_with('.bam')) %>% # remove all old columns 
  mutate(across(where(is.numeric), as.character)) %>%
  mutate(across(where(is.character), as.numeric)) # convert everything to numeric

print(colnames(count_table))
print(ncol(count_table))

# Experimental conditions  are defined and stored in the condition vector. Sample information is stored in a data matrix.
condition <- c("In_vitro","In_vitro", "In_vitro",
               "DPI_10", "DPI_10", "DPI_10",
               "DPI_04", "DPI_04", "DPI_04",
               "DPI_07", "DPI_07", "DPI_07")

# generate count matrix
count_matrix <- as.matrix(count_table)

# Define color palette using RColorBrewer and assign colors to treatments
num_conditions <- length(unique(condition))
color_palette <- brewer.pal(num_conditions, "Set2")
treatment_colors <- setNames(color_palette, unique(condition))


############################################################################
############################## DESeq2 ######################################
############################################################################

# Create DESeqDataSet object and extract normalized counts
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = data.frame(condition = factor(condition)),
                              design = ~ condition)

dds = estimateSizeFactors(dds)
normalized_dds = counts(dds, normalized = T)

# Regularized-logarithm transformation and sample distances for normalized counts
rld <- rlog(dds)
rlog_dds <- assay(rld)
colnames(rlog_dds)

# Boxplot for normalized counts
boxplot(rlog_dds, cex.axis=1, main="rlog transformed data", outline=TRUE, col=treatment_colors[condition])
legend("topright", legend=unique(condition), fill=treatment_colors[unique(condition)], title="Treatment")
dev.copy2pdf(file="boxplot_rlog_distribution_final.pdf")

# Calculate Euclidean distance matrix for normalized counts and plot distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )

# As heatmap
hm <- pheatmap(sampleDistMatrix, 
               clustering_distance_rows = "euclidean", 
               clustering_distance_cols = "euclidean",
               clustering_method = "complete", 
               fontsize = 10, 
               main = "Normalized counts",
               labels_col = colnames(count_table), # Add sample names on the right
               labels_row = colnames(count_table), # Add row names on the right
               color = colorRampPalette(c("purple4", "plum3", "lavenderblush"))(100), # Customize colors
               show_colnames = TRUE, # Show sample names on the heatmap
               show_rownames = TRUE, # Hide row names
               cellwidth = 20, # Adjust cell width
               cellheight = 20, # Adjust cell height
               border_color = NA, # Remove grey line separating samples
               angle_col = 45) # Rotate column names by -45 degrees
dev.copy2pdf(file="count_matrix_heatmap_final.pdf")

# As PCA
colData(rld)$condition <- recode(colData(rld)$condition, "DPI_04" = "4 dpi", "DPI_07" = "7 dpi", "DPI_10" = "10 dpi", "In_vitro" = "in vitro")
plotPCA(rld)
dev.copy2pdf(file="plotPCA_final.pdf")

### Set factor levels and reference and perform DE analysis ###
dds$condition <- factor(dds$condition, levels = c("In_vitro","DPI_10", "DPI_04", "DPI_07"))
dds$condition <- relevel(dds$condition, ref = "In_vitro")
dds$condition 

####################################################################################################

### Differential expression test - LRT(Likelihood ratio test) ###
dds <- DESeq(dds, test="LRT", reduced=~1)
resultsNames(dds) 

# Inspect the estimated dispersions
plotDispEsts(dds)
dev.copy2pdf(file="plotDispEsts_final.pdf")

# Filtering low count genes  
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
colnames(dds)

#write.table(normalized_dds, "counts_Normalized_Zpa796_invitro_inplanta.csv", sep = "\t")

# Calculate TMP:
# Normalize for gene length (reads per kilobase - RPK)
gene_length <- as.matrix(gene_length)
x_count_table <- count_table/c(gene_length)

# Normalize to read depth (here depth = colSums of the count matrix):
tpm <- t(t(x_count_table)*1e6/colSums(x_count_table))
write.table(tpm, "tpm_counts_Zpa796_invitro_inplanta_DEseq2.csv", sep = "\t", col.names = NA)

########################################### CONTRASTS ###############################################

####################################
######## DPI_10 vs In_vitro ########
####################################
### DPI_10

# Test the differentially expressed genes between contrasted samples.
res_DPI_10_vs_In_vitro <- results(dds, contrast = c("condition", "DPI_10", "In_vitro"), test = "Wald", alpha = 0.05, lfcThreshold = 2)
summary(res_DPI_10_vs_In_vitro)

# Differentially expressed genes are visualized using a MA plot
plotMA(res_DPI_10_vs_In_vitro, ylim=c(-8,8), main= "condition_DPI_10_vs_In_vitro", alpha=0.05)
abline(h=c(-4,4),col='green3')
dev.copy2pdf(file='plotMA_DPI_10_005.pdf')

# Results of the analysis, including differentially expressed genes and diagnostic plots, are written to output files.
write.csv(as.data.frame(res_DPI_10_vs_In_vitro), file = "Expressed_Zpa796_DPI_10-In_Vitro_DEseq2_005.csv")

####################################
######## DPI_04 vs In_vitro ########
####################################
### DPI_04

# Test the differentially expressed genes between contrasted samples.
res_DPI_04_vs_In_vitro <- results(dds, contrast = c("condition", "DPI_04", "In_vitro"), test = "Wald", alpha = 0.05, lfcThreshold = 2)
summary(res_DPI_04_vs_In_vitro)

# Differentially expressed genes are visualized using a MA plot
plotMA(res_DPI_04_vs_In_vitro, ylim=c(-8,8), main= "condition_DPI_04_vs_In_vitro", alpha=0.05)
abline(h=c(-4,4),col='green3')
dev.copy2pdf(file='plotMA_DPI_04_005.pdf')

# Results of the analysis, including differentially expressed genes and diagnostic plots, are written to output files.
write.csv(as.data.frame(res_DPI_04_vs_In_vitro), file = "Expressed_Zpa796_DPI_04-In_Vitro_DEseq2_005.csv")

####################################
######## DPI_07 vs In_vitro ########
####################################
### DPI_07

# Test the differentially expressed genes between contrasted samples.
res_DPI_07_vs_In_vitro <- results(dds, contrast = c("condition", "DPI_07", "In_vitro"), test = "Wald", alpha = 0.05, lfcThreshold = 2)
summary(res_DPI_07_vs_In_vitro)

# Differentially expressed genes are visualized using a MA plot
plotMA(res_DPI_07_vs_In_vitro, ylim=c(-8,8), main= "condition_DPI_07_vs_In_vitro", alpha=0.05)
abline(h=c(-4,4),col='green3')
dev.copy2pdf(file='plotMA_DPI_07_005.pdf')

# Results of the analysis, including differentially expressed genes and diagnostic plots, are written to output files.
write.csv(as.data.frame(res_DPI_07_vs_In_vitro), file = "Expressed_Zpa796_DPI_07-In_Vitro_DEseq2_005.csv")

summary(res_DPI_04_vs_In_vitro)
summary(res_DPI_07_vs_In_vitro)
summary(res_DPI_10_vs_In_vitro)


