#PART1:-------------------DATA ANALYSIS STRATEGY-------------------------

# 1. LOADING DATA AND LIBRARIES

library ("DESeq2")
library("SummarizedExperiment")
library("apeglm")
library ("EnhancedVolcano")
library("tidyverse")
library("pheatmap")

setwd("C:/BI/1st Sem/Next Generation Sequencing/Project2")
getwd()

counts <- read.delim("counts.txt")
metadata <- read.delim("metadata.txt")

# 2. VIEWING DATA

View(counts)
View(metadata)

dim(counts)
dim(metadata)

# Check for missing values
sum(is.na(counts))
sum(is.na(metadata))

# 3. FILTER LOW-COUNT GENES

# Keep genes with at least 10 counts in at least 10 samples
keep <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep, ]
dim(counts_filtered)

# 4. CREATE DESeq2 OBJECT

dds = DESeqDataSetFromMatrix(counts_filtered,metadata,~Category)

# 5. VARIANCE STABILIZING TRANSFORMATION

vsd <- vst(dds)
#View Data
head(assay(vsd))

# 6. PCA FOR OUTLIER DETECTION
plotPCA(vsd, intgroup = c("Category"))

#CODE BY CHATGPT
pcaData <- plotPCA(vsd, intgroup="Category", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color=Category)) +
  geom_point(size=3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Crohn Disease VS NonInflammtory Bowel Disease Control") +
  theme_bw()

colData(dds)
str(colData(dds))
#CODE BY CHATGPT ENDED




# Enhanced version with custom colors and shapes
pca_plot_enhanced <- ggplot(pcaData, aes(PC1, PC2, color=Category, shape=Category)) +
  geom_point(size=4, alpha = 0.6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis: Crohn's Disease vs Controls") +
  scale_color_manual(values = c("red", "blue")) +  # Custom colors
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

print(pca_plot_enhanced)
ggsave("PCA_plot_enhanced.png", width = 9, height = 7, dpi = 300, bg = "white")



# Look at the data
head(pcaData)

# 7. IDENTIFY OUTLIERS STATISTICALLY
library(dplyr)

# Calculate group centroids
centroids <- pcaData %>%
  group_by(Category) %>%
  summarise(center_PC1 = mean(PC1), center_PC2 = mean(PC2))

# Calculate distance from centroid for each sample
pcaData <- pcaData %>%
  left_join(centroids, by = "Category") %>%
  mutate(distance = sqrt((PC1 - center_PC1)^2 + (PC2 - center_PC2)^2))

# Identify outliers (>3 SD from group center)
distance_stats <- pcaData %>%
  group_by(Category) %>%
  summarise(
    mean_dist = mean(distance),
    sd_dist = sd(distance),
    threshold = mean_dist + 3 * sd_dist
  )

pcaData <- pcaData %>%
  left_join(distance_stats, by = "Category") %>%
  mutate(is_outlier = distance > threshold)

# Get outlier names
outliers <- pcaData %>% filter(is_outlier) %>% pull(name)

cat("Outliers identified:", length(outliers), "\n")
print(outliers)


#----------------PART2: DESEQ ANALYSIS---------------------------

# 1. Relevel to make "normal" appear first in the data
dds$Category <-relevel (dds$Category, ref = "non inflammatory bowel disease control")

# 2. Run DESeq2
dds <- DESeq(dds)

#Sanity Check to ensure success of above code
resultsNames (dds)


# 3. Build result tables
res <- results(dds)

#mcols = matrix columns #Observe imp info from results
mcols(res, use.names = TRUE)


# 4.  First, create the extreme_sig object
extreme_sig <- subset(res, pvalue < 10e-8 & abs(log2FoldChange) > 3)

# 5. Count genes
extreme_up <- sum(extreme_sig$log2FoldChange > 0, na.rm = TRUE)
extreme_down <- sum(extreme_sig$log2FoldChange < 0, na.rm = TRUE)

cat("=== EXTREMELY STRINGENT RESULTS ===\n")
cat("Significant genes:", nrow(extreme_sig), "\n")
cat("Up-regulated:", extreme_up, "\n")
cat("Down-regulated:", extreme_down, "\n")

# Now create the volcano plot with the count in caption
volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                pCutoff = 10e-8,        # Extremely stringent
                                FCcutoff = 3,           # Extremely stringent
                                pointSize = 3.0,        # Larger points
                                labSize = 4.0,          # Larger labels
                                labFace = 'bold',       # Bold labels
                                boxedLabels = TRUE,     # Nice label boxes
                                drawConnectors = TRUE,  # Connector lines
                                widthConnectors = 0.75,
                                colConnectors = 'grey30',
                                max.overlaps = 20,
                                
                                # Proper titles and labels
                                title = 'Differential Gene Expression in Crohn\'s Disease',
                                subtitle = 'Extremely Stringent Analysis: Pediatric Ileal Biopsies',
                                caption = paste0('Cutoffs: |log2FC| > 3 (≥8-fold change), p-value < 1e-8\n',
                                                 'Total significant genes: ', nrow(extreme_sig)),
                                xlab = bquote(~Log[2]~ 'Fold Change'),
                                ylab = bquote(~-Log[10]~ 'P-value'),
                                
                                # Better colors
                                col = c('grey70', 'blue', 'green', 'red2'),
                                colAlpha = 0.8,
                                legendPosition = 'bottom',
                                legendLabSize = 12,
                                legendIconSize = 3.0) +
  
  # Additional ggplot customization
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)
  )

# Display the plot
print(volcano_plot)

# Save as high-quality image
ggsave("Final_Volcano_Plot.png", 
       plot = volcano_plot,
       width = 10, 
       height = 8, 
       dpi = 300,
       bg = "white")
