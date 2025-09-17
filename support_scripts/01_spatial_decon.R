#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SpatialDecon")
#BiocManager::install(version = "3.19")
#BiocManager::install("NanoStringGeoMxTools")
#install.packages("NanoStringGeoMxTools")
library(SpatialDecon)
browseVignettes("SpatialDecon")
library(SeuratObject)
library(NanoStringGeoMxTools)
library(DESeq2)
#===============================================================================
#                           CONTINUING FROM DESEQ2
#===============================================================================
rld <- rlog(dds, blind = TRUE)
vst <- vst(dds, blind = TRUE)
rlog_counts <- assay(rld)
vst_counts <- assay(vst)
class(vst_counts)
dim(rlog_counts)
nrow(rlog_counts)
head(rlog_counts)
#checking for negprobe
grep("NegProbe", rownames(rlog_counts))
bg_pro = derive_GeoMx_background(norm=rlog_counts, probepool = rep(1, nrow(rlog_counts)), negnames = "NegProbe-WTX")
#bg_pro = derive_GeoMx_background(norm=vst, probepool = rep(1, nrow(vst_counts)), negnames = "NegProbe-WTX")
#vst
# Check the structure of prostate_data
#-------- loading prostate data
load("Prostate_Henry.RData")
write.csv(x = profile_matrix, file = "Prostate_matrix.csv", 
          row.names = TRUE, quote = FALSE)
profile_matrix <- as.matrix(profile_matrix)
str(profile_matrix)

# perform spatial decon
decon_results <- spatialdecon(norm = rlog_counts,
                              bg = bg_pro,
                              X = profile_matrix,
                              align_genes = TRUE)
str(decon_results)
heatmap(decon_results$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

# Trying to remove the NAs....
# Handle missing values (replace NA with row/column mean, for example)
rlog_counts[is.na(rlog_counts)] <- rowMeans(rlog_counts, na.rm = TRUE)
rlog_counts[is.na(rlog_counts)] <- colMeans(rlog_counts, na.rm = TRUE)
decon_results <- spatialdecon(norm = rlog_counts, bg = bg_pro, X = safeTME, align_genes = TRUE)
heatmap(decon_results$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))
# does not change results.
# Check for NA or zero values in the beta matrix
any(is.na(decon_results$beta))
any(decon_results$beta == 0)
missing_values <- decon_results$beta[c("mDCs", "B.memory", "monocytes.C"), ]
print(missing_values)
library(pheatmap)

# Define a color palette
color_palette <- colorRampPalette(c("blue", "green", "yellow"))(50)

# Plot the heatmap with a custom color palette
pheatmap(decon_results$beta, 
         labels_row = rownames(decon_results$beta), 
         main = "Spatial Deconvolution Heatmap",
         fontsize_row = 10, 
         fontsize_col = 8,
         color = color_palette)
#---- clustering by high vs low
colData(dds)$U6ATAC_RNAscopeClass <- colData$U6ATAC_RNAscopeClass
colData(dds)

# Get cluster labels
cluster_labels <- cutree(cluster_results, k = 2)  # Assuming you want 2 clusters (high vs low U6atac)

# Step 4: Visualization
# Visualize the clustered results, e.g., heatmap
# You can use pheatmap or other heatmap functions for visualization
library(pheatmap)

# Plot heatmap
pheatmap(cluster_data, 
         cluster_cols = FALSE,  # Do not cluster columns (genes)
         cluster_rows = FALSE,  # Do not cluster rows (samples)
         annotation_col = colData(dds)$U6ATAC_RNAscopeClass,  # Color samples by U6atac status
         annotation_colors = list(U6ATAC_RNAscopeClass = c("low" = "blue", "high" = "red")),  # Define colors for U6atac status
         main = "Clustered Heatmap of DESeq Results"
)


# checking alignment ------------
# Align genes between vst_counts and safeTME
aligned_genes <- intersect(rownames(vst_counts), rownames(safeTME$genes))
if (length(aligned_genes) < length(rownames(safeTME$genes))) {
  warning("Not all genes from safeTME are present in vst_counts")
}

