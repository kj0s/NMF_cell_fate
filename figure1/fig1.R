# === 1. Prepare RNA matrix for NMF ===
library(Seurat)
library(NMF)
library(dplyr)
library(ggplot2)

# Use the RNA assay and extract the data slot
rna_matrix <- GetAssayData(ST227, assay = "RNA", slot = "data")

# NMF requires non-negative values
# If log-normalized, exponentiate and subtract 1
rna_matrix <- exp(as.matrix(rna_matrix)) - 1
rna_matrix[rna_matrix < 0] <- 0

# Optionally, subset to highly variable genes
hvg <- VariableFeatures(ST227)
rna_matrix <- rna_matrix[hvg, ]

# === 2. Run NMF ===
# Choose rank = number of latent programs/factors (e.g., 10)
nmf_res <- nmf(rna_matrix, rank = 10, method = "brunet", nrun = 30, seed = 1234)

# Extract factor matrices
W <- basis(nmf_res)   # genes × factors
H <- coef(nmf_res)    # factors × cells

# === 3. Add NMF factors as a DimReduc object ===
ST227[["NMF"]] <- CreateDimReducObject(
  embeddings = t(H),  # cells × factors
  key = "NMF_",
  assay = "RNA"
)

# === 4. Run UMAP and clustering based on NMF factors ===
ST227 <- RunUMAP(ST227, reduction = "NMF", dims = 1:10)
ST227 <- FindNeighbors(ST227, reduction = "NMF", dims = 1:10)
ST227 <- FindClusters(ST227, resolution = 0.5)

# === 5. Visualize UMAP with clusters ===
DimPlot(ST227, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Clusters from NMF factors")
DimPlot(ST227, reduction = "umap", group.by = "celltype") + ggtitle("Annotated cell types")

# === 6. Link NMF factors to barcodes / fate ===
# H matrix: factors × cells
H_df <- as.data.frame(t(H))  # cells × factors
H_df$cell <- rownames(H_df)

# Merge with barcode metadata (example if ST227@meta.data contains barcodes)
H_df <- left_join(H_df, ST227@meta.data %>% rownames_to_column("cell"), by = "cell")

# Visualize factor usage across donors or experimental conditions
library(reshape2)
H_long <- melt(H_df, id.vars = c("cell", "donor_id", "celltype"),
               variable.name = "Factor", value.name = "Score")

ggplot(H_long, aes(x = Factor, y = Score, fill = donor_id)) +
  geom_boxplot() +
  facet_wrap(~celltype) +
  theme_minimal() +
  ggtitle("NMF factor usage by donor and cell type")

# === 7. Optional: Heatmap of top genes per factor ===
library(ComplexHeatmap)
top_genes <- apply(W, 2, function(x) names(sort(x, decreasing = TRUE))[1:20])
top_genes_flat <- unique(as.vector(top_genes))
Heatmap(W[top_genes_flat, ], name = "W", cluster_rows = TRUE, cluster_columns = TRUE)
