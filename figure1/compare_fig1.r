# =======================
# 1. Load libraries
# =======================
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(NMF)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(plyr)

# =======================
# 2. Load data
# =======================
setwd("/vast/projects/ST227/ST227/R/")
ST227 <- readRDS("ST227_Day3_21.rds")

# =======================
# 3. Preprocess & multimodal neighbors
# =======================
ST227 <- FindMultiModalNeighbors(
  ST227, reduction.list = list("pca", "sketch_apca"), 
  dims.list = list(1:27, 1:19), modality.weight.name = list("RNA.weight", "ADT.weight")
)

ST227 <- RunUMAP(ST227, nn.name = "weighted.nn", n.neighbors = 30L, 
                 reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ST227 <- FindClusters(ST227, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

# =======================
# 4. Remove contaminating clusters and rename
# =======================
ST227 <- subset(ST227, idents = c("19", "29"), invert = TRUE)

new.cluster.ids <- c("Ery", "cDC2", "Ery", "Ery", "CD14+ Mono", "T cells", "Pre-B cells",
                     "Mast Cells", "Neutrophils", "Megakaryocytes", "Megakaryocytes", "DC3",
                     "Eosinophils", "HSC", "NK cells", "Basophils", "Mast Cells", "cDC1",
                     "CD16+ Mono", "mregDC", "GMP", "cDC2", "NK cells", "Macrophages",
                     "T cells", "cDC2", "Megakaryocytes", "Eosinophils")

ST227$celltype <- plyr::mapvalues(Idents(ST227), from = levels(Idents(ST227)), to = new.cluster.ids)
Idents(ST227) <- ST227$celltype

# =======================
# 5. SingleR annotation
# =======================
ST227.sce.fig <- as.SingleCellExperiment(ST227, assay = "sketch_RNA")
ref <- fetchReference("novershtern_hematopoietic", "2024-02-26")

pred.main <- SingleR(test = ST227.sce.fig, ref = ref, labels = ref$label.main)
ST227[["SingleR.labels"]] <- pred.main$labels

# =======================
# 6. NMF analysis
# =======================
# Use RNA assay
rna_matrix <- GetAssayData(ST227, assay = "RNA", slot = "data")
rna_matrix <- exp(as.matrix(rna_matrix)) - 1
rna_matrix[rna_matrix < 0] <- 0

# Optional: subset to HVG
hvg <- VariableFeatures(ST227)
rna_matrix <- rna_matrix[hvg, ]

# Run NMF (rank = 10)
nmf_res <- nmf(rna_matrix, rank = 10, method = "brunet", nrun = 30, seed = 1234)

W <- basis(nmf_res)   # genes × factors
H <- coef(nmf_res)    # factors × cells

# Add NMF factors to Seurat object
ST227[["NMF"]] <- CreateDimReducObject(
  embeddings = t(H),
  key = "NMF_",
  assay = "RNA"
)

# UMAP & clustering on NMF
ST227 <- RunUMAP(ST227, reduction = "NMF", dims = 1:10, reduction.name = "nmf.umap")
ST227 <- FindNeighbors(ST227, reduction = "NMF", dims = 1:10)
ST227 <- FindClusters(ST227, resolution = 0.5)

# =======================
# 7. Visualization
# =======================
# Compare UMAPs side by side
library(patchwork)

p1 <- DimPlot(ST227, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + 
  ggtitle("Annotated clusters (Seurat + SingleR)") + NoLegend()

p2 <- DimPlot(ST227, reduction = 'nmf.umap', label = TRUE, repel = TRUE, label.size = 2.5) + 
  ggtitle("Clusters from NMF factors") + NoLegend()

# Combine plots
p1 + p2

# =======================
# 8. Optional: Heatmap of top genes per NMF factor
# =======================
top_genes <- apply(W, 2, function(x) names(sort(x, decreasing = TRUE))[1:20])
top_genes_flat <- unique(as.vector(top_genes))
Heatmap(W[top_genes_flat, ], name = "NMF W matrix", cluster_rows = TRUE, cluster_columns = TRUE)
