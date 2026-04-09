
library(Seurat)
library(ggplot2)
library(patchwork)
obj <- readRDS("ST223.annotated.rds")


# 3. Check reductions
Reductions(obj)
# optional
p1 <- DimPlot(obj, reduction = "wnn.umap",
              label = TRUE, repel = TRUE) + NoLegend()

p2 <- DimPlot(obj, reduction = "wnn.umap",
              group.by = "seurat_clusters")

p1 + p2
# =========================
# Set assays
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2) %>%
  ScaleData() %>%
  RunPCA(reduction.name = "apca")

# =========================
# WNN graph
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:18),
  modality.weight.name = "RNA.weight"
)

# =========================
# UMAP
obj <- RunUMAP(obj,
               nn.name = "weighted.nn",
               reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_")



# =========================
# Clustering
obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.8)

# visualising
DimPlot(obj, reduction = "wnn.umap",
        label = TRUE, repel = TRUE) + NoLegend()
# --
DimPlot(obj, reduction = "wnn.umap",
        group.by = "seurat_clusters")
# -- which clusters are rna driven, which are protein driven
VlnPlot(obj,
        features = "RNA.weight",
        group.by = "seurat_clusters",
        pt.size = 0.1) + NoLegend()
FeaturePlot(obj,
            reduction = "wnn.umap",
            features = c("CD34", "CD38", "MPO", "GATA1"))
FeaturePlot(obj,
            reduction = "wnn.umap",
            features = "RNA.weight")
p_rna <- DimPlot(obj, reduction = "rna.umap")
p_adt <- DimPlot(obj, reduction = "adt.umap")

p_rna + p_adt
p1 <- DimPlot(obj, reduction = "wnn.umap")
p2 <- DimPlot(obj, reduction = "rna.umap")

p1 + p2
ggsave("WNN_UMAP_clusters.png", width = 6, height = 5)
