# NMF_py

* **Why use Non-negative Matrix Factorization (NMF) for this data?**
  NMF is well-suited to single-cell RNA data because gene expression is inherently non-negative and often reflects additive biological processes (e.g. pathways, activation programs). Unlike PCA, which mixes positive and negative contributions, NMF forces each cell to be represented as a combination of only positive components. This makes the resulting factors easier to interpret biologically as “programs” or “modules” of co-expressed genes.

This is valuable because we are trying to relate transcriptional state to lineage barcodes. NMF provides a more mechanistic representation of state (which programs are active) rather than just grouping cells into clusters. This allows finer interpretation of whether a lineage corresponds to a single program or a mixture of programs.

---

* **What steps does NMF involve that are different from the previous clustering method?**
  The current workflow uses PCA followed by graph-based clustering. PCA works by projecting the data into orthogonal components that maximise variance, and it allows both positive and negative values. This is then used to construct a nearest-neighbour graph for clustering.

  NMF instead factorises the expression matrix into two matrices: one representing gene weights per factor and the other representing factor usage per cell. This is done iteratively using optimisation to minimize reconstruction error under a non-negativity constraint. There is no concept of orthogonality, and the components are not ordered by variance but by how well they reconstruct the data collectively.

  Practically, this means we skip scaling (to avoid negative values), run NMF on normalized data, and then use the resulting factor matrix in place of PCA embeddings when building the neighbour graph and performing clustering.

---

* **What outputs will we get if we use NMF here?**
  NMF produces two main outputs. The first is a gene-by-factor matrix (W), which tells you which genes define each latent program. Each column can be interpreted as a transcriptional program or pathway. The second is a factor-by-cell matrix (H), which tells you how strongly each program is expressed in each cell.

  From this, we obtain a cell embedding (derived from H) that can be used for UMAP and clustering, just like PCA embeddings. The advantage here is one can directly examine which genes drive each factor and how cells mix these factors. This allows you to describe cells as both belonging to a cluster, and as having graded contributions from multiple biological programs.

---

* **How would we compare the outputs of PCA-based and NMF-based methods?**
structural: We will compare UMAP visualisations generated from PCA versus NMF embeddings to see whether the same cell populations separate similarly or whether new structure emerges. We will also compare clustering results (e.g. overlap between cluster assignments) to assess whether NMF identifies more refined or biologically coherent groups.


you can compare how well each method aligns with barcode-defined lineages: measure whether cells sharing the same barcode are more tightly grouped in PCA space or show more consistent factor usage in NMF space. If NMF better captures lineage-specific programs or reveals divergence within lineages, it provides stronger insight into fate decisions.

