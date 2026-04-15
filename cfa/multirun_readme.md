script function is to integrate 3 data types i.e. barcode lineage data, stored in the tensor, rna expression, and udt/protein markers. 

it uses multi view factor analysis, leans shared latent structures from the data. 

*shared latent structure* is an underlying, observed factor that explains correllations across multiple diff data sources. it maps heterogenous imputs into s common, lower dimensional space to identify shared patterns.

It runs many models, builds consensus factors, aligns and compares models (optimal transport, Hungarian matching). 

**Outputs:**
1. latent factors
2. umap, heatmap graphs
3. cluster structure

**code itself**
he three data types are organised into a tensor; looks like:

```X[i, j, k] =
    abundance of barcode i 
    at timepoint j
    in celltype k
```
then, the data is normalised using ```log1p = log(1+x)```

` size is (5034, 5, 8) `
 
` 5034 barcodes, 5 timepoints, 8 cell types. `

` The numbers (e.g., 626.93) are likely normalised gene expression or UMI counts. If a cell has a value in D7_Ery but not in others, it was captured or detected as an ery cell specifically at Day 7. `

__If the number is 0.0, it means no cell with that specific barcode was found to have turned into that cell type. if you have a value X in D7 for some barcode Y, it means that at day 7, a normalised count of X cells with that barcode turned into that cell type.__
