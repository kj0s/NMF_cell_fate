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
then, the data is normalised using ```log1p = log(1+x)
