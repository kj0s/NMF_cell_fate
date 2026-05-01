# similarity to mofa 2, check if we can replicate findigs w mofa if not similar

clone barcodes done at bulk level, nmf needs 3 modalities at the same resolution. 
meaqn used for barcodes here again ,same as nmf. 
-- 
mofa provides a Variance Explained table or heatmap. It explicitly tells you: "Factor 1 explains 20% of the RNA variation but 0% of the Methylation variation,"
 
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

`Xe` `Xm` are two matrices : abundance x timepoint and celltype x timepoint  most probably. `Xb` is X[i,j,k]. 
U is latent representation
V is how to reconstruct real data from U. 

basically asking if we can compress the data into a form s.t.
```Xb = complex version of U + Vb0 + Vb1``` 
from where you can reconstruct the data.

the fit() function just guesses some U, V, uses the current guess:
`Xe_hat = U @ Ve.T`
`Xm_hat = U @ Vm.T`

then it measures the error by checking `||Xe - Xe_hat||^2`

then U and V are updated. this means they update multiplicatively; the numerator pulls U towards explaining data, denominator prevents explosion of size. 

# end result

at the end, U is the main output, each row is now a learned feature vector. 

V tells you how each hidden factor influences each dataset

and the aim is to rebuild the data well. 
__If the number is 0.0, it means no cell with that specific barcode was found to have turned into that cell type. if you have a value X in D7 for some barcode Y, it means that at day 7, a normalised count of X cells with that barcode turned into that cell type.__

## understand how the multiplicative updates work:

the code performs the following:
for some value of `U[i,j]` it asks, *if i increase this val, will it make a closer reconstruction of my data?* aka does it minimise loss term. 

for `U *= (lamda_e*(Xe @ Ve) + lamda_m*(Xm @ Vm) + ...) / (...)`, we see the term `Xe @ Ve`. 
1. Xe is the data (n x features)
2. Ve is factor matrix (features x k)
each entry is hence (n x k), therefore each i,j tells us how well does factor j align with sample i in the data.

-> if `(n x k)[i,j]` is large, this is explained well w the data, and since this is in the numerator, it increases U. 
denominator has the loss and redundancy terms.
