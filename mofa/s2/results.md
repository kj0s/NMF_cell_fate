Overall: issue why we cant compare top n genes is one factor in NMF has multiple, connected best-matches in mofa. hence, taking the top n from any factors would not lead you to accurately represent the diversity in both MOFA na NMF. After hyperparameter validation, the model now has .. efficiency. this is trained on the same daa with the same no. of barcodes and variance of the previous model. 

Provided the model works well, I believe we must find a better technique to evaluate the genes given importance by each model. Moreover, when I saved my nmf gene list, I found only 1714 as compared to the total 4674 we started off with. There was aggressive pruning somewhere in the model, or some things naturally tended to zero. this is worth exploring. 

Another thing of note is not all of the genes deemed important in the nmf model are present in the mofa model. only around 60% of those genes are present. This means the zero variance genes i removed from mofa are part of the nmf model? 4674 genes for mofa, 4787 original genes. retrying the set overlap without removing the zero variance genes - 1052 overlap. increase of 35. still missing 78 genes? original variance used is **0.40** this was found using the following code block:

```
# done to try find the 113 missing genes in total mofa count. 
df_rna_copy = pd.read_csv(
    "/vast/projects/Sisseq/human-haematopoiesis-sis-seq/data/preprocessed_filtered_hsc/X_rna_variance_0.40.csv",
    index_col=0
)
feat_rna_with_var  = df_rna_copy.columns.tolist()
len(feat_rna_with_var)
len(feat_rna)
common_genes = sorted(set(nmf_genes) & set(feat_rna_with_var))
print(f"Common genes: {len(common_genes)} (NMF had {len(nmf_genes)}, MOFA RNA had {len(feat_rna_with_var)})")
```

*remaining genes found* the final gene list seems to be coming from a non hsc dataset. mofa seems to be run on the overall data here instead of just hscs. thats why the gene list showed such difference. we will rerun nmf wth the hsc dataset and resave, retry. All genes found! run at 0.0, with preprocessed, and that  is the entire list. The NMF needs to be rerun with just hsc for a more accurate answer

Some better techniques to evaluate feature importance based on each model are 
