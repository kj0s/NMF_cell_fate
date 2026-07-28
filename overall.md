Consensus nmf or modified ?
Standard cnmf does not imply it can support tensor. The code has been modified for tensor support. 
Prediction capacity not part of original cnmf/nmf inputs. 
Cnmf clusters only Ve , this clusters Ve,Vm,Vb0, Vb1 as one vector. Can not find this anywhere. 


Overall:
Across every fate-view normalization strategy tested — CLR transform, proportion-normalization + log1p, and plain log1p — the Fate_Tracking view explains close to 0% of variance in every MOFA factor, with the partial exception of plain log1p, which recovers a small amount of signal (1.4% / 0.2% / 0.5% / 0.1% across Factors 1–4, then 0.0% elsewhere).
The central pattern worth interrogating: the two transforms explicitly designed to strip out clone-size (magnitude) effects CLR and proportion-normalization= are exactly the two that drove Fate_Tracking's R2 to zero everywhere. The one transform that preserves magnitude (plain log1p) is the only one that recovers any shared signal at all. Taken together, this is more consistent with MOFA's fate-related factors being driven primarily by total clone size, not differentiation pattern, than with a normalization implementation problem. 


this model is run on rna,adt,and fate data. the fate data has had the least effect on the r2 value. the sample level mismatch is cleared as we are avg expression values per barcode. the rna has been library normalised, then log1p. adt seems to have undergone clr transform. through inspection, I have learnt Looking closer at the fate view, MOFA seems to be differentiating clones based on how many cells they made overall, rather than on their actual fate pattern. Log1p this view doesn't fix that problem as it shrinks the size difference between clones, but doesn't remove it. So two clones with the exact same differentiation pattern (e.g. splitting 50/50 into the same two cell types) can still end up looking very different to MOFA just because one clone expanded into far more cells than the other [explain what kind of visualissation/strategy could have got this output].

The mofa parameters that work best, even if its counterintuitive to mofa documentation ,are 
``` ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_matrix(
    [[mat_rna], [mat_adt], [mat_fate]],
    likelihoods=["gaussian", "gaussian", "gaussian"],
    views_names=["RNA", "ADT", "Fate_Tracking"],
    groups_names=["single_group"],
    samples_names=[barcodes],
    features_names=[feat_rna, feat_adt, fate_features],
)
ent.set_model_options(factors=15, spikeslab_weights=False, ard_weights=False)
ent.set_train_options(convergence_mode="medium", iter=1000, verbose=True, seed=43) ```

some methods were used to try normalize fate more effectively, they produce ‘better’ figures sometimes but the overlap doesn’t get much better overall. Methods tried- clr transform; made fate completely 0 r2 value all over. Method 2: 
```X_raw = np.dstack([df.to_numpy() for df in df_all])
clone_totals = X_raw.reshape(X_raw.shape[0], -1).sum(axis=1, keepdims=True)
clone_totals = np.where(clone_totals == 0, 1, clone_totals)
X_prop = X_raw / clone_totals[:, :, None]
X_norm = np.log1p(X_prop * X_raw.reshape(X_raw.shape[0], -1).mean())

mat_fate_norm = X_norm.reshape(X_norm.shape[0], -1).astype(np.float64)
mat_fate=mat_fate_norm
```

also made fate have zero r2 in every factor. Just log1p the fate matrix was tried, gave the best output so far. 1.40.2/0.5/0.1 and then 0.0 % variance explained everywhere. Expression per day graphs still not as great, not very clear for this. 
I believe hyperperamaters have been optimized best they can, there is a possibility these new normalisations are better for the fate and it isn’t showing up because mofa says they bias larger modalities and fate is much smaller and hence is not explaining much of data. 

