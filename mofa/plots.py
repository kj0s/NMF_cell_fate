Xb_flat = Xb.reshape(n, -1)

Xs = {
    "b": Xb_flat,
    "e": Xe,
    "m": Xm
}

params = {
    "lamda_b": 1.0,
    "lamda_e": 1.0,
    "lamda_m": 1.0,
    "mu_1": 0.1,
    "mu_2": 0.1,
    "mu_lap": 0.0
}
## REPLACE W THE ONES IN THE DOCUMENTATION
model = MOFAStyleModel(k=10, params=params)
model.fit(Xs)

Wb = model.get_weights("b")

plt.figure(figsize=(10, 7.5))

for i in range(model.k):
    plt.subplot(4, 4, i+1)

    mat = Wb[:, i].reshape(len(celltypes), len(timepoints))

    plt.imshow(mat.T)
    plt.colorbar()
    plt.title(f"Factor {i}")

plt.tight_layout()
plt.show()


n = X_bc.shape[0]
X_bc_flat = X_bc.reshape(n, -1)
"b": X_bc_flat 

cfm.fit({
    "b": X_bc_flat,
    "e": X_rna.to_numpy() + 1e-9,
    "m": X_adt.to_numpy() + 1e-9
})

lamda_b = 1 / np.linalg.norm(X_bc_flat)**2
lamda_e = 1 / np.linalg.norm(X_rna)**2
lamda_m = 1 / np.linalg.norm(X_adt)**2
