"""
Per-factor fate trajectory heatmaps for a MOFA model.

Two paths are provided — pick the one matching how you trained MOFA:

  PATH A: fate was included as a training VIEW.
          -> pull loadings (W) for the fate view directly, one column per factor.
          -> this is the correct approach if "fate is a feature of each factor",
             i.e. fate already lives inside the trained model.

  PATH B: fate was held out, MOFA trained on RNA/markers only.
          -> correlate each factor's score (Z) against the held-out fate tensor.
          -> this is the correct approach if you're testing whether an
             unsupervised factor happens to align with real fate
             (the same logic as testing WNN clusters against fate).

Both paths:
  - mask to positive associations only (as requested)
  - use an EXPLICIT timepoint order (never alphabetically sort "d7","d10",...)
  - plot one subplot per factor: celltype on y-axis, timepoint on x-axis

NOTE ON MASKING TO POSITIVE ONLY: this hides any factor that is *negatively*
associated with a lineage (a real, common pattern — e.g. a myeloid-priming
factor that actively predicts *against* erythroid output). Consider a
diverging colormap (RdBu_r, vmin=-max|r|, vmax=max|r|) as a second version
before you commit only to the positive-masked one for the actual figure.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests  # pip install statsmodels

# ---- set this to your real, chronological timepoint order ----
TIMEPOINTS_ORDER = ["d7", "d10", "d14", "d21"]  # EDIT to match your data


# ==========================================================
# PATH A: fate as a trained view -> use loadings directly
# ==========================================================
def plot_factor_fate_from_loadings(model, fate_view_name="fate"):
    """
    model: your trained mofapy2 BayesNet-style model object (e.g. ent.model)
    Assumes fate feature names encode "<celltype>_<timepoint>", e.g. "HSC_d7".
    If your naming differs, adjust `parse_name` below.
    """
    view_names = model.data_opts["view_names"]  # confirm this attribute name in your build
    fate_idx = view_names.index(fate_view_name)

    W = model.getExpectations()["W"]["E"][fate_idx]          # (n_fate_features, n_factors)
    feature_names = model.data_opts["features_names"][fate_idx]

    n_factors = W.shape[1]

    def parse_name(name):
        m = re.match(r"(.+)_(d\d+)$", name)
        if not m:
            raise ValueError(f"Feature name '{name}' doesn't match '<celltype>_<timepoint>' pattern")
        return m.group(1), m.group(2)

    celltypes = sorted({parse_name(f)[0] for f in feature_names})

    fig, axes = plt.subplots(1, n_factors, figsize=(4 * n_factors, 4), sharey=True)
    axes = np.atleast_1d(axes)

    for k in range(n_factors):
        mat = pd.DataFrame(np.nan, index=celltypes, columns=TIMEPOINTS_ORDER)
        for i, name in enumerate(feature_names):
            ct, tp = parse_name(name)
            if tp in TIMEPOINTS_ORDER:
                mat.loc[ct, tp] = W[i, k]

        mat_pos = mat.where(mat > 0)  # keep only positive loadings
        sns.heatmap(
            mat_pos, cmap="Reds", vmin=0, ax=axes[k],
            cbar=True, mask=mat_pos.isna(),
        )
        axes[k].set_title(f"Factor {k + 1}")
        axes[k].set_xlabel("Timepoint")
        if k == 0:
            axes[k].set_ylabel("Cell type")

    plt.tight_layout()
    plt.show()


# ==========================================================
# PATH B: fate held out -> correlate Z against fate tensor
# ==========================================================
def plot_factor_fate_from_correlation(Z, X_lineage, celltype_names, fdr_alpha=0.05):
    """
    Z:          (n_clones, n_factors) — MUST be at the SAME resolution as X_lineage
                (i.e. per clone, not per cell — pseudobulk first if needed).
    X_lineage:  (n_clones, n_timepoints, n_celltypes)
    celltype_names: list matching the last axis of X_lineage, in order.
    """
    n_clones_z, n_factors = Z.shape
    n_clones_x, n_timepoints, n_celltypes = X_lineage.shape
    assert n_clones_z == n_clones_x, (
        f"Resolution mismatch: Z has {n_clones_z} rows, X_lineage has {n_clones_x}. "
        "Z and the fate tensor must be indexed by the same clones."
    )
    assert n_timepoints == len(TIMEPOINTS_ORDER), "TIMEPOINTS_ORDER length must match X_lineage's timepoint axis"

    fig, axes = plt.subplots(1, n_factors, figsize=(4 * n_factors, 4), sharey=True)
    axes = np.atleast_1d(axes)

    all_pvals = []  # for FDR correction across the whole grid, done once at the end

    corr_mats, pval_mats = [], []
    for k in range(n_factors):
        corr_mat = np.zeros((n_celltypes, n_timepoints))
        pval_mat = np.ones((n_celltypes, n_timepoints))
        for t in range(n_timepoints):
            for c in range(n_celltypes):
                r, p = pearsonr(Z[:, k], X_lineage[:, t, c])
                corr_mat[c, t] = r
                pval_mat[c, t] = p
        corr_mats.append(corr_mat)
        pval_mats.append(pval_mat)
        all_pvals.extend(pval_mat.ravel())

    # FDR correction across ALL factor x celltype x timepoint tests jointly
    _, qvals_flat, _, _ = multipletests(all_pvals, alpha=fdr_alpha, method="fdr_bh")
    qvals_iter = iter(qvals_flat)

    for k in range(n_factors):
        corr_mat = corr_mats[k]
        qmat = np.array([next(qvals_iter) for _ in range(corr_mat.size)]).reshape(corr_mat.shape)

        # keep only positive AND significant after FDR correction
        sig_pos = (corr_mat > 0) & (qmat < fdr_alpha)
        masked = np.where(sig_pos, corr_mat, np.nan)

        sns.heatmap(
            masked, cmap="Reds", vmin=0, ax=axes[k],
            xticklabels=TIMEPOINTS_ORDER, yticklabels=celltype_names,
            cbar=True, mask=np.isnan(masked),
        )
        axes[k].set_title(f"Factor {k + 1}")
        axes[k].set_xlabel("Timepoint")
        if k == 0:
            axes[k].set_ylabel("Cell type")

    plt.tight_layout()
    plt.show()


# ==========================================================
# USAGE
# ==========================================================
# Path A (fate trained as a view):
# plot_factor_fate_from_loadings(ent.model, fate_view_name="fate")

# Path B (fate held out, testing correlation):
# Z = ent.model.getExpectations()["Z"]["E"]        # MUST be clone-resolution, not cell-resolution
# plot_factor_fate_from_correlation(Z, x_bc_3d, celltype_names=lineages)
