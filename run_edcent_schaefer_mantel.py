"""
Mantel test: pairwise ED_centroid (mm) vs pairwise beta_r similarity (Schaefer-400).
Reuses cached cope_schaefer400.npy, lesion_array_4mm.npz (centroids), and slim PKL.
"""
import sys
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
import config
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial.distance import cdist
from scipy.stats import rankdata
import glob

INT_X     = Path(config.INT_DIR) / 'x_distance_analysis'
LESION_SRC = INT_X / 'lesion_masks' / 'NiftiLesions_original'
COPE_DIR   = Path(config.INT_DIR) / '2_Correlates' / 'allsessions_glm_zscored' / 'acute'

# ── 1. Reconstruct common_ids ─────────────────────────────────────────────────
lesion_files = sorted(glob.glob(str(LESION_SRC / 'sub-PAT*_A_lesion_ifhfix.nii.gz')))
lesion_ids   = [Path(f).name.replace('.nii.gz','').split('_A_')[0].replace('sub-','')
                for f in lesion_files]
cope_files   = sorted(glob.glob(str(COPE_DIR / 'PAT*_cope_fe.nii.gz')))
cope_set     = {Path(f).name.replace('.nii.gz','').replace('_cope_fe','') for f in cope_files}
common_ids   = [pid for pid in lesion_ids if pid in cope_set]
N            = len(common_ids)
print(f'N = {N}')

# ── 2. Load centroids → ED_centroid matrix ────────────────────────────────────
d          = np.load(str(INT_X / 'lesion_array_4mm.npz'), allow_pickle=True)
centroids  = d['centroids']   # (N, 3) MNI mm
valid      = ~np.any(np.isnan(centroids), axis=1)
print(f'Valid centroids: {valid.sum()}/{N}')

ED_centroid = cdist(centroids, centroids, 'euclidean').astype(np.float32)
ED_centroid[~valid, :] = np.nan
ED_centroid[:, ~valid] = np.nan
np.fill_diagonal(ED_centroid, np.nan)

# ── 3. Load cached Schaefer cope → pairwise Pearson r ─────────────────────────
cope_schaef = np.load(str(INT_X / 'cope_schaefer400.npy'))   # (N, 400)
print(f'cope_schaefer shape: {cope_schaef.shape}')

Xc    = cope_schaef - np.nanmean(cope_schaef, axis=1, keepdims=True)
norms = np.linalg.norm(Xc, axis=1, keepdims=True)
norms[norms == 0] = 1.0
Xn    = Xc / norms
beta_r_schaef = (Xn @ Xn.T).astype(np.float32)
np.fill_diagonal(beta_r_schaef, np.nan)

tri = np.tril_indices(N, -1)
print(f'beta_r_schaefer mean={np.nanmean(beta_r_schaef[tri]):.4f}  std={np.nanstd(beta_r_schaef[tri]):.4f}')
print(f'ED_centroid     mean={np.nanmean(ED_centroid[tri]):.1f} mm  std={np.nanstd(ED_centroid[tri]):.1f} mm')

# ── 4. Log-volume distance (partial Mantel covariate) ─────────────────────────
meta     = pd.read_pickle(config.PKL_SLIM)
acute    = meta[meta['Session'] == 'acute'].drop_duplicates('ID').set_index('ID')
log_vols = np.array([
    np.log(float(acute.loc[pid, 'lesion_volume_mm3'])) if pid in acute.index else np.nan
    for pid in common_ids
], dtype=np.float32)
vol_dist = np.abs(log_vols[:, None] - log_vols[None, :]).astype(np.float32)
np.fill_diagonal(vol_dist, np.nan)

# ── 5. Mantel functions ───────────────────────────────────────────────────────
def mantel_spearman(rdm_x, rdm_y, n_perm=5000, seed=42):
    tri_   = np.tril_indices(rdm_x.shape[0], -1)
    x_vals = rdm_x[tri_]
    y_vals = rdm_y[tri_]
    valid  = np.isfinite(x_vals) & np.isfinite(y_vals)
    xv, yv = x_vals[valid], y_vals[valid]
    obs_r  = float(np.corrcoef(rankdata(xv), rankdata(yv))[0, 1])
    rng    = np.random.default_rng(seed)
    perm_r = np.array([
        float(np.corrcoef(rankdata(rng.permutation(xv)), rankdata(yv))[0, 1])
        for _ in range(n_perm)
    ])
    p = ((np.abs(perm_r) >= np.abs(obs_r)).sum() + 1) / (n_perm + 1)
    return obs_r, p, len(xv)

def partial_mantel(rdm_x, rdm_y, rdm_cov, n_perm=5000, seed=42):
    tri_   = np.tril_indices(rdm_x.shape[0], -1)
    x, y, z = rdm_x[tri_], rdm_y[tri_], rdm_cov[tri_]
    valid   = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = rankdata(x[valid]), rankdata(y[valid]), rankdata(z[valid])
    def partial_r(a, b, c):
        rab = np.corrcoef(a, b)[0, 1]
        rac = np.corrcoef(a, c)[0, 1]
        rbc = np.corrcoef(b, c)[0, 1]
        denom = np.sqrt((1 - rac**2) * (1 - rbc**2))
        return (rab - rac * rbc) / denom if denom > 0 else 0.0
    obs_r  = partial_r(x, y, z)
    rng    = np.random.default_rng(seed)
    perm_r = np.array([partial_r(rankdata(rng.permutation(x)), y, z) for _ in range(n_perm)])
    p = ((np.abs(perm_r) >= np.abs(obs_r)).sum() + 1) / (n_perm + 1)
    return obs_r, p

# ── 6. Run ────────────────────────────────────────────────────────────────────
print('\nRunning Mantel tests (5000 permutations)...')
r_raw,  p_raw,  n_pairs = mantel_spearman(ED_centroid, beta_r_schaef)
r_part, p_part          = partial_mantel(ED_centroid, beta_r_schaef, vol_dist)

print(f'ED_centroid × beta_r_schaefer  raw:     r={r_raw:+.4f}  p={p_raw:.4f}  (n_pairs={n_pairs})')
print(f'ED_centroid × beta_r_schaefer  partial: r={r_part:+.4f}  p={p_part:.4f}')
